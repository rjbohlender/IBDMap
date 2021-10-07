//
// Created by Bohlender,Ryan James on 11/12/20.
//

#include "parser.hpp"
#include "inputvalidator.hpp"

void Parser::count_breakpoints(std::istream &is) {
    nbreakpoints = 0;
    std::string line;
    while (std::getline(is, line)) {
        if (boost::starts_with(line, "#") || boost::starts_with(line, "chr"))
            continue;
        nbreakpoints++;
    }
    if (params.verbose) {
        std::cerr << "total breakpoints: " << nbreakpoints << std::endl;
    }
}

void Parser::parse_input(std::istream &is) {
    std::string line;
    unsigned long submitted = 0;
    double cur_dist = 0;
    double last_dist = 0;

    // Initialize ThreadPool
    ThreadPool<Statistic> threadpool(params);

    // Maintain a single vector that we just update with each line
    arma::sp_vec last(samples->size() * (samples->size() - 1) / 2.);
    arma::sp_vec data(samples->size() * (samples->size() - 1) / 2.);

    arma::wall_clock timer;
    long lineno = -1;
    std::map<std::string, int> indices;
    long long total_delta = 0;
    while (std::getline(is, line)) {
        if (params.dash && boost::starts_with(line, "##")) {
            // Processing Header
            // Example:
            // ####cluster: clusterID:IIDs:hapIDs
            // ##singleton: IID1-IID2:hapID1-hapID2:cM
            size_t header_chars = 0;
            for (const auto c : line) {
                if (c == '#') {
                    header_chars++;
                }
            }
            line = line.substr(header_chars, line.size());// Strip the leading header characters.

            RJBUtil::Splitter<std::string> hsplit(line, ":");

            // There are two breakpoint classes to handle. cluster and singleton.
            // IDs that are expected to be present:
            // clusterID, IIDs, hapIDs, IID1-IID2, hapID1-hapID2, cM
            for (auto it = hsplit.begin() + 1; it != hsplit.end(); it++) {
                indices[*it] = std::distance(hsplit.begin() + 1, it);
            }

            continue;
        }

        // Legacy format handling
        lineno++;
        if (lineno == 0) {// Skip the header
            if (indices.find("ID1-ID2") == indices.end()) {
                indices["ID1-ID2"] = 1;
                indices["cM"] = 0;
            }
            continue;
        }

        RJBUtil::Splitter<std::string_view> splitter(line, "\t", true);

        std::string chrom;

        if (!boost::starts_with(splitter[0], "chr")) {
            std::stringstream chrss;
            chrss << "chr" << splitter[0];
            chrom = chrss.str();
        }

        int pos = std::stoi(splitter[1]);

        // Initialize breakpoint
        Breakpoint bp{};
        bp.breakpoint = std::make_pair(chrom, splitter[1]);

        if (params.dash) {
            // chr pos n.cluster n.cluster.haplotype n.cluster.pair n.singleton.pair cluster.add cluster.del singleton.add singleton.del
            RJBUtil::Splitter<std::string> cluster_add(splitter[splitter.size() - 4], " ");
            RJBUtil::Splitter<std::string> cluster_del(splitter[splitter.size() - 3], " ");
            RJBUtil::Splitter<std::string> singleton_add(splitter[splitter.size() - 2], " ");
            RJBUtil::Splitter<std::string> singleton_del(splitter[splitter.size() - 1], " ");

            update_data(data, indices, cluster_del, bp, -1, true);
            update_data(data, indices, singleton_del, bp, -1, false);
            update_data(data, indices, cluster_add, bp, 1, true);
            update_data(data, indices, singleton_add, bp, 1, false);
        } else {
            RJBUtil::Splitter<std::string> additions(splitter[splitter.size() - 2], " ");
            RJBUtil::Splitter<std::string> deletions(splitter[splitter.size() - 1], " ");

            update_data(data, indices, additions, bp, 1, false);
            update_data(data, indices, deletions, bp, -1, false);
        }

        // Must follow data update -- Data format is just change from previous breakpoint so skipping update is impossible
        std::pair<std::pair<int, double>, std::pair<int, double>> nearest = gmap.find_nearest(chrom, pos);
        if (nearest.first.first != nearest.second.first) {
            cur_dist = (pos - nearest.first.first) * (nearest.second.second - nearest.first.second) / (nearest.second.first - nearest.first.first) + nearest.first.second;
        } else {
            cur_dist = nearest.first.second;
        }

        if (cur_dist - last_dist < params.min_dist) {
            continue;
        }

        if (params.range) {
            if (pos < (*params.range)[0] || pos > (*params.range)[1]) {
                continue;
            }
        }
        if (params.exclude) {
            bool exclude_region = false;
            for (const auto &v : *params.exclude) {
                exclude_region |= pos >= v.first && pos <= v.second;
            }
            if (exclude_region) {
                fmt::print(std::cerr, "{}\n", pos);
                continue;
            }
        }

        if (params.rsquared) {
            double r2 = cor(data, last);
            if (params.verbose) {
                std::cerr << "r2: " << r2 << std::endl;
            }
            if (r2 > *params.rsquared) {
                continue;
            }
        }

#ifndef NDEBUG
        if (params.verbose) {
            std::cerr << "dist: " << cur_dist - last_dist << std::endl;
            std::cerr << "fill rate: " << arma::accu(data) / data.n_elem << std::endl;
        }
#endif
        last_dist = cur_dist;
        last = data;

        Statistic stat(data,
                       bp,
                       indexer,
                       reporter,
                       params,
                       groups);
        threadpool.submit(std::move(stat));
        submitted++;
    }
    while (threadpool.ntasks > 0) {
        std::this_thread::sleep_for(std::chrono::nanoseconds(100000000));
    }
}

bool Parser::check_sample_list(const std::string &sample_pair) {
    RJBUtil::Splitter<std::string_view> vals(sample_pair, ":");
    return params.sample_list->find(vals[0]) == params.sample_list->end() && params.sample_list->find(vals[1]) == params.sample_list->end();
}

void Parser::update_data(arma::sp_vec &data,
                         std::map<std::string, int> &indices,
                         RJBUtil::Splitter<std::string> &changes,
                         Breakpoint &bp,
                         int value,
                         bool cluster) {
    std::string iid_key;
    if (params.dash) {
        if (cluster) {
            iid_key = "IIDs";
        } else {
            iid_key = "IID1-IID2";
        }
    } else {
        iid_key = "ID1-ID2";
    }
    for (auto &entry : changes) {
        if (entry == "NA") {
            break;
        }
        if (params.sample_list) {
            if (check_sample_list(entry)) {
                continue;
            }
        }


        if (params.dash && cluster) {
            RJBUtil::Splitter<std::string> vals(entry, ":");
            RJBUtil::Splitter<std::string> iids(vals[indices[iid_key]], ",");

            std::sort(iids.begin(), iids.end());

            for (auto it1 = iids.begin(); it1 != iids.end(); it1++) {
                for (auto it2 = it1 + 1; it2 != iids.end(); it2++) {
                    arma::sword row_idx = (*indexer)[0].translate(*it1, *it2);

                    int id_idx = indices["clusterID"];
                    auto ids = fmt::format("{},{}", *it1, *it2);
                    if (row_idx < 0) {
                        continue;
                    } else {
                        if (info) {
                            if ((*info).filter_segment(vals[indices["clusterID"]], params)) {
                                continue;
                            }
                        }
                        if (value > 0) {
                            bp.ibd_pairs.emplace_back(std::make_pair(*it1, *it2));
                        }
                    }
                    data(row_idx) += value;
                }
            }

        } else {
            RJBUtil::Splitter<std::string_view> vals(entry, ":");
            RJBUtil::Splitter<std::string> pairs(vals[indices[iid_key]], "-");

            std::sort(pairs.begin(), pairs.end());
            arma::sword row_idx = (*indexer)[0].translate(pairs[0], pairs[1]);

            auto ids = fmt::format("{},{}", pairs[0], pairs[1]);
            if (row_idx < 0) {
                continue;
            } else {
                if (params.cM) {
                    try {
                        double test_val = std::stod(vals[indices["cM"]]);
                        if (test_val < *params.cM) {
                            continue;
                        }
                    } catch (std::invalid_argument &e) {
                        throw(std::runtime_error("Attempted to filter on segment length, in data that lacks segment lengths."));
                    }
                }
                if (value > 0) {
                    try {
                        bp.segment_lengths.push_back(std::stod(vals[indices["cM"]]));
                    } catch (std::invalid_argument &e) {
                        bp.segment_lengths.push_back(nan("1"));
                    }
                    bp.ibd_pairs.emplace_back(std::make_pair(pairs[0], pairs[1]));
                }
            }
            data(row_idx) += value;
        }
    }
}

void Parser::parse_pheno(std::istream &is) {
    int iid = 0;
    int phe = 1;
    std::string line;
    arma::uword lineno = 0;
    std::map<std::vector<bool>, std::vector<arma::uword>> fill_patterns;

    while (std::getline(is, line)) {
        if (boost::starts_with(line, "#") || lineno == 0) {// Skip the header
            lineno++;
            continue;
        }
        RJBUtil::Splitter<std::string_view> splitter(line, " \t");

        samples->push_back(splitter[iid]);
        std::vector<bool> pattern;
        for (int i = 1; i < splitter.size(); i++) {
            if (phenotypes.size() < i) {
                phenotypes.emplace_back();
            }
            if (splitter[i] == "NA") {
                pattern.push_back(false);
                phenotypes[i - 1].push_back(-1);
            } else {
                pattern.push_back(true);
                phenotypes[i - 1].push_back(std::stoi(splitter[i]));
                if (phenotypes[i - 1].back() != 0 && phenotypes[i - 1].back() != 1) {
                    if (params.verbose) {
                        fmt::print(std::cerr, "{} {}\n", splitter[0], splitter[i]);
                    }
                }
                if (params.swap) {// Swap case-control status
                    switch (phenotypes[i - 1].back()) {
                        case 1:
                            phenotypes[i - 1].back() = 0;
                            break;
                        case 0:
                            phenotypes[i - 1].back() = 1;
                            break;
                        default:
                            break;
                    }
                }
            }
        }
        if (fill_patterns.count(pattern) == 0) {
            fill_patterns.emplace(std::make_pair(pattern, std::vector<arma::uword>({lineno - 1})));
        } else {
            fill_patterns[pattern].push_back(lineno - 1);
        }
        lineno++;
    }
    for (const auto &v : phenotypes) {
        int case_count = 0;
        int control_count = 0;
        for (const auto &p : v) {
            switch (p) {
                case 1:
                    case_count++;
                    break;
                case 0:
                    control_count++;
                    break;
                default:
                    break;
            }
        }

        indexer->emplace_back(Indexer(case_count, control_count, (*samples), v));
    }
    if (fill_patterns.size() > 1) {
        groups = std::vector<std::vector<arma::uword>>();
        for (const auto &v : fill_patterns) {
            groups->push_back(v.second);
        }
        if (params.verbose) {
            std::cerr << "Groups: " << groups->size() << std::endl;
            std::cerr << "Group sizes: ";
        }
        for (const auto &v : fill_patterns) {
            if (params.verbose) {
                std::cerr << v.second.size() << " ";
            }
        }
        if (params.verbose) {
            std::cerr << std::endl;
        }
    }
}

Parser::Parser(const std::string &input_path, const std::string &pheno_path, Parameters params_, std::shared_ptr<Reporter> reporter_, GeneticMap &gmap_)
    : nbreakpoints(0), params(std::move(params_)), reporter(std::move(reporter_)), gmap(std::move(gmap_)) {
    // Default construct shared ptrs
    samples = std::make_shared<std::vector<std::string>>();
    phenotypes = std::vector<std::vector<int>>();
    indexer = std::make_shared<std::vector<Indexer>>();
    if (params.info) {
        Source info_source(*params.info);
        std::istream info_stream(&(*info_source.streambuf));
        info = Info(info_stream, params);
    }

    Source bp_source(input_path);
    std::istream bp_is(&(*bp_source.streambuf));
    std::ifstream pheno_ifs(pheno_path);

    if (params.verbose) {
        std::cerr << "Counting breakpoints\n";
    }
    count_breakpoints(bp_is);

    Source input_source(input_path);
    std::istream input_s(&(*input_source.streambuf));

    if (params.verbose) {
        std::cerr << "Parsing phenotypes\n";
    }
    parse_pheno(pheno_ifs);

    for (const auto &idx : (*indexer)) {
        if (params.verbose) {
            std::cerr << "ncases: " << idx.case_count << " ncontrols: " << idx.cont_count << std::endl;
        }
    }

    if (params.verbose) {
        std::cerr << "Parsing data\n";
    }
    parse_input(input_s);
}
