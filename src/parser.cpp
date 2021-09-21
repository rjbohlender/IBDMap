//
// Created by Bohlender,Ryan James on 11/12/20.
//

#include "parser.hpp"

void Parser::count_breakpoints(std::istream &is) {
    nbreakpoints = 0;
    std::string line;
    while (std::getline(is, line)) {
        if (boost::starts_with(line, "#"))
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

    // Generate permutations if we have covariates
    Permute permute(params.seed);
    std::optional<std::shared_ptr<std::vector<std::vector<int32_t>>>> o_perms;
    if (covariates && phenotypes.size() == 1) {
        generate_cov_adj_perms(permute, o_perms);
    }

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
                       phenotypes,
                       reporter,
                       params,
                       groups,
                       o_perms);
        threadpool.submit(std::move(stat));
        submitted++;
    }
    while (threadpool.ntasks > 0) {
        std::this_thread::sleep_for(std::chrono::nanoseconds(100000000));
    }
}

void Parser::generate_cov_adj_perms(Permute &permute,
                                    std::optional<std::shared_ptr<std::vector<std::vector<int32_t>>>> &o_perms) {
    auto permutation_ptr = std::make_shared<std::vector<std::vector<int32_t>>>();

    arma::vec Y = arma::conv_to<arma::vec>::from(phenotypes[0]);
    Binomial link;
    GLM<Binomial> fit(*covariates, Y, link, params);
    arma::vec odds = fit.mu_ / (1 - fit.mu_);

    std::string lr_path = params.output_path.empty() ? "lr.txt" : params.output_path + ".lr.txt";

    std::ofstream lr_out(lr_path);
    lr_out << "Sample\tProb\tOdds" << std::endl;
    for (int i = 0; i < odds.n_elem; i++) {
        lr_out << cov_samples[i] << "\t" << fit.mu_(i) << "\t" << odds(i) << std::endl;
    }

    if (params.verbose) {
        std::cerr << "LR Output: max pr: " << arma::max(fit.mu_) << " min pr: " << arma::min(fit.mu_) << std::endl;
        std::cerr << "LR Output: max odds: " << arma::max(odds) << " min odds: " << arma::min(odds) << std::endl;
    }
    lr_out.close();

    permute.get_permutations(permutation_ptr, odds, (*indexer)[0].case_count, params.nperms, params.nthreads - 2);
    o_perms = permutation_ptr;
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

    bool notified = false;

    while (std::getline(is, line)) {
        if (boost::starts_with(line, "#") || lineno == 0) {// Skip the header
            lineno++;
            continue;
        }
        RJBUtil::Splitter<std::string_view> splitter(line, " \t");
        if (splitter.size() > 2 && !notified) {
            std::cerr << "Multiple phenotypes provided. Covariates will be ignored." << std::endl;
            notified = true;
        }

        if (!splitter.empty() && skip.find(splitter[0]) != skip.end()) {
            // Skip samples with missing cov values; don't increment lineno because we're treating them as if they don't exist
            continue;
        }

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
                        fmt::print(std::cerr, "{} {}\n", splitter[0], splitter[1]);
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

        indexer->emplace_back(Indexer(case_count, control_count, (*samples), phenotypes.back()));
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
    if (covariates) {
        int i = 0;
        arma::uvec idx(samples->size());
        for (const auto &v : (*samples)) {
            if (cov_samples[i] != v) {
                auto swap_val = std::find(cov_samples.begin(), cov_samples.end(), v);
                if (swap_val == cov_samples.end())
                    throw(std::runtime_error("ERROR: Sample not present in covariate file."));
                int j = static_cast<int>(std::distance(cov_samples.begin(), swap_val));
                idx(i) = j;
                std::swap(cov_samples[i], cov_samples[j]);
            } else {
                idx(i) = i;
            }
            i++;
        }
        covariates = (*covariates).rows(idx);
    }
}

void Parser::parse_cov(std::istream &is) {
    std::string line;
    unsigned long lineno = 0;
    unsigned long nfields = 0;

    std::map<std::string, std::vector<double>> data;
    std::vector<std::vector<std::string>> unconvertible;

    while (std::getline(is, line)) {
        RJBUtil::Splitter<std::string> splitter(line, " \t");
        if (lineno == 0) {// Skip the header
            if (splitter.size() > 0) {
                nfields = splitter.size() - 1;// Not counting the sample field so we don't have to subtract all the time.
                unconvertible.resize(nfields);
                lineno++;
            } else {
                throw(std::runtime_error("Header line of covariate file is empty. Please include a header line."));
            }
            continue;
        }

        if (splitter.empty()) {
            continue;
        }
        std::string sampleid = splitter[0];
        if (splitter.size() < nfields + 1) {// Skip samples where we have a missing column or two
            skip.emplace(sampleid);
            continue;
        }

        cov_samples.push_back(sampleid);
        data[sampleid] = std::vector<double>(nfields, 0);

        for (int i = 0; i < nfields; i++) {
            try {
                data[sampleid][i] = std::stod(splitter[i + 1]);
            } catch (...) {
                unconvertible[i].push_back(splitter[i + 1]);
            }
        }
        lineno++;
    }

    // Handle unconvertible fields by treating them as factors with levels -- convert to dummy variables
    int fieldno = 0;
    int offset = 0;
    for (const auto &field : unconvertible) {
        if (!field.empty()) {
            std::set<std::string> unique(field.begin(), field.end());
            std::map<std::string, int> levels;

            if (params.verbose) {
                std::cerr << "In reading covariates, could not convert column " << fieldno + 1 << " to double." << std::endl;
                std::cerr << "Levels: ";
            }

            for (auto it = unique.begin(); it != unique.end(); it++) {
                levels.emplace(std::make_pair(*it, std::distance(unique.begin(), it)));
                if (params.verbose) {
                    std::cerr << *it << " : " << std::distance(unique.begin(), it) << " ";
                }
            }
            if (params.verbose) {
                std::cerr << std::endl;
            }

            int sampleno = 0;
            int nlevels = levels.size() - 1;
            for (const auto &v : field) {// Convert to dummy variable
                for (int j = 0; j < nlevels; j++) {
                    if (j == 0) {
                        if (j == levels[v]) {
                            data[cov_samples[sampleno]][fieldno + offset] = 1.0;
                        } else {
                            data[cov_samples[sampleno]][fieldno + offset] = 0.0;
                        }
                    } else {
                        if (j == levels[v]) {
                            data[cov_samples[sampleno]].insert(data[cov_samples[sampleno]].begin() + fieldno + offset + j, 1.0);
                        } else {
                            data[cov_samples[sampleno]].insert(data[cov_samples[sampleno]].begin() + fieldno + offset + j, 0.0);
                        }
                    }
                }
                sampleno++;
            }
            offset += nlevels - 1;
            nfields += nlevels - 1;
        }
        fieldno++;
    }

    arma::mat design(cov_samples.size(), nfields + 1);
    if (params.verbose) {
        std::cerr << "Design.n_rows: " << design.n_rows << std::endl;
        std::cerr << "Design.n_cols: " << design.n_cols << std::endl;
    }
    int i = 0;
    for (const auto &s : cov_samples) {
        design(i, 0) = 1;
        int j = 1;
        for (const auto &v : data[s]) {
            design(i, j) = v;
            j++;
        }
        i++;
    }

    covariates = design;
}

Parser::Parser(const std::string &input_path,
               const std::string &pheno_path,
               std::optional<std::string> cov_path,
               Parameters params_,
               std::shared_ptr<Reporter> reporter_,
               GeneticMap &gmap_)
    : nbreakpoints(0), params(std::move(params_)), reporter(std::move(reporter_)), gmap(std::move(gmap_)) {
    // Default construct shared ptrs
    samples = std::make_shared<std::vector<std::string>>();
    phenotypes = std::vector<std::vector<int>>();
    indexer = std::make_shared<std::vector<Indexer>>();
    if (params.info) {
        Source info_source(*params.info);
        std::istream info_stream(&(*info_source.streambuf));
        info = Info(info_stream);
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

    if (cov_path) {
        if (params.verbose) {
            std::cerr << "Parsing covariates\n";
        }

        Source cov_source(*cov_path);
        std::istream cov_is(&(*cov_source.streambuf));

        parse_cov(cov_is);
    }

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
