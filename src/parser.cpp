//
// Created by Bohlender,Ryan James on 11/12/20.
//

#include "parser.hpp"
#include "inputvalidator.hpp"
#include <unordered_set>

template <typename T>
void Parser<T>::parse_input(std::istream &is) {
    std::string line;
    unsigned long submitted = 0;
    double cur_dist = 0;
    double last_dist = 0;

    // Initialize ThreadPool
    ThreadPool<Statistic<T>> threadpool(params);

    // Maintain a single vector that we just update with each line
    arma::SpCol<int32_t> last(pheno.samples->size() * (pheno.samples->size() - 1) / 2.);
    arma::SpCol<int32_t> data(pheno.samples->size() * (pheno.samples->size() - 1) / 2.);

    arma::wall_clock timer;
    long lineno = -1;
    std::map<std::string, std::map<std::string, int>> indices;
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

            RJBUtil::Splitter<std::string> hsplit(line, ": ");

            // There are two breakpoint classes to handle. cluster and singleton.
            // IDs that are expected to be present:
            // clusterID, IIDs, hapIDs, IID1-IID2, hapID1-hapID2, cM
            for (auto it = hsplit.begin() + 1; it != hsplit.end(); it++) {
                indices[hsplit[0]][*it] = std::distance(hsplit.begin() + 1, it);
            }

            continue;
        }

        // Legacy format handling
        lineno++;
        if (!params.dash) {
            if (indices.find("singleton") == indices.end()) {
                indices["singleton"] = std::map<std::string, int>();
            }
            if (lineno == 0) {// Skip the header
                if (indices["singleton"].find("ID1-ID2") == indices["singleton"].end()) {
                    indices["singleton"]["ID1-ID2"] = 1;
                }
                if (indices["singleton"].find("cM") == indices["singleton"].end()) {
                    indices["singleton"]["cM"] = 0;
                }
                continue;
            }
        } else {
            if(lineno == 0) {
                continue;
            }
        }

        boost::char_separator<char> sep{"\t"};
        boost::tokenizer<boost::char_separator<char>> tok{line, sep};
        auto tok_it = tok.begin();

        std::string chrom = *(tok_it);
        tok_it++;

        if (!boost::starts_with(chrom, "chr")) {
            chrom = std::string("chr") + chrom;
        }

        int pos = std::stoi(*tok_it);

        // Initialize breakpoint
        Breakpoint bp{};
        bp.breakpoint = std::make_pair(chrom, *tok_it);
        tok_it++;

        if (params.dash) {
            // chr pos n.cluster n.cluster.haplotype n.cluster.pair n.singleton.pair cluster.add cluster.del singleton.add singleton.del
            for (int i = 0; i < 4; i++) {
                tok_it++;
            }
            boost::char_separator<char> inner_sep {" "};
            RJBUtil::Splitter<std::string> cluster_add(*tok_it, " ");
            tok_it++;
            RJBUtil::Splitter<std::string> cluster_del(*tok_it, " ");
            tok_it++;
            RJBUtil::Splitter<std::string> singleton_add(*tok_it, " ");
            tok_it++;
            RJBUtil::Splitter<std::string> singleton_del(*tok_it, " ");

            update_data(data, indices["cluster"], cluster_del, bp, -1, true);
            update_data(data, indices["singleton"], singleton_del, bp, -1, false);
            update_data(data, indices["cluster"], cluster_add, bp, 1, true);
            update_data(data, indices["singleton"], singleton_add, bp, 1, false);
        } else {
            // chr     pos  npairs  nsegments   add     del
            if (params.oldformat) {
                tok_it++;
                tok_it++;
            }
            // chr     pos     add     del
            boost::char_separator<char> inner_sep {" "};
            boost::tokenizer<boost::char_separator<char>> additions {*tok_it, inner_sep};
            update_data(data, indices["singleton"], additions, bp, 1, false);
            tok_it++;
            boost::tokenizer<boost::char_separator<char>> deletions {*tok_it, inner_sep};
            update_data(data, indices["singleton"], deletions, bp, -1, false);
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
            if (check_range(pos)) {
                continue;
            }
        }
        if (params.exclude) {
            if (check_exclude(pos)) {
                fmt::print(std::cerr, "Excluded breakpoint: {}\n", pos);
                continue;
            }
        }

        if (params.rsquared) {
            if(check_r2(data, last)) {
                continue;
            }
        }
        if (arma::accu(data) < params.threshold) {
            continue;
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
                       pheno.indexer,
                       reporter,
                       params,
                       pheno.phenotypes);
        threadpool.submit(std::move(stat));
        submitted++;
    }
    while (threadpool.ntasks > 0) {
        std::this_thread::sleep_for(std::chrono::nanoseconds(100000000));
    }
}

template <typename T>
bool Parser<T>::check_sample_list(const std::string &sample_pair) {
    RJBUtil::Splitter<std::string_view> vals(sample_pair, ":");
    return params.sample_list->find(vals[0]) == params.sample_list->end() && params.sample_list->find(vals[1]) == params.sample_list->end();
}

template <typename T>
bool Parser<T>::check_range(int pos) {
    return pos < (*params.range)[0] || pos > (*params.range)[1];
}

template <typename T>
bool Parser<T>::check_exclude(int pos) {
    bool exclude_region = false;
    for (const auto &v : *params.exclude) {
        exclude_region |= pos >= v.first && pos <= v.second;
    }
    return exclude_region;
}

template <typename T>
bool Parser<T>::check_r2(const arma::SpCol<int32_t> &data, const arma::SpCol<int32_t> &last) {
    double r2 = cor(data, last);
    if (params.verbose) {
        std::cerr << "r2: " << r2 << std::endl;
    }
    return r2 > *params.rsquared;
}

template <typename T>
void Parser<T>::update_data(arma::SpCol<int32_t> &data,
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
            if (info) {
                if ((*info).filter_segment(vals[indices["clusterID"]], params)) {
                    continue;
                }
            }

            RJBUtil::Splitter<std::string> iids(vals[indices[iid_key]], ",");

            std::sort(iids.begin(), iids.end());

            for (auto it1 = iids.begin(); it1 != iids.end(); it1++) {
                for (auto it2 = it1 + 1; it2 != iids.end(); it2++) {
                    arma::sword row_idx = (*pheno.indexer).translate(*it1, *it2);

                    int id_idx = indices["clusterID"];
                    auto ids = fmt::format("{},{}", *it1, *it2);
                    if (row_idx < 0) {
                        continue;
                    } else {
                        if (value > 0) {
                            bp.ibd_pairs.emplace_back(*it1, *it2);
                        }
                    }
                    data(row_idx) += value;
                }
            }
        } else {
            RJBUtil::Splitter<std::string_view> vals(entry, ":");
            RJBUtil::Splitter<std::string> pairs(vals[indices[iid_key]], "-");

            std::sort(pairs.begin(), pairs.end());
            arma::sword row_idx = (*pheno.indexer).translate(pairs[0], pairs[1]);

            auto ids = fmt::format("{},{}", pairs[0], pairs[1]);
            if (row_idx < 0) {
                continue;
            } else {
                if (params.cM > 0) {
                    try {
                        double test_val = std::stod(vals[indices["cM"]]);
                        if (test_val < params.cM) {
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
                    bp.ibd_pairs.emplace_back(pairs[0], pairs[1]);
                }
            }
            data(row_idx) += value;
        }
    }
}

template <typename T>
void Parser<T>::update_data(arma::SpCol<int32_t> &data,
                         std::map<std::string, int> &indices,
                         boost::tokenizer<boost::char_separator<char>> &changes,
                         Breakpoint &bp,
                         int value,
                         bool cluster) {
    std::string iid_key;
    std::unordered_set<std::string> seen;
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
                    arma::sword row_idx = (*pheno.indexer).translate(*it1, *it2);

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
                            bp.ibd_pairs.emplace_back(*it1, *it2);
                        }
                    }
                    data(row_idx) += value;
                }
            }
        } else {
            RJBUtil::Splitter<std::string_view> vals(entry, ":");
            RJBUtil::Splitter<std::string> pairs(vals[indices[iid_key]], "-");

            std::sort(pairs.begin(), pairs.end());
            if (pheno.lookup.find(pairs[0]) == pheno.lookup.end() || pheno.lookup.find(pairs[1]) == pheno.lookup.end()) {
                continue;
            }
            arma::sword row_idx = (*pheno.indexer).translate(pairs[0], pairs[1]);

            auto ids = fmt::format("{},{}", pairs[0], pairs[1]);
            if (seen.find(ids) == seen.end()) {
                seen.insert(ids);
            } else {
                continue;
            }
            if (row_idx < 0) {
                continue;
            } else {
                auto test = (*pheno.indexer).back_translate(row_idx);
                if ((*pheno.samples)[test.first] != pairs[0] || (*pheno.samples)[test.second] != pairs[1]) {
                    fmt::print(std::cerr, "Failed to map back. {} == {} ; {} == {}\n", (*pheno.samples)[test.first], pairs[0], (*pheno.samples)[test.second], pairs[1]);
                }
                if (params.cM) {
                    try {
                        double test_val = std::stod(vals[indices["cM"]]);
                        if (test_val < params.cM) {
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
                    bp.ibd_pairs.emplace_back(pairs[0], pairs[1]);
                }
            }
            data(row_idx) += value;
        }
    }
}

template <typename T>
Parser<T>::Parser(Parameters params_, std::shared_ptr<Reporter> reporter_, GeneticMap &gmap_, Phenotypes<T> &pheno_)
    : params(std::move(params_)), reporter(std::move(reporter_)), gmap(std::move(gmap_)), pheno(pheno_) {
    // Default construct shared ptrs
    if (params.info) {
        Source info_source(*params.info);
        std::istream info_stream(&(*info_source.streambuf));
        info = Info(info_stream, params);
    }

    Source input_source(params.input);
    std::istream input_s(&(*input_source.streambuf));

    if (params.verbose) {
        std::cerr << "Parsing data\n";
    }
    parse_input(input_s);
}

template class Parser<pheno_vector>;
template class Parser<compressed_pheno_vector>;
