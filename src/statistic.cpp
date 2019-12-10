//
// Created by Bohlender,Ryan James on 9/23/19.
//

#include "statistic.hpp"

Statistic::Statistic(arma::sp_colvec &&data_,
                     Breakpoint &bp_,
                     std::vector<Indexer> &indexer_,
                     Parser<std::string> &parser_,
                     std::shared_ptr<Reporter> reporter_,
                     Parameters &params_) :
        data(std::move(data_)), indexer(indexer_), parser(parser_), params(params_), bp(bp_),
        reporter(std::move(reporter_)) {
}

double Statistic::calculate(std::vector<int> &phenotypes_, double case_case, double case_cont, double cont_cont) {
    std::vector<std::string> &samples = parser.samples;
    double statistic;

    double ibd_pairs = arma::accu(data);

    double case_case_rate = 0;
    double case_cont_rate = 0;
    double cont_cont_rate = 0;

    // Save time allocating the full amount
    if (pairs.empty()) {
        rows.reserve(ibd_pairs);
        pairs.reserve(ibd_pairs);
        for (auto it = data.begin(); it != data.end(); it++) {
#ifdef LOWMEM
            auto[p1, p2] = indexer[0].back_translate_alt(it.row());
#else
            auto [p1, p2] = indexer.back_translate(it.row());
#endif
            pairs.emplace_back(std::make_pair(p1, p2));
            rows.push_back(it.row()); // Store rows for later lookup
            switch (phenotypes_[p1]) {
                case 1:
                    switch (phenotypes_[p2]) {
                        case 1:
                            case_case_rate += 1. / case_case;
                            break;
                        case 0:
                            case_cont_rate += 1. / case_cont;
                            break;
                        case -1:
                            break;
                        default:
                            throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
                    }
                    break;
                case 0:
                    switch (phenotypes_[p2]) {
                        case 1:
                            case_cont_rate += 1. / case_cont;
                            break;
                        case 0:
                            cont_cont_rate += 1. / cont_cont;
                            break;
                        case -1:
                            break;
                        default:
                            throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
                    }
                    break;
                case -1:
                    break;
                default:
                    throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
            }
        }
        permuted_cscs.push_back(case_case_rate);
        permuted_cscn.push_back(case_cont_rate);
    } else {
        for (int i = 0; i < pairs.size(); i++) {
            const auto &p = pairs[i];
            const auto &r = rows[i];
            switch (phenotypes_[p.first]) {
                case 1:
                    switch (phenotypes_[p.second]) {
                        case 1:
                            case_case_rate += 1. / case_case;
                            break;
                        case 0:
                            case_cont_rate += 1. / case_cont;
                            break;
                        case -1:
                            break;
                        default:
                            throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
                    }
                    break;
                case 0:
                    switch (phenotypes_[p.second]) {
                        case 1:
                            case_cont_rate += 1. / case_cont;
                            break;
                        case 0:
                            cont_cont_rate += 1. / cont_cont;
                            break;
                        case -1:
                            break;
                        default:
                            throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
                    }
                    break;
                case -1:
                    break;
                default:
                    throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
            }
        }
        permuted_cscs.push_back(case_case_rate);
        permuted_cscn.push_back(case_cont_rate);
    }
    statistic = case_case_rate - case_cont_rate;
    return statistic;
}

void Statistic::run() {
    unsigned long k = 0;
    for (; k < indexer.size(); k++) {
        original.push_back(
                calculate(indexer[k].phenotypes, indexer[k].case_case, indexer[k].case_cont, indexer[k].cont_cont));
        successes.push_back(0);
        permutations.push_back(0);
    }

    std::mt19937 gen(params.seed);

    arma::wall_clock timer;
    auto &phenotypes = parser.phenotypes;
    while (permutations[k - 1] < params.nperms) {
        for (unsigned long i = phenotypes[0].size() - 1; i > 0; i--) {
            std::uniform_int_distribution<> dis(0, i);
            int j = dis(gen);
            for (k = phenotypes.size() - 1; k > 0; k--) {
                int tmp = phenotypes[k][j];
                phenotypes[k][j] = phenotypes[k][i];
                phenotypes[k][i] = tmp;
            }
        }

        k = 0;
        for (auto &v : phenotypes) {
            double val = calculate(v, indexer[k].case_case, indexer[k].case_cont, indexer[k].cont_cont);
            if (val > original[k])
                successes[k]++;
            permutations[k]++;
            permuted[k].push_back(val);
            k++;
        }
    }

    // Output
    std::stringstream iss;
    k = 0;
    for (const auto &v : permuted) {
        iss << bp.breakpoint.first << "\t" << bp.breakpoint.second << "\t" << original[k];
        for (int i = 0; i < params.nperms; i++) {
            iss << "\t" << v[i];
        }
        iss << std::endl;
        k++;
    }
    reporter->submit(iss.str());
    done = true;
}

