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

double
Statistic::calculate(std::vector<int> &phenotypes_, double cscs_count, double cscn_count, double cncn_count, int k) {
    std::vector<std::string> &samples = parser.samples;
    double statistic;

    double ibd_pairs = arma::accu(data);

    double cscs = 0;
    double cscn = 0;
    double cncn = 0;

    // Save time allocating the full amount
    if (pairs.empty()) {
        rows.reserve(ibd_pairs);
        pairs.reserve(ibd_pairs);
        for (auto it = data.begin(); it != data.end(); it++) {
            auto [p1, p2] = indexer[0].back_translate(it.row());
            pairs.emplace_back(std::make_pair(p1, p2));
            rows.push_back(it.row()); // Store rows for later lookup
            switch (phenotypes_[p1]) {
                case 1:
                    switch (phenotypes_[p2]) {
                        case 1:
                            cscs += 1.;
                            break;
                        case 0:
                            cscn += 1.;
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
                            cscn += 1.;
                            break;
                        case 0:
                            cncn += 1.;
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
    } else {
        for (int i = 0; i < pairs.size(); i++) {
            const auto &p = pairs[i];
            const auto &r = rows[i];
            switch (phenotypes_[p.first]) {
                case 1:
                    switch (phenotypes_[p.second]) {
                        case 1:
                            cscs += 1.;
                            break;
                        case 0:
                            cscn += 1.;
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
                            cscn += 1.;
                            break;
                        case 0:
                            cncn += 1.;
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
    }
    permuted_cscs[k].push_back(cscs);
    permuted_cscn[k].push_back(cscn);
    permuted_cncn[k].push_back(cncn);

    statistic = cscs / cscs_count - cscn / cscn_count;

    return statistic;
}

void Statistic::run() {
    unsigned long k = 0;
    for (; k < indexer.size(); k++) {
        permuted.emplace_back();
        permuted_cscs.emplace_back();
        permuted_cscn.emplace_back();
        permuted_cncn.emplace_back();
        original.push_back(
                calculate(indexer[k].phenotypes, indexer[k].case_case, indexer[k].case_cont, indexer[k].cont_cont, 0));
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
            for (k = 0; k < phenotypes.size(); k++) { // Permute all phenotypes the same way.
                int tmp = phenotypes[k][j];
                phenotypes[k][j] = phenotypes[k][i];
                phenotypes[k][i] = tmp;
            }
        }

        k = 0;
        for (auto &v : phenotypes) {
            double val = calculate(v, indexer[k].case_case, indexer[k].case_cont, indexer[k].cont_cont, k);
            if (val > original[k])
                successes[k]++;
            permutations[k]++;
            permuted[k].push_back(val);
            k++;
        }
    }

    // Output
    std::stringstream iss;
    for (int k = 0; k < permuted_cscs.size(); k++) {
        iss << bp.breakpoint.first << "\t" << bp.breakpoint.second;
        for (int i = 0; i < params.nperms + 1; i++) { // nperms + 1 because the original values are also in the permuted set
            iss << "\t" << permuted_cscs[k][i] << "\t" << permuted_cscn[k][i] << "\t" << permuted_cncn[k][i];
        }
        iss << std::endl;
    }
    reporter->submit(iss.str());
    done = true;
}

