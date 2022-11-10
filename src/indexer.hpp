//
// Created by Bohlender,Ryan James on 12/13/19.
//

#ifndef CARVAIBD_SRC_INDEXER_HPP
#define CARVAIBD_SRC_INDEXER_HPP

#define ARMA_DONT_USE_WRAPPER

#include "types.hpp"
#include "indexsort.hpp"
#include <armadillo>
#include <unordered_map>
#include <utility>
#include <vector>

/**
 * @brief Class to handle converting individual indices into row index
 *
 * Invertible mapping function between sample pairs and rows. Pairs are given
 * as a single index, while pairs of individual indices correspond to the coordinates
 * in the cartesian product of the samples, excluding matches to self. Autozygous
 * samples are dropped currently.
 */
struct Indexer {
    // Class counts and categories
    arma::uword case_count;
    arma::uword cont_count;
    arma::uword sz;
    std::vector<std::string> samples;
    pheno_vector phenotypes;
    std::vector<arma::sword> transitions;// Transition points between pairing sets
    std::unordered_map<std::string, int> ordered_positions;
    std::unordered_map<std::string, int> ordered_cc_positions;

    // Pairs
    arma::uword case_case;
    arma::uword case_cont;
    arma::uword cont_cont;

    Indexer() = delete;

    Indexer(arma::uword case_count_,
            arma::uword cont_count_,
            std::vector<std::string> samples_,
            pheno_vector phenotypes_)
        : case_count(case_count_), cont_count(cont_count_), sz(samples_.size()), phenotypes(std::move(phenotypes_)),
          samples(std::move(samples_)), case_case(0), case_cont(0), cont_cont(0) {
        setup(case_count_, cont_count_);
#ifndef NDEBUG
        test_run();
#endif
    }

    void setup(arma::uword case_count_,
               arma::uword cont_count_) {
        IndexSort indexSort(samples);
        indexSort.sort_vector(samples);
        indexSort.sort_vector(phenotypes);// Both must be sorted
        case_case = case_count * (case_count - 1.) / 2.;
        case_cont = case_count * cont_count;
        cont_cont = cont_count * (cont_count - 1.) / 2.;

        int case_idx = 0;
        int cont_idx = 0;
        int current = 0;

        arma::sword start = 0;
        arma::sword k = samples.size() - 1;
        for (; k > 0; k--) {
            start += k;
            transitions.push_back(start);
        }

        for (const auto &s : samples) {
            switch (phenotypes[current]) {
                case 1:
                    ordered_positions[s] = current;
                    ordered_cc_positions[s] = case_idx;
                    case_idx++;
                    break;
                case 0:
                    ordered_positions[s] = current;
                    ordered_cc_positions[s] = cont_idx;
                    cont_idx++;
                    break;
                case -1:
                    break;
                default:
                    throw(std::runtime_error("Invalid phenotype value."));
            }
            current++;
        }
    }

    arma::sword translate(const std::string &a, const std::string &b) {
        auto finda = std::lower_bound(samples.begin(), samples.end(), a);
        auto findb = std::lower_bound(samples.begin(), samples.end(), b);
        if (*finda != a || *findb != b || a == b) {
            return -1;
        }
        int i = std::distance(samples.begin(), finda);
        int j = std::distance(samples.begin(), findb);
        if (i == j) {
            throw(std::runtime_error("Searching for two of the same pair."));
        }
        if (i > j) {
            std::swap(i, j);
        }

        arma::sword start = 0;
        // Replace loop with function for the value -- 1 - 2 orders of magnitude faster
        long double sz_ = sz;
        start = (2. * sz_ - 1. - i) * i / 2.;

        return start + j - i - 1;
    }

    std::pair<int32_t, int32_t> back_translate(arma::uword row) {
        arma::sword i = -1;
        arma::sword j = -1;

        auto bound = std::lower_bound(transitions.begin(), transitions.end(), row);

        i = std::distance(transitions.begin(), bound);
        if (row == *bound) {
            i += 1;
            j = i + 1;
        } else {
            j = i > 0 ? row - *(bound - 1) + i + 1 : row + 1;
        }

        return std::make_pair(i, j);
    }

    static void test_run() {
        std::vector<std::string> test_samples{
                "A", "B", "C", "D", "E", "F", "G", "H", "I", "J"};

        std::vector<int> test_phenotypes{
                1, 1, 1, 1, 1, 0, 0, 0, 0, 0};

        arma::sword N = 10;
        std::vector<int> test_transitions;
        int start = 0;
        for (int k = N - 1; k > 0; k--) {
            start += k;
            test_transitions.push_back(start);
        }

        auto test_translate = [&](std::string &a, std::string &b) {
            auto finda = std::lower_bound(test_samples.begin(), test_samples.end(), a);
            auto findb = std::lower_bound(test_samples.begin(), test_samples.end(), b);
            if (*finda != a || *findb != b) {
                return -1ll;
            }
            int i = std::distance(test_samples.begin(), finda);
            int j = std::distance(test_samples.begin(), findb);
            std::cerr << "i: " << i << " j: " << j << std::endl;
            if (i == j) {
                throw(std::runtime_error("Searching for two of the same pair."));
            }
            if (i > j) {
                auto tmp = i;
                i = j;
                j = tmp;
            }

            arma::sword start = 0;
            arma::sword k = test_samples.size() - 1;
            for (; k > test_samples.size() - i - 1; k--)
                start += k;
            return start + j - i - 1;
        };

        auto test_back_translate = [&](arma::sword row) {
            arma::sword i = -1;
            arma::sword j = -1;

            auto bound = std::lower_bound(test_transitions.begin(), test_transitions.end(), row);

            i = std::distance(test_transitions.begin(), bound);
            if (row == *bound) {
                i += 1;
                j = i + 1;
            } else {
                j = i > 0 ? row - *(bound - 1) + i + 1 : row + 1;
            }
            std::cerr << "bound: " << *bound << std::endl;

            return std::make_pair(i, j);
        };

        // Testing A:B; Correct answer: 0
        auto row = test_translate(test_samples[0], test_samples[1]);
        std::cerr << "Row returned: " << row << " Correct row: " << 0 << std::endl;
        row = test_translate(test_samples[0], test_samples[2]);
        std::cerr << "Row returned: " << row << " Correct row: " << 1 << std::endl;
        row = test_translate(test_samples[0], test_samples[3]);
        std::cerr << "Row returned: " << row << " Correct row: " << 2 << std::endl;
        row = test_translate(test_samples[1], test_samples[2]);
        std::cerr << "Row returned: " << row << " Correct row: " << 9 << std::endl;
        row = test_translate(test_samples[8], test_samples[9]);
        std::cerr << "Row returned: " << row << " Correct row: " << 44 << std::endl;

        auto p = test_back_translate(0);
        std::cerr << "Pair returned: " << p.first << "," << p.second << " Correct pair: "
                  << "0,1" << std::endl;
        p = test_back_translate(44);
        std::cerr << "Pair returned: " << p.first << "," << p.second << " Correct pair: "
                  << "8,9" << std::endl;
        p = test_back_translate(9);
        std::cerr << "Pair returned: " << p.first << "," << p.second << " Correct pair: "
                  << "1,2" << std::endl;
        p = test_back_translate(10);
        std::cerr << "Pair returned: " << p.first << "," << p.second << " Correct pair: "
                  << "1,3" << std::endl;
        p = test_back_translate(5);
        std::cerr << "Pair returned: " << p.first << "," << p.second << " Correct pair: "
                  << "0,6" << std::endl;
    }
};

#endif//CARVAIBD_SRC_INDEXER_HPP
