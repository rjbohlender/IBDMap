//
// Created by Bohlender,Ryan James on 12/13/19.
//

#ifndef CARVAIBD_SRC_INDEXER_HPP
#define CARVAIBD_SRC_INDEXER_HPP

#include <vector>
#include <armadillo>
#include "absl/container/flat_hash_map.h"
#include "IndexSort.hpp"

/**
 * @brief Class to handle converting individual indices into row index
 *
 * Invertible mapping function between sample pairs and rows.
 */
struct Indexer {
  // Class counts and categories
  arma::uword case_count;
  arma::uword cont_count;
  std::vector<int> phenotypes;
  std::vector<std::string> samples;
  std::vector<arma::sword> transitions; // Transition points between pairing sets
  absl::flat_hash_map<std::string, int> ordered_positions;
  absl::flat_hash_map<std::string, int> ordered_cc_positions;
  // std::unordered_map<std::string, int> ordered_positions;
  // std::unordered_map<std::string, int> ordered_cc_positions;


  // Pairs
  arma::uword case_case;
  arma::uword case_cont;
  arma::uword cont_cont;

  Indexer() = default;

  Indexer(arma::uword case_count_,
          arma::uword cont_count_,
          std::vector<std::string> samples_,
          std::vector<int> phenotypes_)
      : case_count(case_count_), cont_count(cont_count_), phenotypes(std::move(phenotypes_)),
        samples(std::move(samples_)) {
    setup(case_count_, cont_count_);
#ifndef NDEBUG
    test_run();
#endif
  }

  void setup(arma::uword case_count_,
             arma::uword cont_count_) {
    IndexSort indexSort(samples);
    indexSort.sort_vector(samples);
    indexSort.sort_vector(phenotypes); // Both must be sorted
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
      case 1:ordered_positions[s] = current;
        ordered_cc_positions[s] = case_idx;
        case_idx++;
        break;
      case 0:ordered_positions[s] = current;
        ordered_cc_positions[s] = cont_idx;
        cont_idx++;
        break;
      default:throw (std::runtime_error("Invalid phenotype value."));
      }
      current++;
    }
  }

  arma::sword translate(std::string a, std::string b) {
    auto finda = std::lower_bound(samples.begin(), samples.end(), a);
    auto findb = std::lower_bound(samples.begin(), samples.end(), b);
    if (*finda != a || *findb != b) {
      return -1;
    }
    int i = std::distance(samples.begin(), finda);
    int j = std::distance(samples.begin(), findb);
    if (i == j) {
      throw (std::runtime_error("Searching for two of the same pair."));
    }
    if (i > j) {
      auto tmp = i;
      i = j;
      j = tmp;
    }

    arma::sword start = 0;
    arma::sword k = samples.size() - 1;
    for (; k > samples.size() - i - 1; k--)
      start += k;
    return start + j - i - 1;
  }

  std::pair<arma::sword, arma::sword> back_translate(arma::uword row) {
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

  void test_run() {
    std::vector<std::string> test_samples{
        "A", "B", "C", "D", "E", "F", "G", "H", "I", "J"
    };

    std::vector<int> test_phenotypes{
        1, 1, 1, 1, 1, 0, 0, 0, 0, 0
    };

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
        throw (std::runtime_error("Searching for two of the same pair."));
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
    std::cerr << "Pair returned: " << p.first << "," << p.second << " Correct pair: " << "0,1" << std::endl;
    p = test_back_translate(44);
    std::cerr << "Pair returned: " << p.first << "," << p.second << " Correct pair: " << "8,9" << std::endl;
    p = test_back_translate(9);
    std::cerr << "Pair returned: " << p.first << "," << p.second << " Correct pair: " << "1,2" << std::endl;
    p = test_back_translate(10);
    std::cerr << "Pair returned: " << p.first << "," << p.second << " Correct pair: " << "1,3" << std::endl;
    p = test_back_translate(5);
    std::cerr << "Pair returned: " << p.first << "," << p.second << " Correct pair: " << "0,6" << std::endl;
  }
};


#endif //CARVAIBD_SRC_INDEXER_HPP
