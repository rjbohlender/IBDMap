//
// Created by Bohlender,Ryan James on 12/13/19.
//

#ifndef CARVAIBD_SRC_INDEXER_HPP
#define CARVAIBD_SRC_INDEXER_HPP

#define ARMA_DONT_USE_WRAPPER

#include <vector>
#include <armadillo>
#include <utility>
#include <unordered_map>
#include "indexsort.hpp"

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
  std::vector<int> phenotypes;
  std::vector<arma::sword> transitions; // Transition points between pairing sets
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
		  std::vector<int> phenotypes_)
	  : case_count(case_count_), cont_count(cont_count_), sz(samples_.size()), phenotypes(std::move(phenotypes_)),
		samples(std::move(samples_)), case_case(0), case_cont(0), cont_cont(0) {
	setup();
#ifndef NDEBUG
	test_run();
#endif
  }

  /**
   * @brief Sort samples and phenotypes and map the individuals to the correct positions.
   *
   * For n samples there are n^2 positions that we need to map. We exclude the
   * autozygous pairs. There is not an absolute order on pairs, for example does
   * (1, 2) occur before or after (2, 1)? We impose our own order on the samples
   * so that they occur in a specific way. We sort the samples lexicographically
   * and order the samples so that case-case pairs occur first, then case-control
   * pairs, then control-control pairs, in descending order.
   *
   * Envisioning a 2d plane with samples on both the x and y axes, we have
   * a point corresponding to every pair of samples. In the output space, case-case
   * pairs occur before case-control pairs, followed by control-control pairs.
   */
  void setup() {
	IndexSort indexSort(samples); // Ensure that both samples and phenotypes are sorted in the same order.
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
	for (; k > 0; k--) { // How many rows correspond to each sample.
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
	  case -1: break;
	  default:throw (std::runtime_error("Invalid phenotype value."));
	  }
	  current++;
	}
  }

  /**
   * @brief Convert a pair of sample ids to a position in an ordered array
   * @param a The first sample
   * @param b The second sample
   * @return index of the pair in an array that corresponds to all the distinct pairs.
   */
  arma::sword translate(const std::string &a, const std::string &b) {
	auto finda = std::lower_bound(samples.begin(), samples.end(), a);
	auto findb = std::lower_bound(samples.begin(), samples.end(), b);
	if (*finda != a || *findb != b || a == b) {
	  return -1;
	}
	long i = std::distance(samples.begin(), finda);
	long j = std::distance(samples.begin(), findb);
	if (i == j) {
	  throw (std::runtime_error("Searching for two of the same pair."));
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

  /**
   * @brief The inverse of the translate operation, map from pair space to pairs of individuals
   * @param row The row of the pairs.
   * @return The pair of samples corresponding to the given row.
   *
   * Our translate function is bijective and thus has an inverse. In other words,
   * all samples have a distinct output for a distinct input, so we can "undo"
   * the translate operation by enforcing a specific relationship between the
   * order the cases and controls appear in, and the order that they are mapped
   * into IBD pairs.
   */
  std::pair<arma::sword, arma::sword> back_translate(arma::uword row) {
	arma::sword i = -1;
	arma::sword j = -1;

	auto bound = std::lower_bound(transitions.begin(), transitions.end(), row);

	i = std::distance(transitions.begin(), bound);
	if (row == *bound) {
	  i += 1;
	  j = i + 1;
	} else {
	  if (i > 0) {
		j = row - *(bound - 1) + i + 1;
	  } else {
		j = row + 1;
	  }
	}

	return std::make_pair(i, j);
  }

  /**
   * @brief Logic test for the implementation of indexing the pairs.
   */
  static void test_run() {
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

	/**
	 * @brief Test function for mapping from the ordered rows to the product coordinates
	 */
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
