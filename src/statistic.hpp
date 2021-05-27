//
// Created by Bohlender,Ryan James on 9/23/19.
//

#ifndef CARVAIBD_STATISTIC_HPP
#define CARVAIBD_STATISTIC_HPP

#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <future>
#include "reporter.hpp"
#include "parameters.hpp"
#include "indexer.hpp"
#include "breakpoint.hpp"
#include <boost/optional.hpp>

/**
 * @brief IBD Statistic Calculation Class
 */
class Statistic {
  arma::sp_mat data;
  std::shared_ptr<std::vector<Indexer>> indexer;
  std::vector<std::vector<int>> phenotypes;
  Parameters params;
  arma::uword keep_index;  // Index of the retained statistic.

  std::shared_ptr<Reporter> reporter;

  std::vector<std::vector<std::pair<arma::sword, arma::sword>>> pairs;

  std::mt19937 gen;

  Breakpoint bp;

  std::vector<double> original;
  std::vector<double> successes;
  std::vector<size_t> permutations;

  std::optional<std::vector<std::vector<arma::uword>>> groups;
  std::optional<std::shared_ptr<std::vector<std::vector<int32_t>>>> perms;

  arma::vec cscs;
  arma::vec cscn;
  arma::vec cncn;

  std::vector<std::vector<double>> permuted;
  std::vector<std::vector<double>> permuted_cscs;
  std::vector<std::vector<double>> permuted_cscn;
  std::vector<std::vector<double>> permuted_cncn;

  static void x1(int y, double &cscs, double &cscn);
  static void x0(int y, double &cscn, double &cncn);
  void initialize();
  static void joint_shuffle(std::vector<std::vector<int>> &phen, std::mt19937 &gen);
  void joint_permute();
  static void group_unpack(std::vector<std::vector<int>> &p_original,
						   const std::vector<std::vector<int>> &p_tmp,
						   const std::vector<arma::uword> &group_indices);
  static void group_pack(const std::vector<std::vector<int>> &p_original,
						 std::vector<std::vector<int>> &p_tmp,
						 const std::vector<arma::uword> &groupIndices);
  void build_output(std::stringstream &ss);

public:
  bool done = false;

  Statistic(arma::sp_mat data_,
			Breakpoint bp_,
			std::shared_ptr<std::vector<Indexer>> indexer_,
			std::vector<std::vector<int>> phenotypes_,
			std::shared_ptr<Reporter> reporter_,
			Parameters params_,
			std::optional<std::vector<std::vector<arma::uword>>> groups_,
			std::optional<std::shared_ptr<std::vector<std::vector<int32_t>>>> perms_);

  arma::vec calculate(std::vector<int> &phenotypes_, double cscs_count, double cscn_count, double cncn_count, size_t k);

  void run();
  void cleanup();
  static void test_statistic();
  static void test_group_permutation(unsigned int seed);
};

#endif //CARVAIBD_STATISTIC_HPP
