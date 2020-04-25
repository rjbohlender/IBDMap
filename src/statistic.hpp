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
  arma::sp_colvec data;
  std::vector<Indexer> &indexer;
  std::vector<std::string> &samples;
  std::vector<std::vector<int>> phenotypes;
  Parameters &params;

  std::shared_ptr<Reporter> reporter;

  std::vector<std::pair<arma::sword, arma::sword>> pairs;

  void test_group_permutation();

public:
  bool done = false;

  Breakpoint bp;

  boost::optional<std::vector<std::vector<arma::uword>>> groups;
  boost::optional<std::shared_ptr<std::vector<std::vector<int32_t>>>> perms;

  std::vector<double> original;
  std::vector<double> successes;
  std::vector<double> permutations;

  std::vector<std::vector<double>> permuted;
  std::vector<std::vector<double>> permuted_cscs;
  std::vector<std::vector<double>> permuted_cscn;
  std::vector<std::vector<double>> permuted_cncn;

  std::vector<arma::sword> rows;

  Statistic(arma::sp_colvec data_,
            Breakpoint bp_,
            std::vector<Indexer> &indexer_,
            std::vector<std::string> &samples_,
            std::vector<std::vector<int>> &phenotypes_,
            std::shared_ptr<Reporter> reporter_,
            Parameters &params_,
            boost::optional<std::vector<std::vector<arma::uword>>> groups_,
            boost::optional<std::shared_ptr<std::vector<std::vector<int32_t>>>> perms_);

  double calculate(std::vector<int> &phenotypes_, double cscs_count, double cscn_count, double cncn_count, int k);

  void run();
  void cleanup();
};

#endif //CARVAIBD_STATISTIC_HPP
