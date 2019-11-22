//
// Created by Bohlender,Ryan James on 9/23/19.
//

#ifndef CARVAIBD_STATISTIC_HPP
#define CARVAIBD_STATISTIC_HPP

#include <armadillo>
#include <future>
#include "parser.hpp"
#include "reporter.hpp"
#include "parameters.hpp"

/**
 * @brief IBD Statistic Calculation Class
 */
class Statistic {
  enum class State {
    avg,
    no_avg
  };
  enum class Range {
    case_control,
    control_control
  };
  arma::sp_colvec data;
  Indexer &indexer;
  Parser<std::string> &parser;
  Parameters &params;

  State current_state;
  Range current_range;

  std::shared_ptr<Reporter> reporter;

  std::vector<std::pair<arma::sword, arma::sword>> pairs;
public:
  bool done = false;

  Breakpoint &bp;

  double original;
  double successes;
  double permutations;

  double case_case_ibd;
  double case_cont_ibd;
  double cont_cont_ibd;

  std::vector<double> permuted;
  std::vector<double> permuted_cscs;
  std::vector<double> permuted_cscn;

  std::vector<arma::sword> rows;
  std::vector<double> avgs;
  std::vector<double> stats;

  Statistic(arma::sp_colvec &&data_,
			Breakpoint &bp_,
			Indexer &indexer_,
			Parser<std::string> &parser_,
			std::shared_ptr<Reporter> reporter_,
			Parameters &params_);

  double calculate();
  double calculate(std::vector<int> &phenotypes_);
  void calc_avg(const std::vector<int> &phenotypes_);

  std::pair<double, double> calc_original_rates();
  std::pair<double, double> calc_rates(std::vector<int> &phenotypes_);

  void run();
};

#endif //CARVAIBD_STATISTIC_HPP
