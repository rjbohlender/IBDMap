//
// Created by Bohlender,Ryan James on 9/23/19.
//

#ifndef CARVAIBD_STATISTIC_HPP
#define CARVAIBD_STATISTIC_HPP

#include <armadillo>

/**
 * @brief IBD Statistic Calculation Class
 */
class Statistic {
  enum class State {
    no_avg_no_cont_cont,
    no_avg_cont_cont,
    avg_no_cont_cont,
    avg_cont_cont
  };
  arma::vec &data;
  double case_case_avg;
  double case_control_avg;
  double control_control_avg;

  // Indices
  arma::uvec case_case;
  arma::uvec case_control;
  arma::uvec control_control;

  State current_state;

  double calculate();

public:
  Statistic(arma::vec &data_,
			arma::uvec case_case_,
			arma::uvec case_control_,
			double case_case_avg_,
			double case_control_avg_);
  Statistic(arma::vec &data_,
			arma::uvec case_case_,
			arma::uvec case_control_,
			arma::uvec control_control_,
			double case_case_avg_,
			double case_control_avg_,
			double control_control_avg_);
  Statistic(arma::vec &data_,
			arma::uvec case_case_,
			arma::uvec case_control_);
  Statistic(arma::vec &data_,
			arma::uvec case_case_,
			arma::uvec case_control_,
			arma::uvec control_control_);
};

#endif //CARVAIBD_STATISTIC_HPP
