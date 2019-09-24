//
// Created by Bohlender,Ryan James on 9/23/19.
//

#include "statistic.hpp"

double Statistic::calculate() {
  double statistic;
  double case_case_rate = arma::sum(data(case_case)) / static_cast<double>(case_case.n_rows);
  double case_control_rate;
  double control_control_rate;

  switch (current_state) {
  case (State::no_avg_no_cont_cont):
    case_control_rate = arma::sum(data(case_control)) / static_cast<double>(case_control.n_rows);
	statistic = case_case_rate - case_control_rate;
    break;
  case (State::no_avg_cont_cont):
	control_control_rate = arma::sum(data(control_control)) / static_cast<double>(control_control.n_rows);
	statistic = case_case_rate - control_control_rate;
    break;
  case (State::avg_no_cont_cont):
	case_control_rate = arma::sum(data(case_control)) / static_cast<double>(case_control.n_rows);
	statistic = (case_case_rate - case_case_avg) - (case_control_rate - case_control_avg);
	break;
  case (State::avg_cont_cont):
	control_control_rate = arma::sum(data(control_control)) / static_cast<double>(control_control.n_rows);
	statistic = (case_case_rate - case_case_avg) - (control_control_rate - control_control_avg);
	break;
  default: throw (std::runtime_error("Not one of the allowed states in Statistic."));
  }

  return statistic;
}
