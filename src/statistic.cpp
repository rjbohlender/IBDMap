//
// Created by Bohlender,Ryan James on 9/23/19.
//

#include "statistic.hpp"

Statistic::Statistic(arma::sp_colvec &&data_,
					 Breakpoint &bp_,
					 Indexer &indexer_,
					 Parser<std::string> &parser_,
					 std::shared_ptr<Reporter> reporter_,
					 Parameters &params_) :
	data(std::move(data_)), indexer(indexer_), parser(parser_), params(params_), bp(bp_),
	reporter(std::move(reporter_)) {
  original = 0;
  successes = 0;
  permutations = 0;
  current_state = State::avg;
  current_range = Range::case_control;
}

double Statistic::calculate() {
  double statistic;
  double case_case_rate = arma::accu(data.rows(0, indexer.case_case - 1)) / static_cast<double>(indexer.case_case);
  double other_rate;

  case_case_ibd = arma::accu(data.rows(0, indexer.case_case - 1));
  case_cont_ibd = arma::accu(data.rows(indexer.case_case, indexer.case_case + indexer.case_cont - 1));
  cont_cont_ibd = arma::accu(data.rows(indexer.case_case + indexer.case_cont, data.n_rows - 1));

  arma::uword first_idx;
  arma::uword secnd_idx;
  double n;

  switch (current_range) {
  case Range::case_control: first_idx = indexer.case_case;
	secnd_idx = indexer.case_case + indexer.case_cont - 1;
	n = indexer.case_cont;
	break;
  case Range::control_control: first_idx = indexer.case_case + indexer.case_cont;
	secnd_idx = data.n_rows - 1;
	n = indexer.cont_cont;
	break;
  default:throw (std::runtime_error("Not one of the allowed ranged in Statistic."));
  }
  other_rate = arma::accu(data.rows(first_idx, secnd_idx)) / n;
  switch (current_state) {
  case (State::no_avg): statistic = case_case_rate - other_rate;
	break;
  case (State::avg): statistic = (case_case_rate - other_rate);
	break;
  default: throw (std::runtime_error("Not one of the allowed states in Statistic."));
  }

  return statistic;
}

double Statistic::calculate(std::vector<int> &phenotypes_) {
  std::vector<std::string> &samples = parser.samples;
  double statistic;

  double ibd_pairs = arma::accu(data);

  double case_case_rate = 0;
  double case_cont_rate = 0;
  double cont_cont_rate = 0;

  // Save time allocating the full amount
  avgs.push_back(0);
  if (pairs.empty()) {
	rows.reserve(ibd_pairs);
	pairs.reserve(ibd_pairs);
	for (auto it = data.begin(); it != data.end(); it++) {
#ifdef LOWMEM
	  auto[p1, p2] = indexer.back_translate_alt(it.row());
#else
	  auto [p1, p2] = indexer.back_translate(it.row());
#endif
	  pairs.emplace_back(std::make_pair(p1, p2));
	  rows.push_back(it.row()); // Store rows for later lookup
	  switch (phenotypes_[p1]) {
	  case 1:
		switch (phenotypes_[p2]) {
		case 1: case_case_rate += 1. / indexer.case_case;
		  avgs.back() += parser.row_pair_avg[it.row()] / params.total_breakpoints;
		  break;
		case 0: case_cont_rate += 1. / indexer.case_cont;
		  avgs.back() -= parser.row_pair_avg[it.row()] / params.total_breakpoints;
		  break;
		default: throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
		}
		break;
	  case 0:
		switch (phenotypes_[p2]) {
		case 1: case_cont_rate += 1. / indexer.case_cont;
		  avgs.back() -= parser.row_pair_avg[it.row()] / params.total_breakpoints;
		  break;
		case 0: cont_cont_rate += 1. / indexer.cont_cont;
		  break;
		default: throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
		}
		break;
	  default: throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
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
		case 1: case_case_rate += 1. / indexer.case_case;
		  avgs.back() += parser.row_pair_avg[r] / params.total_breakpoints;
		  break;
		case 0: case_cont_rate += 1. / indexer.case_cont;
		  avgs.back() -= parser.row_pair_avg[r] / params.total_breakpoints;
		  break;
		default: throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
		}
		break;
	  case 0:
		switch (phenotypes_[p.second]) {
		case 1: case_cont_rate += 1. / indexer.case_cont;
		  avgs.back() -= parser.row_pair_avg[r] / params.total_breakpoints;
		  break;
		case 0: cont_cont_rate += 1. / indexer.cont_cont;
		  break;
		default: throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
		}
		break;
	  default: throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
	  }
	}
	permuted_cscs.push_back(case_case_rate);
	permuted_cscn.push_back(case_cont_rate);
  }
  statistic = case_case_rate - case_cont_rate;
  stats.push_back(statistic);
  return statistic;
}

std::pair<double, double> Statistic::calc_rates(std::vector<int> &phenotypes_) {
  std::vector<std::string> &samples = parser.samples;

  double ibd_pairs = arma::accu(data);

  double case_case_rate = 0;
  double case_cont_rate = 0;
  double cont_cont_rate = 0;
  unsigned long first = 0;
  arma::sword row;
  if (pairs.empty()) {
	pairs.reserve(ibd_pairs);
	for (auto it = data.begin(); it != data.end(); it++) {
#ifdef LOWMEM
	  auto[p1, p2] = indexer.back_translate_alt(it.row());
#else
	  auto [p1, p2] = indexer.back_translate(it.row());
#endif
	  pairs.emplace_back(std::make_pair(p1, p2));
	  switch (phenotypes_[p1]) {
	  case 1:
		switch (phenotypes_[p2]) {
		case 1: case_case_rate++;
		  break;
		case 0: case_cont_rate++;
		  break;
		default: throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
		}
		break;
	  case 0:
		switch (phenotypes_[p2]) {
		case 1: case_cont_rate++;
		  break;
		case 0: cont_cont_rate++;
		  break;
		default: throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
		}
		break;
	  default: throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
	  }
	}
	permuted_cscs.push_back(case_case_rate);
	permuted_cscn.push_back(case_cont_rate);
  } else {
	for (const auto &p : pairs) {
	  switch (phenotypes_[p.first]) {
	  case 1:
		switch (phenotypes_[p.second]) {
		case 1: case_case_rate++;
		  break;
		case 0: case_cont_rate++;
		  break;
		default: throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
		}
		break;
	  case 0:
		switch (phenotypes_[p.second]) {
		case 1: case_cont_rate++;
		  break;
		case 0: cont_cont_rate++;
		  break;
		default: throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
		}
		break;
	  default: throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
	  }
	}
	permuted_cscs.push_back(case_case_rate);
	permuted_cscn.push_back(case_cont_rate);
  }
  return std::make_pair(case_case_rate, case_cont_rate);
}

void Statistic::calc_avg(const std::vector<int> &phenotypes_) {
  avgs.push_back(0);
  for (auto it = data.begin(); it != data.end(); it++) {
#ifdef LOWMEM
	auto[p1, p2] = indexer.back_translate_alt(it.row());
#else
	auto [p1, p2] = indexer.back_translate(it.row());
#endif
	// pairs.emplace_back(std::make_pair(p1, p2));
	switch (phenotypes_[p1]) {
	case 1:
	  switch (phenotypes_[p2]) {
	  case 1:
		avgs.back() += parser.row_pair_avg[it.row()] / params.total_breakpoints;
		break;
	  case 0:
		avgs.back() -= parser.row_pair_avg[it.row()] / params.total_breakpoints;
		break;
	  default: throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
	  }
	  break;
	case 0:
	  switch (phenotypes_[p2]) {
	  case 1:
		avgs.back() -= parser.row_pair_avg[it.row()] / params.total_breakpoints;
		break;
	  case 0:
		break;
	  default: throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
	  }
	  break;
	default: throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
	}
  }
}

void Statistic::run() {
  original = calculate(indexer.phenotypes);
  calc_avg(indexer.phenotypes);
  successes = 0;
  permutations = 0;

  std::vector<int> phenotypes(indexer.case_count + indexer.cont_count, 0);
  for (int i = 0; i < indexer.case_count; i++) {
	phenotypes[i] = 1;
  }

  std::mt19937 gen(params.seed);

  arma::wall_clock timer;
  while (permutations < params.nperms) {
	for (unsigned long i = phenotypes.size() - 1; i > 0; i--) {
	  std::uniform_int_distribution<> dis(0, i);
	  int j = dis(gen);
	  int tmp = phenotypes[j];
	  phenotypes[j] = phenotypes[i];
	  phenotypes[i] = tmp;
	}

	// timer.tic();
	double val = calculate(phenotypes);
	// std::cerr << "Time to calculate: " << timer.toc() << std::endl;
	if (val > original)
	  successes++;
	permutations++;
	permuted.push_back(val);
  }
  std::stringstream iss;
#if 0
  iss << bp.breakpoint.first << "\t" << bp.breakpoint.second << "\tP-value:\t"
	  << (successes + 1) / (permutations + 1) << "\tSuccesses:\t" << successes << "\tPerms:\t" << permutations
	  << "\tDelta:\t" << original << "\tcscs:\t" << case_case_ibd << "\tcscn:\t" << case_cont_ibd << "\tcncn:\t" << cont_cont_ibd << std::endl;
#else
  iss << bp.breakpoint.first << "\t" << bp.breakpoint.second;
  for(int i = 0; i < params.nperms; i++) {
	iss << "\t" << stats[i];
  }
  iss << std::endl;
  iss << bp.breakpoint.first << "\t" << bp.breakpoint.second;
  for(int i = 0; i < params.nperms; i++) {
    iss << "\t" << avgs[i];
  }
  iss << std::endl;
#endif

  reporter->submit(iss.str());
  done = true;
}

std::pair<double, double> Statistic::calc_original_rates() {
  double case_case_rate = arma::accu(data.rows(0, indexer.case_case - 1));
  double other_rate;

  case_case_ibd = arma::accu(data.rows(0, indexer.case_case - 1));
  case_cont_ibd = arma::accu(data.rows(indexer.case_case, indexer.case_case + indexer.case_cont - 1));
  cont_cont_ibd = arma::accu(data.rows(indexer.case_case + indexer.case_cont, data.n_rows - 1));

  arma::uword first_idx;
  arma::uword secnd_idx;
  double n;

  switch (current_range) {
  case Range::case_control: first_idx = indexer.case_case;
	secnd_idx = indexer.case_case + indexer.case_cont - 1;
	n = indexer.case_cont;
	break;
  case Range::control_control: first_idx = indexer.case_case + indexer.case_cont;
	secnd_idx = data.n_rows - 1;
	n = indexer.cont_cont;
	break;
  default:throw (std::runtime_error("Not one of the allowed ranged in Statistic."));
  }
  other_rate = arma::accu(data.rows(first_idx, secnd_idx));

  return std::make_pair(case_case_rate, other_rate);
}
