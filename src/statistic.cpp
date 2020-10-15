//
// Created by Bohlender,Ryan James on 9/23/19.
//

#include <utility>
#include <boost/algorithm/string/predicate.hpp>
#include "statistic.hpp"
#include "split.hpp"
#include <Range/Enumerate.hpp>

Statistic::Statistic(arma::sp_colvec data_,
					 Breakpoint bp_,
					 std::shared_ptr<std::vector<Indexer>> indexer_,
					 std::vector<std::vector<int>> phenotypes_,
					 std::shared_ptr<Reporter> reporter_,
					 Parameters params_,
					 std::optional<std::vector<std::vector<arma::uword>>> groups_,
					 std::optional<std::shared_ptr<std::vector<std::vector<int32_t>>>> perms_) :
	data(std::move(data_)), indexer(std::move(indexer_)), phenotypes(std::move(phenotypes_)), params(std::move(params_)),
	bp(std::move(bp_)),
	reporter(std::move(reporter_)), groups(std::move(groups_)), perms(std::move(perms_)), gen(params.seed) {
#ifndef NDEBUG
  if (params.enable_testing) {
	test_group_permutation();
  }
#endif
}

double
Statistic::calculate(std::vector<int> &phenotypes_, double cscs_count, double cscn_count, double cncn_count, size_t k) {
  double statistic;

  double cscs = 0;
  double cscn = 0;
  double cncn = 0;

  if (pairs.empty()) {
	for (auto it = data.begin(); it != data.end(); ++it) {
	  auto p = (*indexer)[k].back_translate(it.row());
	  try {
		pairs.emplace_back(p);
	  } catch (std::length_error &e) {
		std::cerr << "Failed to emplace or push at " << __LINE__ << std::endl;
		throw (e);
	  }
	}
  }

  for (auto &p : pairs) {
	auto& [p1, p2] = p;
	int x = phenotypes_[p1];
	int y = phenotypes_[p2];
	if (x == 1) {
	  x1(y, cscs, cscn);
	} else if(x == 0) {
	  x0(y, cscn, cncn);
	} else if(x == -1) {
	} else {
	  std::cerr << "Phenotype: " << x << " p1: " << p1 << std::endl;
	  throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
	}
  }

  permuted_cscs[k].push_back(cscs);
  permuted_cscn[k].push_back(cscn);
  permuted_cncn[k].push_back(cncn);

  if (params.contcont) {
	statistic = cscs / cscs_count - cscn / cscn_count - cncn / cncn_count;
  } else {
	statistic = cscs / cscs_count - cscn / cscn_count;
  }

  return statistic;
}

void Statistic::x1(int y, double &cscs, double &cscn) {
  if(y == 1) {
	cscs += 1.;
  } else if(y == 0) {
	cscn += 0.;
  } else if(y == -1) {
  } else {
	std::cerr << "Phenotype: " << y << std::endl;
	throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
  }
}

void Statistic::x0(int y, double &cscn, double &cncn) {
  if(y == 1) {
	cscn += 1.;
  } else if(y == 0) {
	cncn += 0.;
  } else if(y == -1) {
  } else {
	std::cerr << "Phenotype: " << y << std::endl;
	throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
  }
}

void Statistic::test_group_permutation() const {
  std::map<std::vector<bool>, std::vector<arma::uword>> tfill_patterns;
  boost::optional<std::vector<std::vector<arma::uword>>> tgroups;

  std::stringstream test_data;

  test_data << "#header\n";
  test_data << "test1\tNA\t1\t1\tNA\n";
  test_data << "test2\tNA\t0\t0\tNA\n";
  test_data << "test3\tNA\t1\t0\tNA\n";
  test_data << "test4\tNA\t0\t1\tNA\n";
  test_data << "test5\t1\tNA\tNA\t1\n";
  test_data << "test6\t0\tNA\tNA\t0\n";
  test_data << "test7\t1\tNA\tNA\t0\n";
  test_data << "test8\t0\tNA\tNA\t1\n";

  std::string line;
  arma::uword lineno = 0;

  std::vector<std::vector<int>> tphenotypes;

  while (std::getline(test_data, line)) {
	if (boost::starts_with(line, "#") || lineno == 0) { // Skip the header
	  lineno++;
	  continue;
	}
	RJBUtil::Splitter<std::string_view> splitter(line, " \t");
	std::vector<bool> pattern;
	for (int i = 1; i < splitter.size(); i++) {
	  if (tphenotypes.size() < i) {
		tphenotypes.emplace_back();
	  }
	  if (splitter[i] == "NA") {
		pattern.push_back(false);
		tphenotypes[i - 1].push_back(-1);
	  } else {
		pattern.push_back(true);
		tphenotypes[i - 1].push_back(std::stoi(splitter[i]));
		if (tphenotypes[i - 1].back() != 0 && tphenotypes[i - 1].back() != 1) {
		  std::cerr << splitter[0] << " " << splitter[1] << std::endl;
		}
	  }
	}
	if (tfill_patterns.count(pattern) == 0) {
	  tfill_patterns.emplace(std::make_pair(pattern, std::vector<arma::uword>({lineno - 1})));
	} else {
	  tfill_patterns[pattern].push_back(lineno - 1);
	}
	lineno++;
  }
  for (const auto &tp : tphenotypes) {
	int case_count = 0;
	int control_count = 0;
	for (const auto &v : tp) {
	  switch (v) {
	  case 1:case_count++;
		break;
	  case 0:control_count++;
		break;
	  default:break;
	  }
	}
  }
  if (tfill_patterns.size() > 1) {
	tgroups = std::vector<std::vector<arma::uword>>();
	for (const auto &v : tfill_patterns) {
	  tgroups->push_back(v.second);
	}
	std::cerr << "Groups: " << tgroups->size() << std::endl;
	std::cerr << "Group sizes: ";
	for (const auto &v : tfill_patterns) {
	  std::cerr << v.second.size() << " ";
	  for (const auto &k : v.second) {
		std::cerr << "idx: " << k << " ";
	  }
	}
	std::cerr << std::endl;
  }

  std::cerr << "Test phenotypes\n";
  unsigned long k;
  for (k = 0; k < tphenotypes.size(); k++) {
	for (const auto &v : tphenotypes[k]) {
	  std::cerr << v << " ";
	}
	std::cerr << std::endl;
  }

  std::mt19937 tgen(params.seed);
  for (const auto &v : *tgroups) { // For each of the set of group indices
	std::vector<std::vector<int>> tphenotypes_tmp;
	group_pack(tphenotypes, tphenotypes_tmp, v);
	joint_shuffle(tphenotypes_tmp, tgen);
	group_unpack(tphenotypes, tphenotypes_tmp, v);
  }
  std::cerr << "Test phenotypes after shuffle\n";
  for (k = 0; k < tphenotypes.size(); k++) {
	for (const auto &v : tphenotypes[k]) {
	  std::cerr << v << " ";
	}
	std::cerr << std::endl;
  }
}

void Statistic::run() {
  initialize();

  arma::vec odds;

  unsigned long k = permutations.size();
  while (permutations[k - 1] < params.nperms) {
	if (perms && (*indexer).size() == 1) {
	  phenotypes = *(*perms);
	} else {
	  if (groups) {
		for (const auto &groupIndices : *groups) {
		  std::vector<std::vector<int>> phenotypes_tmp;
		  group_pack(phenotypes, phenotypes_tmp, groupIndices);
		  joint_shuffle(phenotypes_tmp, gen);
		  group_unpack(phenotypes, phenotypes_tmp, groupIndices);
		}
	  } else {
		joint_shuffle(phenotypes, gen);
	  }
	}
	joint_permute();
  }

  if (params.verbose) {
	std::cerr << bp.breakpoint.first << " " << bp.breakpoint.second << " cscs: " << permuted_cscs[0][0] << " cscn: "
			  << permuted_cscn[0][0] << std::endl;
  }

  // Output
  std::stringstream ss;
  build_output(ss);

  reporter->submit(ss.str());
  cleanup();
  done = true;
}

void Statistic::build_output(std::stringstream &ss) {
  for (auto [i, o] : Enumerate(original)) {
	ss << bp.breakpoint.first << "\t" << bp.breakpoint.second << "\t" << o;
	for (int j = 0; j < params.nperms; j++) {
	  ss << "\t" << permuted[i][j];
	}
	ss << std::endl;
  }
}

void Statistic::group_pack(const std::vector<std::vector<int>> &p_original,
						   std::vector<std::vector<int>> &p_tmp,
						   const std::vector<arma::uword> &groupIndices) {
  for(auto [i, p] : Enumerate(p_original)) {
	try {
	  p_tmp.emplace_back(std::vector<int>(groupIndices.size(), 0));
	} catch (std::length_error &e) {
	  std::cerr << "failed to emplace at line " << __LINE__ << std::endl;
	}
	int x = 0;
	for (const auto &j : groupIndices) {
	  p_tmp[i][x] = p[j];
	  x++;
	}
  }
}

void Statistic::group_unpack(std::vector<std::vector<int>> &p_original,
							 const std::vector<std::vector<int>> &p_tmp,
							 const std::vector<arma::uword> &group_indices) {
  for(auto [i, p] : Enumerate(p_tmp)) {
	int x = 0;
	for (const auto &j : group_indices) {
	  p_original[i][j] = p[x];
	  x++;
	}
  }
}

void Statistic::joint_permute() {
  double val;
  for (auto [i, p] : Enumerate(phenotypes)) {
    if (perms) {
	  val = calculate(p, (*indexer)[0].case_case, (*indexer)[0].case_cont, (*indexer)[0].cont_cont, 0);
    } else {
	  val = calculate(p, (*indexer)[i].case_case, (*indexer)[i].case_cont, (*indexer)[i].cont_cont, i);
    }
	if (val > original[0]) {
	  successes[0]++;
	}
	permutations[0]++;
	permuted[0].push_back(val);
  }
}

void Statistic::joint_shuffle(std::vector<std::vector<int>> &phen, std::mt19937 &gen) {
  for (size_t i = phen[0].size() - 1; i > 0; i--) {
	std::uniform_int_distribution<> dis(0, i);
	int j = dis(gen);
	for (auto &p : phen) { // Permute all phenotypes the same way.
	  std::swap(p[i], p[j]);
	}
  }
}

void Statistic::initialize() {
  for (auto [i, idx] : Enumerate(*indexer)) {
	permuted.emplace_back();
	permuted_cscs.emplace_back();
	permuted_cscn.emplace_back();
	permuted_cncn.emplace_back();
	original.push_back(
		calculate(idx.phenotypes, idx.case_case, idx.case_cont, idx.cont_cont, i));
	successes.push_back(0);
	permutations.push_back(0);
  }
}

void Statistic::cleanup() {
  data.reset();
  permuted.clear();
  permuted.shrink_to_fit();
  permuted_cscs.clear();
  permuted_cscs.shrink_to_fit();
  permuted_cscn.clear();
  permuted_cscn.shrink_to_fit();
  permuted_cncn.clear();
  permuted_cncn.shrink_to_fit();
}

