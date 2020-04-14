//
// Created by Bohlender,Ryan James on 9/23/19.
//

#include <utility>
#include <boost/algorithm/string/predicate.hpp>
#include "statistic.hpp"
#include "split.hpp"

Statistic::Statistic(arma::sp_colvec data_,
                     Breakpoint bp_,
                     std::vector<Indexer> &indexer_,
                     std::vector<std::string> &samples_,
                     std::vector<std::vector<int>> &phenotypes_,
                     std::shared_ptr<Reporter> reporter_,
                     Parameters &params_,
                     boost::optional<std::vector<std::vector<arma::uword>>> groups_,
                     boost::optional<std::shared_ptr<std::vector<std::vector<int32_t>>>> perms_) :
    data(std::move(data_)), indexer(indexer_), samples(samples_), phenotypes(phenotypes_), params(params_), bp(std::move(bp_)),
    reporter(std::move(reporter_)), groups(std::move(groups_)), perms(perms_) {
#ifndef NDEBUG
  test_group_permutation();
#endif
}

double
Statistic::calculate(std::vector<int> &phenotypes_, double cscs_count, double cscn_count, double cncn_count, int k) {
  double statistic;

  double cscs = 0;
  double cscn = 0;
  double cncn = 0;

  if (pairs.empty()) {
    // Save time allocating the full amount
    double ibd_pairs = arma::accu(data);

    rows.reserve(ibd_pairs);
    pairs.reserve(ibd_pairs);
    for (auto it = data.begin(); it != data.end(); it++) {
      auto[p1, p2] = indexer[0].back_translate(it.row());
      pairs.emplace_back(std::make_pair(p1, p2));
      rows.push_back(it.row()); // Store rows for later lookup
      switch (phenotypes_[p1]) {
      case 1:
        switch (phenotypes_[p2]) {
        case 1:cscs += 1.;
          break;
        case 0:cscn += 1.;
          break;
        case -1:break;
        default:
          std::cerr << "Phenotype: " << phenotypes_[p2] << " p2: " << p2 << std::endl;
          throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
        }
        break;
      case 0:
        switch (phenotypes_[p2]) {
        case 1:cscn += 1.;
          break;
        case 0:cncn += 1.;
          break;
        case -1:break;
        default:
          std::cerr << "Phenotype: " << phenotypes_[p2] << " p2: " << p2 << std::endl;
          throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
        }
        break;
      case -1:break;
      default:
        std::cerr << "Phenotype: " << phenotypes_[p1] << " p1: " << p1 << std::endl;
        throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
      }
    }
  } else {
    for (int i = 0; i < pairs.size(); i++) {
      const auto &p = pairs[i];
      const auto &r = rows[i];
      switch (phenotypes_[p.first]) {
      case 1:
        switch (phenotypes_[p.second]) {
        case 1:cscs += 1.;
          break;
        case 0:cscn += 1.;
          break;
        case -1:break;
        default:
          std::cerr << "Phenotype: " << phenotypes_[p.second] << " p.second: " << p.second << std::endl;
          throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
        }
        break;
      case 0:
        switch (phenotypes_[p.second]) {
        case 1:cscn += 1.;
          break;
        case 0:cncn += 1.;
          break;
        case -1:break;
        default:
          std::cerr << "Phenotype: " << phenotypes_[p.second] << " p.second: " << p.second << std::endl;
          throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
        }
        break;
      case -1:break;
      default:
        std::cerr << "Phenotype: " << phenotypes_[p.first] << " p.first: " << p.first << std::endl;
        throw (std::runtime_error("ERROR: invalid phenotype in calculate."));
      }
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

void Statistic::test_group_permutation() {
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
        if (tphenotypes[i-1].back() != 0 && tphenotypes[i-1].back() != 1) {
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
  for (unsigned long k = 0; k < tphenotypes.size(); k++) {
    int case_count = 0;
    int control_count = 0;
    for (const auto &v : tphenotypes[k]) {
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
  unsigned long k = 0;
  for (k = 0; k < tphenotypes.size(); k++) {
    for (const auto &v : tphenotypes[k]) {
      std::cerr << v << " ";
    }
    std::cerr << std::endl;
  }

  std::mt19937 gen(params.seed);
  for (const auto &v : *tgroups) { // For each of the set of group indices
    std::vector<std::vector<int>> tphenotypes_tmp;
    for (k = 0; k < tphenotypes.size(); k++) {
      tphenotypes_tmp.emplace_back(std::vector<int>(v.size(), 0));
      int x = 0;
      for (const auto &i : v) {
        tphenotypes_tmp[k][x] = tphenotypes[k][i];
        x++;
      }
    }
    for (int i = tphenotypes_tmp[0].size() - 1; i > 0; i--) { // Shuffle the indices
      std::uniform_int_distribution<> dis(0, i);
      int j = dis(gen);
      for (k = 0; k < tphenotypes.size(); k++) {
        std::swap(tphenotypes_tmp[k][i], tphenotypes_tmp[k][j]);
      }
    }
    for (k = 0; k < tphenotypes.size(); k++) {
      int x = 0;
      for (const auto &i : v) {
        tphenotypes[k][i] = tphenotypes_tmp[k][x];
        x++;
      }
    }
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
  unsigned long k = 0;
  arma::vec odds;
  for (; k < indexer.size(); k++) {
    permuted.emplace_back();
    permuted_cscs.emplace_back();
    permuted_cscn.emplace_back();
    permuted_cncn.emplace_back();
    original.push_back(
        calculate(indexer[k].phenotypes, indexer[k].case_case, indexer[k].case_cont, indexer[k].cont_cont, 0));
    successes.push_back(0);
    permutations.push_back(0);
  }

  std::mt19937 gen(params.seed);

  while (permutations[k - 1] < params.nperms) {
    if (perms && indexer.size() == 1) {
      phenotypes[0] = (*(*perms))[permutations[0]];
      double val = calculate(phenotypes[0], indexer[0].case_case, indexer[0].case_cont, indexer[0].cont_cont, 0);
      if (val > original[0])
        successes[0]++;
      permutations[0]++;
      permuted[0].push_back(val);
    } else {
      if (groups) { // If we need to do grouped permutation
        for (const auto &v : *groups) { // For each of the set of group indices
          std::vector<std::vector<int>> phenotypes_tmp;
          for (k = 0; k < phenotypes.size(); k++) {
            phenotypes_tmp.emplace_back(std::vector<int>(v.size(), 0));
            int x = 0;
            for (const auto &i : v) {
              phenotypes_tmp[k][x] = phenotypes[k][i];
              x++;
            }
          }
          for (int i = phenotypes_tmp[0].size() - 1; i > 0; i--) { // Shuffle the indices
            std::uniform_int_distribution<> dis(0, i);
            int j = dis(gen);
            for (k = 0; k < phenotypes.size(); k++) {
              std::swap(phenotypes_tmp[k][i], phenotypes_tmp[k][j]);
            }
          }
          for (k = 0; k < phenotypes.size(); k++) {
            int x = 0;
            for (const auto &i : v) {
              phenotypes[k][i] = phenotypes_tmp[k][x];
              x++;
            }
          }
        }
      } else {
        for (unsigned long i = phenotypes[0].size() - 1; i > 0; i--) {
          std::uniform_int_distribution<> dis(0, i);
          int j = dis(gen);
          for (k = 0; k < phenotypes.size(); k++) { // Permute all phenotypes the same way.
            std::swap(phenotypes[k][i], phenotypes[k][j]);
          }
        }
      }
    }
    k = 0;
    for (auto &v : phenotypes) {
      double val = calculate(v, indexer[k].case_case, indexer[k].case_cont, indexer[k].cont_cont, k);
      if (val > original[k])
        successes[k]++;
      permutations[k]++;
      permuted[k].push_back(val);
      k++;
    }
  }

  std::cerr << bp.breakpoint.first << " " << bp.breakpoint.second << " cscs: " << permuted_cscs[0][0] << " cscn: " << permuted_cscn[0][0] << std::endl;

  // Output
  std::stringstream iss;
  if(bp.breakpoint.second == "8560476") {
    std::cerr << bp.breakpoint.first << "\t" << bp.breakpoint.second << "\t" << original[0];
    for(const auto &v : permuted[0]) {
      std::cerr << "\t" << v;
    }
    std::cerr << std::endl;
  }
  for (k = 0; k < original.size(); k++) {
    iss << bp.breakpoint.first << "\t" << bp.breakpoint.second << "\t" << original[k];
    for (int i = 0; i < params.nperms; i++) {
      iss << "\t" << permuted[k][i];
    }
    iss << std::endl;
  }

  reporter->submit(iss.str());
  data.reset();
  done = true;
}

