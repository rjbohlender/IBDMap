//
// Created by Bohlender,Ryan James on 9/23/19.
//

#include <utility>
#include "statistic.hpp"

Statistic::Statistic(arma::sp_colvec &&data_,
                     Breakpoint bp_,
                     std::vector<Indexer> &indexer_,
                     std::vector<std::string> &samples_,
                     std::vector<std::vector<int>> &phenotypes_,
                     std::shared_ptr<Reporter> reporter_,
                     Parameters &params_,
                     boost::optional<std::vector<std::vector<arma::uword>>> groups_) :
    data(std::move(data_)), indexer(indexer_), samples(samples_), phenotypes(phenotypes_), params(params_), bp(std::move(bp_)),
    reporter(std::move(reporter_)), groups(std::move(groups_)) {
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

  statistic = cscs / cscs_count - cscn / cscn_count;

  return statistic;
}

void Statistic::run() {
  unsigned long k = 0;
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
    if (groups) { // If we need to do grouped permutation
      for (const auto &v : *groups) { // For each of the set of group indices
        std::vector<arma::uword> p = v; // Group indices
        for(int i = p.size() - 1; i > 0; i--) { // Shuffle the indices
          std::uniform_int_distribution<> dis(0, i);
          int j = dis(gen);
          std::swap(p[i], p[j]);
        }
        for(k = 0; k < phenotypes.size(); k++) { // For each of the phenotypes
          for (int i = 0; i < v.size(); i++) { // For each group index
            int m = v[i];
            int n = p[i];
            std::swap(phenotypes[k][n], phenotypes[k][m]); // Swap the values of the phenotypes at the indices
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

  // Output
  std::stringstream iss;
  for (k = 0; k < permuted_cscs.size(); k++) {
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

