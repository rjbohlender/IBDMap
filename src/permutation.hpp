//
// Created by Bohlender,Ryan James on 8/4/18.
//

#ifndef PERMUTE_ASSOCIATE_PERMUTATION_HPP
#define PERMUTE_ASSOCIATE_PERMUTATION_HPP

#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <stocc/randomc.h>
#include <stocc/stocc.h>
#include <memory>

struct Permute {
  Permute();
  explicit Permute(int seed);

  struct PermuteData {
	std::shared_ptr<std::vector<std::vector<int32_t>>> permutations;
	int32_t *m;
	double *odds;
	int ncases;
	int ngroups;
	int offset;
	int nperm;
	int seed;
  };

  void get_permutations(std::shared_ptr<std::vector<std::vector<int32_t>>> permutations,
						arma::vec &odds,
						int ncases,
						int nperm,
						int nthreads);

  static void permute_thread(PermuteData data);

  StochasticLib3 sto;
};
#endif //PERMUTE_ASSOCIATE_PERMUTATION_HPP
