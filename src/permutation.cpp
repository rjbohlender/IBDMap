//
// Created by Bohlender,Ryan James on 8/4/18.
//

#include <stocc/stocc.h>
#include <ctime>
#include <thread>
#include <cassert>
#include <fmt/include/fmt/ostream.h>
#include "permutation.hpp"

Permute::Permute()
	: sto(std::random_device{}()) {}

Permute::Permute(int seed)
	: sto(seed) {}

void Permute::get_permutations(std::shared_ptr<std::vector<std::vector<int32_t>>> permutations,
							   arma::vec &odds,
							   int ncases,
							   int nperm,
							   int nthreads) {
  std::vector<double> odds_ = arma::conv_to<std::vector<double>>::from(odds);
  std::vector<int32_t> m(odds.n_rows, 1);

  // Initialize permutations
  permutations->resize(nperm);
  for (int i = 0; i < nperm; i++) {
	(*permutations)[i].resize(odds.n_rows);
  }

  int step = nperm / nthreads;
  int remaining = nperm;
  std::vector<std::thread> threads;
  for (int i = 0; i < nthreads; i++) {
	int seed = sto.IRandom(0, std::numeric_limits<int>::max());
	int offset = i * step;
	if (remaining < 0) {
	  fmt::print(std::cerr, "Failed during distribution of covariate adjusted permutation jobs.\n");
	  std::exit(-1);
	}

	PermuteData data{
		permutations,
		&m[0],
		&odds_[0],
		static_cast<int>(ncases),
		static_cast<int>(odds.n_rows),
		offset,
		i == nthreads - 1 ? remaining : step,
		seed
	};
	threads.push_back(std::thread(&Permute::permute_thread, data));
	remaining -= step;
  }

  for (auto &t : threads) {
	t.join();
  }
}

void Permute::permute_thread(PermuteData data) {
  StochasticLib3 rng(data.seed);

  for (int i = 0; i < data.nperm; i++) {
	rng.MultiFishersNCHyp(&(*data.permutations)[data.offset + i][0], data.m, data.odds, data.ncases, data.ngroups);
  }

}

