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
  Permute(int seed);
  Permute(const Permute &other);
  Permute &operator=(const Permute &rhs);

  void get_permutations(std::shared_ptr<std::vector<std::vector<int32_t>>> permutations,
						arma::colvec &odds,
						arma::uword ncases,
						arma::uword nperm,
						arma::uword nthreads);

  void permute_thread(std::shared_ptr<std::vector<std::vector<int32_t>>> p,
						int32_t *m,
						double *odds,
						int ncases,
						int ngroups,
						int offset,
						int nperm,
						int seed);
  std::vector<std::vector<int32_t>> permutations_maj_bin(int nperm,
														 arma::vec &odds,
														 arma::uword ncases,
														 arma::uvec &mac_indices,
														 arma::uvec &maj_indices,
														 const std::string &transcript);
  std::vector<std::vector<int32_t>> permutations_mac_bin(int nperm,
														 arma::vec &odds,
														 arma::uword ncases,
														 arma::uvec &mac_indices,
														 arma::uvec &maj_indices,
														 arma::uword &approximate,
														 const std::string &transcript);
  std::vector<std::vector<int32_t>> permutations_bin(int nperm,
													 arma::vec &odds,
													 arma::uword ncases,
													 arma::uvec &mac_indices,
													 arma::uvec &maj_indices,
													 arma::uword &approximate,
													 const std::string &transcript);
  arma::vec calculate_fisher_mean(int32_t n, arma::vec &odds);
  std::vector<std::vector<int32_t>> cases_in_bins(int nperm,
												  arma::colvec &odds,
												  int ncases,
												  std::vector<int32_t> &bin_counts);
  std::vector<int32_t> random_case_count(int nperm,
										 arma::uvec &mac_indices,
										 arma::uvec &maj_indices,
										 arma::vec &prob,
										 int ncases);
  std::vector<int32_t> random_case_count(int nperm,
										 arma::uvec &mac_indices,
										 arma::uvec &maj_indices,
										 arma::vec &prob,
										 int ncases,
										 int n_maj_bins);

  auto unpack(int x, int y, bool shuffle) -> arma::uvec;
  auto reset() -> void;

  StochasticLib3 sto;
  // Preserve group info for transcript
  std::map<std::string, bool> bins_built;
  std::map<std::string, std::vector<double>> odds_;
  std::map<std::string, std::vector<int32_t>> m;
  std::map<std::string, double> mac_bins;
  std::map<std::string, double> maj_bins;
  std::map<std::string, std::vector<std::vector<int32_t>>> ret;
  std::map<std::string, arma::uvec> sort_mac_idx;
  std::map<std::string, arma::uvec> sort_maj_idx;
};
#endif //PERMUTE_ASSOCIATE_PERMUTATION_HPP
