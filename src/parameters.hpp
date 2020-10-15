//
// Created by Bohlender,Ryan James on 9/24/19.
//

#ifndef CARVAIBD_PARAMETERS_HPP
#define CARVAIBD_PARAMETERS_HPP

#include <optional>
#include <armadillo>

/**
 * @brief Runtime parameters
 */
struct Parameters {
  std::string input;
  std::string pheno;
  std::string gmap;
  std::optional<std::string> cov;
  std::optional<std::string> info;
  size_t nperms;
  size_t nthreads;
  std::string output_path;
  unsigned seed;
  std::optional<arma::uword> lower_bound;
  bool swap;
  bool contcont;
  double min_dist;
  std::optional<double> rsquared;
  bool verbose;
  bool enable_testing;
  std::optional<std::vector<int>> range;
  std::optional<double> cM;
  std::optional<double> AF;
};
#endif //CARVAIBD_PARAMETERS_HPP
