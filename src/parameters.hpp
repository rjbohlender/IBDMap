//
// Created by Bohlender,Ryan James on 9/24/19.
//

#ifndef CARVAIBD_PARAMETERS_HPP
#define CARVAIBD_PARAMETERS_HPP

#include <boost/optional.hpp>
#include <armadillo>

/**
 * @brief Runtime parameters
 */
struct Parameters {
  unsigned long nperms;
  unsigned long nthreads;
  std::string output_path;
  unsigned int seed;
  boost::optional<arma::uword> lower_bound;
  bool swap;
  bool contcont;
  double min_dist;
  boost::optional<double> rsquared;
  bool verbose;
  bool enable_testing;
};

#endif //CARVAIBD_PARAMETERS_HPP
