//
// Created by Bohlender,Ryan James on 9/24/19.
//

#ifndef CARVAIBD_PARAMETERS_HPP
#define CARVAIBD_PARAMETERS_HPP

/**
 * @brief Runtime parameters
 */
struct Parameters {
  unsigned long nperms;
  unsigned long nthreads;
  std::string output_path;
  unsigned int seed;
};

#endif //CARVAIBD_PARAMETERS_HPP
