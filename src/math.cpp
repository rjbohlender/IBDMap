//
// Created by Bohlender,Ryan James on 11/12/20.
//

#include "math.hpp"

double cor(const arma::sp_vec &X, const arma::sp_vec &Y) {
  auto N = static_cast<double>(X.n_elem);
  double mxy = arma::as_scalar(arma::mean(X % Y)); // as_scalar requrired for compatibility with older versions
  double mx = arma::mean(X);
  double my = arma::mean(Y);
  double sdx = std::sqrt(arma::var(X));
  double sdy = std::sqrt(arma::var(Y));

  double r = N / (N - 1) * (mxy - mx * my) / (sdx * sdy);
  return std::pow(r, 2);
};
