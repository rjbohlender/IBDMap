//
// Created by Bohlender,Ryan James on 10/15/18.
//

#ifndef PERMUTE_ASSOCIATE_LINK_HPP
#define PERMUTE_ASSOCIATE_LINK_HPP

#define ARMA_DONT_USE_WRAPPER

#include <armadillo>

struct Family {
  virtual arma::vec link(const arma::mat &X, const arma::vec &beta) noexcept = 0;
  virtual arma::vec link(const arma::vec &mu) noexcept = 0;
  virtual arma::vec linkinv(const arma::mat &X, const arma::vec &beta) noexcept = 0;
  virtual arma::vec linkinv(const arma::vec &eta) noexcept = 0;
  virtual arma::vec variance(const arma::vec &mu) noexcept = 0;
  virtual arma::vec mueta(const arma::vec &eta) noexcept = 0;
  virtual arma::vec dev_resids(const arma::vec &y, const arma::vec &mu, const arma::vec &weight) noexcept = 0;
  virtual void initialize(arma::vec &y, arma::vec &n, arma::vec &mu, arma::vec &weight) noexcept = 0;
  virtual double aic(arma::vec &y, arma::vec &n, arma::vec &mu, arma::vec &weight, double dev, double rank) noexcept = 0;
};

#endif //PERMUTE_ASSOCIATE_LINK_HPP
