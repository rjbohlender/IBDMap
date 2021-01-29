//
// Created by Bohlender,Ryan James on 10/15/18.
//

#ifndef PERMUTE_ASSOCIATE_BINOMIAL_HPP
#define PERMUTE_ASSOCIATE_BINOMIAL_HPP

#include "family.hpp"

struct Binomial : Family {
  enum class LinkID {
    Logit,
    Probit,
    cloglog,
    Cauchit,
    Log
  };

  static const std::vector<std::string> links;
  static const std::string family;
  const LinkID linkid;
  const std::string linkname;

  explicit Binomial(const std::string &link="logit");
  arma::vec link(const arma::mat &X, const arma::vec &beta) noexcept override;
  arma::vec link(const arma::vec &mu) noexcept override;
  arma::vec linkinv(const arma::mat &X, const arma::vec &beta) noexcept override;
  arma::vec linkinv(const arma::vec &eta) noexcept override;
  arma::vec variance(const arma::vec &mu) noexcept override;
  arma::vec mueta(const arma::vec &eta) noexcept override;
  arma::vec dev_resids(const arma::vec &y, const arma::vec &mu, const arma::vec &weight) noexcept override;
  void initialize(arma::vec &y, arma::vec &n, arma::vec &mu, arma::vec &weight) noexcept override;
  double aic(arma::vec &y, arma::vec &n, arma::vec &mu, arma::vec &weight, double dev, double rank) noexcept override;

  LinkID check_linkid(const std::string &link);
};

#endif //PERMUTE_ASSOCIATE_BINOMIAL_HPP
