//
// Created by Bohlender,Ryan James on 10/15/18.
//

#include <algorithm>

#include "gaussian.hpp"

const std::vector<std::string> Gaussian::links{"identity", "inverse", "log"};

const std::string Gaussian::family{"gaussian"};

Gaussian::Gaussian(const std::string &link)
    : linkid(check_linkid(link)), linkname(link) {}

arma::vec Gaussian::link(const arma::mat &X, const arma::vec &beta) noexcept {
  switch (linkid) {
  case Gaussian::LinkID::Identity:return X * beta;
  case Gaussian::LinkID::Inverse:return 1. / (X * beta);
  case Gaussian::LinkID::Log:return arma::log(1. / (X * beta));
  default:return arma::vec();
  }
}

arma::vec Gaussian::link(const arma::vec &mu) noexcept {
  switch (linkid) {
  case Gaussian::LinkID::Identity:return mu;
  case Gaussian::LinkID::Inverse:return 1. / (mu);
  case Gaussian::LinkID::Log:return arma::log(1. / (mu));
  default:return arma::vec();
  }
}

arma::vec Gaussian::linkinv(const arma::mat &X, const arma::vec &beta) noexcept {
  switch (linkid) {
  case Gaussian::LinkID::Identity:return X * beta;
  case Gaussian::LinkID::Inverse:return 1. / (X * beta);
  case Gaussian::LinkID::Log:return arma::exp(X * beta);
  default:return arma::vec();
  }
}

arma::vec Gaussian::linkinv(const arma::vec &eta) noexcept {
  switch (linkid) {
  case Gaussian::LinkID::Identity:return eta;
  case Gaussian::LinkID::Inverse:return 1. / (eta);
  case Gaussian::LinkID::Log:return arma::exp(eta);
  default:return arma::vec();
  }
}

arma::vec Gaussian::variance(const arma::vec &mu) noexcept {
  // All branches are identical. Can be replaced.
  switch (linkid) {
  case Gaussian::LinkID::Identity:return arma::vec(arma::size(mu), arma::fill::ones);
  case Gaussian::LinkID::Inverse:return arma::vec(arma::size(mu), arma::fill::ones);
  case Gaussian::LinkID::Log:return arma::vec(arma::size(mu), arma::fill::ones);
  default:return arma::vec();
  }
}

arma::vec Gaussian::mueta(const arma::vec &eta) noexcept {
  switch (linkid) {
  case Gaussian::LinkID::Identity:return arma::vec(arma::size(eta), arma::fill::ones);
  case Gaussian::LinkID::Inverse:return -1. / arma::pow(eta, 2);
  case Gaussian::LinkID::Log:return arma::exp(eta);
  default:return arma::vec();
  }
}

Gaussian::LinkID Gaussian::check_linkid(const std::string &link) {
  auto ok = std::find(links.cbegin(), links.cend(), link);
  if (ok == links.cend())
    throw (std::logic_error("Wrong link argument to Gaussian."));

  if (link == "identity") {
    return Gaussian::LinkID::Identity;
  } else if (link == "log") {
    return Gaussian::LinkID::Log;
  } else {
    return Gaussian::LinkID::Inverse;
  }
}

arma::vec Gaussian::dev_resids(const arma::vec &y, const arma::vec &mu, const arma::vec &weight) noexcept {
  return weight % arma::pow(y - mu, 2);
}



