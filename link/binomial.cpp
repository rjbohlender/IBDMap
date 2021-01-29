//
// Created by Bohlender,Ryan James on 10/15/18.
//

#include <algorithm>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/binomial.hpp>

#include "binomial.hpp"

const std::vector<std::string> Binomial::links {"logit", "probit", "cloglog", "Cauchit", "Log"};

const std::string Binomial::family {"binomial"};

Binomial::Binomial(const std::string &link)
: linkid(check_linkid(link)), linkname(link) {}

arma::vec Binomial::link(const arma::mat &X, const arma::vec &beta) noexcept {
  switch(linkid) {
	case Binomial::LinkID::Logit:
	  return arma::log(X * beta / (1. - X * beta));
	case Binomial::LinkID::Probit: {
	  arma::vec ret(beta.n_elem);
	  auto vit = ret.begin();
	  boost::math::normal dist(0.0, 1.0);
	  for (auto &v : arma::vec(X * beta)) {
		*vit = boost::math::quantile(dist, v);
		vit++;
	  }
	  return ret;
	}
	case Binomial::LinkID::cloglog:
	  return arma::log(-arma::log((1. - X * beta)));
	case Binomial::LinkID::Cauchit: {
	  arma::vec ret(beta.n_elem);
	  auto vit = ret.begin();
	  boost::math::cauchy dist(0.0, 1.0);
	  for (auto &v : arma::vec(X * beta)) {
		*vit = boost::math::quantile(dist, v);
		vit++;
	  }
	  return ret;
	}
	case Binomial::LinkID::Log:
	  return arma::log(X * beta);
  }
}

arma::vec Binomial::link(const arma::vec &mu) noexcept {
  switch(linkid) {
  case Binomial::LinkID::Logit:
	return arma::log(mu / (1. - mu));
  case Binomial::LinkID::Probit: {
	arma::vec ret(mu.n_elem);
	auto vit = ret.begin();
	boost::math::normal dist(0.0, 1.0);
	for (auto &v : arma::vec(mu)) {
	  *vit = boost::math::quantile(dist, v);
	  vit++;
	}
	return ret;
  }
  case Binomial::LinkID::cloglog:
	return arma::log(-arma::log((1. - mu)));
  case Binomial::LinkID::Cauchit: {
	arma::vec ret(mu.n_elem);
	auto vit = ret.begin();
	boost::math::cauchy dist(0.0, 1.0);
	for (auto &v : arma::vec(mu)) {
	  *vit = boost::math::quantile(dist, v);
	  vit++;
	}
	return ret;
  }
  case Binomial::LinkID::Log:
	return arma::log(mu);
  }
}

arma::vec Binomial::linkinv(const arma::mat &X, const arma::vec &beta) noexcept {
  switch(linkid) {
  case Binomial::LinkID::Logit:
	return 1. / (1. + arma::exp(-X * beta));
  case Binomial::LinkID::Probit:
    return arma::normcdf(X * beta);
  case Binomial::LinkID::cloglog: {
    arma::vec ret(beta.n_elem);
	auto vit = ret.begin();
	for (auto &v : arma::vec(X * beta)) {
	  *vit = -std::expm1(-std::exp(v));
	  vit++;
	}
	ret = arma::clamp(ret, std::numeric_limits<double>::min(), 1. - std::numeric_limits<double>::min());
	return ret;
  }
  case Binomial::LinkID::Cauchit: {
	arma::vec ret(beta.n_elem);
	auto vit = ret.begin();
	boost::math::cauchy dist(0.0, 1.0);
	for (auto &v : arma::vec(X * beta)) {
	  *vit = boost::math::cdf(dist, v);
	  vit++;
	}
	return ret;
  }
  case Binomial::LinkID::Log:
    return arma::exp(X * beta);
  }
}

arma::vec Binomial::linkinv(const arma::vec &eta) noexcept {
  switch(linkid) {
  case Binomial::LinkID::Logit:
	return 1. / (1. + arma::exp(-eta));
  case Binomial::LinkID::Probit:
	return arma::normcdf(eta);
  case Binomial::LinkID::cloglog: {
	arma::vec ret(eta.n_elem);
	auto vit = ret.begin();
	for (auto &v : arma::vec(eta)) {
	  *vit = -std::expm1(-std::exp(v));
	  vit++;
	}
	ret = arma::clamp(ret, std::numeric_limits<double>::min(), 1. - std::numeric_limits<double>::min());
	return ret;
  }
  case Binomial::LinkID::Cauchit: {
	arma::vec ret(eta.n_elem);
	auto vit = ret.begin();
	boost::math::cauchy dist(0.0, 1.0);
	for (auto &v : arma::vec(eta)) {
	  *vit = boost::math::cdf(dist, v);
	  vit++;
	}
	return ret;
  }
  case Binomial::LinkID::Log:
	return arma::exp(eta);
  }
}

arma::vec Binomial::variance(const arma::vec &mu) noexcept {
  return mu % (1. - mu);
}

arma::vec Binomial::mueta(const arma::vec &eta) noexcept {
  switch(linkid) {
  case LinkID::Logit:
	return arma::exp(eta) / arma::pow((arma::exp(eta) + 1.), 2);
  case LinkID::Probit:
    return arma::normpdf(eta);
  case LinkID::cloglog: {
    arma::vec re_eta = arma::clamp(eta, arma::min(eta), 700);
    arma::vec ret = arma::exp(re_eta) * arma::exp(-arma::exp(re_eta));
	return arma::clamp(ret, std::numeric_limits<double>::min(), arma::max(ret));
  }
  case LinkID::Cauchit: {
    boost::math::cauchy dist(0., 1.);
	arma::vec ret(eta.n_elem);
	auto vit = ret.begin();
	for (auto &v : eta) {
	  *vit = boost::math::pdf(dist, v);
	  vit++;
	}
	ret = arma::clamp(ret, std::numeric_limits<double>::min(), 1. - std::numeric_limits<double>::min());
  }
  case LinkID::Log:
    return arma::clamp(arma::exp(eta), std::numeric_limits<double>::min(), 1. - std::numeric_limits<double>::min());
  }
}

Binomial::LinkID Binomial::check_linkid(const std::string &link) {
  auto ok = std::find(links.cbegin(), links.cend(), link);
  if(ok == links.cend())
    throw(std::logic_error("Wrong link argument to Binomial."));

  if(link == "logit") {
    return Binomial::LinkID::Logit;
  } else if(link == "probit") {
    return Binomial::LinkID::Probit;
  } else if(link == "cloglog") {
    return Binomial::LinkID::cloglog;
  } else if(link == "cauchit") {
    return Binomial::LinkID::Cauchit;
  } else {
    return Binomial::LinkID::Log;
  }
}

arma::vec Binomial::dev_resids(const arma::vec &y, const arma::vec &mu, const arma::vec &weight) noexcept {
  auto y_log_y = [](const arma::vec &y, const arma::vec &mu) {
    arma::uvec idx = arma::find(y > 0);
    arma::vec ret(y.n_elem, arma::fill::zeros);
    ret(idx) = y(idx) % arma::log(y(idx) / mu(idx));
    return ret;
  };
  arma::vec ret = 2. * weight % (y_log_y(y, mu) + y_log_y((1. - y), (1. - mu)));
  return ret;
}

void Binomial::initialize(arma::vec &y, arma::vec &n, arma::vec &mu, arma::vec &weight) noexcept {
  n = arma::vec(y.n_elem, arma::fill::ones);
  mu = (weight % y + 0.5) / (weight + 1);
}

double Binomial::aic(arma::vec &y, arma::vec &n, arma::vec &mu, arma::vec &weight, double dev, double rank) noexcept {
  // aic <- function(y, n, mu, wt, dev) {
  //        m <- if(any(n > 1)) n else wt
  //	-2*sum(ifelse(m > 0, (wt/m), 0)*
  //               dbinom(round(m*y), round(m), mu, log=TRUE))
  //    }
  arma::vec m;
  if(arma::any(n > 1)) {
    m = n;
  } else {
    m = weight;
  }
  double aic_val = 0;
  for(int i = 0; i < y.n_elem; i++) {
    boost::math::binomial_distribution<> binom(boost::math::round(m(i)), mu(i));
    aic_val += -2 * std::log(boost::math::pdf(binom, boost::math::round(m(i) * y(i))));
  }
  return aic_val + 2 * rank;
}
