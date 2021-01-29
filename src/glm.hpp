//
// Created by Bohlender,Ryan James on 10/16/18.
//

#ifndef PERMUTE_ASSOCIATE_GLM_HPP
#define PERMUTE_ASSOCIATE_GLM_HPP

#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <boost/math/tools/toms748_solve.hpp>

#include "../link/family.hpp"
#include "parameters.hpp"

template <typename LinkT>
struct GLM {
  LinkT link;
  arma::vec beta_; // coefficients
  arma::vec mu_;   // fitted.values
  arma::vec eta_;  // linear.predictors
  arma::vec var_; // variance
  arma::vec weights_;
  arma::vec n_;
  double dev_;
  arma::uvec indices_;
  bool success_;

  GLM(arma::mat &X, arma::vec &Y, LinkT &link, Parameters &params) {
    indices_ = arma::regspace<arma::uvec>(0, X.n_cols - 1);
    weights_ = arma::vec(X.n_rows, arma::fill::ones);

    if (X.n_cols == 0) {
      eta_ = arma::vec(X.n_rows, arma::fill::zeros);
      mu_ = link.linkinv(eta_);
      dev_ = arma::accu(link.dev_resids(Y, mu_, arma::vec(X.n_rows, arma::fill::ones)));
      return;
    } else {
      link.initialize(Y, n_, mu_, weights_);
      eta_ = link.linkinv(mu_);
      dev_ = arma::accu(link.dev_resids(Y, mu_, arma::vec(X.n_rows, arma::fill::ones)));
    }

	if (params.optimizer == "irls") {
	  beta_ = irls(X, Y);
	} else if(params.optimizer == "irls_svdnewton") {
	  beta_ = irls_svdnewton(X, Y);
	} else if(params.optimizer == "irls_qr") {
	  beta_ = irls_qr(X, Y);
	} else if(params.optimizer == "irls_qr_R") {
	  beta_ = irls_qr_R(X, Y);
	} else if(params.optimizer == "gradient_descent") {
	  beta_ = gradient_descent(X, Y);
	}

	mu_ = link.linkinv(X.cols(indices_), beta_);
    eta_ = link.link(mu_);

	success_ = !(mu_.has_nan() || eta_.has_nan());
  }

  GLM(arma::mat &X, arma::vec &Y, LinkT &link, arma::vec initial_beta, Parameters params) {
    indices_ = arma::regspace<arma::uvec>(0, X.n_cols - 1);
    weights_ = arma::vec(X.n_rows, arma::fill::ones);

    if (initial_beta.n_rows < X.n_cols) {
      initial_beta.insert_rows(initial_beta.n_rows, X.n_cols - initial_beta.n_rows); // Zeroed by default
    }

    if (X.n_cols == 0) {
      eta_ = arma::vec(X.n_rows, arma::fill::zeros);
      mu_ = link.linkinv(eta_);
      dev_ = arma::accu(link.dev_resids(Y, mu_, arma::vec(X.n_rows, arma::fill::ones)));
      return;
    } else {
      mu_ = link.linkinv(X, initial_beta);
      eta_ = link.link(mu_);
      dev_ = arma::accu(link.dev_resids(Y, mu_, arma::vec(X.n_rows, arma::fill::ones)));
    }

    if (params.optimizer == "irls") {
	  beta_ = irls(X, Y);
    } else if(params.optimizer == "irls_svdnewton") {
      beta_ = irls_svdnewton(X, Y);
    } else if(params.optimizer == "irls_qr") {
	  beta_ = irls_qr(X, Y);
	} else if(params.optimizer == "irls_qr_R") {
	  beta_ = irls_qr_R(X, Y);
	} else if(params.optimizer == "gradient_descent") {
	  beta_ = gradient_descent(X, Y);
	}

	mu_ = link.linkinv(X.cols(indices_), beta_);
    eta_ = link.link(mu_);

	success_ = !(mu_.has_nan() || eta_.has_nan());
  }
  // Algorithms for finding the optimum
  auto gradient_descent(arma::mat &X, arma::colvec &Y) -> arma::vec;
  auto irls_svdnewton(arma::mat &X, arma::colvec &Y) -> arma::vec;
  auto irls_qr(arma::mat &X, arma::colvec &Y) -> arma::vec;
  auto irls_qr_R(arma::mat &X, arma::colvec &Y) -> arma::vec;
  auto irls(arma::mat &X, arma::colvec &Y) -> arma::vec;
  auto svd_subset(arma::mat &X) -> arma::uvec;
  auto newton_rhapson(arma::mat &X, arma::vec &Y) -> arma::vec;
};

template<typename LinkT>
auto GLM<LinkT>::gradient_descent(arma::mat &X, arma::colvec &Y) -> arma::vec {
  // std::cerr << "Running gradient descent.\n";
  auto iterations = 0ull;
  const auto max_iter = 1000000ull;
  auto alpha = 0.0005; // Learning rate
  auto tol = 1e-10;
  auto m = static_cast<double>(X.n_cols);
  arma::mat A = X;
  auto b = arma::vec(A.n_cols, arma::fill::randn);
  auto grad = b;
  do {
    // Vectorized update
    grad = alpha * (A.t() * (link.linkinv(A, b) - Y)) / m;
    b -= grad;

    iterations++;
  } while(iterations < max_iter && arma::norm(grad) > tol);

  dev_ = arma::sum(link.dev_resids(Y, link.linkinv(A, b), arma::vec(A.n_rows, arma::fill::ones)));

  return b;
}

template<typename LinkT>
auto GLM<LinkT>::irls_svdnewton(arma::mat &X, arma::colvec &Y) -> arma::vec {
  const auto tol = 1e-25;
  const auto max_iter = 50;
  auto iter = 0;

  arma::uword m = X.n_rows;
  arma::uword n = X.n_cols;

  // SVD
  arma::mat U, V;
  arma::vec S;

  bool success = arma::svd_econ(U, S, V, X, "both", "dc");
  if(!success) {
    arma::svd_econ(U, S, V, X, "both", "std");
  }

  // Matrices and Vectors
  eta_ = arma::vec (m, arma::fill::randn);
  arma::vec s(n, arma::fill::randn);
  arma::vec weights(m, arma::fill::ones);
  arma::vec s_old = s;
  arma::uvec good;
  arma::mat Ugood;
  arma::mat Wgood;

  // Handle rank deficiency
  arma::uvec tiny = (S / S(0) < tol);
  if(arma::sum(tiny) > 0) {
    std::cerr << "Dropping columns to correct for rank deficiency" << std::endl;
	indices_ = svd_subset(X);
	std::cerr << "indices: " << indices_.t();
	arma::svd_econ(U, S, V, X.cols(indices_));
  }

  good = arma::find(weights > arma::datum::eps * 2);
  arma::vec diff = s - s_old;

  do {
    arma::vec g = link.linkinv(eta_(good));
    arma::vec varg = link.variance(g);
    arma::vec gprime = link.mueta(eta_(good));

    arma::vec z(m);
    arma::vec W(m);

    z(good) = eta_(good) + (Y(good) - g) / gprime;
    W(good) = weights(good) % arma::pow(gprime, 2) / varg;

	good = arma::find(W > arma::datum::eps * 2);
	if(sum(W > arma::datum::eps * 2) < m)
	  std::cerr << "Warning: Encountered tiny weights in IRLS.\n";

	arma::mat C;
	Ugood = U.rows(good);
	Wgood = W(good);
    arma::mat UWU = Ugood.t() * (Ugood.each_col() % Wgood);
    success = arma::chol(C, UWU);
    if (!success) {
      throw(std::runtime_error("Cholesky decomposition of UWU matrix failed."));
    }

    s_old = s;
    s = solve(arma::trimatl(C.t()), Ugood.t() * (Wgood % z));
	s = solve(arma::trimatu(C), s);

    eta_ = Ugood * s;
    mu_ = link.linkinv(eta_);

    iter++;

    dev_ = arma::sum(link.dev_resids(Y, g, weights));

  // } while(iter < max_iter && std::abs(dev_ - devold) / (0.1 + std::abs(dev_)) > tol);
  } while(iter < max_iter && arma::accu(arma::sqrt(diff.t() * diff)) > tol);

  if(arma::accu(arma::sqrt(diff.t() * diff)) > tol) {
    std::cerr << "IRLS failed to converge." << std::endl;
  }

  dev_ = arma::accu(link.dev_resids(Y, mu_, weights_));
  beta_ = (V * (arma::diagmat(1. / S) * (Ugood.t() * eta_(good))));
  return beta_;
}
template<typename LinkT>

auto GLM<LinkT>::svd_subset(arma::mat &X) -> arma::uvec {
  arma::mat U, Vt;
  arma::vec S;
  arma::mat A = X;
  // Rescale columns of X to unit variance
  arma::rowvec sds = arma::sqrt(arma::pow(arma::sum(A.each_row() - arma::mean(A), 0), 2));

  arma::svd_econ(U, S, Vt, A);
  arma::uvec n = arma::find(S < 2 * std::numeric_limits<double>::epsilon());

  arma::uword k = X.n_cols;
  if(n.n_elem > 0 and n.n_elem < k) {
    k = n.n_elem - 1;
  }
  // Find the ordering
  arma::vec lens(Vt.n_cols, arma::fill::zeros);
  for(arma::uword i = 0; i < Vt.n_cols; i++) {
	lens(i) = arma::norm(Vt.col(i), 1);
  }
  arma::uvec idx = arma::sort_index(lens, "descend");
  return idx(arma::span(0, k-1));
}

template<typename LinkT>
auto GLM<LinkT>::irls_qr(arma::mat &X, arma::colvec &Y) -> arma::vec {
  const auto tol = 1e-8;
  const auto max_iter = 25;
  auto iter = 0;

  arma::uword m = X.n_rows;
  arma::uword n = X.n_cols;

  // A = QR
  arma::mat Q, R;
  bool qr_success = arma::qr_econ(Q, R, X);
  if (!qr_success) {
    qr_success = arma::qr(Q, R, X);
    if (!qr_success) {
      throw(std::runtime_error("Failed QR decomposition of X."));
    }
  }

  // Matrices and Vectors
  arma::vec s(n, arma::fill::zeros);
  double devold;

  do {
	mu_ = link.linkinv(eta_);
	var_ = link.variance(mu_);
    if(!arma::is_finite(var_)) {
      throw(std::runtime_error("Non-finite values in V(mu)."));
    }
	if(arma::any(var_ == 0)) {
	  throw(std::runtime_error("0s in V(mu)."));
	}
	arma::vec mueta = link.mueta(eta_);

	arma::vec z = eta_ + (Y - mu_) / mueta;
	arma::vec W = arma::sqrt(weights_ % arma::pow(mueta, 2) / var_);

	arma::mat Lt;
	arma::mat QWQ = Q.t() * arma::diagmat(W) * Q;
    bool chol_success = arma::chol(Lt, QWQ);
	double jitter = 1e-9;
	while(!chol_success && jitter < 1.) {
	  QWQ.diag() += jitter;
      chol_success = arma::chol(Lt, QWQ);
	  jitter *= 10;
	}
	if (!chol_success) {
	  throw(std::runtime_error("Cholesky decomposition of QWQ matrix failed."));
	}

	s = arma::solve(arma::trimatl(Lt.t()), Q.t() * (W % z));
	s = arma::solve(arma::trimatu(Lt), s);

    eta_ = Q * s;
	mu_ = link.linkinv(eta_);

	iter++;

	devold = dev_;
	dev_ = arma::accu(link.dev_resids(Y, mu_, weights_));

  } while(iter < max_iter && std::abs(dev_ - devold) / (0.1 + std::abs(dev_)) > tol);

  if(std::abs(dev_ - devold) / (0.1 + std::abs(dev_)) > tol) {
	std::cerr << "IRLS failed to converge." << std::endl;
  }

  return arma::solve(R.t(), Q.t() * eta_);
}

template<typename LinkT>
auto GLM<LinkT>::irls_qr_R(arma::mat &X, arma::colvec &Y) -> arma::vec {
  const auto tol = 1e-8;
  const auto max_iter = 25;
  auto iter = 0;

  arma::uword m = X.n_rows;
  arma::uword n = X.n_cols;

  // Matrices and Vectors
  eta_ = link.link(mu_);
  arma::vec s(n, arma::fill::zeros);

  dev_ = arma::accu(link.dev_resids(Y, mu_, weights_));
  double devold = dev_;

  do {
    mu_ = link.linkinv(eta_);
    var_ = link.variance(mu_);
    if(!arma::is_finite(var_)) {
      throw(std::runtime_error("Non-finite values in V(mu)."));
    }
    if(arma::any(var_ == 0)) {
      throw(std::runtime_error("0s in V(mu)."));
    }
    arma::vec mueta = link.mueta(eta_);
    if(arma::any(mueta == 0)) {
      throw(std::runtime_error("0s in mueta."));
    }

    arma::vec z = eta_ + (Y - mu_) / mueta;
    arma::vec W = arma::sqrt(weights_ % arma::pow(mueta, 2) / var_);

    // A = QR
    arma::mat Q, R;
    bool qr_success = arma::qr_econ(Q, R, X.each_col() % W);
    if (!qr_success) {
      qr_success = arma::qr(Q, R, X.each_col() % W);
      if (!qr_success) {
        throw(std::runtime_error("Failed QR decomposition of X."));
      }
    }

    arma::vec zw = z % W;
    beta_ = arma::solve(arma::trimatu(R), Q.t() * zw);

    eta_ = X * beta_;
    mu_ = link.linkinv(eta_);

    devold = dev_;
    dev_ = arma::accu(link.dev_resids(Y, mu_, weights_));

    iter++;
  } while(iter < max_iter && std::abs(dev_ - devold) / (0.1 + std::abs(dev_)) > tol);

  if(std::abs(dev_ - devold) / (0.1 + std::abs(dev_)) > tol) {
    std::cerr << "IRLS failed to converge." << std::endl;
  }

  return beta_;
}

template<typename LinkT>
auto GLM<LinkT>::irls(arma::mat &X, arma::colvec &Y) -> arma::vec {
  const auto tol = 1e-8;
  const auto max_iter = 25;
  auto iter = 0;

  beta_ = arma::vec(X.n_cols, arma::fill::zeros);
  arma::vec xold;
  double devold;
  for(iter = 0; iter < max_iter; iter++) {
    eta_ = X * beta_;
    mu_ = link.linkinv(eta_);
	arma::vec muprime = link.mueta(eta_);
	arma::vec z = eta_ + (Y - mu_) / muprime;
	arma::vec W = arma::pow(muprime, 2) / link.variance(mu_);
	xold = beta_;
	beta_ = arma::solve(X.t() * (X.each_col() % W), X.t() * (W % z));
	devold = dev_;
    dev_ = arma::sum(link.dev_resids(Y, link.linkinv(X * beta_), arma::vec(X.n_rows, arma::fill::ones)));
    if(!std::isfinite(dev_)) {
      int i = 0;
      while (!std::isfinite(dev_)) {
        if(i >= max_iter) {
          throw(std::runtime_error("Unable to correct step size for divergent IRLS step."));
        }
		beta_ = (beta_ + xold) / 2.;
		eta_ = X * beta_;
		mu_ = link.linkinv(eta_);
        dev_ = arma::sum(link.dev_resids(Y, mu_, arma::vec(X.n_rows, arma::fill::ones)));
        i++;
      }
    }
	//if(arma::norm(x - xold) < tol) {
    if(std::abs(dev_ - devold) / (0.1 + std::abs(dev_)) < tol) {
	  break;
	}
  }
  return beta_;
}
template<typename LinkT>
auto GLM<LinkT>::newton_rhapson(arma::mat &X, arma::vec &Y) -> arma::vec {
  int max_iter = 100;
  int iter = 0;
  double last_dev = -99999;
  beta_ = arma::vec(X.n_cols, arma::fill::zeros);
  while (iter < max_iter) {
    eta_ = X * beta_;
    mu_ = link.linkinv(eta_);
    var_ = link.variance(mu_);

    arma::mat D = X.t() * arma::diagmat(var_) * X;
    arma::vec r = X.t() * (Y - mu_);
    arma::vec delta_beta = arma::solve(D, r);

    beta_ += delta_beta;
    dev_ = arma::accu(link.dev_resids(Y, mu_, arma::vec(X.n_rows, arma::fill::ones)));
    if(std::abs(dev_ - last_dev) < 1e-3) {
      iter = 0;
      break;
    }

    last_dev = dev_;
    iter++;
  }
  return beta_;
}

#endif //PERMUTE_ASSOCIATE_GLM_HPP
