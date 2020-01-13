//
// Created by Bohlender,Ryan James on 10/16/18.
//

#ifndef PERMUTE_ASSOCIATE_GLM_HPP
#define PERMUTE_ASSOCIATE_GLM_HPP

#define ARMA_DONT_USE_WRAPPER

#include <armadillo>

#include "../link/family.hpp"

template <typename LinkT>
struct GLM {
  LinkT link;
  arma::vec beta_; // coefficients
  arma::vec mu_;   // fitted.values
  arma::vec eta_;  // linear.predictors
  double dev_;
  arma::uvec indices_;
  bool success;

  GLM(arma::mat &X, arma::vec &Y, LinkT &link) {
    indices_ = arma::regspace<arma::uvec>(0, X.n_cols - 1);
	try{
      // beta_ = gradient_descent(X, Y);
      // beta_ = irls_svdnewton(X, Y);
      // beta_ = irls_qr(X, Y);
      beta_ = irls(X, Y);
	} catch(std::exception &e) {
      std::cerr << "IRLS failed; Using gradient descent." << std::endl;
      //beta_ = gradient_descent(X, Y);
      //std::cerr << e.what();
      try {
		beta_ = gradient_descent(X, Y);
		std::cerr << beta_;
	  } catch(std::exception &e) {
		beta_ = arma::vec(X.n_cols);
		mu_ = arma::vec(X.n_rows);
		eta_ = arma::vec(X.n_rows);
		beta_.fill(arma::datum::nan);
		mu_.fill(arma::datum::nan);
		eta_.fill(arma::datum::nan);
	  }
    }

    mu_ = link.linkinv(X.cols(indices_), beta_);
    eta_ = link.link(mu_);

	success = !(mu_.has_nan() || eta_.has_nan());
  }
  // Algorithms for finding the optimum
  auto gradient_descent(arma::mat &X, arma::colvec &Y) -> arma::vec;
  auto irls_svdnewton(arma::mat &X, arma::colvec &Y) -> arma::vec;
  auto irls_qr(arma::mat &X, arma::colvec &Y) -> arma::vec;
  auto irls(arma::mat &X, arma::colvec &Y) -> arma::vec;
  auto svd_subset(arma::mat &X) -> arma::uvec;
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
  arma::vec eta(m, arma::fill::randn);
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

    arma::vec g = link.linkinv(eta(good));
    arma::vec varg = link.variance(g);
    arma::vec gprime = link.mueta(eta(good));

    arma::vec z(m);
    arma::vec W(m);

    z(good) = eta(good) + (Y(good) - g) / gprime;
    W(good) = weights(good) % arma::pow(gprime, 2) / varg;

	good = arma::find(W > arma::datum::eps * 2);
	if(sum(W > arma::datum::eps * 2) < m)
	  std::cerr << "Warning: Encountered tiny weights in IRLS.\n";

	arma::mat C;
	Ugood = U.rows(good);
	Wgood = W(good);
    arma::mat UWU = Ugood.t() * (Ugood.each_col() % Wgood);
    success = arma::chol(C, UWU);
#if 0
    double jitter = 1e-9;
    while(!success && jitter < 1.) {
      UWU.diag() += jitter;
      success = arma::chol(C, UWU);
      jitter *= 10;
    }
#endif
    if (!success) {
      throw(std::runtime_error("Cholesky decomposition of UWU matrix failed."));
    }

    s_old = s;
    s = solve(arma::trimatl(C.t()), Ugood.t() * (Wgood % z));
	s = solve(arma::trimatu(C), s);

    eta = Ugood * s;

    iter++;

    dev_ = arma::sum(link.dev_resids(Y, g, weights));

  // } while(iter < max_iter && std::abs(dev_ - devold) / (0.1 + std::abs(dev_)) > tol);
  } while(iter < max_iter && arma::accu(arma::sqrt(diff.t() * diff)) > tol);

  if(arma::accu(arma::sqrt(diff.t() * diff)) > tol) {
    std::cerr << "IRLS failed to converge." << std::endl;
  }

  return (V * (arma::diagmat(1. / S) * (Ugood.t() * eta(good))));
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

  // SVD
  arma::mat Q, R;

  bool success = arma::qr_econ(Q, R, X);

  // Matrices and Vectors
  arma::vec eta(m, arma::fill::randn);
  arma::vec s(n, arma::fill::randn);
  arma::vec weights(m, arma::fill::ones);

  dev_ = 0;
  double devold;

  do {
	devold = dev_;

	arma::vec g = link.linkinv(eta);
	arma::vec varg = link.variance(g);
	arma::vec gprime = link.mueta(eta);

	arma::vec z(m);
	arma::vec W(m);

	z = eta + (Y - g) / gprime;
	W = weights % arma::pow(gprime, 2) / varg;

	arma::mat C;
	arma::mat QWQ = Q.t() * (Q.each_col() % W);
	success = arma::chol(C, QWQ);
	double jitter = 1e-9;
	while(!success && jitter < 1.) {
	  QWQ.diag() += jitter;
	  success = arma::chol(C, QWQ);
	  jitter *= 10;
	}
	if (!success) {
	  throw(std::runtime_error("Cholesky decomposition of UWU matrix failed."));
	}

	s = arma::solve(arma::trimatl(C.t()), Q.t() * (W % z));
	s = arma::solve(arma::trimatu(C), s);

	eta = Q * s;

	iter++;

	dev_ = arma::sum(link.dev_resids(Y, g, weights));

  } while(iter < max_iter && std::abs(dev_ - devold) / (0.1 + std::abs(dev_)) > tol);

  if(std::abs(dev_ - devold) / (0.1 + std::abs(dev_)) <= tol) {
	std::cerr << "IRLS failed to converge." << std::endl;
  }

  return arma::solve(R.t(), Q.t() * eta);
}

template<typename LinkT>
auto GLM<LinkT>::irls(arma::mat &X, arma::colvec &Y) -> arma::vec {
  const auto tol = 1e-8;
  const auto max_iter = 25;
  auto iter = 0;

  arma::vec x(X.n_cols, arma::fill::zeros);
  arma::vec xold;
  for(iter = 0; iter < max_iter; iter++) {
    arma::vec eta = X * x;
    arma::vec g = link.linkinv(eta);
	arma::vec gprime = link.mueta(eta);
	arma::vec z = eta + (Y - g) / gprime;
	arma::vec W = arma::pow(gprime, 2) / link.variance(g);
	xold = x;
	x = arma::solve(X.t() * (X.each_col() % W), X.t() * (W % z));
	if(arma::norm(x - xold) < tol) {
	  break;
	}
  }
  dev_ = arma::sum(link.dev_resids(Y, link.linkinv(X * x), arma::vec(X.n_rows, arma::fill::ones)));
  return x;
}

#endif //PERMUTE_ASSOCIATE_GLM_HPP
