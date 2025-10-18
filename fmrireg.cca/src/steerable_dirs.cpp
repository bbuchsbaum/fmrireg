// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat steerable_dirs_3d() {
  const double a = 2.0 / std::sqrt(10.0 + 2.0*std::sqrt(5.0));
  const double b = (1.0 + std::sqrt(5.0)) / std::sqrt(10.0 + 2.0*std::sqrt(5.0));
  arma::mat N(3,6);
  N.col(0) = arma::vec({  a, 0.0,  b});
  N.col(1) = arma::vec({ -a, 0.0,  b});
  N.col(2) = arma::vec({  b,   a, 0.0});
  N.col(3) = arma::vec({  b,  -a, 0.0});
  N.col(4) = arma::vec({ 0.0,  b,   a});
  N.col(5) = arma::vec({ 0.0,  b,  -a});
  return N;
}

