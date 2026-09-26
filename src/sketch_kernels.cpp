// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

static inline bool inv_sympd_safe(mat& out, const mat& A) {
  try {
    out = inv_sympd(A);
    return true;
  } catch (...) {
    try {
      out = pinv(A);
      return true;
    } catch (...) {
      return false;
    }
  }
}

// In-place Walsh-Hadamard transform on each column; n must be power of 2
static void fwht_cols(mat& X) {
  int n = X.n_rows;
  for (int len = 1; len < n; len <<= 1) {
    int step = len << 1;
    for (int i = 0; i < n; i += step) {
      for (int j = 0; j < len; ++j) {
        rowvec a = X.row(i + j);
        rowvec b = X.row(i + j + len);
        X.row(i + j)        = a + b;
        X.row(i + j + len)  = a - b;
      }
    }
  }
}

// SRHT apply: M (T x k), rows (m), signs (T), perm (T), scale scalar
// [[Rcpp::export]]
arma::mat cpp_srht_apply(const arma::mat& M,
                         const arma::uvec& rows,
                         const arma::vec& signs,
                         const arma::uvec& perm,
                         const double scale) {
  int T = M.n_rows, K = M.n_cols;
  if ((int)signs.n_elem != T) {
    Rcpp::stop("cpp_srht_apply: length(signs) must equal nrow(M).");
  }
  if ((int)perm.n_elem != T) {
    Rcpp::stop("cpp_srht_apply: length(perm) must equal nrow(M).");
  }
  // 1) D (random signs)
  mat X = M.each_col() % signs;
  // 2) H (Hadamard) on power-of-two padded length
  int T2 = 1; while (T2 < T) T2 <<= 1;
  if (perm.n_elem > 0 && perm.max() >= (uword)T2) {
    Rcpp::stop("cpp_srht_apply: perm contains out-of-bounds indices.");
  }
  if (rows.n_elem > 0 && rows.max() >= (uword)T) {
    Rcpp::stop("cpp_srht_apply: rows contains out-of-bounds indices.");
  }
  mat Xpad(T2, K, fill::zeros);
  Xpad.rows(0, T-1) = X;
  fwht_cols(Xpad);
  // 3) P (permute)
  mat XP(T, K);
  for (int i = 0; i < T; ++i) XP.row(i) = Xpad.row(perm(i));
  // 4) R (row sample) + scale
  mat out(rows.n_elem, K);
  for (uword i = 0; i < rows.n_elem; ++i) out.row(i) = XP.row(rows(i)) * scale;
  return out;
}

// One IHS iteration helper
static void ihs_iter(const mat& X, const mat& Z, int m, mat& M, mat& Ginv_out) {
  int T = X.n_rows;
  if (m <= 0 || m > T) {
    Rcpp::stop("ihs_iter: sketch size m must satisfy 0 < m <= nrow(X).");
  }
  // Build SRHT plan
  arma::vec signs = 2.0 * randu<vec>(T) - 1.0; signs.transform( [](double v){ return v>=0 ? 1.0 : -1.0; } );
  arma::uvec perm = randperm(T);
  arma::uvec order = sort_index(randu<vec>(T));
  arma::uvec rows = order.subvec(0, m - 1);
  double scale = std::sqrt( (double)T / (double)m );
  mat Xs = cpp_srht_apply(X, rows, signs, perm, scale);
  // cpp_srht_apply uses an unnormalised Hadamard transform, so
  // G = Xs' Xs estimates T * X' X. The returned Ginv keeps that sketched
  // scale because callers pair it with residual variances computed from
  // SRHT-sketched residuals (which carry the same factor T).
  mat G = Xs.t() * Xs;
  mat Ginv;
  if (!inv_sympd_safe(Ginv, G)) {
    Rcpp::stop("ihs_iter: unable to invert sketched Gram matrix.");
  }
  // Iterative Hessian sketch (Pilanci & Wainwright, 2016): only the Hessian
  // is sketched; the gradient X'(Z - XM) uses the full data. Sketching the
  // gradient as well makes every iteration re-solve an independent sketched
  // problem, so the iterates never converge to the least-squares solution.
  mat E = Z - X * M;
  mat dM = ((double)T * Ginv) * (X.t() * E);
  // dM is a descent direction for ||Z - XM||_F^2 (Ginv is positive
  // definite), but a poorly conditioned sketch can overshoot. Halve the step
  // until the full-data residual sum of squares does not increase, so the
  // iterates are monotone and more iterations never move away from the
  // least-squares solution.
  const double rss0 = accu(square(E));
  double step = 1.0;
  for (int k = 0; k < 30; ++k) {
    if (accu(square(Z - X * (M + step * dM))) <= rss0) break;
    step *= 0.5;
  }
  M += step * dM;
  Ginv_out = Ginv;
}

// IHS latent solve: returns M and Ginv after iters
// [[Rcpp::export]]
Rcpp::List cpp_ihs_latent(const arma::mat& X, const arma::mat& Z,
                          const int m, const int iters) {
  int T = X.n_rows, p = X.n_cols;
  if (iters < 1) {
    Rcpp::stop("cpp_ihs_latent: iters must be >= 1.");
  }
  if (m <= 0 || m > T) {
    Rcpp::stop("cpp_ihs_latent: sketch size m must satisfy 0 < m <= nrow(X).");
  }
  if (!X.is_finite() || !Z.is_finite()) {
    Rcpp::stop("cpp_ihs_latent: X and Z must be finite.");
  }
  // Warm start from the classical sketch-and-solve solution. IHS contracts
  // the error relative to its starting point, and fMRI data carry a large
  // baseline that a zero start would leave almost entirely unfitted after the
  // default few iterations; sketch-and-solve fits any part of Z in span(X)
  // exactly, so iterations only refine the noise-driven error.
  arma::vec signs = 2.0 * randu<vec>(T) - 1.0;
  signs.transform( [](double v){ return v>=0 ? 1.0 : -1.0; } );
  arma::uvec perm = randperm(T);
  arma::uvec order = sort_index(randu<vec>(T));
  arma::uvec rows = order.subvec(0, m - 1);
  double scale = std::sqrt( (double)T / (double)m );
  mat Xs = cpp_srht_apply(X, rows, signs, perm, scale);
  mat Zs = cpp_srht_apply(Z, rows, signs, perm, scale);
  mat Ginv;
  if (!inv_sympd_safe(Ginv, Xs.t() * Xs)) {
    Rcpp::stop("cpp_ihs_latent: unable to invert sketched Gram matrix.");
  }
  mat M = Ginv * (Xs.t() * Zs);
  for (int t = 0; t < iters; ++t) {
    ihs_iter(X, Z, m, M, Ginv);
  }
  return Rcpp::List::create(
    Rcpp::Named("M")    = M,
    Rcpp::Named("Ginv") = Ginv
  );
}
