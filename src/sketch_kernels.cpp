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

// In-place Walsh-Hadamard transform on each column; n must be power of 2.
// The transform is unnormalised (entries +-1), so H H' = n I and H = H'.
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

static int next_pow2(int T) {
  int T2 = 1;
  while (T2 < T) T2 <<= 1;
  return T2;
}

static void srht_check_plan(int T, const uvec& rows, const vec& signs,
                            const uvec& perm, const char* who) {
  if ((int)signs.n_elem != T) {
    Rcpp::stop("%s: length(signs) must equal the series length.", who);
  }
  if ((int)perm.n_elem != T) {
    Rcpp::stop("%s: length(perm) must equal the series length.", who);
  }
  if (perm.n_elem > 0 && perm.max() >= (uword)T) {
    Rcpp::stop("%s: perm contains out-of-bounds indices.", who);
  }
  if (rows.n_elem > 0 && rows.max() >= (uword)T) {
    Rcpp::stop("%s: rows contains out-of-bounds indices.", who);
  }
}

// SRHT apply: S M for M (T x k), with S = scale * R P H D (m x T).
//   D: random signs (T), H: unnormalised Walsh-Hadamard transform on the
//   zero-padded power-of-two length, P: permutation of the first T outputs,
//   R: selection of m of them.
// Every row of H D has squared norm T over the unpadded coordinates and
// E_D |(H D r)_i|^2 = ||r||^2, so scale = 1 / sqrt(m) gives the normalised
// sketch E ||S r||^2 = ||r||^2 (make_srht_plan() uses that scale). When T is
// a power of two, S S' = (T / m) I_m exactly.
// [[Rcpp::export]]
arma::mat cpp_srht_apply(const arma::mat& M,
                         const arma::uvec& rows,
                         const arma::vec& signs,
                         const arma::uvec& perm,
                         const double scale) {
  int T = M.n_rows, K = M.n_cols;
  srht_check_plan(T, rows, signs, perm, "cpp_srht_apply");
  mat X = M.each_col() % signs;
  int T2 = next_pow2(T);
  mat Xpad(T2, K, fill::zeros);
  Xpad.rows(0, T-1) = X;
  fwht_cols(Xpad);
  mat out(rows.n_elem, K);
  for (uword i = 0; i < rows.n_elem; ++i) {
    out.row(i) = Xpad.row(perm(rows(i))) * scale;
  }
  return out;
}

// Adjoint of cpp_srht_apply: S' B for B (m x k), returning T x k. With
// cpp_srht_apply it gives S S' = S (S' I_m), the Gram matrix of the sketch
// rows needed for the conditional variance of sketch-and-solve estimates.
// [[Rcpp::export]]
arma::mat cpp_srht_adjoint(const arma::mat& B,
                           const arma::uvec& rows,
                           const arma::vec& signs,
                           const arma::uvec& perm,
                           const double scale) {
  int T = signs.n_elem, K = B.n_cols;
  srht_check_plan(T, rows, signs, perm, "cpp_srht_adjoint");
  if (B.n_rows != rows.n_elem) {
    Rcpp::stop("cpp_srht_adjoint: nrow(B) must equal length(rows).");
  }
  int T2 = next_pow2(T);
  mat Xpad(T2, K, fill::zeros);
  for (uword i = 0; i < rows.n_elem; ++i) {
    Xpad.row(perm(rows(i))) += B.row(i) * scale;
  }
  fwht_cols(Xpad);  // H is symmetric
  mat out = Xpad.rows(0, T - 1);
  out.each_col() %= signs;
  return out;
}

// Draw a normalised SRHT sketch of X (m rows).
static mat srht_sketch_random(const mat& X, int m, const mat* Z, mat* Zs) {
  int T = X.n_rows;
  arma::vec signs = 2.0 * randu<vec>(T) - 1.0;
  signs.transform( [](double v){ return v>=0 ? 1.0 : -1.0; } );
  arma::uvec perm = randperm(T);
  arma::uvec order = sort_index(randu<vec>(T));
  arma::uvec rows = order.subvec(0, m - 1);
  const double scale = 1.0 / std::sqrt((double)m);
  if (Z != nullptr && Zs != nullptr) {
    *Zs = cpp_srht_apply(*Z, rows, signs, perm, scale);
  }
  return cpp_srht_apply(X, rows, signs, perm, scale);
}

// Iterative Hessian sketch (Pilanci & Wainwright, 2016) for the multi-RHS
// least-squares problem min ||Z - X M||_F^2.
//
// Only the Hessian is sketched: each iteration draws a fresh normalised SRHT
// S (E S'S = I), so G = X'S'SX estimates X'X, and takes the step
// G^{-1} X'(Z - XM) with the exact full-data gradient. The iterates therefore
// converge to the least-squares solution itself, not to the solution of one
// sketched problem. The step is halved until the full-data residual sum of
// squares does not increase, so iterations are monotone.
//
// Iteration stops after `iters` iterations or, when tol > 0, as soon as every
// coefficient's step is below tol times its OLS standard error
// sqrt([X'X]^{-1}_ii * rss_j / (T - p)) (plus a relative floor of sqrt(eps)
// |M_ij| for columns fitted exactly).
//
// The exact gradient already touches the full data, so the exact (X'X)^{-1}
// (p x p) and the exact residuals Z - XM are returned at no extra cost in
// order: callers report OLS-exact covariance and residual variance with
// T - p degrees of freedom instead of sketch-based approximations.
// [[Rcpp::export]]
Rcpp::List cpp_ihs_latent(const arma::mat& X, const arma::mat& Z,
                          const int m, const int iters,
                          const double tol = 0.0) {
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
  if (!std::isfinite(tol) || tol < 0) {
    Rcpp::stop("cpp_ihs_latent: tol must be a finite non-negative number.");
  }
  if (Z.n_rows != X.n_rows) {
    Rcpp::stop("cpp_ihs_latent: X and Z must have the same number of rows.");
  }
  mat XtXinv;
  if (!inv_sympd_safe(XtXinv, X.t() * X)) {
    Rcpp::stop("cpp_ihs_latent: unable to invert X'X.");
  }
  // Warm start from the classical sketch-and-solve solution. IHS contracts
  // the error relative to its starting point, and fMRI data carry a large
  // baseline that a zero start would leave almost entirely unfitted;
  // sketch-and-solve fits any part of Z in span(X) exactly, so iterations
  // only refine the noise-driven error.
  mat Zs;
  mat Xs = srht_sketch_random(X, m, &Z, &Zs);
  mat Ginv;
  if (!inv_sympd_safe(Ginv, Xs.t() * Xs)) {
    Rcpp::stop("cpp_ihs_latent: unable to invert sketched Gram matrix.");
  }
  mat M = Ginv * (Xs.t() * Zs);
  mat E = Z - X * M;
  const double df = std::max(1, T - p);
  const double rel_floor = std::sqrt(arma::datum::eps);
  vec xdiag = sqrt(clamp(XtXinv.diag(), 0.0, arma::datum::inf));
  int used = 0;
  bool converged = false;
  for (int t = 0; t < iters; ++t) {
    mat Xh = srht_sketch_random(X, m, nullptr, nullptr);
    if (!inv_sympd_safe(Ginv, Xh.t() * Xh)) {
      Rcpp::stop("cpp_ihs_latent: unable to invert sketched Gram matrix.");
    }
    mat dM = Ginv * (X.t() * E);
    const double rss0 = accu(square(E));
    double step = 1.0;
    mat Enew = Z - X * (M + dM);
    for (int k = 0; k < 30; ++k) {
      if (accu(square(Enew)) <= rss0) break;
      step *= 0.5;
      Enew = Z - X * (M + step * dM);
    }
    M += step * dM;
    ++used;
    rowvec rss_cols = sum(square(E), 0);
    E = Enew;
    if (tol > 0) {
      // se_ij = sqrt([X'X]^{-1}_ii) * sqrt(rss_j / df)
      mat se = xdiag * sqrt(rss_cols / df);
      mat bound = tol * se + rel_floor * abs(M);
      if (all(vectorise(abs(step * dM) <= bound))) {
        converged = true;
        break;
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("M")         = M,
    Rcpp::Named("Ginv")      = XtXinv,
    Rcpp::Named("residuals") = E,
    Rcpp::Named("iters")     = used,
    Rcpp::Named("converged") = converged
  );
}
