// src/meta_kernels.cpp
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

inline bool solve_sympd_safe(mat& out, const mat& A, const vec& b) {
  // solve A * x = b, with A symmetric PD if possible; fall back to pinv
  bool ok = false;
  try {
    out = solve(A, b, solve_opts::likely_sympd);
    ok = true;
  } catch (...) {
    try {
      out = pinv(A) * b;
      ok = true;
    } catch (...) {
      ok = false;
    }
  }
  return ok;
}

inline bool inv_sympd_safe(mat& out, const mat& A) {
  bool ok = false;
  try {
    out = inv_sympd(A);
    ok = true;
  } catch (...) {
    try {
      out = pinv(A);
      ok = true;
    } catch (...) {
      ok = false;
    }
  }
  return ok;
}

inline double huber_weight(double rstd, double c) {
  double a = std::fabs(rstd);
  if (a <= c) return 1.0;
  return (c / a);
}

// Compute Q(tau2) for meta-regression: Q = r' W r, with W = diag(1/(v + tau2))
// Also returns beta and (optionally) XtWX^{-1} for reuse.
// X (m×k), y (m), v (m)
// returns Q, fills beta_hat, invXtWX (if non-null)
inline double Q_tau2(const mat& X, const vec& y, const vec& v, double tau2,
                     vec& beta_hat, mat* invXtWX_ptr = nullptr) {
  vec w = 1.0 / (v + tau2);
  // Weighted LS
  mat Xw = X; Xw.each_col() %= sqrt(w);
  mat XtWX = Xw.t() * Xw;
  vec Xty = X.t() * (w % y);

  mat invXtWX;
  if (!inv_sympd_safe(invXtWX, XtWX)) {
    beta_hat.set_size(X.n_cols); beta_hat.fill(datum::nan);
    if (invXtWX_ptr) invXtWX_ptr->reset();
    return datum::nan;
  }
  beta_hat = invXtWX * Xty;

  vec r = y - X * beta_hat;
  double Q = as_scalar( (w % r) .t() * r ); // sum_i w_i r_i^2
  if (invXtWX_ptr) *invXtWX_ptr = invXtWX;
  return Q;
}

// DL estimator for tau^2, generalized to meta-regression.
// Uses weights w0 = 1/v, FE fit (tau2 = 0), Q0 = sum w0 r^2, df = m - p,
// denominator C = sum(w0) - sum(w0^2 * t_i), where t_i = x_i' (X'W0X)^{-1} x_i
inline double tau2_DL(const mat& X, const vec& y, const vec& v) {
  const uword m = y.n_elem;
  const uword p = X.n_cols;
  vec w0 = 1.0 / v;

  // FE (tau2 = 0)
  mat Xw = X; Xw.each_col() %= sqrt(w0);
  mat XtWX = Xw.t() * Xw;
  vec Xty = X.t() * (w0 % y);

  mat invXtWX;
  if (!inv_sympd_safe(invXtWX, XtWX)) return 0.0;

  vec beta_fe = invXtWX * Xty;
  vec r = y - X * beta_fe;

  double Q0 = as_scalar( (w0 % r).t() * r );
  double df = (double)m - (double)p;

  // Compute t_i = x_i' invXtWX x_i; then denom C = sum(w0) - sum(w0^2 * t_i)
  vec ti(m, fill::zeros);
  for (uword i = 0; i < m; ++i) {
    rowvec xi = X.row(i);
    double ti_i = as_scalar( xi * invXtWX * xi.t() );
    ti(i) = ti_i;
  }
  double sumw = sum(w0);
  double denom = sumw - sum( (w0 % w0) % ti );

  double num = Q0 - df;
  double tau2 = (denom > 0.0) ? std::max(0.0, num / denom) : 0.0;
  if (!arma::is_finite(tau2)) tau2 = 0.0;
  return tau2;
}

// Paule–Mandel via simple bisection on f(tau2) = Q(tau2) - df.
// Monotone decreasing in tau2, bracket [0, tau_hi] until sign change.
inline double tau2_PM(const mat& X, const vec& y, const vec& v,
                      double tol = 1e-6, int max_iter = 50,
                      bool* warned_nonconverged = nullptr) {
  const uword m = y.n_elem;
  const uword p = X.n_cols;
  const double df = (double)m - (double)p;

  vec beta_hat;
  double f0 = Q_tau2(X, y, v, 0.0, beta_hat) - df;
  if (!arma::is_finite(f0) || f0 <= 0.0) return 0.0;

  // Find upper bound
  double tau_lo = 0.0, tau_hi = std::max(1.0, 1000.0 * median(v));
  for (int k = 0; k < 20; ++k) {
    double fhi = Q_tau2(X, y, v, tau_hi, beta_hat) - df;
    if (!arma::is_finite(fhi)) tau_hi *= 2.0; // try larger
    else if (fhi <= 0.0) break;
    else tau_hi *= 2.0;
    if (tau_hi > 1e12) break;
  }

  // Bisection
  double a = tau_lo, b = tau_hi, fb = Q_tau2(X, y, v, b, beta_hat) - df;
  if (!arma::is_finite(fb)) fb = -1.0; // force progress
  double mid = a;
  bool converged = false;
  for (int it = 0; it < max_iter; ++it) {
    mid = 0.5 * (a + b);
    double fm = Q_tau2(X, y, v, mid, beta_hat) - df;
    if (!arma::is_finite(fm)) { a = mid; continue; }
    if (std::fabs(fm) < tol) { converged = true; break; }
    if (fm > 0.0) { a = mid; } else { b = mid; fb = fm; }
  }
  if (!arma::is_finite(mid)) mid = 0.0;
  if (mid < 0.0) mid = 0.0;
  if (!converged && warned_nonconverged != nullptr) {
    *warned_nonconverged = true;
  }
  return mid;
}

// One feature (column) meta-regression with method and optional robust Huber IRLS.
// Returns: beta (k), se (k), z (k), tau2, Q_fe, I2_fe, df
struct MetaResult {
  vec beta, se, z;
  double tau2;
  double Q_fe;
  double I2_fe;
  double df;
  bool ok;
  bool pm_nonconverged;
  mat invXtWX;
  MetaResult() : tau2(0.0), Q_fe(NA_REAL), I2_fe(NA_REAL), df(NA_REAL), ok(false), pm_nonconverged(false) {}
};

enum Method { FE=0, DL=1, PM=2, REML_ALIAS=3 };

MetaResult fit_one(const mat& Xfull, const vec& yfull, const vec& vfull,
                   int method_code, bool robust, double huber_c, int robust_iter) {

  MetaResult out;
  // Filter finite rows - find indices where all conditions are true
  uvec finite_y = find_finite(yfull);
  uvec finite_v = find_finite(vfull);
  uvec positive_v = find(vfull > 0);
  
  // Find intersection of all conditions
  uvec keep = intersect(intersect(finite_y, finite_v), positive_v);
  if (keep.n_elem < Xfull.n_cols) return out; // not enough data

  mat X = Xfull.rows(keep);
  vec y = yfull.elem(keep);
  vec v = vfull.elem(keep);

  const uword m = y.n_elem;
  const uword p = X.n_cols;
  if (m <= p) return out;

  // FE Q and I2 (using w0 = 1/v)
  vec w0 = 1.0 / v;
  mat Xw0 = X; Xw0.each_col() %= sqrt(w0);
  mat XtWX0 = Xw0.t() * Xw0;
  vec Xty0 = X.t() * (w0 % y);

  mat invXtWX0;
  if (!inv_sympd_safe(invXtWX0, XtWX0)) return out;
  vec beta_fe = invXtWX0 * Xty0;
  vec r0 = y - X * beta_fe;
  double Q0 = as_scalar( (w0 % r0).t() * r0 );
  double df = (double)m - (double)p;
  double I2 = NA_REAL;
  if (Q0 > 0.0) {
    double i2 = (Q0 - df) / Q0;
    I2 = (i2 < 0.0) ? 0.0 : i2;
  }
  out.Q_fe = Q0;
  out.I2_fe = I2;
  out.df = df;

  // tau2 by chosen method
  double tau2 = 0.0;
  Method M = (method_code==0 ? FE : method_code==1 ? DL : method_code==2 ? PM : REML_ALIAS);
  if (M == FE) {
    tau2 = 0.0;
  } else if (M == DL) {
    tau2 = tau2_DL(X, y, v);
  } else { // PM or REML_ALIAS -> use PM solver (robust & fast)
    bool pm_warn = false;
    tau2 = tau2_PM(X, y, v, 1e-6, 50, &pm_warn);
    out.pm_nonconverged = pm_warn;
  }
  if (!arma::is_finite(tau2) || tau2 < 0.0) tau2 = 0.0;

  // Base weights
  vec w = 1.0 / (v + tau2);
  vec aw = ones<vec>(m);

  // Optional robust Huber IRLS (few sweeps)
  vec beta_hat(p, fill::zeros);
  mat invXtWX;
  for (int it = 0; it < (robust ? robust_iter : 1); ++it) {
    vec weff = w % aw;
    mat XwE = X; XwE.each_col() %= sqrt(weff);
    mat XtWX = XwE.t() * XwE;
    vec Xty  = X.t() * (weff % y);
    if (!inv_sympd_safe(invXtWX, XtWX)) { out.ok = false; return out; }
    beta_hat = invXtWX * Xty;
    if (!robust) break;
    // Update robust weights
    vec r = y - X * beta_hat;
    for (uword i = 0; i < m; ++i) {
      double rstd = r(i) / std::sqrt(v(i) + tau2);
      aw(i) = huber_weight(rstd, huber_c);
    }
  }

  // Final variance and z
  vec se = sqrt( invXtWX.diag() );
  vec z = beta_hat / se;

  out.beta = beta_hat;
  out.se   = se;
  out.z    = z;
  out.tau2 = tau2;
  out.ok   = true;
  out.invXtWX = invXtWX;
  return out;
}


// [[Rcpp::export]]
Rcpp::List meta_fit_cpp(const arma::mat& Y,      // subjects x features
                        const arma::mat& V,      // subjects x features (variances)
                        const arma::mat& X,      // subjects x predictors (incl. intercept)
                        const std::string method, // "fe","dl","pm","reml"
                        const std::string robust, // "none","huber"
                        const double huber_c = 1.345,
                        const int robust_iter = 2,
                        const int n_threads = 0) {

  const uword S = Y.n_rows;
  const uword P = Y.n_cols;
  const uword Sp = V.n_rows;
  const uword Pp = V.n_cols;
  if (S != Sp || P != Pp) stop("Y and V must have the same shape.");
  if (X.n_rows != S)      stop("X must have the same number of rows as Y.");

  int mcode = 0;
  if (method == "fe") mcode = 0; else if (method=="dl") mcode = 1;
  else if (method=="pm") mcode = 2; else if (method=="reml") mcode = 3;
  else stop("Unknown method: ", method);

  bool do_robust = (robust == "huber");

  arma::mat B(X.n_cols, P); B.fill(datum::nan);
  arma::mat SE(X.n_cols, P); SE.fill(datum::nan);
  arma::mat Z(X.n_cols, P); Z.fill(datum::nan);
  arma::rowvec TAU2(P); TAU2.fill(datum::nan);
  arma::rowvec Qfe(P); Qfe.fill(datum::nan);
  arma::rowvec I2(P); I2.fill(datum::nan);
  arma::rowvec DF(P); DF.fill(datum::nan);
  arma::uvec   OK(P,   fill::zeros);
  arma::uvec PM_WARN(P, fill::zeros);

  // Parallel over features
#ifdef _OPENMP
  int nthreads = n_threads;
  if (nthreads <= 0) {
    nthreads = omp_get_max_threads();
  }
#pragma omp parallel for num_threads(nthreads) schedule(static)
#else
  (void)n_threads; // Suppress unused parameter warning
#endif
  for (int j = 0; j < (int)P; ++j) {
    vec y = Y.col(j);
    vec v = V.col(j);

    MetaResult res = fit_one(X, y, v, mcode, do_robust, huber_c, robust_iter);
    if (!res.ok) continue;

    B.col(j)   = res.beta;
    SE.col(j)  = res.se;
    Z.col(j)   = res.z;
    TAU2(j)    = res.tau2;
    Qfe(j)     = res.Q_fe;
    I2(j)      = res.I2_fe;
    DF(j)      = res.df;
    OK(j)      = 1u;
    PM_WARN(j) = res.pm_nonconverged ? 1u : 0u;
  }

  if (arma::accu(PM_WARN) > 0) {
    Rcpp::warning("Paule-Mandel solver did not fully converge for some features; using last iterate");
  }

  return List::create(
    _["beta"] = B,
    _["se"]   = SE,
    _["z"]    = Z,
    _["tau2"] = TAU2,
    _["Q_fe"] = Qfe,
    _["I2_fe"]= I2,
    _["df"]   = DF,
    _["ok"]   = OK
  );
}

// [[Rcpp::export]]
Rcpp::List meta_fit_contrasts_cpp(const arma::mat& Y,      // S x P
                                  const arma::mat& V,      // S x P (variances)
                                  const arma::mat& X,      // S x K
                                  const arma::mat& Cmat,   // K x J (rows are predictors, cols are contrasts)
                                  const std::string method,
                                  const std::string robust,
                                  const double huber_c = 1.345,
                                  const int robust_iter = 2,
                                  const int n_threads = 0) {
  using namespace arma;

  const uword S = Y.n_rows, P = Y.n_cols;
  if (V.n_rows != S || V.n_cols != P) stop("meta_fit_contrasts_cpp: Y and V must match.");
  if (X.n_rows != S) stop("meta_fit_contrasts_cpp: X rows must equal nrow(Y).");
  if (Cmat.n_rows != X.n_cols) stop("meta_fit_contrasts_cpp: nrow(C) must equal ncol(X).");

  int mcode = 0;
  if (method == "fe") mcode = 0; else if (method=="dl") mcode = 1; else if (method=="pm") mcode = 2; else if (method=="reml") mcode = 3; else stop("Unknown method");
  bool do_robust = (robust == "huber");

  const uword K = X.n_cols;
  const uword J = Cmat.n_cols;

  mat B(K, P); B.fill(datum::nan);
  mat SE(K, P); SE.fill(datum::nan);
  mat Z(K, P); Z.fill(datum::nan);
  rowvec TAU2(P); TAU2.fill(datum::nan);
  rowvec Qfe(P); Qfe.fill(datum::nan);
  rowvec I2(P); I2.fill(datum::nan);
  rowvec DF(P); DF.fill(datum::nan);
  uvec   OK(P,   fill::zeros);
  uvec PM_WARN(P, fill::zeros);

  mat CB(J, P); CB.fill(datum::nan);
  mat CSE(J, P); CSE.fill(datum::nan);
  mat CZ(J, P); CZ.fill(datum::nan);

#ifdef _OPENMP
  int nthreads = n_threads;
  if (nthreads <= 0) nthreads = omp_get_max_threads();
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
  for (int j = 0; j < (int)P; ++j) {
    vec y = Y.col(j);
    vec v = V.col(j);
    MetaResult res = fit_one(X, y, v, mcode, do_robust, huber_c, robust_iter);
    if (!res.ok) continue;

    B.col(j)   = res.beta;
    SE.col(j)  = res.se;
    Z.col(j)   = res.z;
    TAU2(j)    = res.tau2;
    Qfe(j)     = res.Q_fe;
    I2(j)      = res.I2_fe;
    DF(j)      = res.df;
    OK(j)      = 1u;
    PM_WARN(j) = res.pm_nonconverged ? 1u : 0u;

    // Compute contrasts
    for (uword c = 0; c < J; ++c) {
      vec cw = Cmat.col(c);
      double est = dot(cw, res.beta);
      double varc = as_scalar(cw.t() * res.invXtWX * cw);
      double sec = (varc > 0.0 && arma::is_finite(varc)) ? std::sqrt(varc) : datum::nan;
      double zc = (std::isfinite(sec) && sec > 0.0) ? (est / sec) : datum::nan;
      CB(c, j)  = est;
      CSE(c, j) = sec;
      CZ(c, j)  = zc;
    }
  }

  if (arma::accu(PM_WARN) > 0) {
    Rcpp::warning("Paule-Mandel solver did not fully converge for some features; using last iterate");
  }

  return List::create(
    _["beta"] = B,
    _["se"]   = SE,
    _["z"]    = Z,
    _["tau2"] = TAU2,
    _["Q_fe"] = Qfe,
    _["I2_fe"]= I2,
    _["df"]   = DF,
    _["ok"]   = OK,
    _["c_beta"] = CB,
    _["c_se"] = CSE,
    _["c_z"] = CZ
  );
}

// [[Rcpp::export]]
Rcpp::List meta_fit_cov_cpp(const arma::mat& Y,      // S x P betas
                            const arma::mat& V,      // S x P variances
                            const arma::mat& X,      // S x K
                            const std::string method,
                            const std::string robust,
                            const double huber_c = 1.345,
                            const int robust_iter = 2,
                            const int n_threads = 0) {
  using namespace arma;

  const uword S = Y.n_rows, P = Y.n_cols;
  if (V.n_rows != S || V.n_cols != P) Rcpp::stop("meta_fit_cov_cpp: Y/V shape mismatch.");
  if (X.n_rows != S) Rcpp::stop("meta_fit_cov_cpp: X rows must equal nrow(Y).");

  int mcode = 0;
  if (method == "fe") mcode = 0;
  else if (method == "dl") mcode = 1;
  else if (method == "pm") mcode = 2;
  else if (method == "reml") mcode = 3;
  else Rcpp::stop("Unknown method: %s", method.c_str());

  bool do_robust = (robust == "huber");

  const uword K = X.n_cols;
  const uword TSZ = K * (K + 1) / 2;
  mat B(K, P); B.fill(datum::nan);
  mat SE(K, P); SE.fill(datum::nan);
  mat Z(K, P); Z.fill(datum::nan);
  rowvec TAU2(P); TAU2.fill(datum::nan);
  rowvec Qfe(P); Qfe.fill(datum::nan);
  rowvec I2(P); I2.fill(datum::nan);
  rowvec DF(P); DF.fill(datum::nan);
  uvec OK(P, fill::zeros);
  uvec PM_WARN(P, fill::zeros);
  mat COVTRI(TSZ, P); COVTRI.fill(datum::nan);

#ifdef _OPENMP
  int nthreads = n_threads;
  if (nthreads <= 0) nthreads = omp_get_max_threads();
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
  for (int j = 0; j < (int)P; ++j) {
    vec y = Y.col(j);
    vec v = V.col(j);
    MetaResult res = fit_one(X, y, v, mcode, do_robust, huber_c, robust_iter);
    if (!res.ok) continue;
    B.col(j)   = res.beta;
    SE.col(j)  = res.se;
    Z.col(j)   = res.z;
    TAU2(j)    = res.tau2;
    Qfe(j)     = res.Q_fe;
    I2(j)      = res.I2_fe;
    DF(j)      = res.df;
    OK(j)      = 1u;
    PM_WARN(j) = res.pm_nonconverged ? 1u : 0u;

    // pack upper triangle
    uword idx = 0;
    for (uword a = 0; a < K; ++a) {
      for (uword b = a; b < K; ++b) {
        COVTRI(idx++, j) = res.invXtWX(a, b);
      }
    }
  }

  if (arma::accu(PM_WARN) > 0) {
    Rcpp::warning("Paule-Mandel solver did not fully converge for some features; using last iterate");
  }

  return Rcpp::List::create(
    Rcpp::Named("beta") = B,
    Rcpp::Named("se") = SE,
    Rcpp::Named("z") = Z,
    Rcpp::Named("tau2") = TAU2,
    Rcpp::Named("Q_fe") = Qfe,
    Rcpp::Named("I2_fe") = I2,
    Rcpp::Named("df") = DF,
    Rcpp::Named("ok") = OK,
    Rcpp::Named("cov_tri") = COVTRI
  );
}

//' Meta-regression with ONE voxelwise covariate
//' 
//' @param Y S x P matrix of effect sizes
//' @param V S x P matrix of variances
//' @param X S x K design matrix
//' @param C S x P matrix of voxelwise covariates
//' @param method Meta-analysis method
//' @param robust Robust estimation method
//' @param huber_c Huber tuning constant
//' @param robust_iter Number of IRLS iterations
//' @param n_threads Number of OpenMP threads
//' @return List with results
//' @export
// [[Rcpp::export]]
Rcpp::List meta_fit_vcov_cpp(const arma::mat& Y,  // S x P betas
                             const arma::mat& V,  // S x P variances
                             const arma::mat& X,  // S x K
                             const arma::mat& C,  // S x P (voxelwise covariate)
                             const std::string method,
                             const std::string robust,
                             const double huber_c = 1.345,
                             const int robust_iter = 2,
                             const int n_threads = 0) {
  using namespace arma;
  
  const uword S = Y.n_rows, P = Y.n_cols;
  
  if (V.n_rows != S || V.n_cols != P) {
    Rcpp::stop("meta_fit_vcov_cpp: Y/V shape mismatch.");
  }
  if (X.n_rows != S) {
    Rcpp::stop("meta_fit_vcov_cpp: X rows must equal nrow(Y).");
  }
  if (C.n_rows != S || C.n_cols != P) {
    Rcpp::stop("meta_fit_vcov_cpp: C must be S x P.");
  }
  
  int mcode = 0;
  if (method == "fe") mcode = 0;
  else if (method == "dl") mcode = 1;
  else if (method == "pm") mcode = 2;
  else if (method == "reml") mcode = 3;
  else Rcpp::stop("Unknown method: %s", method.c_str());
  
  bool do_robust = (robust == "huber");
  
  const uword K = X.n_cols + 1;  // augmented dimension
  mat B(K, P); B.fill(datum::nan);
  mat SE(K, P); SE.fill(datum::nan);
  mat Z(K, P); Z.fill(datum::nan);
  rowvec TAU2(P); TAU2.fill(datum::nan);
  rowvec Qfe(P); Qfe.fill(datum::nan);
  rowvec I2(P); I2.fill(datum::nan);
  rowvec DF(P); DF.fill(datum::nan);
  uvec OK(P, fill::zeros);
  uvec PM_WARN(P, fill::zeros);
  
  int nth = n_threads;
#ifdef _OPENMP
  if (nth <= 0) nth = omp_get_max_threads();
#pragma omp parallel for num_threads(nth) schedule(static)
#endif
  for (int j = 0; j < (int)P; ++j) {
    // Build augmented design matrix for feature j
    mat Xj(X.n_rows, X.n_cols + 1);
    Xj.cols(0, X.n_cols - 1) = X;
    Xj.col(X.n_cols) = C.col(j);
    
    vec y = Y.col(j);
    vec v = V.col(j);
    
    // Use existing fit_one function
    MetaResult res = fit_one(Xj, y, v, mcode, do_robust, huber_c, robust_iter);
    
    if (!res.ok) continue;
    
    B.col(j) = res.beta;
    SE.col(j) = res.se;
    Z.col(j) = res.z;
    TAU2(j) = res.tau2;
    Qfe(j) = res.Q_fe;
    I2(j) = res.I2_fe;
    DF(j) = res.df;
    OK(j) = 1u;
    PM_WARN(j) = res.pm_nonconverged ? 1u : 0u;
  }

  if (arma::accu(PM_WARN) > 0) {
    Rcpp::warning("Paule-Mandel solver did not fully converge for some features; using last iterate");
  }
  
  return Rcpp::List::create(
    Rcpp::Named("beta") = B,
    Rcpp::Named("se") = SE,
    Rcpp::Named("z") = Z,
    Rcpp::Named("tau2") = TAU2,
    Rcpp::Named("Q_fe") = Qfe,
    Rcpp::Named("I2_fe") = I2,
    Rcpp::Named("df") = DF,
    Rcpp::Named("ok") = OK,
    Rcpp::Named("method") = method,
    Rcpp::Named("robust") = robust
  );
}
