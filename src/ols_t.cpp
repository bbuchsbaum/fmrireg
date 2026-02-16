// OLS t-test and ANCOVA kernels for fmri_ttest
// Fast Student t-tests and Welch tests across features

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Safe inverse of XtX (symmetric); pinv fallback
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

//' OLS t-test / ANCOVA across features
//' 
//' @param Y S x P matrix (subjects x features)
//' @param X S x K design matrix with intercept if desired
//' @return List with beta (K x P), se (K x P), t (K x P), df (scalar), ok (P)
//' @export
// [[Rcpp::export]]
List ols_t_cpp(const arma::mat& Y, const arma::mat& X) {
  const uword S = Y.n_rows, P = Y.n_cols, K = X.n_cols;
  
  if (X.n_rows != S) {
    stop("ols_t_cpp: nrow(X) must equal nrow(Y).");
  }
  
  mat XtX = X.t() * X;
  mat invXtX;
  
  if (!inv_sympd_safe(invXtX, XtX)) {
    stop("ols_t_cpp: design matrix not invertible.");
  }
  
  // Compute coefficients
  mat B = invXtX * (X.t() * Y);  // K x P
  
  // Compute residuals
  mat R = Y - X * B;  // S x P
  
  // Residual sum of squares
  rowvec RSS = sum(square(R), 0);  // 1 x P
  
  // Degrees of freedom
  double df = (double)S - (double)K;
  if (df <= 0) {
    stop("ols_t_cpp: non-positive residual degrees of freedom.");
  }
  
  // Residual variance
  rowvec s2 = RSS / df;  // 1 x P
  
  // Standard errors: SE_kj = sqrt(s2_j * invXtX_kk)
  vec diagInv = invXtX.diag();  // K
  mat SE = sqrt(diagInv * s2);  // K x P (outer product)
  
  // T-statistics
  mat T = B / SE;  // K x P
  
  // Identify valid columns (finite t-values)
  uvec ok = conv_to<uvec>::from(find_finite(T.row(0).t()));
  
  return List::create(
    Named("beta") = B,
    Named("se") = SE,
    Named("t") = T,
    Named("df") = df,
    Named("ok") = ok
  );
}

//' Welch two-sample t-test across features
//' 
//' @param Y S x P matrix
//' @param g_in Length S vector of group indicators (1/2 or 0/1)
//' @return List with muA, muB, t, df (Welch), nA, nB
//' @export
// [[Rcpp::export]]
List welch_t_cpp(const arma::mat& Y, const IntegerVector& g_in) {
  const uword S = Y.n_rows, P = Y.n_cols;
  
  if ((uword)g_in.size() != S) {
    stop("welch_t_cpp: group length must match nrow(Y).");
  }
  
  // Map to {0, 1}
  std::vector<int> gS(S);
  for (uword i = 0; i < S; ++i) {
    gS[i] = (g_in[i] == 2 || g_in[i] == 1) ? (g_in[i] - 1) : g_in[i];
  }
  
  // Separate indices for each group
  std::vector<uword> ia, ib;
  ia.reserve(S); 
  ib.reserve(S);
  
  for (uword i = 0; i < S; ++i) {
    if (gS[i] == 0) {
      ia.push_back(i);
    } else {
      ib.push_back(i);
    }
  }
  
  const uword nA = ia.size(), nB = ib.size();
  
  if (nA < 2 || nB < 2) {
    stop("welch_t_cpp: need at least 2 subjects per group.");
  }
  
  // Compute means
  rowvec muA(P, fill::zeros), muB(P, fill::zeros);
  
  for (uword j = 0; j < P; ++j) {
    double sA = 0, sB = 0;
    for (uword a = 0; a < nA; ++a) sA += Y(ia[a], j);
    for (uword b = 0; b < nB; ++b) sB += Y(ib[b], j);
    muA(j) = sA / (double)nA;
    muB(j) = sB / (double)nB;
  }
  
  // Compute variances
  rowvec vA(P, fill::zeros), vB(P, fill::zeros);
  
  for (uword j = 0; j < P; ++j) {
    double sA = 0, sB = 0;
    for (uword a = 0; a < nA; ++a) {
      double d = Y(ia[a], j) - muA(j);
      sA += d * d;
    }
    for (uword b = 0; b < nB; ++b) {
      double d = Y(ib[b], j) - muB(j);
      sB += d * d;
    }
    vA(j) = sA / (double)(nA - 1);
    vB(j) = sB / (double)(nB - 1);
  }
  
  // Standard error and t-statistic
  rowvec se2 = vA / (double)nA + vB / (double)nB;  // 1 x P
  rowvec t = (muA - muB) / sqrt(se2 + 1e-300);
  
  // Welch-Satterthwaite degrees of freedom
  rowvec num = square(se2);
  rowvec den = (square(vA) / ((double)nA * (double)nA * (double)(nA - 1))) +
               (square(vB) / ((double)nB * (double)nB * (double)(nB - 1)));
  rowvec df = num / (den + 1e-300);
  
  return List::create(
    Named("muA") = muA,
    Named("muB") = muB,
    Named("t") = t,
    Named("df") = df,
    Named("nA") = (int)nA,
    Named("nB") = (int)nB
  );
}

//' OLS with ONE voxelwise covariate
//' 
//' @param Y S x P matrix of outcomes
//' @param X S x K design matrix
//' @param C S x P matrix of voxelwise covariates
//' @return List with beta ((K+1) x P), se, t, df
//' @export
// [[Rcpp::export]]
List ols_t_vcov_cpp(const arma::mat& Y, const arma::mat& X, const arma::mat& C) {
  const uword S = Y.n_rows, P = Y.n_cols, K = X.n_cols;
  
  if (X.n_rows != S || C.n_rows != S || C.n_cols != P) {
    stop("ols_t_vcov_cpp: dimension mismatch.");
  }
  
  if (S <= (K + 1)) {
    stop("ols_t_vcov_cpp: not enough rows (subjects) for K+1 parameters.");
  }
  
  // Precompute X'X and its inverse
  mat Xt = X.t();
  mat XtX = Xt * X;
  mat Ainv;
  
  if (!inv_sympd_safe(Ainv, XtX)) {
    stop("ols_t_vcov_cpp: X'X not invertible.");
  }
  
  double df = (double)S - (double)(K + 1);
  mat B(K + 1, P, fill::value(datum::nan));
  mat SE(K + 1, P, fill::value(datum::nan));
  mat T(K + 1, P, fill::value(datum::nan));
  
  // Process each feature
  for (uword j = 0; j < P; ++j) {
    vec c = C.col(j);
    vec y = Y.col(j);
    
    // Build augmented inverse via block inversion
    vec Xtc = Xt * c;  // K x 1
    double ssc = dot(c, c);  // scalar
    vec AinvB = Ainv * Xtc;  // K x 1
    double Ssch = ssc - as_scalar(Xtc.t() * AinvB);
    
    mat Minv;
    Minv.set_size(K + 1, K + 1);
    
    if (!std::isfinite(Ssch) || Ssch <= 1e-12) {
      // Fallback: full inverse of augmented matrix
      mat M(K + 1, K + 1, fill::zeros);
      M(span(0, K - 1), span(0, K - 1)) = XtX;
      M(span(0, K - 1), K) = Xtc;
      M(K, span(0, K - 1)) = Xtc.t();
      M(K, K) = ssc;
      
      if (!inv_sympd_safe(Minv, M)) {
        Minv = pinv(M);
      }
    } else {
      // Block inversion formula
      double Sinv = 1.0 / Ssch;
      mat TL = Ainv + (AinvB * AinvB.t()) * Sinv;  // K x K
      vec TR = -AinvB * Sinv;  // K x 1
      rowvec BL = -AinvB.t() * Sinv;  // 1 x K
      double BR = Sinv;  // 1 x 1
      
      Minv(span(0, K - 1), span(0, K - 1)) = TL;
      Minv(span(0, K - 1), K) = TR;
      Minv(K, span(0, K - 1)) = BL;
      Minv(K, K) = BR;
    }
    
    // Compute coefficients: RHS = [X'y; c'y]
    vec rhs(K + 1);
    rhs(span(0, K - 1)) = Xt * y;
    rhs(K) = dot(c, y);
    vec bj = Minv * rhs;
    
    // Compute residuals and variance
    vec r = y - X * bj.head(K) - c * bj(K);
    double s2 = dot(r, r) / df;
    
    // Standard errors and t-statistics
    vec sej = sqrt(s2 * Minv.diag());
    vec tj = bj / sej;
    
    B.col(j) = bj;
    SE.col(j) = sej;
    T.col(j) = tj;
  }
  
  return List::create(
    Named("beta") = B,
    Named("se") = SE,
    Named("t") = T,
    Named("df") = df
  );
}