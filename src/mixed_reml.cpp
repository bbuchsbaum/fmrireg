// src/mixed_reml.cpp
//
// Batched REML / BLUP kernel for estimate_betas(method = "mixed").
//
// fmrilss::mixed_solve() fits y = X b + Z u + e, u ~ N(0, Vu I), e ~ N(0, Ve I)
// for one response and redoes the spectral decomposition of the design (two
// n x n decompositions) on every call. In estimate_betas() every voxel shares
// the design, so here (EMMA / rrBLUP):
//
//   * mixed_reml_spectral_cpp() forms the REML rotation (theta, Q) once per
//     design, with the same arithmetic as fmrilss::mixed_solve (K = NULL),
//     including its choice between the SVD/QR ("cholesky") and eigen routes.
//   * mixed_reml_blup_cpp() then does, per voxel, omega = Q'y, the 1-D REML
//     search over lambda = Ve / Vu with fmrilss's objective, optimiser
//     (roptim L-BFGS-B, numerical gradient, maxit = 1000), start value and
//     bounds, and Vu / Ve exactly as fmrilss computes them. The REML surface is
//     often flat in lambda (weak trial signal), where L-BFGS-B stops wherever
//     rounding noise in the finite-difference gradient leads it; reproducing
//     fmrilss's arithmetic keeps those voxels on the same path.
//   * GLS fixed effects and BLUPs use the thin SVD Z = U D V'
//     (lambda H^-1 = I - U U' + U diag(h) U', h = lambda / (d^2 + lambda)),
//     so no n x n matrix is formed per voxel.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(roptim)]]

#include <RcppArmadillo.h>
#include <roptim.h>
#include <limits>
#include <string>

using namespace Rcpp;

namespace {

// Identical to REML_obj_fun in fmrilss/src/mixed_solve.cpp.
struct RemlObjective : public roptim::Functor {
  int n_p;
  arma::vec theta;
  arma::vec omega_sq;

  RemlObjective(int n_p_, const arma::vec& theta_)
    : n_p(n_p_), theta(theta_) {}

  double operator()(const arma::vec& lambda_vec) override {
    double lambda = lambda_vec[0];
    arma::vec theta_lambda = theta + lambda;

    if (arma::any(theta_lambda <= 0)) {
      return std::numeric_limits<double>::infinity();
    }

    double sum_omega_sq_div_theta_lambda = arma::sum(omega_sq / theta_lambda);
    double obj = n_p * std::log(sum_omega_sq_div_theta_lambda) +
      arma::sum(arma::log(theta_lambda));
    return obj;
  }
};

} // namespace

// REML rotation of a mixed-model design, as in fmrilss::mixed_solve (K = NULL):
// theta (n - p eigenvalues of S Z Z' S on col(S)) and Q (n x (n - p)).
// [[Rcpp::export]]
List mixed_reml_spectral_cpp(const arma::mat& X, const arma::mat& Z) {
  const int n_filtered = X.n_rows;
  const int p = X.n_cols;
  const int m = Z.n_cols;
  if (static_cast<int>(Z.n_rows) != n_filtered) stop("nrow(Z) != nrow(X)");

  arma::mat XtX = X.t() * X;
  int rank_X = arma::rank(XtX);
  if (rank_X < p) stop("X not full rank");

  arma::mat XtXinv = arma::inv_sympd(XtX);
  arma::mat S = arma::eye<arma::mat>(n_filtered, n_filtered) - X * XtXinv * X.t();

  std::string spectral_method = (n_filtered <= m + p) ? "eigen" : "cholesky";

  arma::vec theta;
  arma::mat Q;

  if (spectral_method == "cholesky") {
    const arma::mat& ZBt = Z;
    arma::mat SZBt = S * ZBt;
    arma::mat U_svd_SZBt;
    arma::vec s_svd_SZBt;
    arma::mat V_svd_SZBt;
    bool svd_success = arma::svd(U_svd_SZBt, s_svd_SZBt, V_svd_SZBt, SZBt);

    if (!svd_success) {
      SZBt += arma::randn<arma::mat>(SZBt.n_rows, SZBt.n_cols) * 1e-10;
      arma::svd(U_svd_SZBt, s_svd_SZBt, V_svd_SZBt, SZBt);
    }

    arma::mat QR_mat = arma::join_rows(X, U_svd_SZBt);
    arma::mat Q_full;
    arma::mat R_full;
    arma::qr(Q_full, R_full, QR_mat);

    Q = Q_full.cols(p, Q_full.n_cols - 1);

    int R_start = p;
    int R_end = std::min(p + m - 1, static_cast<int>(R_full.n_cols) - 1);
    arma::mat R = R_full.submat(R_start, R_start, R_end, R_end);

    arma::mat R_squared = arma::square(R);
    arma::mat R_squared_t = R_squared.t();

    arma::vec ans;
    bool solve_success = arma::solve(ans, R_squared_t, arma::square(s_svd_SZBt));

    if (!solve_success) {
      spectral_method = "eigen";
    } else {
      theta = arma::vec(n_filtered - p);
      if (ans.n_elem > 0) {
        arma::uword end_idx = std::min<arma::uword>(ans.n_elem - 1, n_filtered - p - 1);
        theta.subvec(0, end_idx) = ans;
      }
      arma::uword start_idx = std::min<arma::uword>(ans.n_elem, n_filtered - p);
      if (start_idx <= static_cast<arma::uword>(n_filtered - p - 1)) {
        theta.subvec(start_idx, n_filtered - p - 1).zeros();
      }
    }
  }

  if (spectral_method == "eigen") {
    double offset = std::sqrt(static_cast<double>(n_filtered));
    arma::mat Hb = Z * Z.t() + offset * arma::eye<arma::mat>(n_filtered, n_filtered);

    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, Hb);
    eigval = arma::reverse(eigval);
    arma::vec phi = eigval - offset;
    if (phi.min() < -1e-6) stop("K not positive semi-definite.");

    arma::mat SHbS = S * Hb * S;
    arma::vec eigval_SHbS;
    arma::mat eigvec_SHbS;
    arma::eig_sym(eigval_SHbS, eigvec_SHbS, SHbS);

    eigval_SHbS = arma::reverse(eigval_SHbS);
    eigvec_SHbS = arma::fliplr(eigvec_SHbS);

    theta = eigval_SHbS.subvec(0, n_filtered - p - 1) - offset;
    Q = eigvec_SHbS.cols(0, n_filtered - p - 1);
  }

  return List::create(
    Named("theta") = NumericVector(theta.begin(), theta.end()),
    Named("Q") = Q,
    Named("method") = spectral_method
  );
}

// Per-voxel REML search and BLUPs for a design prepared once in R.
// Y: n x V data; Q, theta: from mixed_reml_spectral_cpp(); d, UtX, Vz: thin
// SVD pieces of Z (Z = U diag(d) Vz', UtX = U'X); G = X'(I - U U')X;
// C = X'(I - U U')Y (p x V); B = U'Y (k x V).
// [[Rcpp::export]]
List mixed_reml_blup_cpp(const arma::mat& Y,
                         const arma::mat& Q,
                         const arma::vec& theta,
                         double lower,
                         double upper,
                         const arma::vec& d,
                         const arma::mat& UtX,
                         const arma::mat& G,
                         const arma::mat& C,
                         const arma::mat& B,
                         const arma::mat& Vz) {
  const arma::uword V = Y.n_cols;
  const arma::uword p = UtX.n_cols;
  const arma::uword m = Vz.n_rows;
  const int df = static_cast<int>(theta.n_elem);

  if (C.n_cols != V || B.n_cols != V) {
    stop("mixed_reml_blup_cpp: inconsistent number of responses");
  }
  if (Q.n_rows != Y.n_rows || Q.n_cols != theta.n_elem ||
      B.n_rows != d.n_elem || UtX.n_rows != d.n_elem ||
      Vz.n_cols != d.n_elem || G.n_rows != p || G.n_cols != p ||
      C.n_rows != p) {
    stop("mixed_reml_blup_cpp: inconsistent design dimensions");
  }

  arma::vec lambda_out(V), Vu_out(V), Ve_out(V);
  arma::uvec ok(V, arma::fill::ones);
  arma::mat beta_out(p, V, arma::fill::zeros);
  arma::mat u_out(m, V, arma::fill::zeros);
  const arma::vec d2 = arma::square(d);

  RemlObjective obj(df, theta);
  const double lambda_start = std::max(lower, std::min(arma::mean(theta), upper));

  for (arma::uword v = 0; v < V; ++v) {
    if ((v & 255u) == 0u) Rcpp::checkUserInterrupt();
    arma::vec y = Y.col(v);
    arma::vec omega = Q.t() * y;
    obj.omega_sq = arma::square(omega);

    roptim::Roptim<RemlObjective> opt("L-BFGS-B");
    opt.control.trace = 0;
    opt.control.maxit = 1000;
    opt.control.abstol = 1e-8;
    opt.control.reltol = 1e-8;
    arma::vec lo(1), hi(1);
    lo[0] = lower;
    hi[0] = upper;
    opt.set_lower(lo);
    opt.set_upper(hi);

    arma::vec lam(1);
    lam[0] = lambda_start;
    opt.minimize(obj, lam);
    const double lambda = lam[0];

    const double Vu = arma::sum(obj.omega_sq / (theta + lambda)) / df;
    lambda_out[v] = lambda;
    Vu_out[v] = Vu;
    Ve_out[v] = lambda * Vu;

    // GLS fixed effects: (lambda X'H^-1 X) beta = lambda X'H^-1 y.
    const arma::vec h = lambda / (d2 + lambda);
    arma::mat hA = UtX;
    hA.each_col() %= h;
    const arma::mat M = G + UtX.t() * hA;
    const arma::vec rhs = C.col(v) + hA.t() * B.col(v);
    arma::vec beta;
    bool solved = arma::solve(beta, M, rhs, arma::solve_opts::no_approx);
    if (!solved || !beta.is_finite()) {
      ok[v] = 0;
      beta_out.col(v).fill(NA_REAL);
      u_out.col(v).fill(NA_REAL);
      continue;
    }
    beta_out.col(v) = beta;

    // BLUP: u = Z'H^-1 (y - X beta) = Vz diag(d / (d^2 + lambda)) (U'y - U'X beta).
    const arma::vec r = B.col(v) - UtX * beta;
    u_out.col(v) = Vz * ((d / (d2 + lambda)) % r);
  }

  return List::create(
    Named("lambda") = NumericVector(lambda_out.begin(), lambda_out.end()),
    Named("Vu") = NumericVector(Vu_out.begin(), Vu_out.end()),
    Named("Ve") = NumericVector(Ve_out.begin(), Ve_out.end()),
    Named("beta") = beta_out,
    Named("u") = u_out,
    Named("ok") = LogicalVector(ok.begin(), ok.end())
  );
}
