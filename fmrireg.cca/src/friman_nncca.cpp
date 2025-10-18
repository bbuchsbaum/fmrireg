// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// Unconstrained CCA via whitening + SVD
static inline void cca_unconstrained(const arma::mat& Cxx, const arma::mat& Cyy,
                                     const arma::mat& Cxy, double& rho,
                                     arma::vec& wx, arma::vec& wy, double shrink = 0.0) {
  int p = Cxx.n_rows, q = Cyy.n_rows;
  arma::mat Cxxs = Cxx, Cyys = Cyy;
  if (shrink > 0.0) {
    Cxxs.diag() += shrink * arma::trace(Cxx) / p;
    Cyys.diag() += shrink * arma::trace(Cyy) / q;
  }
  arma::mat Lx = arma::chol(Cxxs, "lower");
  arma::mat Ly = arma::chol(Cyys, "lower");
  arma::mat B  = arma::solve(Lx, Cxy);
  B            = arma::solve(Ly, B.t()).t();
  arma::mat Ux, Uy; arma::vec s;
  arma::svd(Ux, s, Uy, B);
  arma::vec ux = Ux.col(0), uy = Uy.col(0);
  wx = arma::solve(Lx.t(), ux);
  wy = arma::solve(Ly.t(), uy);
  double nx = std::sqrt( arma::as_scalar(wx.t()*Cxx*wx) );
  double ny = std::sqrt( arma::as_scalar(wy.t()*Cyy*wy) );
  if (nx > 0) wx /= nx; if (ny > 0) wy /= ny; rho = s(0);
}

// Exact non-negative CCA by subset enumeration
static inline void nncca_subset_exact(const arma::mat& Cxx, const arma::mat& Cyy,
                                      const arma::mat& Cxy, double& rho,
                                      arma::vec& wx_full, arma::vec& wy_full,
                                      bool simplex_x = false, bool simplex_y = false,
                                      double shrink = 0.0) {
  int p = Cxx.n_rows, q = Cyy.n_rows; rho = -1.0;
  wx_full.zeros(p); wy_full.zeros(q);
  std::vector<std::vector<int>> subsX, subsY;
  subsX.reserve((1u<<p) - 1u); subsY.reserve((1u<<q) - 1u);
  for (int mask = 1; mask < (1<<p); ++mask) {
    std::vector<int> idx; idx.reserve(p);
    for (int i=0;i<p;++i) if (mask&(1<<i)) idx.push_back(i);
    subsX.push_back(std::move(idx));
  }
  for (int mask = 1; mask < (1<<q); ++mask) {
    std::vector<int> idx; idx.reserve(q);
    for (int j=0;j<q;++j) if (mask&(1<<j)) idx.push_back(j);
    subsY.push_back(std::move(idx));
  }
  for (const auto& ix : subsX) {
    arma::uvec I = arma::conv_to<arma::uvec>::from(ix); arma::mat Cxxs = Cxx.submat(I,I);
    for (const auto& iy : subsY) {
      arma::uvec J = arma::conv_to<arma::uvec>::from(iy);
      arma::mat Cyys = Cyy.submat(J,J); arma::mat Cxys = Cxy.submat(I,J);
      double r; arma::vec wx, wy; cca_unconstrained(Cxxs, Cyys, Cxys, r, wx, wy, shrink);
      if (!arma::all(wx >= 0.0) || !arma::all(wy >= 0.0)) continue;
      if (simplex_x) { double sx = arma::accu(wx); if (sx <= 0) continue; wx /= sx; }
      if (simplex_y) { double sy = arma::accu(wy); if (sy <= 0) continue; wy /= sy; }
      if (r > rho) { rho = r; wx_full.zeros(); wy_full.zeros();
        for (size_t a=0; a<ix.size(); ++a) wx_full(ix[a]) = wx(a);
        for (size_t b=0; b<iy.size(); ++b) wy_full(iy[b]) = wy(b);
      }
    }
  }
  if (rho < 0) {
    double r; arma::vec wx, wy; cca_unconstrained(Cxx, Cyy, Cxy, r, wx, wy, shrink);
    wx_full = arma::clamp(wx, 0.0, arma::datum::inf);
    wy_full = arma::clamp(wy, 0.0, arma::datum::inf);
    rho = std::max(0.0, r);
  }
}

// [[Rcpp::export]]
Rcpp::List friman_nncca2xk(const arma::mat& Cxx, const arma::mat& Cyy,
                           const arma::mat& Cxy, double shrink=0.0,
                           bool simplex_x=false, bool simplex_y=false) {
  double rho; arma::vec wx, wy;
  nncca_subset_exact(Cxx, Cyy, Cxy, rho, wx, wy, simplex_x, simplex_y, shrink);
  return Rcpp::List::create(
    Rcpp::Named("rho") = rho,
    Rcpp::Named("wx")  = wx,
    Rcpp::Named("wy")  = wy
  );
}

