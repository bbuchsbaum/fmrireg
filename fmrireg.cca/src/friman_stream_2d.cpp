// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

Rcpp::List friman_base_responses_2d(const arma::mat& slice, arma::vec spacing_mm, double fwhm_mm);
arma::mat steerable_dirs_2d();

struct AB2 { double A, B; };

// Simple calibration for 2D oriented low-pass (fallback constants)
static inline AB2 calibrate_AB2d_simple() { return AB2{1.0, -1.0}; }

// Unconstrained CCA via whitening + SVD
static inline void cca_unconstrained(const arma::mat& Cxx, const arma::mat& Cyy,
                                     const arma::mat& Cxy, double& rho,
                                     arma::vec& wx, arma::vec& wy, double shrink = 0.0) {
  int p = Cxx.n_rows, q = Cyy.n_rows;
  arma::mat Cxxs = Cxx, Cyys = Cyy;
  if (shrink > 0.0) { Cxxs.diag() += shrink * arma::trace(Cxx) / p; Cyys.diag() += shrink * arma::trace(Cyy) / q; }
  arma::mat Lx = arma::chol(Cxxs, "lower"); arma::mat Ly = arma::chol(Cyys, "lower");
  arma::mat B  = arma::solve(Lx, Cxy); B = arma::solve(Ly, B.t()).t();
  arma::mat Ux, Uy; arma::vec s; arma::svd(Ux, s, Uy, B);
  arma::vec ux = Ux.col(0), uy = Uy.col(0);
  wx = arma::solve(Lx.t(), ux); wy = arma::solve(Ly.t(), uy);
  double nx = std::sqrt( arma::as_scalar(wx.t()*Cxx*wx) ); double ny = std::sqrt( arma::as_scalar(wy.t()*Cyy*wy) );
  if (nx > 0) wx /= nx; if (ny > 0) wy /= ny; rho = s(0);
}

// Exact non-negative CCA by subset enumeration (tiny problems)
static inline void nncca_subset_exact(const arma::mat& Cxx, const arma::mat& Cyy,
                                      const arma::mat& Cxy, double& rho,
                                      arma::vec& wx_full, arma::vec& wy_full,
                                      bool simplex_x = false, bool simplex_y = false,
                                      double shrink = 0.0) {
  int p = Cxx.n_rows, q = Cyy.n_rows; rho = -1.0; wx_full.zeros(p); wy_full.zeros(q);
  std::vector<std::vector<int>> subsX, subsY; subsX.reserve((1u<<p)-1u); subsY.reserve((1u<<q)-1u);
  for (int mask = 1; mask < (1<<p); ++mask) { std::vector<int> idx; for (int i=0;i<p;++i) if (mask&(1<<i)) idx.push_back(i); subsX.push_back(std::move(idx)); }
  for (int mask = 1; mask < (1<<q); ++mask) { std::vector<int> idx; for (int j=0;j<q;++j) if (mask&(1<<j)) idx.push_back(j); subsY.push_back(std::move(idx)); }
  for (const auto& ix : subsX) {
    arma::uvec I = arma::conv_to<arma::uvec>::from(ix); arma::mat Cxxs = Cxx.submat(I,I);
    for (const auto& iy : subsY) {
      arma::uvec J = arma::conv_to<arma::uvec>::from(iy); arma::mat Cyys = Cyy.submat(J,J); arma::mat Cxys = Cxy.submat(I,J);
      double r; arma::vec wx, wy; cca_unconstrained(Cxxs, Cyys, Cxys, r, wx, wy, shrink);
      if (!arma::all(wx >= 0.0) || !arma::all(wy >= 0.0)) continue;
      if (simplex_x) { double sx = arma::accu(wx); if (sx <= 0) continue; wx /= sx; }
      if (simplex_y) { double sy = arma::accu(wy); if (sy <= 0) continue; wy /= sy; }
      if (r > rho) { rho = r; wx_full.zeros(); wy_full.zeros(); for (size_t a=0; a<ix.size(); ++a) wx_full(ix[a]) = wx(a); for (size_t b=0; b<iy.size(); ++b) wy_full(iy[b]) = wy(b); }
    }
  }
  if (rho < 0) { double r; arma::vec wx, wy; cca_unconstrained(Cxx, Cyy, Cxy, r, wx, wy, shrink);
    wx_full = arma::clamp(wx, 0.0, arma::datum::inf); wy_full = arma::clamp(wy, 0.0, arma::datum::inf); rho = std::max(0.0, r); }
}

// [[Rcpp::export]]
Rcpp::List friman_stream_weights_2d(const Rcpp::List& vols,
                                    arma::vec spacing_mm2d,
                                    double fwhm_mm,
                                    const arma::mat& Ytilde,   // T×2
                                    const arma::uvec& mask_lin,
                                    const arma::uvec& mask_z,
                                    double shrink = 0.0,
                                    bool simplex = true) {
  int T = vols.size();
  int V = mask_lin.n_elem;
  arma::mat N = steerable_dirs_2d();
  int k = N.n_cols;
  AB2 ab;
  double c = 1.0 / (2.0*std::sqrt(2.0*std::log(2.0)));
  double sigma_vox_avg = (fwhm_mm * c) / ((spacing_mm2d(0)+spacing_mm2d(1))/2.0);
  try {
    Rcpp::List abList = friman_calibrate_AB_2d(sigma_vox_avg, 1.0, 0.0, -1);
    ab.A = Rcpp::as<double>(abList["A"]);
    ab.B = Rcpp::as<double>(abList["B"]);
  } catch(...) {
    ab = calibrate_AB2d_simple();
  }

  arma::mat Cyy = Ytilde.t() * Ytilde; // 2×2

  std::vector<arma::mat> Cxx_dir(V, arma::mat(k,k,arma::fill::zeros));
  std::vector<arma::mat> Cxy_dir(V, arma::mat(k,2,arma::fill::zeros));
  arma::vec v_center(V, arma::fill::zeros);
  std::vector<arma::vec> c_cross(V, arma::vec(k, arma::fill::zeros));
  arma::mat x_cY(V, 2, arma::fill::zeros);

  for (int t=0; t<T; ++t) {
    arma::cube vol = Rcpp::as<arma::cube>(vols[t]); // [nx,ny,nz]
    int nx=vol.n_rows, ny=vol.n_cols, nz=vol.n_slices;
    // Compute 2D base responses per slice
    std::vector<Rcpp::List> base(nz);
    for (int z=0; z<nz; ++z) base[z] = friman_base_responses_2d(vol.slice(z), spacing_mm2d, fwhm_mm);
    double y1 = Ytilde(t,0), y2 = Ytilde(t,1);
    for (int v=0; v<V; ++v) {
      unsigned idx = mask_lin(v); int z = mask_z(v);
      int xy = idx - z*(nx*ny); int y = xy / nx; int x = xy - y*nx;
      arma::mat G0  = base[z]["G0"]; arma::mat Gxx = base[z]["Gxx"]; arma::mat Gyy = base[z]["Gyy"]; arma::mat Gxy = base[z]["Gxy"];
      double x_c = G0(x,y);
      arma::vec xdir(k, arma::fill::zeros);
      for (int i=0; i<k; ++i) {
        double nx_i=N(0,i), ny_i=N(1,i);
        double D2 = nx_i*nx_i*Gxx(x,y) + ny_i*ny_i*Gyy(x,y) + 2.0*nx_i*ny_i*Gxy(x,y);
        double Ori = ab.A * x_c + ab.B * D2;
        xdir(i) = Ori;
      }
      v_center(v) += x_c * x_c;
      c_cross[v]  += x_c * xdir;
      x_cY(v,0)   += x_c * y1; x_cY(v,1) += x_c * y2;
      Cxx_dir[v]  += xdir * xdir.t();
      Cxy_dir[v].col(0) += xdir * y1; Cxy_dir[v].col(1) += xdir * y2;
    }
  }

  arma::mat w_dir(V, k, arma::fill::zeros);
  arma::mat w_step2(V, 2, arma::fill::zeros);
  arma::vec rho1(V, arma::fill::zeros), rho2(V, arma::fill::zeros);
  for (int v=0; v<V; ++v) {
    double r1; arma::vec wx(k), wy(2); nncca_subset_exact(Cxx_dir[v], Cyy, Cxy_dir[v], r1, wx, wy, simplex, true, shrink);
    rho1(v) = r1; w_dir.row(v) = wx.t();
    arma::mat Cxx2(2,2), Cxy2(2,2);
    double v_oriented = arma::as_scalar( wx.t() * Cxx_dir[v] * wx );
    double cov_co     = arma::as_scalar( wx.t() * c_cross[v] );
    Cxx2(0,0) = v_center(v);
    Cxx2(0,1) = v_center(v) - cov_co; Cxx2(1,0) = Cxx2(0,1);
    Cxx2(1,1) = v_center(v) - 2.0*cov_co + v_oriented;
    arma::rowvec coY = (wx.t() * Cxy_dir[v]); // 1×2
    Cxy2.row(0) = x_cY.row(v); Cxy2.row(1) = x_cY.row(v) - coY;
    double r2; arma::vec w2x(2), w2y(2); nncca_subset_exact(Cxx2, Cyy, Cxy2, r2, w2x, w2y, true, true, shrink);
    rho2(v) = r2; w_step2(v,0) = w2x(0); w_step2(v,1) = w2x(1);
  }

  return Rcpp::List::create(_["w_dir"]=w_dir, _["w_step2"]=w_step2, _["rho1"]=rho1, _["rho2"]=rho2);
}
