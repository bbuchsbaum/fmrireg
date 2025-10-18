// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// Forward decls
Rcpp::List friman_base_responses_3d(const arma::cube& vol, arma::vec spacing_mm, double fwhm_mm);
arma::mat steerable_dirs_3d();

struct AB { double A, B; };
// Use exported calibrator; fallback to simple constants on error
static inline AB calibrate_AB_3d_avg(const arma::vec& spacing_mm, double fwhm_mm) {
  double c = 1.0 / (2.0*std::sqrt(2.0*std::log(2.0)));
  double sigma_vox_avg = (fwhm_mm * c) / arma::mean(spacing_mm);
  AB ab;
  try {
    Rcpp::Function calib = Rcpp::Environment::namespace_env("fmrireg.cca")["friman_calibrate_AB_3d"];
    Rcpp::NumericVector n = Rcpp::NumericVector::create(1.0, 0.0, 0.0);
    Rcpp::List abList = calib(sigma_vox_avg, n, -1);
    ab.A = Rcpp::as<double>(abList["A"]);
    ab.B = Rcpp::as<double>(abList["B"]);
  } catch(...) {
    ab.A = 1.0; ab.B = -1.0;
  }
  return ab;
}

// [[Rcpp::export]]
Rcpp::List friman_pass2_xts_sts_3d(const Rcpp::List& vols,
                                   arma::vec spacing_mm,
                                   double fwhm_mm,
                                   const arma::mat& X,
                                   const arma::mat& w_dir,
                                   const arma::mat& w_step2,
                                   const arma::uvec& mask_linear_idx) {
  int T = vols.size();
  int p = X.n_cols;
  int V = mask_linear_idx.n_elem;
  arma::mat XtS(p, V, arma::fill::zeros);
  arma::rowvec StS(V, arma::fill::zeros);

  AB ab = calibrate_AB_3d_avg(spacing_mm, fwhm_mm);
  arma::mat N = steerable_dirs_3d();

  for (int t=0; t<T; ++t) {
    arma::cube vol = Rcpp::as<arma::cube>(vols[t]);
    Rcpp::List resp = friman_base_responses_3d(vol, spacing_mm, fwhm_mm);
    arma::cube G0  = resp["G0"], Gxx=resp["Gxx"], Gyy=resp["Gyy"], Gzz=resp["Gzz"],
               Gxy = resp["Gxy"], Gxz=resp["Gxz"], Gyz=resp["Gyz"];

    int nx=G0.n_rows, ny=G0.n_cols, nz=G0.n_slices;
    for (int v=0; v<V; ++v) {
      unsigned idx = mask_linear_idx(v);
      int z = idx / (nx*ny);
      int rem = idx - z*(nx*ny);
      int y = rem / nx;
      int x = rem - y*nx;

      double x_c = G0(x,y,z);

      double x_o = 0.0;
      for (int i=0; i<N.n_cols; ++i) {
        double nx_i=N(0,i), ny_i=N(1,i), nz_i=N(2,i);
        double D2 = (nx_i*nx_i)*Gxx(x,y,z) + (ny_i*ny_i)*Gyy(x,y,z) + (nz_i*nz_i)*Gzz(x,y,z)
                  + 2.0*nx_i*ny_i*Gxy(x,y,z) + 2.0*nx_i*nz_i*Gxz(x,y,z) + 2.0*ny_i*nz_i*Gyz(x,y,z);
        double Ori = ab.A * x_c + ab.B * D2;
        x_o += w_dir(v,i) * Ori;
      }

      double s = w_step2(v,0) * x_c + w_step2(v,1) * (x_c - x_o);
      for (int j=0; j<p; ++j) XtS(j,v) += X(t,j) * s;
      StS(v) += s*s;
    }
  }
  return Rcpp::List::create(_["XtS"]=XtS, _["StS"]=StS);
}

// [[Rcpp::export]]
arma::mat friman_apply_series_3d(const Rcpp::List& vols,
                                 arma::vec spacing_mm,
                                 double fwhm_mm,
                                 const arma::mat& w_dir,
                                 const arma::mat& w_step2,
                                 const arma::uvec& mask_linear_idx) {
  int T = vols.size();
  int V = mask_linear_idx.n_elem;
  arma::mat S(T, V, arma::fill::zeros);

  AB ab = calibrate_AB_3d(spacing_mm, fwhm_mm);
  arma::mat N = steerable_dirs_3d();

  for (int t=0; t<T; ++t) {
    arma::cube vol = Rcpp::as<arma::cube>(vols[t]);
    Rcpp::List resp = friman_base_responses_3d(vol, spacing_mm, fwhm_mm);
    arma::cube G0  = resp["G0"], Gxx=resp["Gxx"], Gyy=resp["Gyy"], Gzz=resp["Gzz"],
               Gxy = resp["Gxy"], Gxz=resp["Gxz"], Gyz=resp["Gyz"];

    int nx=G0.n_rows, ny=G0.n_cols, nz=G0.n_slices;
    for (int v=0; v<V; ++v) {
      unsigned idx = mask_linear_idx(v);
      int z = idx / (nx*ny); int rem = idx - z*(nx*ny); int y = rem / nx; int x = rem - y*nx;
      double x_c = G0(x,y,z);
      double x_o = 0.0;
      for (int i=0; i<N.n_cols; ++i) {
        double nx_i=N(0,i), ny_i=N(1,i), nz_i=N(2,i);
        double D2 = (nx_i*nx_i)*Gxx(x,y,z) + (ny_i*ny_i)*Gyy(x,y,z) + (nz_i*nz_i)*Gzz(x,y,z)
                  + 2.0*nx_i*ny_i*Gxy(x,y,z) + 2.0*nx_i*nz_i*Gxz(x,y,z) + 2.0*ny_i*nz_i*Gyz(x,y,z);
        double Ori = ab.A * x_c + ab.B * D2;
        x_o += w_dir(v,i) * Ori;
      }
      double s = w_step2(v,0) * x_c + w_step2(v,1) * (x_c - x_o);
      S(t, v) = s;
    }
  }
  return S;
}
