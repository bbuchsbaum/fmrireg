// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

Rcpp::List friman_base_responses_2d(const arma::mat& slice, arma::vec spacing_mm, double fwhm_mm);
arma::mat steerable_dirs_2d();

struct AB2 { double A, B; };
static inline AB2 calibrate_AB2d_simple() { return AB2{1.0, -1.0}; }

// [[Rcpp::export]]
Rcpp::List friman_pass2_xts_sts_2d(const Rcpp::List& vols,
                                   arma::vec spacing_mm2d,
                                   double fwhm_mm,
                                   const arma::mat& X,
                                   const arma::mat& w_dir,
                                   const arma::mat& w_step2,
                                   const arma::uvec& mask_lin,
                                   const arma::uvec& mask_z) {
  int T = vols.size(), p = X.n_cols, V = mask_lin.n_elem;
  arma::mat XtS(p, V, arma::fill::zeros);
  arma::rowvec StS(V, arma::fill::zeros);

  AB2 ab;
  try {
    Rcpp::List abL = Rcpp::Environment::namespace_env("fmrireg.cca")["friman_calibrate_AB_2d"](1.0);
    // Above call placeholder; we will calibrate below using spacing/fwhm avg sigma
  } catch(...) {}
  // Calibrate using average sigma in vox (fallback simple if issues)
  double c = 1.0 / (2.0*std::sqrt(2.0*std::log(2.0)));
  double sigma_vox_avg = (fwhm_mm * c) / ((spacing_mm2d(0)+spacing_mm2d(1))/2.0);
  try {
    Rcpp::List abList = friman_calibrate_AB_2d(sigma_vox_avg, 1.0, 0.0, -1);
    ab.A = Rcpp::as<double>(abList["A"]);
    ab.B = Rcpp::as<double>(abList["B"]);
  } catch(...) {
    ab = calibrate_AB2d_simple();
  }
  arma::mat N = steerable_dirs_2d();

  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    arma::mat XtS_local(p, V, arma::fill::zeros);
    arma::rowvec StS_local(V, arma::fill::zeros);

    #ifdef _OPENMP
    #pragma omp for schedule(dynamic)
    #endif
    for (int t=0; t<T; ++t) {
      arma::cube vol = Rcpp::as<arma::cube>(vols[t]);
      int nx=vol.n_rows, ny=vol.n_cols, nz=vol.n_slices;
      std::vector<Rcpp::List> base(nz);
      for (int z=0; z<nz; ++z) base[z] = friman_base_responses_2d(vol.slice(z), spacing_mm2d, fwhm_mm);

      for (int v=0; v<V; ++v) {
        unsigned idx = mask_lin(v); int z = mask_z(v);
        int xy = idx - z*(nx*ny); int y = xy / nx; int x = xy - y*nx;
        arma::mat G0  = base[z]["G0"]; arma::mat Gxx = base[z]["Gxx"]; arma::mat Gyy = base[z]["Gyy"]; arma::mat Gxy = base[z]["Gxy"];
        double x_c = G0(x,y);
        double x_o = 0.0;
        for (int i=0; i<3; ++i) {
          double nx_i=N(0,i), ny_i=N(1,i);
          double D2 = nx_i*nx_i*Gxx(x,y) + ny_i*ny_i*Gyy(x,y) + 2.0*nx_i*ny_i*Gxy(x,y);
          double Ori = ab.A * x_c + ab.B * D2;
          x_o += w_dir(v,i) * Ori;
        }
        double s = w_step2(v,0) * x_c + w_step2(v,1) * (x_c - x_o);
        for (int j=0; j<p; ++j) XtS_local(j,v) += X(t,j) * s;
        StS_local(v) += s*s;
      }
    }

    #ifdef _OPENMP
    #pragma omp critical
    #endif
    {
      XtS += XtS_local; StS += StS_local;
    }
  }

  return Rcpp::List::create(_["XtS"]=XtS, _["StS"]=StS);
}

// [[Rcpp::export]]
arma::mat friman_apply_series_2d(const Rcpp::List& vols,
                                 arma::vec spacing_mm2d,
                                 double fwhm_mm,
                                 const arma::mat& w_dir,
                                 const arma::mat& w_step2,
                                 const arma::uvec& mask_lin,
                                 const arma::uvec& mask_z) {
  int T = vols.size(); int V = mask_lin.n_elem;
  arma::mat S(T, V, arma::fill::zeros);
  AB2 ab = calibrate_AB2d_simple();
  arma::mat N = steerable_dirs_2d();
  for (int t=0; t<T; ++t) {
    arma::cube vol = Rcpp::as<arma::cube>(vols[t]);
    int nx=vol.n_rows, ny=vol.n_cols, nz=vol.n_slices;
    std::vector<Rcpp::List> base(nz);
    for (int z=0; z<nz; ++z) base[z] = friman_base_responses_2d(vol.slice(z), spacing_mm2d, fwhm_mm);
    for (int v=0; v<V; ++v) {
      unsigned idx = mask_lin(v); int z = mask_z(v);
      int xy = idx - z*(nx*ny); int y = xy / nx; int x = xy - y*nx;
      arma::mat G0  = base[z]["G0"]; arma::mat Gxx = base[z]["Gxx"]; arma::mat Gyy = base[z]["Gyy"]; arma::mat Gxy = base[z]["Gxy"];
      double x_c = G0(x,y);
      double x_o = 0.0;
      for (int i=0; i<3; ++i) {
        double nx_i=N(0,i), ny_i=N(1,i);
        double D2 = nx_i*nx_i*Gxx(x,y) + ny_i*ny_i*Gyy(x,y) + 2.0*nx_i*ny_i*Gxy(x,y);
        double Ori = ab.A * x_c + ab.B * D2;
        x_o += w_dir(v,i) * Ori;
      }
      double s = w_step2(v,0) * x_c + w_step2(v,1) * (x_c - x_o);
      S(t, v) = s;
    }
  }
  return S;
}
