// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// 1D Gaussian and derivatives
static inline arma::vec gauss1d(double sigma, int radius) {
  int L = 2*radius + 1; arma::vec k(L);
  double s2 = sigma*sigma, norm = 0.0;
  for (int i=-radius, j=0; i<=radius; ++i, ++j) { double v = std::exp(-0.5 * (i*i)/s2); k(j)=v; norm+=v; }
  return k / norm;
}
static inline arma::vec gauss1d_d1(double sigma, int radius) {
  int L = 2*radius + 1; arma::vec g = gauss1d(sigma, radius), k(L); double s2=sigma*sigma;
  for (int i=-radius, j=0; i<=radius; ++i, ++j) k(j) = -(i/s2) * g(j); return k;
}
static inline arma::vec gauss1d_d2(double sigma, int radius) {
  int L = 2*radius + 1; arma::vec g = gauss1d(sigma, radius), k(L); double s2=sigma*sigma, s4=s2*s2;
  for (int i=-radius, j=0; i<=radius; ++i, ++j) k(j) = ((i*i - s2)/s4) * g(j); return k;
}

// 2D separable convolution (reflect boundaries)
static inline arma::mat separable_conv2(const arma::mat& img,
                                        const arma::vec& kx,
                                        const arma::vec& ky) {
  int nx = img.n_rows, ny = img.n_cols, Rx = (kx.n_elem-1)/2, Ry = (ky.n_elem-1)/2;
  arma::mat tmp(nx, ny, arma::fill::zeros), out(nx, ny, arma::fill::zeros);
  for (int y=0; y<ny; ++y) {
    for (int x=0; x<nx; ++x) {
      double acc=0.0; for (int r=-Rx; r<=Rx; ++r) { int xx=x+r; if (xx<0) xx=-xx; if (xx>=nx) xx=2*nx-xx-2; acc += kx(r+Rx) * img(xx,y);} tmp(x,y)=acc;
    }
  }
  for (int x=0; x<nx; ++x) {
    for (int y=0; y<ny; ++y) {
      double acc=0.0; for (int r=-Ry; r<=Ry; ++r) { int yy=y+r; if (yy<0) yy=-yy; if (yy>=ny) yy=2*ny-yy-2; acc += ky(r+Ry) * tmp(x,yy);} out(x,y)=acc;
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat steerable_dirs_2d() {
  arma::mat N(2,3);
  N.col(0) = arma::vec({1.0, 0.0});
  const double r3_2 = std::sqrt(3.0)/2.0;
  N.col(1) = arma::vec({0.5,  r3_2});
  N.col(2) = arma::vec({-0.5, r3_2});
  return N;
}

// [[Rcpp::export]]
Rcpp::List friman_base_responses_2d(const arma::mat& slice, arma::vec spacing_mm, double fwhm_mm) {
  double c = 1.0 / (2.0*std::sqrt(2.0*std::log(2.0)));
  arma::vec sigma_vox = (fwhm_mm * c) / spacing_mm; // (sx, sy)
  auto k0x = gauss1d(sigma_vox(0), std::ceil(3*sigma_vox(0)));
  auto k0y = gauss1d(sigma_vox(1), std::ceil(3*sigma_vox(1)));
  auto k1x = gauss1d_d1(sigma_vox(0), k0x.n_elem/2);
  auto k1y = gauss1d_d1(sigma_vox(1), k0y.n_elem/2);
  auto k2x = gauss1d_d2(sigma_vox(0), k0x.n_elem/2);
  auto k2y = gauss1d_d2(sigma_vox(1), k0y.n_elem/2);

  arma::mat G0  = separable_conv2(slice, k0x, k0y);
  arma::mat Gxx = separable_conv2(slice, k2x, k0y);
  arma::mat Gyy = separable_conv2(slice, k0x, k2y);
  arma::mat Gxy = separable_conv2(slice, k1x, k1y);

  return Rcpp::List::create(_["G0"]=G0, _["Gxx"]=Gxx, _["Gyy"]=Gyy, _["Gxy"]=Gxy);
}

