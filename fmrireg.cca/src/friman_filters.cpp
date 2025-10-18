// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

static inline arma::vec gauss1d(double sigma, int radius) {
  int L = 2*radius + 1; arma::vec k(L);
  double s2 = sigma*sigma, norm = 0.0;
  for (int i=-radius, j=0; i<=radius; ++i, ++j) {
    double v = std::exp(-0.5 * (i*i)/s2); k(j) = v; norm += v;
  }
  return k / norm;
}
static inline arma::vec gauss1d_d1(double sigma, int radius) {
  int L = 2*radius + 1; arma::vec g = gauss1d(sigma, radius), k(L);
  double s2 = sigma*sigma;
  for (int i=-radius, j=0; i<=radius; ++i, ++j) k(j) = -(i/s2) * g(j);
  return k;
}
static inline arma::vec gauss1d_d2(double sigma, int radius) {
  int L = 2*radius + 1; arma::vec g = gauss1d(sigma, radius), k(L);
  double s2 = sigma*sigma, s4 = s2*s2;
  for (int i=-radius, j=0; i<=radius; ++i, ++j) k(j) = ((i*i - s2)/s4) * g(j);
  return k;
}

static inline void convolve_along_x(const arma::cube& in, arma::cube& out, const arma::vec& k) {
  int nx = in.n_rows, ny = in.n_cols, nz = in.n_slices, R = (k.n_elem-1)/2;
  out.zeros(nx,ny,nz);
  for (int z=0; z<nz; ++z)
    for (int y=0; y<ny; ++y)
      for (int x=0; x<nx; ++x) {
        double acc=0.0;
        for (int r=-R; r<=R; ++r) {
          int xx = x+r; if (xx<0) xx=-xx; if (xx>=nx) xx = 2*nx-xx-2;
          acc += k(r+R) * in(xx,y,z);
        }
        out(x,y,z) = acc;
      }
}

static inline arma::cube separable_conv3(const arma::cube& vol,
                                         const arma::vec& kx,
                                         const arma::vec& ky,
                                         const arma::vec& kz) {
  arma::cube tmp1(vol.n_rows, vol.n_cols, vol.n_slices);
  arma::cube tmp2(vol.n_rows, vol.n_cols, vol.n_slices);
  convolve_along_x(vol, tmp1, kx);
  int nx=vol.n_rows, ny=vol.n_cols, nz=vol.n_slices, R=(ky.n_elem-1)/2;
  for (int z=0; z<nz; ++z)
    for (int x=0; x<nx; ++x)
      for (int y=0; y<ny; ++y) {
        double acc=0.0;
        for (int r=-R; r<=R; ++r) {
          int yy=y+r; if (yy<0) yy=-yy; if (yy>=ny) yy=2*ny-yy-2;
          acc += ky(r+R) * tmp1(x,yy,z);
        }
        tmp2(x,y,z) = acc;
      }
  arma::cube out(nx,ny,nz);
  R = (kz.n_elem-1)/2;
  for (int y=0; y<ny; ++y)
    for (int x=0; x<nx; ++x)
      for (int z=0; z<nz; ++z) {
        double acc=0.0;
        for (int r=-R; r<=R; ++r) {
          int zz=z+r; if (zz<0) zz=-zz; if (zz>=nz) zz=2*nz-zz-2;
          acc += kz(r+R) * tmp2(x,y,zz);
        }
        out(x,y,z) = acc;
      }
  return out;
}

// [[Rcpp::export]]
Rcpp::List friman_base_responses_3d(const arma::cube& vol,
                                    arma::vec spacing_mm,
                                    double fwhm_mm) {
  double c = 1.0 / (2.0*std::sqrt(2.0*std::log(2.0)));
  arma::vec sigma_vox = (fwhm_mm * c) / spacing_mm; // per-axis
  auto k0x = gauss1d(sigma_vox(0), std::ceil(3*sigma_vox(0)));
  auto k0y = gauss1d(sigma_vox(1), std::ceil(3*sigma_vox(1)));
  auto k0z = gauss1d(sigma_vox(2), std::ceil(3*sigma_vox(2)));
  auto k1x = gauss1d_d1(sigma_vox(0), k0x.n_elem/2);
  auto k1y = gauss1d_d1(sigma_vox(1), k0y.n_elem/2);
  auto k1z = gauss1d_d1(sigma_vox(2), k0z.n_elem/2);
  auto k2x = gauss1d_d2(sigma_vox(0), k0x.n_elem/2);
  auto k2y = gauss1d_d2(sigma_vox(1), k0y.n_elem/2);
  auto k2z = gauss1d_d2(sigma_vox(2), k0z.n_elem/2);

  arma::cube G0 = separable_conv3(vol, k0x, k0y, k0z);
  arma::cube Gxx = separable_conv3(vol, k2x, k0y, k0z);
  arma::cube Gyy = separable_conv3(vol, k0x, k2y, k0z);
  arma::cube Gzz = separable_conv3(vol, k0x, k0y, k2z);
  arma::cube Gxy = separable_conv3(vol, k1x, k1y, k0z);
  arma::cube Gxz = separable_conv3(vol, k1x, k0y, k1z);
  arma::cube Gyz = separable_conv3(vol, k0x, k1y, k1z);

  return Rcpp::List::create(
    _["G0"]=G0, _["Gxx"]=Gxx, _["Gyy"]=Gyy, _["Gzz"]=Gzz,
    _["Gxy"]=Gxy, _["Gxz"]=Gxz, _["Gyz"]=Gyz
  );
}

