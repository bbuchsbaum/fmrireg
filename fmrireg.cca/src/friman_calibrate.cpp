// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

static inline double gauss_iso_2d(int x, int y, double sigma) {
  double r2 = double(x)*x + double(y)*y;
  return std::exp(-0.5 * r2 / (sigma*sigma));
}

static inline double gauss_iso_3d(int x, int y, int z, double sigma) {
  double r2 = double(x)*x + double(y)*y + double(z)*z;
  return std::exp(-0.5 * r2 / (sigma*sigma));
}

static inline double giso_weight_2d(int x, int y, double sigma_half) {
  double r2 = double(x)*x + double(y)*y;
  return std::exp(-0.5 * r2 / (sigma_half*sigma_half));
}

static inline double giso_weight_3d(int x, int y, int z, double sigma_half) {
  double r2 = double(x)*x + double(y)*y + double(z)*z;
  return std::exp(-0.5 * r2 / (sigma_half*sigma_half));
}

// Directional 2nd derivative of isotropic Gaussian in continuous space (2D)
static inline double D2n_gauss_2d(int x, int y, double sigma, double nx, double ny) {
  double proj = nx*double(x) + ny*double(y);
  double G = gauss_iso_2d(x,y,sigma);
  return ((proj*proj - sigma*sigma) / (sigma*sigma*sigma*sigma)) * G;
}

// Directional 2nd derivative of isotropic Gaussian in continuous space (3D)
static inline double D2n_gauss_3d(int x, int y, int z, double sigma, double nx, double ny, double nz) {
  double proj = nx*double(x) + ny*double(y) + nz*double(z);
  double G = gauss_iso_3d(x,y,z,sigma);
  return ((proj*proj - sigma*sigma) / (sigma*sigma*sigma*sigma)) * G;
}

// [[Rcpp::export]]
Rcpp::List friman_calibrate_AB_2d(double sigma, double nx = 1.0, double ny = 0.0, int radius = -1) {
  if (radius < 0) radius = std::max(3, int(std::ceil(3.0*sigma)));
  int R = radius;
  int N = (2*R+1) * (2*R+1);
  arma::vec t(N); arma::mat M(N, 2);
  int row = 0;
  for (int y = -R; y <= R; ++y) {
    for (int x = -R; x <= R; ++x) {
      double G0 = gauss_iso_2d(x,y,sigma);
      double giso = giso_weight_2d(x,y, sigma/2.0);
      double r = std::sqrt(double(x)*x + double(y)*y);
      double cos2 = (r > 0.0) ? std::pow((nx*double(x) + ny*double(y))/r, 2.0) : 1.0;
      // Eq. (11): (4/3) * (1 - giso) * (cos^2 - 1/4)
      double gi = (4.0/3.0) * (1.0 - giso) * (cos2 - 0.25);
      double target = gi * G0;
      double D2 = D2n_gauss_2d(x,y,sigma, nx, ny);
      t(row) = target; M(row,0) = G0; M(row,1) = D2; ++row;
    }
  }
  arma::vec ab = arma::solve(M, t);
  return Rcpp::List::create(_["A"] = ab(0), _["B"] = ab(1));
}

// [[Rcpp::export]]
Rcpp::List friman_calibrate_AB_3d(double sigma, Rcpp::NumericVector n, int radius = -1) {
  if (n.size() < 3) Rcpp::stop("n must have length 3");
  double nx = n[0], ny = n[1], nz = n[2];
  if (radius < 0) radius = std::max(3, int(std::ceil(3.0*sigma)));
  int R = radius;
  int N = (2*R+1) * (2*R+1) * (2*R+1);
  arma::vec t(N); arma::mat M(N, 2);
  int row = 0;
  for (int z = -R; z <= R; ++z) {
    for (int y = -R; y <= R; ++y) {
      for (int x = -R; x <= R; ++x) {
        double G0 = gauss_iso_3d(x,y,z,sigma);
        double giso = giso_weight_3d(x,y,z, sigma/2.0);
        double r = std::sqrt(double(x)*x + double(y)*y + double(z)*z);
        double cos2 = (r > 0.0) ? std::pow((nx*double(x) + ny*double(y) + nz*double(z))/r, 2.0) : 1.0;
        // Eq. (14): (1 - giso) * (cos^2 - 1/6)
        double gi = (1.0 - giso) * (cos2 - 1.0/6.0);
        double target = gi * G0;
        double D2 = D2n_gauss_3d(x,y,z,sigma, nx, ny, nz);
        t(row) = target; M(row,0) = G0; M(row,1) = D2; ++row;
      }
    }
  }
  arma::vec ab = arma::solve(M, t);
  return Rcpp::List::create(_["A"] = ab(0), _["B"] = ab(1));
}

