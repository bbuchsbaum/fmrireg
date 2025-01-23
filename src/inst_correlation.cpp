#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector instantaneous_correlation_rcpp(NumericVector x, NumericVector y, double eta, int offset = 0) {
  int N = x.size();
  
  double a = exp(-eta);
  
  NumericVector x_mean(N);
  NumericVector y_mean(N);
  NumericVector s_xy(N);
  NumericVector s_xx(N);
  NumericVector s_yy(N);
  NumericVector rho(N);
  
  // Initial values
  x_mean[0] = x[0];
  int y_index = 0 - offset;
  double y_value = (y_index >= 0 && y_index < N) ? y[y_index] : 0.0;
  y_mean[0] = y_value;
  s_xy[0] = x[0] * y_value;
  s_xx[0] = x[0] * x[0];
  s_yy[0] = y_value * y_value;
  rho[0] = 0.0;  // Correlation is undefined at first point
  
  // Recursive computation
  for (int k = 1; k < N; k++) {
    // Update indices and values with offset
    double xk = x[k];
    int y_index = k - offset;
    double yk = (y_index >= 0 && y_index < N) ? y[y_index] : 0.0;
    
    // Update means
    x_mean[k] = a * x_mean[k - 1] + (1 - a) * xk;
    y_mean[k] = a * y_mean[k - 1] + (1 - a) * yk;
    
    // Update covariances
    s_xy[k] = a * s_xy[k - 1] + (1 - a) * xk * yk;
    s_xx[k] = a * s_xx[k - 1] + (1 - a) * xk * xk;
    s_yy[k] = a * s_yy[k - 1] + (1 - a) * yk * yk;
    
    // Compute covariance and variances
    double cov_xy = s_xy[k] - x_mean[k] * y_mean[k];
    double var_xx = s_xx[k] - x_mean[k] * x_mean[k];
    double var_yy = s_yy[k] - y_mean[k] * y_mean[k];
    
    // Compute instantaneous correlation
    double denom = sqrt(var_xx * var_yy);
    rho[k] = (denom > 0) ? cov_xy / denom : 0.0;
  }
  
  return rho;
}
