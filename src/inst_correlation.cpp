#include <Rcpp.h>
using namespace Rcpp;

// helper ------------------------------------------------------------------
inline double y_at(const NumericVector &y,
                   int                  j,          // index we want
                   int                  N,
                   const std::string   &fill,
                   double               last_left,  // y[0]
                   double               last_right) // y[N-1]
{
  if (j >= 0 && j < N)          return y[j];
  if (fill == "zero")           return 0.0;
  if (fill == "na")             return NA_REAL;
  /* fill == "last" */
  return (j < 0) ? last_left : last_right;
}

// [[Rcpp::export]]
NumericVector instantaneous_correlation_rcpp(NumericVector  x,
                                             NumericVector  y,
                                             double         eta      = NA_REAL,
                                             double         tau_half = NA_REAL,
                                             int            offset   = 0,
                                             int            warmup   = -1,
                                             std::string    fill     = "zero")
{
  // ---- argument checking -------------------------------------------------
  if (x.size() != y.size())
    stop("x and y must have the same length");
  if (NumericVector::is_na(eta) && NumericVector::is_na(tau_half))
    stop("Provide either 'eta' or 'tau_half'");
  if (!NumericVector::is_na(tau_half))          // tau_half wins
    eta = std::log(2.0) / tau_half;
  if (eta <= 0.0)
    stop("eta must be > 0");

  const int    N   = x.size();
  const double a   = std::exp(-eta);
  const double eps = std::numeric_limits<double>::epsilon();

  if (warmup < 0)                      // default warm‑up
    warmup = std::ceil(4.0 / eta);
  warmup = std::min(warmup, N);

  NumericVector mu_x(N), mu_y(N),
                s_xy(N), s_xx(N), s_yy(N),
                rho (N);

  // ---- initial point -----------------------------------------------------
  double y0 = y_at(y, -offset, N, fill, y[0], y[N-1]);

  mu_x[0] = x[0];
  mu_y[0] = y0;
  s_xy[0] = x[0] * y0;
  s_xx[0] = x[0] * x[0];
  s_yy[0] = y0   * y0;
  rho [0] = NA_REAL;                         // undefined at first sample

  // ---- recursion ---------------------------------------------------------
  for (int k = 1; k < N; ++k) {
    const double xk = x[k];
    const double yk = y_at(y, k - offset, N, fill, y[0], y[N-1]);

    // EWM means
    mu_x[k] = a * mu_x[k-1] + (1.0 - a) * xk;
    mu_y[k] = a * mu_y[k-1] + (1.0 - a) * yk;

    // EWM second‑order moments
    s_xy[k] = a * s_xy[k-1] + (1.0 - a) * xk * yk;
    s_xx[k] = a * s_xx[k-1] + (1.0 - a) * xk * xk;
    s_yy[k] = a * s_yy[k-1] + (1.0 - a) * yk * yk;

    const double cov = s_xy[k] - mu_x[k] * mu_y[k];
    const double varx = s_xx[k] - mu_x[k] * mu_x[k];
    const double vary = s_yy[k] - mu_y[k] * mu_y[k];
    const double denom = std::sqrt(std::max(varx * vary, eps));

    rho[k] = (denom > 0.0) ? cov / denom : NA_REAL;
  }

  // ---- warm‑up period set to NA -----------------------------------------
  for (int i = 0; i < warmup; ++i)
    rho[i] = NA_REAL;

  return rho;
}