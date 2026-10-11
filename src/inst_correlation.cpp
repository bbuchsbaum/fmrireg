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
  if (N == 0)                          // nothing to do; y[0] / y[N-1] would be out of bounds
    return NumericVector(0);

  const double a   = std::exp(-eta);
  const double eps = std::numeric_limits<double>::epsilon();

  if (warmup < 0)                      // default warm‑up
    warmup = std::ceil(4.0 / eta);
  warmup = std::min(warmup, N);

  NumericVector rho(N, NA_REAL);

  // EWM state. A sample whose x or y is missing (fill = "na" outside the
  // series, or NA in the data) leaves the state unchanged and gets rho = NA;
  // the recursion starts at the first complete pair. Without this, one NA
  // propagates through the recursion into every later value.
  bool   started = false;
  double mu_x = 0.0, mu_y = 0.0, s_xy = 0.0, s_xx = 0.0, s_yy = 0.0;

  for (int k = 0; k < N; ++k) {
    const double xk = x[k];
    const double yk = y_at(y, k - offset, N, fill, y[0], y[N-1]);
    if (ISNAN(xk) || ISNAN(yk))
      continue;

    if (!started) {                    // initial point: rho undefined
      mu_x = xk;      mu_y = yk;
      s_xy = xk * yk; s_xx = xk * xk; s_yy = yk * yk;
      started = true;
      continue;
    }

    // EWM means
    mu_x = a * mu_x + (1.0 - a) * xk;
    mu_y = a * mu_y + (1.0 - a) * yk;

    // EWM second‑order moments
    s_xy = a * s_xy + (1.0 - a) * xk * yk;
    s_xx = a * s_xx + (1.0 - a) * xk * xk;
    s_yy = a * s_yy + (1.0 - a) * yk * yk;

    const double cov  = s_xy - mu_x * mu_y;
    const double varx = s_xx - mu_x * mu_x;
    const double vary = s_yy - mu_y * mu_y;
    const double denom = std::sqrt(std::max(varx * vary, eps));

    rho[k] = (denom > 0.0) ? cov / denom : NA_REAL;
  }

  // ---- warm‑up period set to NA -----------------------------------------
  for (int i = 0; i < warmup; ++i)
    rho[i] = NA_REAL;

  return rho;
}
