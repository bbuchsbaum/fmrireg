
#' Despike a Univariate Time Series
#'
#' Robustly remove (or flag) isolated spikes using a running median
#' + MAD test (Hoaglin et al. JASA 1986).
#'
#' @param x       numeric vector – the raw time‑series
#' @param k       odd integer – window length of the running median (default 7)
#' @param thresh  positive – threshold in robust z‑scores (default 6)
#' @param method  how to treat detected spikes:
#'                *"median"* (replace by local median – default),
#'                *"interp"* (linear interpolate),
#'                *"na"* (set to NA, just flag)
#'
#' @return A list with two elements  
#'   * **clean** – despiked series (same length as `x`)  
#'   * **spike_idx** – integer indices of points deemed spikes
#' @examples
#' set.seed(1); y <- sin(1:300 / 20) + rnorm(300, 0, .1)
#' y[c(50, 180)] <- 4           # add two huge spikes
#' out <- despike(y, k = 11)
#' matplot(cbind(y, out$clean), type = "l", lty = 1,
#'         col = c("grey50","red"), ylab = "value")
#' points(out$spike_idx, y[out$spike_idx], pch = 19)
#' @export
despike <- function(x, k = 7, thresh = 6, method = c("median","interp","na"))
{
  stopifnot(is.numeric(x), length(k) == 1L, k %% 2 == 1L, k >= 3,
            thresh > 0)
  method <- match.arg(method)

  med_run <- stats::runmed(x, k = k)
  delta   <- x - med_run
  mad_raw <- stats::median(abs(delta))
  mad     <- 1.4826 * mad_raw                      # ~= sigma under normality

  if (mad < .Machine$double.eps)                  # flat signal → nothing to do
    return(list(clean = x, spike_idx = integer(0)))

  z       <- abs(delta / mad)
  idx     <- which(z > thresh)

  y <- x                                     # copy
  if (length(idx)) {
    if (method == "median") {
      y[idx] <- med_run[idx]
    } else if (method == "interp") {
      y[idx] <- approx(seq_along(x)[-idx], y[-idx],
                       xout = idx, rule = 2)$y
    } else {                                # "na"
      y[idx] <- NA_real_
    }
  }
  list(clean = y, spike_idx = idx)
}

#' Exponentially‑Weighted Instantaneous Correlation
#'
#' Computes, online and in **O(N)** time, the correlation between two
#' signals using an EWMA window:  
#' \eqn{\alpha = e^{-η}} or equivalently
#' half‑life `τ½ = log(2)/η`.
#'
#' @param x,y    numeric vectors of equal length
#' @param eta    decay constant (mutually exclusive with `tau_half`)
#' @param tau_half half‑life of the EWMA window (overrides `eta` if given)
#' @param offset integer lag so that \eqn{\rho(t)} compares
#'               \eqn{x[t]} with \eqn{y[t-offset]}
#' @param warmup integer – set the first `warmup` points to `NA`
#'               (`default = ceiling(4/eta)`)
#' @param fill   value for `y` outside valid range (`"zero"`, `"na"`, `"last"`)
#'
#' @return numeric vector `rho` of length `length(x)`
#'
#' @section Use cases:
#' * **Real‑time functional connectivity** in real‑time fMRI
#' * **Streaming quality control** – detect sudden desynchronisation
#' * **Motion artefact flagging** – correlate FD/FD' with global signal
#' * **Econometrics** – high‑frequency lead‑lag correlations
#'
#' @examples
#' x <- rnorm(1000); y <- filter(x, 0.9, "recursive") + rnorm(1000, 0, .3)
#' rho <- instantaneous_correlation(x, y, tau_half = 10)
#' plot(rho, type = "l")
#' @export
instantaneous_correlation <- function(x, y,
                                      eta      = NULL,
                                      tau_half = NULL,
                                      offset   = 0,
                                      warmup   = NULL,
                                      fill     = c("zero","na","last"))
{
  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y))
  if (is.null(eta)) {
    if (is.null(tau_half)) stop("Provide either 'eta' or 'tau_half'.")
    eta <- log(2) / tau_half
  }
  a    <- exp(-eta)
  fill <- match.arg(fill)
  N    <- length(x)

  # helper for y[k-offset]
  y_shifted <- function(k) {
    j <- k - offset
    if (j < 1L || j > N) {
      if (fill == "zero") return(0)
      if (fill == "na")   return(NA_real_)
      return(ifelse(j < 1L, y[1L], y[N]))   # "last"
    }
    y[j]
  }

  mu_x <- mu_y <- s_xy <- s_xx <- s_yy <- rho <- numeric(N)
  eps  <- .Machine$double.eps

  mu_x[1] <- x[1]
  mu_y[1] <- y_shifted(1)
  s_xy[1] <- x[1] * mu_y[1]
  s_xx[1] <- x[1]^2
  s_yy[1] <- mu_y[1]^2

  for (k in 2:N) {
    xk <- x[k]; yk <- y_shifted(k)

    mu_x[k] <- a * mu_x[k-1] + (1-a) * xk
    mu_y[k] <- a * mu_y[k-1] + (1-a) * yk

    s_xy[k] <- a * s_xy[k-1] + (1-a) * xk * yk
    s_xx[k] <- a * s_xx[k-1] + (1-a) * xk^2
    s_yy[k] <- a * s_yy[k-1] + (1-a) * yk^2

    num  <- s_xy[k] - mu_x[k]*mu_y[k]
    den  <- sqrt( (s_xx[k]-mu_x[k]^2) * (s_yy[k]-mu_y[k]^2) + eps )
    rho[k] <- ifelse(is.finite(den) & den > 0, num/den, NA_real_)
  }

  warmup <- warmup %||% ceiling(4/eta)
  rho[seq_len(min(warmup, N))] <- NA_real_

  rho
}

#' Compute Instantaneous Correlation between Two Signals (Rcpp Version)
#'
#' Computes the instantaneous correlation between two signals using an exponential decay weighting function. This version uses Rcpp for improved performance.
#'
#' @param x Numeric vector. First time-series signal.
#' @param y Numeric vector. Second time-series signal.
#' @param eta Positive numeric value. Decay parameter for the exponential weighting function.
#' @param offset Integer. Temporal offset between the two signals. Defaults to 0.
#' @return Numeric vector of instantaneous correlation coefficients.
#' @examples
#' # Example signals
#' x <- sin(seq(0, 2 * pi, length.out = 100))
#' y <- cos(seq(0, 2 * pi, length.out = 100))
#' # Compute instantaneous correlation using Rcpp
#' rho <- instantaneous_correlation_rcpp(x, y, eta = 0.1)
#' # Plot the result
#' plot(rho, type = 'l', ylab = 'Instantaneous Correlation', xlab = 'Time')
#' @export
instantaneous_correlation_rcpp <- function(x, y, eta, offset = 0) {
  # Ensure inputs are numeric vectors of the same length
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric vectors.")
  }
  N <- length(x)
  if (length(y) != N) {
    stop("x and y must be of the same length.")
  }
  #.Call('_instantaneous_correlation_rcpp', x, y, eta, offset)
}


