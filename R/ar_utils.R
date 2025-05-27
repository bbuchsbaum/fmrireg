# Autoregressive (AR) Utility Functions

#' Estimate AR(p) Coefficients via Yule-Walker
#'
#' Estimates autoregressive coefficients from a numeric residual time series
#' using the Yule-Walker equations.
#'
#' @param residuals_vec Numeric vector of residuals.
#' @param p_order Integer specifying the AR order `p`.
#' @return Numeric vector of length `p_order` containing the estimated AR coefficients.
#' @keywords internal
#' @noRd
.estimate_ar <- function(residuals_vec, p_order) {
  stopifnot(is.numeric(residuals_vec))
  stopifnot(length(p_order) == 1L, p_order >= 1L)

  res <- residuals_vec - mean(residuals_vec, na.rm = TRUE)
  acvf <- stats::acf(res, lag.max = p_order, plot = FALSE, type = "covariance")$acf
  # acvf[1] is lag 0, next p_order elements correspond to lags 1:p_order
  gamma <- as.numeric(acvf)
  R <- toeplitz(gamma[1:p_order])
  r <- gamma[2:(p_order + 1)]
  phi <- drop(solve(R, r))
  as.numeric(phi)
}
