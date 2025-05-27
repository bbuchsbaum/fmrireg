#' Autoregressive Modeling Utilities
#'
#' Functions for estimating AR coefficients and applying AR whitening
#' transforms. These are used internally by the linear modelling code.
#'
#' @keywords internal
#' @noRd
NULL

#' Estimate AR parameters
#'
#' Wrapper around `.estimate_ar` that performs basic validation and
#' returns the estimated AR(p) coefficients.
#'
#' @param residuals_vec Numeric vector of residuals.
#' @param p_order Integer autoregressive order.
#'
#' @return Numeric vector of length `p_order` with the estimated
#'   coefficients.
#' @keywords internal
#' @noRd
estimate_ar_parameters <- function(residuals_vec, p_order) {
  stopifnot(is.numeric(residuals_vec))
  stopifnot(length(p_order) == 1L, p_order >= 1L)
  if (anyNA(residuals_vec)) {
    stop("NA values detected in 'residuals_vec' for estimate_ar_parameters")
  }
  .estimate_ar(residuals_vec, p_order)
}

#' Apply AR whitening transform
#'
#' Uses `ar_whiten_inplace()` to apply a causal AR filter defined by
#' `phi` to both `X` and `Y` matrices. New matrices are returned and the
#' originals are left unchanged.
#'
#' @param X Design matrix (time points \eqn{\times} predictors).
#' @param Y Data matrix (time points \eqn{\times} observations).
#' @param phi Numeric vector of AR coefficients.
#' @param exact_first Logical, apply exact scaling of the first sample
#'   for AR(1).
#'
#' @return List with components `X` and `Y` containing the whitened
#'   matrices.
#' @keywords internal
#' @noRd
ar_whiten_transform <- function(X, Y, phi, exact_first = FALSE) {
  if (anyNA(X) || anyNA(Y)) {
    stop("NA values detected in 'X' or 'Y' for ar_whiten_transform")
  }
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)

  Xw <- X
  Yw <- Y
  ar_whiten_inplace(Yw, Xw, phi, exact_first)
  list(X = Xw, Y = Yw)
}

