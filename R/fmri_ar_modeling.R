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
#' `phi` to both `X` and `Y` matrices. Returns whitened copies while
#' leaving the originals unchanged (due to R's pass-by-value semantics).
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

  # Debug
  if (getOption("fmrireg.debug.ar", FALSE)) {
    message("ar_whiten_transform called with phi=", paste(phi, collapse=","), 
            " exact_first=", exact_first)
    message("  X dims: ", nrow(X), "x", ncol(X))
    message("  Y dims: ", nrow(Y), "x", ncol(Y))
  }

  # Make explicit copies to ensure consistent behavior
  # The C++ function behavior varies depending on how the package is loaded
  # (modifies in place with devtools::load_all, doesn't with library())
  # Making copies ensures originals are never modified
  X_copy <- X + 0  # Force copy
  Y_copy <- Y + 0  # Force copy
  
  # Call ar_whiten_inplace with copies
  # Note: ar_whiten_inplace takes (Y, X) and returns list(Y=..., X=...)
  result <- ar_whiten_inplace(Y_copy, X_copy, phi, exact_first)
  list(X = result$X, Y = result$Y)
}

