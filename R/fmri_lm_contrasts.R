#' Contrast Computation for fMRI Linear Models
#'
#' Functions to compute contrasts, standard errors, and statistics
#' for fitted GLM models
#'
#' @keywords internal
#' @noRd
#' @importFrom MASS ginv
NULL

#' Sandwich variance estimator for heteroscedasticity-robust inference
#'
#' @description
#' Computes the heteroscedasticity-consistent (HC) covariance matrix
#' for robust standard errors.
#'
#' @param X Design matrix
#' @param residuals Residual matrix
#' @param type HC type ("HC0", "HC1", "HC2", "HC3")
#'
#' @return Sandwich variance-covariance matrix
#' @keywords internal
#' @noRd
compute_sandwich_variance <- function(X, residuals, type = "HC1") {
  n <- nrow(X)
  p <- ncol(X)
  
  # Bread: (X'X)^-1
  XtX_inv <- solve(crossprod(X))
  
  # Meat: X' diag(e^2) X
  e2 <- residuals^2
  
  # Adjustment factor
  adjustment <- switch(type,
    "HC0" = 1,
    "HC1" = n / (n - p),
    "HC2" = {
      # Would need hat matrix diagonal
      warning("HC2 not implemented, using HC1")
      n / (n - p)
    },
    "HC3" = {
      warning("HC3 not implemented, using HC1")
      n / (n - p)
    }
  )
  
  # For multiple voxels, average the meat matrix
  if (ncol(residuals) > 1) {
    meat <- matrix(0, p, p)
    for (v in 1:ncol(residuals)) {
      meat <- meat + crossprod(X, X * e2[, v])
    }
    meat <- meat / ncol(residuals)
  } else {
    meat <- crossprod(X, X * as.vector(e2))
  }
  
  # Sandwich
  sandwich_vcov <- adjustment * XtX_inv %*% meat %*% XtX_inv
  
  sandwich_vcov
}