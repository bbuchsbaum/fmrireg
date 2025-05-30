#' Effective Degrees of Freedom Calculations
#'
#' Functions to compute effective degrees of freedom for various
#' model types including AR, robust, and mixed models
#'
#' @keywords internal
#' @noRd
NULL

#' Compute effective degrees of freedom
#'
#' @description
#' Main dispatcher that computes effective df based on the model type.
#' Handles AR models, robust models, and combinations.
#'
#' @param fit_result Result from GLM fitting
#' @param X Design matrix
#' @param method Method used ("ols", "ar", "robust", "ar_robust")
#'
#' @return Numeric effective degrees of freedom
#' @keywords internal
#' @noRd
compute_effective_df <- function(fit_result, X, method = "ols") {
  n <- nrow(X)
  p <- ncol(X)
  base_df <- n - p
  
  switch(method,
    "ols" = base_df,
    "ar" = compute_effective_df_ar(n, p, fit_result$ar_coef),
    "robust" = compute_effective_df_robust(X, fit_result$robust_weights, 
                                           fit_result$XtXinv),
    "ar_robust" = compute_effective_df_ar_robust(n, p, fit_result$ar_coef,
                                                 X, fit_result$robust_weights,
                                                 fit_result$XtXinv),
    base_df  # Default
  )
}

#' Effective df for AR models
#'
#' @description
#' Adjusts degrees of freedom for autocorrelation using the
#' effective sample size approach.
#'
#' @param n Sample size
#' @param p Number of parameters
#' @param ar_coef AR coefficients (can be a list for multiple runs)
#'
#' @return Effective degrees of freedom
#' @keywords internal
#' @noRd
compute_effective_df_ar <- function(n, p, ar_coef) {
  if (is.null(ar_coef)) {
    return(n - p)
  }
  
  # If list (multiple runs), pool the coefficients
  if (is.list(ar_coef)) {
    phi <- mean(unlist(ar_coef))
  } else {
    phi <- ar_coef[1]  # Use first coefficient for AR(p)
  }
  
  # Effective sample size for AR(1)
  # n_eff = n * (1 - phi^2) / (1 + phi^2)
  # Simplified approximation: n_eff ≈ n * (1 - phi^2)
  ar_factor <- 1 - phi^2
  ar_factor <- max(ar_factor, 0.1)  # Prevent too small values
  
  n_eff <- n * ar_factor
  max(n_eff - p, 1)
}

#' Effective df for robust models
#'
#' @description
#' Computes effective df using the trace of the weighted hat matrix.
#' This accounts for downweighting of outliers.
#'
#' @param X Design matrix
#' @param weights Robust weights
#' @param XtWXinv Weighted (X'WX)^-1 matrix
#'
#' @return Effective degrees of freedom
#' @keywords internal
#' @noRd
compute_effective_df_robust <- function(X, weights, XtWXinv) {
  n <- nrow(X)
  p <- ncol(X)
  
  if (is.null(weights) || all(weights == 1)) {
    return(n - p)
  }
  
  # For robust regression with weights, the effective df is:
  # df = sum(weights) - p
  # This is because downweighted observations contribute less to df
  # 
  # Alternative formulation using hat matrix:
  # H = X(X'WX)^-1 X'W
  # df = tr(W) - tr(H) = sum(weights) - p
  # But the trace of H for weighted regression equals p when computed correctly
  
  # Use the sum of weights approach
  effective_n <- sum(weights)
  
  # Effective df is effective sample size minus parameters
  max(effective_n - p, 1)
}

#' Effective df for AR + robust models
#'
#' @description
#' Combines adjustments for both autocorrelation and robust weights.
#' Uses a multiplicative adjustment approach.
#'
#' @param n Sample size
#' @param p Number of parameters
#' @param ar_coef AR coefficients
#' @param X Design matrix (after whitening)
#' @param weights Robust weights
#' @param XtWXinv Weighted covariance matrix
#'
#' @return Effective degrees of freedom
#' @keywords internal
#' @noRd
compute_effective_df_ar_robust <- function(n, p, ar_coef, X, weights, XtWXinv) {
  # First get AR adjustment
  df_ar <- compute_effective_df_ar(n, p, ar_coef)
  ar_ratio <- df_ar / (n - p)
  
  # Then get robust adjustment on whitened data
  df_robust <- compute_effective_df_robust(X, weights, XtWXinv)
  robust_ratio <- df_robust / (n - p)
  
  # Combined adjustment
  combined_ratio <- ar_ratio * robust_ratio
  
  max(combined_ratio * (n - p), 1)
}

#' Satterthwaite approximation for df
#'
#' @description
#' Approximates degrees of freedom using Satterthwaite's method
#' for linear combinations of chi-squared variables.
#'
#' @param variance_components Vector of variance components
#' @param df_components Degrees of freedom for each component
#'
#' @return Approximate degrees of freedom
#' @keywords internal
#' @noRd
satterthwaite_df <- function(variance_components, df_components) {
  # df ≈ (sum(v_i))^2 / sum(v_i^2 / df_i)
  num <- sum(variance_components)^2
  denom <- sum(variance_components^2 / df_components)
  num / denom
}

#' Kenward-Roger approximation
#'
#' @description
#' Placeholder for Kenward-Roger df approximation for mixed models.
#' This is complex and would require additional dependencies.
#'
#' @param model Mixed model object
#'
#' @return Approximate degrees of freedom
#' @keywords internal
#' @noRd
kenward_roger_df <- function(model) {
  warning("Kenward-Roger approximation not implemented")
  # Would require implementation of:
  # - Adjusted covariance matrix
  # - Bias correction terms
  # - Complex matrix calculations
  NA
}

#' Residual effective df for diagnostics
#'
#' @description
#' Computes effective residual df for use in diagnostic plots
#' and residual variance estimation.
#'
#' @param fit_result GLM fit result
#' @param method Fitting method used
#'
#' @return Residual degrees of freedom
#' @keywords internal
#' @noRd
residual_effective_df <- function(fit_result, method = "ols") {
  if (!is.null(fit_result$effective_df)) {
    return(fit_result$effective_df)
  }
  
  if (!is.null(fit_result$dfres)) {
    # Adjust for special cases
    if (method %in% c("ar", "robust", "ar_robust") && 
        !is.null(fit_result$rank)) {
      # Further adjust if rank deficient
      return(fit_result$dfres)
    }
    return(fit_result$dfres)
  }
  
  # Fallback
  warning("No df information available")
  NA
}