#' Calculate Effective Degrees of Freedom
#'
#' Calculates effective degrees of freedom for robust regression and/or
#' autoregressive models. This adjustment accounts for the loss of efficiency
#' when using robust estimators or AR error structures.
#'
#' @param n Integer. Number of observations.
#' @param p Integer. Number of parameters.
#' @param robust_weights Numeric vector of robust weights (NULL if not robust).
#' @param ar_order Integer. Order of AR model (0 if no AR).
#' @param method Character. Method for DF adjustment ("satterthwaite", "simple").
#'
#' @return Numeric. Adjusted degrees of freedom.
#'
#' @details
#' For robust regression, the effective DF is reduced based on the 
#' downweighting of outliers. For AR models, DF is reduced based on
#' the autocorrelation structure.
#'
#' The Satterthwaite approximation provides a more accurate adjustment
#' but requires additional computation.
#'
#' @keywords internal
#' @noRd
calculate_effective_df <- function(n, p, robust_weights = NULL, 
                                   ar_order = 0, method = c("simple", "satterthwaite")) {
  method <- match.arg(method)
  
  # Base degrees of freedom
  df_base <- n - p
  
  # Adjustment for robust regression
  df_robust_adjust <- 1
  if (!is.null(robust_weights)) {
    # Effective sample size based on weights
    # Weights close to 0 indicate outliers that don't contribute
    eff_n <- sum(robust_weights)
    df_robust_adjust <- eff_n / n
  }
  
  # Adjustment for AR models
  df_ar_adjust <- 1
  if (ar_order > 0) {
    # Simple adjustment: lose ar_order degrees of freedom
    # More sophisticated: could use spectral density at zero
    if (method == "simple") {
      df_ar_adjust <- (n - ar_order) / n
    } else {
      # Satterthwaite-style adjustment would go here
      # For now, use simple method
      df_ar_adjust <- (n - ar_order) / n
    }
  }
  
  # Combined effective DF
  df_effective <- df_base * df_robust_adjust * df_ar_adjust
  
  # Ensure positive
  max(df_effective, 1)
}

#' Calculate Sandwich Variance Estimator
#'
#' Computes the sandwich (robust) variance estimator for regression
#' coefficients. This is useful when model assumptions are violated.
#'
#' @param X Design matrix.
#' @param residuals Vector of residuals.
#' @param XtXinv Inverse of X'X (can be weighted).
#' @param weights Optional weights (for robust regression).
#'
#' @return Matrix. Sandwich variance-covariance matrix.
#'
#' @details
#' The sandwich estimator is computed as:
#' \deqn{V = (X'X)^{-1} X' diag(e^2) X (X'X)^{-1}}
#' 
#' where e are the residuals. This provides valid standard errors
#' even under heteroscedasticity.
#'
#' @keywords internal
#' @noRd
calculate_sandwich_variance <- function(X, residuals, XtXinv, weights = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Squared residuals (meat of sandwich)
  resid2 <- residuals^2
  
  if (!is.null(weights)) {
    # Weight the squared residuals
    resid2 <- resid2 * weights
  }
  
  # Calculate meat: X' diag(resid2) X
  if (is.matrix(residuals) && ncol(residuals) > 1) {
    # Multiple response case
    # Average over voxels for now (could be improved)
    resid2_avg <- rowMeans(resid2)
    meat <- crossprod(X * sqrt(resid2_avg))
  } else {
    meat <- crossprod(X * sqrt(resid2))
  }
  
  # Sandwich: (X'X)^-1 meat (X'X)^-1
  sandwich <- XtXinv %*% meat %*% XtXinv
  
  # Small sample correction
  correction <- n / (n - p)
  
  sandwich * correction
}

#' Adjust Statistics for Effective Degrees of Freedom
#'
#' Updates t-statistics and p-values using effective degrees of freedom
#' rather than nominal degrees of freedom.
#'
#' @param t_stats Vector of t-statistics.
#' @param df_nominal Nominal degrees of freedom.
#' @param df_effective Effective degrees of freedom.
#'
#' @return List with adjusted t-statistics and p-values.
#'
#' @keywords internal
#' @noRd
adjust_stats_for_effective_df <- function(t_stats, df_nominal, df_effective) {
  # Adjustment factor for t-statistics
  # Based on the relationship between t-distributions with different df
  adjust_factor <- sqrt(df_effective / df_nominal)
  
  # Adjusted t-statistics
  t_adj <- t_stats * adjust_factor
  
  # P-values using effective df
  p_values <- 2 * pt(abs(t_adj), df = df_effective, lower.tail = FALSE)
  
  list(
    t_stats = t_adj,
    p_values = p_values,
    df = df_effective
  )
}