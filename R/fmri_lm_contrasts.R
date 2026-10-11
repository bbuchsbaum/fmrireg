#' Contrast Computation for fMRI Linear Models
#'
#' Functions to compute contrasts, standard errors, and statistics
#' for fitted GLM models
#'
#' @keywords internal
#' @noRd
#' @importFrom MASS ginv
NULL

#' Compute contrast estimates and statistics
#'
#' @description
#' Computes contrast estimates, standard errors, t-statistics and p-values
#' for a given contrast matrix and GLM fit results.
#'
#' @param fit_result Result from solve_integrated_glm or similar
#' @param contrast_matrix Contrast matrix (rows are contrasts, cols are parameters)
#' @param contrast_names Optional names for contrasts
#'
#' @return List with:
#'   - estimate: Contrast estimates (ncontrasts x nvoxels)
#'   - stderr: Standard errors
#'   - tstat: t-statistics
#'   - pvalue: Two-sided p-values
#'   - df: Degrees of freedom used
#' @keywords internal
#' @noRd
compute_contrast <- function(fit_result, contrast_matrix, contrast_names = NULL) {
  # Ensure contrast matrix
  if (!is.matrix(contrast_matrix)) {
    contrast_matrix <- matrix(contrast_matrix, nrow = 1)
  }
  
  # Check dimensions
  if (ncol(contrast_matrix) != nrow(fit_result$betas)) {
    stop("Contrast matrix columns must match number of parameters")
  }
  
  # Compute contrast estimates: C * beta
  estimates <- contrast_matrix %*% fit_result$betas
  
  # Compute variance of contrasts: C * vcov * C'
  if (!is.null(fit_result$XtXinv)) {
    vcov <- fit_result$XtXinv
  } else {
    stop("No variance-covariance matrix available in fit results")
  }
  
  contrast_var <- contrast_matrix %*% vcov %*% t(contrast_matrix)
  
  # Handle multiple voxels
  nvoxels <- ncol(fit_result$betas)
  ncontrasts <- nrow(contrast_matrix)
  
  # Standard errors: se[c, v] = sqrt(Var_c * sigma2_v), for all voxels at once
  sigma2 <- rep_len(as.numeric(fit_result$sigma2), nvoxels)
  stderr <- sqrt(outer(diag(contrast_var), sigma2))
  dim(stderr) <- c(ncontrasts, nvoxels)
  
  # T-statistics
  tstat <- estimates / stderr
  
  # Degrees of freedom
  df <- fit_result$effective_df %||% fit_result$dfres
  
  # P-values (two-sided)
  pvalue <- 2 * pt(-abs(tstat), df = df)
  
  # Add names if provided
  if (!is.null(contrast_names)) {
    rownames(estimates) <- contrast_names
    rownames(stderr) <- contrast_names
    rownames(tstat) <- contrast_names
    rownames(pvalue) <- contrast_names
  }
  
  list(
    estimate = estimates,
    stderr = stderr,
    tstat = tstat,
    pvalue = pvalue,
    df = df,
    contrast_matrix = contrast_matrix
  )
}

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