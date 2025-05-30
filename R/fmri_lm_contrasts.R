#' Contrast Computation for fMRI Linear Models
#'
#' Functions to compute contrasts, standard errors, and statistics
#' for fitted GLM models
#'
#' @keywords internal
#' @noRd
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
  
  # Standard errors
  if (length(fit_result$sigma2) == 1) {
    # Single sigma2
    stderr <- sqrt(diag(contrast_var) * fit_result$sigma2)
    stderr <- matrix(stderr, ncontrasts, nvoxels)
  } else {
    # Different sigma2 per voxel
    stderr <- matrix(NA, ncontrasts, nvoxels)
    for (v in 1:nvoxels) {
      stderr[, v] <- sqrt(diag(contrast_var) * fit_result$sigma2[v])
    }
  }
  
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

#' Compute F-statistic for multi-row contrasts
#'
#' @description
#' Computes F-statistics and p-values for contrasts with multiple rows
#' (testing multiple linear combinations simultaneously).
#'
#' @param fit_result GLM fit results
#' @param contrast_matrix Contrast matrix (can have multiple rows)
#'
#' @return List with:
#'   - fstat: F-statistics (1 x nvoxels)
#'   - pvalue: p-values
#'   - df1: Numerator degrees of freedom
#'   - df2: Denominator degrees of freedom
#' @keywords internal
#' @noRd
compute_f_statistic <- function(fit_result, contrast_matrix) {
  if (!is.matrix(contrast_matrix)) {
    contrast_matrix <- matrix(contrast_matrix, nrow = 1)
  }
  
  # Number of contrasts (restrictions)
  q <- nrow(contrast_matrix)
  
  # Contrast estimates
  Cb <- contrast_matrix %*% fit_result$betas
  
  # Variance of contrasts
  vcov <- fit_result$XtXinv
  CVC <- contrast_matrix %*% vcov %*% t(contrast_matrix)
  
  # Check if CVC is invertible
  CVC_inv <- tryCatch(
    solve(CVC),
    error = function(e) {
      warning("Contrast variance matrix is singular")
      MASS::ginv(CVC)  # Use generalized inverse
    }
  )
  
  nvoxels <- ncol(fit_result$betas)
  fstat <- numeric(nvoxels)
  
  # Compute F-statistic for each voxel
  for (v in 1:nvoxels) {
    cb_v <- Cb[, v, drop = FALSE]
    sigma2_v <- if (length(fit_result$sigma2) == 1) {
      fit_result$sigma2
    } else {
      fit_result$sigma2[v]
    }
    
    # F = (Cb)' (C vcov C')^-1 (Cb) / (q * sigma2)
    fstat[v] <- as.numeric(t(cb_v) %*% CVC_inv %*% cb_v) / (q * sigma2_v)
  }
  
  # Degrees of freedom
  df1 <- q
  df2 <- fit_result$effective_df %||% fit_result$dfres
  
  # P-values
  pvalue <- pf(fstat, df1, df2, lower.tail = FALSE)
  
  list(
    fstat = matrix(fstat, nrow = 1),
    pvalue = matrix(pvalue, nrow = 1),
    df1 = df1,
    df2 = df2,
    contrast_matrix = contrast_matrix
  )
}

#' Apply contrast weights to fMRI model terms
#'
#' @description
#' Helper function to construct contrast matrices from contrast
#' specifications and term information.
#'
#' @param contrast_spec Contrast specification object
#' @param term_info Term information from fMRI model
#' @param beta_names Names of all parameters
#'
#' @return Numeric contrast matrix
#' @keywords internal
#' @noRd
construct_contrast_matrix <- function(contrast_spec, term_info, beta_names) {
  # This would interface with the existing contrast system
  # For now, return a placeholder
  
  if (is.matrix(contrast_spec)) {
    return(contrast_spec)
  }
  
  if (is.numeric(contrast_spec)) {
    return(matrix(contrast_spec, nrow = 1))
  }
  
  # Would handle contrast objects here
  stop("Complex contrast specifications not yet implemented")
}

#' Compute voxelwise contrasts efficiently
#'
#' @description
#' Optimized function for computing many contrasts across many voxels.
#' Uses matrix operations to avoid loops where possible.
#'
#' @param fit_result GLM fit results
#' @param contrast_list List of contrast matrices
#' @param parallel Whether to use parallel processing
#'
#' @return List of contrast results
#' @keywords internal
#' @noRd
compute_voxelwise_contrasts <- function(fit_result, contrast_list, parallel = FALSE) {
  if (parallel && getOption("fmrireg.num_threads", 1) > 1) {
    # Parallel version would go here
    message("Parallel contrast computation not yet implemented")
  }
  
  # Sequential version
  results <- vector("list", length(contrast_list))
  
  for (i in seq_along(contrast_list)) {
    contrast <- contrast_list[[i]]
    
    # Determine if F-test needed
    if (is.matrix(contrast) && nrow(contrast) > 1) {
      results[[i]] <- compute_f_statistic(fit_result, contrast)
    } else {
      results[[i]] <- compute_contrast(fit_result, contrast)
    }
  }
  
  results
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
      Xe <- X * e2[, v]
      meat <- meat + crossprod(Xe)
    }
    meat <- meat / ncol(residuals)
  } else {
    Xe <- X * as.vector(e2)
    meat <- crossprod(Xe)
  }
  
  # Sandwich
  sandwich_vcov <- adjustment * XtX_inv %*% meat %*% XtX_inv
  
  sandwich_vcov
}