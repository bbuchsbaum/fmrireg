#' Integrated GLM Solver with AR and Robust Options
#'
#' Main entry point for solving GLM with optional AR whitening and robust fitting
#'
#' @keywords internal
#' @noRd
NULL

#' Solve GLM with integrated AR and robust options
#'
#' @description
#' Coordinates the complete GLM solving pipeline including:
#' - Initial OLS fit
#' - Optional AR whitening (iterative)
#' - Optional robust fitting (IRLS)
#' - Combined AR+Robust ("whiten then robustly weight")
#'
#' @param X Design matrix
#' @param Y Response matrix
#' @param config fmri_lm_config object with all options
#' @param run_indices List of indices for each run
#'
#' @return List with full results including coefficients, residuals,
#'   standard errors, AR parameters, robust weights, etc.
#' @keywords internal
#' @noRd
solve_integrated_glm <- function(X, Y, config, run_indices = NULL) {
  stopifnot(inherits(config, "fmri_lm_config"))
  
  # Initial projection and context
  proj <- .fast_preproject(X)
  glm_ctx <- glm_context(X = X, Y = Y, proj = proj)
  
  # Determine solving strategy
  method <- config$method %||% "auto"
  if (method == "auto") {
    # Determine method based on options
    # Handle both config structures (ar_options vs ar)
    ar_struct <- config$ar$struct %||% 
                 (config$ar_options$cor_struct %||% "iid")
    robust_type <- config$robust$type %||% config$robust %||% FALSE
    
    if (!isFALSE(robust_type) && ar_struct != "iid" && ar_struct != "none") {
      method <- "ar_robust"
    } else if (!isFALSE(robust_type)) {
      method <- "robust"
    } else if (ar_struct != "iid" && ar_struct != "none") {
      method <- "ar"
    } else {
      method <- "ols"
    }
  }
  
  # Solve based on method
  # Extract ar_options in correct format
  ar_options <- if (!is.null(config$ar)) {
    list(
      cor_struct = config$ar$struct,
      iter = config$ar$iter_gls,
      exact_first = config$ar$exact_first,
      p = config$ar$p
    )
  } else {
    config$ar_options
  }
  
  result <- switch(method,
    "ols" = solve_ols_pipeline(glm_ctx),
    "ar" = solve_ar_pipeline(glm_ctx, ar_options, run_indices),
    "robust" = solve_robust_pipeline(glm_ctx, config$robust),
    "ar_robust" = solve_ar_robust_pipeline(glm_ctx, config, run_indices),
    stop("Unknown method: ", method)
  )
  
  # Add method info
  result$method <- method
  result$config <- config
  
  # Compute standard errors
  result <- compute_standard_errors(result, glm_ctx)
  
  result
}

#' OLS pipeline
#' @keywords internal
#' @noRd
solve_ols_pipeline <- function(glm_ctx) {
  result <- solve_glm_core(glm_ctx, return_fitted = TRUE)
  result$residuals <- glm_ctx$Y - result$fitted
  # Add XtXinv for contrast computation
  result$XtXinv <- glm_ctx$proj$XtXinv
  result
}

#' AR pipeline
#' @keywords internal
#' @noRd
solve_ar_pipeline <- function(glm_ctx, ar_options, run_indices) {
  result <- iterative_ar_solve(glm_ctx, ar_options, run_indices)
  
  # Compute effective df
  if (!is.null(result$ar_coef)) {
    n <- nrow(glm_ctx$X)
    p <- ncol(glm_ctx$X)
    phi_pooled <- mean(unlist(result$ar_coef))
    result$effective_df <- compute_ar_effective_df(n, p, phi_pooled)
  }
  
  # Add XtXinv if not present
  if (is.null(result$XtXinv)) {
    result$XtXinv <- glm_ctx$proj$XtXinv
  }
  
  result
}

#' Robust pipeline
#' @keywords internal
#' @noRd
solve_robust_pipeline <- function(glm_ctx, robust_options) {
  # Parse robust options
  if (is.character(robust_options)) {
    robust_config <- list(
      type = robust_options,
      k_huber = 1.345,
      c_tukey = 4.685,
      max_iter = 10,
      scale_scope = "voxel"
    )
  } else {
    robust_config <- robust_options
  }
  
  # Run robust fitting
  robust_result <- robust_iterative_fitter(
    initial_glm_ctx = glm_ctx,
    cfg_robust_options = robust_config,
    X_orig_for_resid = glm_ctx$X,
    sigma_fixed = NULL
  )
  
  # Format output
  list(
    betas = robust_result$betas_robust,
    sigma2 = robust_result$sigma_robust_scale_final^2,
    robust_weights = robust_result$robust_weights_final,
    XtXinv = robust_result$XtWXi_final,
    dfres = robust_result$dfres,
    fitted = glm_ctx$X %*% robust_result$betas_robust,
    residuals = glm_ctx$Y - glm_ctx$X %*% robust_result$betas_robust
  )
}

#' AR + Robust pipeline ("whiten then robustly weight")
#' @keywords internal
#' @noRd
solve_ar_robust_pipeline <- function(glm_ctx, config, run_indices) {
  # Extract ar_options in correct format
  ar_options <- if (!is.null(config$ar)) {
    list(
      cor_struct = config$ar$struct,
      iter = config$ar$iter_gls,
      exact_first = config$ar$exact_first,
      p = config$ar$p
    )
  } else {
    config$ar_options
  }
  
  # Step 1: Iterative AR whitening
  ar_result <- iterative_ar_solve(glm_ctx, ar_options, run_indices)
  
  # Step 2: Apply final AR whitening
  glm_ctx$residuals <- glm_ctx$Y - ar_result$fitted
  glm_ctx_white <- whiten_glm_context(glm_ctx, ar_options, run_indices)
  
  # Step 3: Robust fitting on whitened data
  robust_config <- if (is.character(config$robust)) {
    list(
      type = config$robust,
      k_huber = 1.345,
      c_tukey = 4.685,
      max_iter = 10,
      scale_scope = "voxel"
    )
  } else {
    config$robust
  }
  
  robust_result <- robust_iterative_fitter(
    initial_glm_ctx = glm_ctx_white,
    cfg_robust_options = robust_config,
    X_orig_for_resid = glm_ctx_white$X,
    sigma_fixed = NULL
  )
  
  # Combine results
  list(
    betas = robust_result$betas_robust,
    sigma2 = robust_result$sigma_robust_scale_final^2,
    robust_weights = robust_result$robust_weights_final,
    ar_coef = ar_result$ar_coef,
    XtXinv = robust_result$XtWXi_final,
    dfres = robust_result$dfres,
    fitted = glm_ctx$X %*% robust_result$betas_robust,
    residuals = glm_ctx$Y - glm_ctx$X %*% robust_result$betas_robust,
    effective_df = compute_ar_effective_df(
      nrow(glm_ctx$X), 
      ncol(glm_ctx$X),
      mean(unlist(ar_result$ar_coef))
    )
  )
}

#' Compute standard errors
#' @keywords internal
#' @noRd
compute_standard_errors <- function(result, glm_ctx) {
  # Extract variance-covariance matrix
  if (!is.null(result$XtXinv)) {
    vcov <- result$XtXinv
  } else if (!is.null(glm_ctx$proj$XtXinv)) {
    vcov <- glm_ctx$proj$XtXinv
  } else {
    stop("Cannot compute standard errors: no variance-covariance matrix")
  }
  
  # Use appropriate sigma2
  if (!is.null(result$sigma2)) {
    if (length(result$sigma2) == 1) {
      # Single sigma2 for all voxels
      result$standard_errors <- sqrt(diag(vcov) * result$sigma2)
    } else {
      # Vector of sigma2 (one per voxel)
      result$standard_errors <- sqrt(diag(vcov) %o% result$sigma2)
    }
  }
  
  result
}

#' Compute effective degrees of freedom for robust models
#'
#' @description
#' Adjusts degrees of freedom based on robust weights using
#' the trace of the hat matrix.
#'
#' @param X Design matrix
#' @param weights Robust weights
#' @param XtXinv Weighted (X'WX)^-1 matrix
#'
#' @return Effective degrees of freedom
#' @keywords internal
#' @noRd
compute_robust_effective_df <- function(X, weights, XtXinv) {
  # Hat matrix trace for weighted regression
  # tr(H) = tr(X(X'WX)^-1 X'W)
  W <- diag(weights)
  H_diag <- diag(X %*% XtXinv %*% t(X) %*% W)
  n <- nrow(X)
  n - sum(H_diag)
}