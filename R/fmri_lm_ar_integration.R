#' AR Integration for GLM
#'
#' Functions to integrate AR modeling with the GLM solver pipeline
#'
#' @keywords internal
#' @noRd
NULL

#' Apply AR whitening to GLM context
#'
#' @description
#' Transforms the design matrix X and response Y in a GLM context
#' using estimated or provided AR coefficients. Handles multi-run
#' data by applying whitening separately to each run.
#'
#' @param glm_ctx A glm_context object
#' @param ar_options List with AR configuration:
#'   - cor_struct: "ar1", "ar2", etc. or "none"
#'   - iter: Number of AR estimation iterations
#'   - phi: Optional pre-specified AR coefficients
#' @param run_indices List of indices for each run
#'
#' @return Updated glm_context with whitened matrices and AR info
#' @keywords internal
#' @noRd
whiten_glm_context <- function(glm_ctx, ar_options, run_indices = NULL) {
  stopifnot(is.glm_context(glm_ctx))
  
  if (is.null(ar_options) || ar_options$cor_struct == "none") {
    return(glm_ctx)
  }
  
  X <- glm_ctx$X
  Y <- glm_ctx$Y
  
  # Determine AR order
  ar_order <- switch(ar_options$cor_struct,
    "ar1" = 1,
    "ar2" = 2,
    "ar3" = 3,
    "ar4" = 4,
    "arp" = ar_options$p %||% stop("p must be specified for cor_struct='arp'"),
    stop("Unknown AR structure: ", ar_options$cor_struct)
  )
  
  # If no run indices, treat as single run
  if (is.null(run_indices)) {
    run_indices <- list(1:nrow(X))
  }
  
  # Store estimated phi for each run
  phi_list <- vector("list", length(run_indices))
  
  # Apply whitening by run
  X_whitened <- X * 0  # Initialize
  Y_whitened <- Y * 0
  
  for (i in seq_along(run_indices)) {
    idx <- run_indices[[i]]
    
    # Get phi for this run
    if (!is.null(ar_options$phi)) {
      # Use provided phi
      phi <- ar_options$phi
    } else {
      # Estimate from residuals
      if (is.null(glm_ctx$residuals)) {
        # Need initial OLS residuals
        proj_temp <- .fast_preproject(X[idx, , drop = FALSE])
        resid_temp <- Y[idx, , drop = FALSE] - 
                      X[idx, , drop = FALSE] %*% 
                      (proj_temp$Pinv %*% Y[idx, , drop = FALSE])
        
        # Pool residuals across voxels for stable estimation
        pooled_resid <- as.vector(resid_temp)
        phi <- estimate_ar_parameters(pooled_resid, ar_order)
      } else {
        # Use provided residuals
        pooled_resid <- as.vector(glm_ctx$residuals[idx, ])
        phi <- estimate_ar_parameters(pooled_resid, ar_order)
      }
    }
    
    phi_list[[i]] <- phi
    
    # Apply whitening
    whitened <- ar_whiten_transform(
      X[idx, , drop = FALSE],
      Y[idx, , drop = FALSE],
      phi,
      exact_first = TRUE
    )
    
    X_whitened[idx, ] <- whitened$X
    Y_whitened[idx, ] <- whitened$Y
  }
  
  # Update projection for whitened X
  proj_new <- .fast_preproject(X_whitened)
  
  # Create new context with whitened data
  glm_context(
    X = X_whitened,
    Y = Y_whitened,
    proj = proj_new,
    phi_hat = phi_list,
    sigma_robust_scale = glm_ctx$sigma_robust_scale,
    robust_weights = glm_ctx$robust_weights
  )
}

#' Iterative AR estimation and whitening
#'
#' @description
#' Performs iterative AR parameter estimation and whitening.
#' Alternates between estimating AR parameters from residuals
#' and re-fitting the whitened model.
#'
#' @param glm_ctx Initial GLM context
#' @param ar_options AR configuration options
#' @param run_indices Run structure
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance for AR parameters
#'
#' @return List with final results including AR coefficients
#' @keywords internal
#' @noRd
iterative_ar_solve <- function(glm_ctx, ar_options, run_indices = NULL, 
                               max_iter = NULL, tol = 1e-4) {
  
  if (is.null(max_iter)) {
    max_iter <- ar_options$iter %||% 2
  }
  
  # Initial OLS fit
  result <- solve_glm_core(glm_ctx, return_fitted = TRUE)
  
  if (ar_options$cor_struct == "none" || max_iter == 0) {
    return(result)
  }
  
  # Store AR coefficients
  phi_prev <- NULL
  
  for (iter in 1:max_iter) {
    # Update context with residuals
    glm_ctx$residuals <- glm_ctx$Y - result$fitted
    
    # Apply AR whitening
    glm_ctx_white <- whiten_glm_context(glm_ctx, ar_options, run_indices)
    
    # Check convergence
    if (!is.null(phi_prev) && !is.null(glm_ctx_white$phi_hat)) {
      phi_current <- unlist(glm_ctx_white$phi_hat)
      if (all(abs(phi_current - phi_prev) < tol)) {
        message("AR parameters converged at iteration ", iter)
        break
      }
      phi_prev <- phi_current
    } else if (!is.null(glm_ctx_white$phi_hat)) {
      phi_prev <- unlist(glm_ctx_white$phi_hat)
    }
    
    # Solve whitened system
    result <- solve_glm_core(glm_ctx_white, return_fitted = TRUE)
    
    # Add AR info to result
    result$ar_coef <- glm_ctx_white$phi_hat
    result$ar_order <- switch(ar_options$cor_struct,
      "ar1" = 1, "ar2" = 2, "ar3" = 3, "ar4" = 4,
      "arp" = ar_options$p
    )
  }
  
  result
}

#' Compute effective degrees of freedom for AR models
#'
#' @description
#' Adjusts degrees of freedom to account for autocorrelation
#' in the residuals. Uses Satterthwaite-type approximation.
#'
#' @param n Sample size
#' @param p Number of parameters
#' @param phi AR coefficients
#'
#' @return Effective degrees of freedom
#' @keywords internal
#' @noRd
compute_ar_effective_df <- function(n, p, phi) {
  if (is.null(phi) || length(phi) == 0) {
    return(n - p)
  }
  
  # For AR(1), effective sample size is approximately n * (1 - phi^2)
  # For higher orders, use sum of squared AR coefficients
  ar_factor <- 1 - sum(phi^2)
  ar_factor <- max(ar_factor, 0.1)  # Prevent too small values
  
  effective_n <- n * ar_factor
  max(effective_n - p, 1)
}