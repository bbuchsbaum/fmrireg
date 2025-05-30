#' Bootstrap Methods for fMRI Linear Models
#'
#' Functions for bootstrap inference including confidence intervals
#' and hypothesis testing
#'
#' @keywords internal
#' @noRd
NULL

#' Bootstrap confidence intervals for GLM
#'
#' @description
#' Computes bootstrap confidence intervals using either
#' residual or case resampling. Handles temporal dependencies
#' via block bootstrap.
#'
#' @param fit_result Result from solve_integrated_glm
#' @param X Design matrix
#' @param Y Response matrix
#' @param config fmri_lm_config object
#' @param contrasts List of contrast matrices
#' @param nboot Number of bootstrap iterations
#' @param block_size Block size for temporal data
#' @param confidence Confidence level (default 0.95)
#' @param method Bootstrap method ("residual", "case", "wild")
#' @param parallel Use parallel processing
#' @param run_indices Run structure for multi-run data
#'
#' @return List with bootstrap distributions and confidence intervals
#' @keywords internal
#' @noRd
bootstrap_glm_inference <- function(fit_result, X, Y, config, 
                                   contrasts = NULL,
                                   nboot = 1000,
                                   block_size = NULL,
                                   confidence = 0.95,
                                   method = "residual",
                                   parallel = FALSE,
                                   run_indices = NULL) {
  
  n <- nrow(X)
  p <- ncol(X)
  nvoxels <- ncol(Y)
  
  # Determine block size if not provided
  if (is.null(block_size)) {
    block_size <- if (!is.null(run_indices)) {
      # Use run structure
      min(sapply(run_indices, length)) / 4
    } else {
      # Default based on sample size
      max(10, n / 20)
    }
    block_size <- round(block_size)
  }
  
  # Original estimates
  orig_betas <- fit_result$betas
  orig_contrasts <- if (!is.null(contrasts)) {
    lapply(contrasts, function(C) compute_contrast(fit_result, C))
  }
  
  # Bootstrap storage
  boot_betas <- array(NA, dim = c(nboot, p, nvoxels))
  boot_contrasts <- if (!is.null(contrasts)) {
    lapply(contrasts, function(C) {
      matrix(NA, nboot, nvoxels)
    })
  }
  
  # Fitted values and residuals
  fitted_vals <- X %*% orig_betas
  residuals <- Y - fitted_vals
  
  # Create blocks
  blocks <- create_bootstrap_blocks(n, block_size, run_indices)
  
  # Bootstrap iterations
  pb <- progress::progress_bar$new(total = nboot)
  
  for (b in 1:nboot) {
    pb$tick()
    
    # Generate bootstrap sample
    boot_data <- switch(method,
      "residual" = bootstrap_residual(X, fitted_vals, residuals, blocks),
      "case" = bootstrap_case(X, Y, blocks),
      "wild" = bootstrap_wild(X, fitted_vals, residuals),
      stop("Unknown bootstrap method: ", method)
    )
    
    # Fit bootstrap model
    boot_fit <- solve_integrated_glm(
      boot_data$X, 
      boot_data$Y, 
      config, 
      run_indices
    )
    
    # Store results
    boot_betas[b, , ] <- boot_fit$betas
    
    if (!is.null(contrasts)) {
      for (i in seq_along(contrasts)) {
        boot_con <- compute_contrast(boot_fit, contrasts[[i]])
        boot_contrasts[[i]][b, ] <- boot_con$estimate
      }
    }
  }
  
  # Compute confidence intervals
  alpha <- 1 - confidence
  probs <- c(alpha/2, 1 - alpha/2)
  
  # Beta CIs
  beta_ci <- apply(boot_betas, c(2, 3), quantile, probs = probs, na.rm = TRUE)
  beta_se_boot <- apply(boot_betas, c(2, 3), sd, na.rm = TRUE)
  
  # Contrast CIs
  contrast_ci <- if (!is.null(contrasts)) {
    lapply(boot_contrasts, function(bc) {
      t(apply(bc, 2, quantile, probs = probs, na.rm = TRUE))
    })
  }
  
  list(
    boot_betas = boot_betas,
    boot_contrasts = boot_contrasts,
    beta_ci = beta_ci,
    beta_se_boot = beta_se_boot,
    contrast_ci = contrast_ci,
    nboot = nboot,
    method = method,
    block_size = block_size,
    confidence = confidence
  )
}

#' Create bootstrap blocks
#'
#' @keywords internal
#' @noRd
create_bootstrap_blocks <- function(n, block_size, run_indices = NULL) {
  if (!is.null(run_indices)) {
    # Respect run structure
    blocks <- list()
    for (run in run_indices) {
      run_blocks <- split(
        run,
        ceiling(seq_along(run) / block_size)
      )
      blocks <- c(blocks, run_blocks)
    }
  } else {
    # Single run
    blocks <- split(
      1:n,
      ceiling(seq_along(1:n) / block_size)
    )
  }
  blocks
}

#' Residual bootstrap
#'
#' @keywords internal
#' @noRd
bootstrap_residual <- function(X, fitted, residuals, blocks) {
  n <- nrow(X)
  
  # Resample blocks of residuals
  sampled_blocks <- sample(length(blocks), replace = TRUE)
  resampled_idx <- unlist(blocks[sampled_blocks])
  
  # Adjust to exact size
  if (length(resampled_idx) > n) {
    resampled_idx <- resampled_idx[1:n]
  } else if (length(resampled_idx) < n) {
    extra <- sample(1:n, n - length(resampled_idx), replace = TRUE)
    resampled_idx <- c(resampled_idx, extra)
  }
  
  # New Y = fitted + resampled residuals
  Y_boot <- fitted + residuals[resampled_idx, , drop = FALSE]
  
  list(X = X, Y = Y_boot)
}

#' Case bootstrap
#'
#' @keywords internal
#' @noRd
bootstrap_case <- function(X, Y, blocks) {
  n <- nrow(X)
  
  # Resample blocks of cases
  sampled_blocks <- sample(length(blocks), replace = TRUE)
  resampled_idx <- unlist(blocks[sampled_blocks])
  
  # Adjust to exact size
  if (length(resampled_idx) > n) {
    resampled_idx <- resampled_idx[1:n]
  } else if (length(resampled_idx) < n) {
    extra <- sample(1:n, n - length(resampled_idx), replace = TRUE)
    resampled_idx <- c(resampled_idx, extra)
  }
  
  list(
    X = X[resampled_idx, , drop = FALSE],
    Y = Y[resampled_idx, , drop = FALSE]
  )
}

#' Wild bootstrap
#'
#' @keywords internal
#' @noRd
bootstrap_wild <- function(X, fitted, residuals) {
  n <- nrow(X)
  
  # Rademacher weights
  weights <- sample(c(-1, 1), n, replace = TRUE)
  
  # New Y = fitted + weighted residuals
  Y_boot <- fitted + residuals * weights
  
  list(X = X, Y = Y_boot)
}

#' Bootstrap hypothesis test
#'
#' @description
#' Performs bootstrap hypothesis test for contrasts using
#' the percentile method.
#'
#' @param boot_result Result from bootstrap_glm_inference
#' @param contrast_idx Index of contrast to test
#' @param null_value Null hypothesis value (default 0)
#'
#' @return P-value for two-sided test
#' @keywords internal
#' @noRd
bootstrap_hypothesis_test <- function(boot_result, contrast_idx, null_value = 0) {
  if (is.null(boot_result$boot_contrasts)) {
    stop("No contrasts in bootstrap results")
  }
  
  boot_dist <- boot_result$boot_contrasts[[contrast_idx]]
  
  # Compute p-value for each voxel
  nvoxels <- ncol(boot_dist)
  pvals <- numeric(nvoxels)
  
  for (v in 1:nvoxels) {
    # Two-sided test
    boot_centered <- boot_dist[, v] - mean(boot_dist[, v], na.rm = TRUE)
    observed <- mean(boot_dist[, v], na.rm = TRUE) - null_value
    
    pvals[v] <- mean(abs(boot_centered) >= abs(observed), na.rm = TRUE)
  }
  
  pvals
}

#' BCa confidence intervals
#'
#' @description
#' Computes bias-corrected and accelerated (BCa) bootstrap
#' confidence intervals for more accurate coverage.
#'
#' @param boot_dist Bootstrap distribution (nboot x nvoxels)
#' @param observed Observed estimates
#' @param confidence Confidence level
#'
#' @return Matrix of confidence intervals
#' @keywords internal
#' @noRd
compute_bca_ci <- function(boot_dist, observed, confidence = 0.95) {
  nboot <- nrow(boot_dist)
  nvoxels <- ncol(boot_dist)
  alpha <- 1 - confidence
  
  ci <- matrix(NA, 2, nvoxels)
  
  for (v in 1:nvoxels) {
    boot_v <- boot_dist[, v]
    obs_v <- observed[v]
    
    # Bias correction
    z0 <- qnorm(mean(boot_v <= obs_v, na.rm = TRUE))
    
    # Acceleration (would need jackknife for full implementation)
    a <- 0  # Simplified version
    
    # Adjusted percentiles
    z_alpha <- qnorm(c(alpha/2, 1 - alpha/2))
    p_adj <- pnorm(z0 + (z0 + z_alpha) / (1 - a * (z0 + z_alpha)))
    
    ci[, v] <- quantile(boot_v, probs = p_adj, na.rm = TRUE)
  }
  
  ci
}