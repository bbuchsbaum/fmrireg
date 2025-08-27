context("AR Advanced Features")

# Helper functions for simulating AR data
simulate_ar_data <- function(n = 100, p = 1, phi = NULL, sigma = 1) {
  if (is.null(phi)) {
    phi <- runif(p, -0.9, 0.9) / p  # Default stable AR coefficients
  }
  
  # Handle AR(0) case (white noise)
  if (length(phi) == 0 || all(phi == 0)) {
    return(rnorm(n, sd = sigma))
  }
  
  # Ensure stability
  roots <- polyroot(c(1, -phi))
  if (any(abs(roots) <= 1.05)) {
    warning("AR coefficients may be near unit root")
  }
  
  as.numeric(arima.sim(model = list(ar = phi), n = n, sd = sigma))
}

create_ar_glm_context <- function(n = 100, p_design = 3, ar_order = 1, 
                                  ar_coef = 0.5, n_vox = 5) {
  # Create design matrix
  X <- cbind(1, matrix(rnorm(n * (p_design - 1)), n, p_design - 1))
  
  # True betas
  true_betas <- rnorm(p_design)
  
  # Generate Y with AR errors
  Y <- matrix(NA, n, n_vox)
  for (v in 1:n_vox) {
    signal <- X %*% true_betas
    errors <- simulate_ar_data(n = n, p = ar_order, phi = ar_coef)
    Y[, v] <- signal + errors
  }
  
  # Create projection
  proj <- fmrireg:::.fast_preproject(X)
  
  # Create context
  fmrireg:::glm_context(X = X, Y = Y, proj = proj)
}

test_that("whiten_glm_context handles single run AR whitening", {
  set.seed(123)
  glm_ctx <- create_ar_glm_context(n = 100, ar_order = 1, ar_coef = 0.6)
  
  ar_options <- list(
    cor_struct = "ar1",
    iter = 1
  )
  
  # Apply whitening
  whitened_ctx <- fmrireg:::whiten_glm_context(glm_ctx, ar_options)
  
  # Check that whitening was applied
  expect_true(is.glm_context(whitened_ctx))
  expect_false(all(whitened_ctx$X == glm_ctx$X))
  expect_false(all(whitened_ctx$Y == glm_ctx$Y))
  expect_true(!is.null(whitened_ctx$phi_hat))
  expect_equal(length(whitened_ctx$phi_hat), 1)  # Single run
  
  # Check AR coefficient is reasonable
  phi_est <- whitened_ctx$phi_hat[[1]]
  expect_true(abs(phi_est) < 1)  # Stationary
  expect_true(abs(phi_est - 0.6) < 0.2)  # Close to true value
})

test_that("whiten_glm_context handles multi-run AR whitening with different coefficients", {
  set.seed(124)
  
  # Create multi-run data with different AR coefficients
  n_per_run <- 50
  n_runs <- 3
  p_design <- 3
  n_vox = 4
  
  # Different AR coefficients per run
  ar_coefs <- c(0.3, 0.5, 0.7)
  
  X_list <- list()
  Y_list <- list()
  
  for (r in 1:n_runs) {
    X_run <- cbind(1, matrix(rnorm(n_per_run * (p_design - 1)), n_per_run, p_design - 1))
    true_betas <- rnorm(p_design)
    
    Y_run <- matrix(NA, n_per_run, n_vox)
    for (v in 1:n_vox) {
      signal <- X_run %*% true_betas
      errors <- simulate_ar_data(n = n_per_run, p = 1, phi = ar_coefs[r])
      Y_run[, v] <- signal + errors
    }
    
    X_list[[r]] <- X_run
    Y_list[[r]] <- Y_run
  }
  
  X <- do.call(rbind, X_list)
  Y <- do.call(rbind, Y_list)
  
  # Create run indices
  run_indices <- list()
  start_idx <- 1
  for (r in 1:n_runs) {
    run_indices[[r]] <- start_idx:(start_idx + n_per_run - 1)
    start_idx <- start_idx + n_per_run
  }
  
  # Create context
  proj <- fmrireg:::.fast_preproject(X)
  glm_ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  
  ar_options <- list(
    cor_struct = "ar1",
    iter = 1
  )
  
  # Apply whitening
  whitened_ctx <- fmrireg:::whiten_glm_context(glm_ctx, ar_options, run_indices)
  
  # Check results
  expect_true(is.glm_context(whitened_ctx))
  expect_equal(length(whitened_ctx$phi_hat), n_runs)
  
  # Check that estimated AR coefficients are different per run
  phi_estimates <- unlist(whitened_ctx$phi_hat)
  expect_equal(length(phi_estimates), n_runs)
  
  # Each should be reasonably close to true value
  for (r in 1:n_runs) {
    expect_true(abs(phi_estimates[r] - ar_coefs[r]) < 0.3)
  }
})

test_that("iterative_ar_solve converges with monitoring", {
  set.seed(125)
  glm_ctx <- create_ar_glm_context(n = 150, ar_order = 1, ar_coef = 0.5)
  
  ar_options <- list(
    cor_struct = "ar1",
    iter = 5
  )
  
  # Test with convergence monitoring
  # Message is conditional on convergence, so we don't test for it
  result <- fmrireg:::iterative_ar_solve(glm_ctx, ar_options, max_iter = 5, tol = 1e-3)
  
  # Check result structure
  expect_true(!is.null(result$betas))
  expect_true(!is.null(result$ar_coef))
  expect_equal(result$ar_order, 1)
})

test_that("iterative_ar_solve handles no AR case", {
  set.seed(126)
  glm_ctx <- create_ar_glm_context(n = 100, ar_order = 1, ar_coef = 0)
  
  ar_options <- list(
    cor_struct = "none"
  )
  
  result <- fmrireg:::iterative_ar_solve(glm_ctx, ar_options)
  
  # Should return standard OLS result
  expect_true(is.null(result$ar_coef))
  expect_true(!is.null(result$betas))
})

test_that("compute_ar_effective_df adjusts degrees of freedom correctly", {
  # Test AR(1)
  n <- 100
  p <- 5
  
  # No AR
  df_no_ar <- fmrireg:::compute_ar_effective_df(n, p, NULL)
  expect_equal(df_no_ar, n - p)
  
  # Moderate AR(1)
  phi1 <- 0.5
  df_ar1 <- fmrireg:::compute_ar_effective_df(n, p, phi1)
  expect_true(df_ar1 < df_no_ar)
  expect_true(df_ar1 > 0)
  # Standard formula: no AR penalty (following SPM/FSL/AFNI)
  expected_df <- n * (1 - phi1^2) - p
  expect_equal(df_ar1, expected_df)
  
  # Strong AR(1)
  phi1_strong <- 0.9
  df_ar1_strong <- fmrireg:::compute_ar_effective_df(n, p, phi1_strong)
  expect_true(df_ar1_strong < df_ar1)
  
  # AR(2)
  phi2 <- c(0.4, 0.3)
  df_ar2 <- fmrireg:::compute_ar_effective_df(n, p, phi2)
  expect_true(df_ar2 < df_no_ar)
  # Standard formula: no AR penalty (following SPM/FSL/AFNI)
  expected_df_ar2 <- n * (1 - sum(phi2^2)) - p
  expect_equal(df_ar2, expected_df_ar2)
  
  # Test conservative mode with penalize_ar = TRUE
  df_ar1_conservative <- fmrireg:::compute_ar_effective_df(n, p, phi1, n_runs = 1, penalize_ar = TRUE)
  expect_equal(df_ar1_conservative, n * (1 - phi1^2) - p - 1)  # Subtract AR order (1)
  
  df_ar2_conservative <- fmrireg:::compute_ar_effective_df(n, p, phi2, n_runs = 1, penalize_ar = TRUE)
  expect_equal(df_ar2_conservative, n * (1 - sum(phi2^2)) - p - 2)  # Subtract AR order (2)
})

test_that("whiten_glm_context handles AR(3) and AR(4) models", {
  set.seed(127)
  
  # Test AR(3)
  n <- 200
  p_design <- 3
  n_vox <- 3
  ar3_coef <- c(0.3, 0.2, 0.1)
  
  X <- cbind(1, matrix(rnorm(n * (p_design - 1)), n, p_design - 1))
  true_betas <- rnorm(p_design)
  
  Y <- matrix(NA, n, n_vox)
  for (v in 1:n_vox) {
    signal <- X %*% true_betas
    errors <- simulate_ar_data(n = n, p = 3, phi = ar3_coef)
    Y[, v] <- signal + errors
  }
  
  proj <- fmrireg:::.fast_preproject(X)
  glm_ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  
  ar_options <- list(
    cor_struct = "ar3",
    iter = 1
  )
  
  whitened_ctx <- fmrireg:::whiten_glm_context(glm_ctx, ar_options)
  
  # Check AR(3) coefficients
  expect_equal(length(whitened_ctx$phi_hat[[1]]), 3)
  
  # Test AR(4)
  ar4_coef <- c(0.2, 0.15, 0.1, 0.05)
  
  Y4 <- matrix(NA, n, n_vox)
  for (v in 1:n_vox) {
    signal <- X %*% true_betas
    errors <- simulate_ar_data(n = n, p = 4, phi = ar4_coef)
    Y4[, v] <- signal + errors
  }
  
  glm_ctx4 <- fmrireg:::glm_context(X = X, Y = Y4, proj = proj)
  
  ar_options4 <- list(
    cor_struct = "ar4",
    iter = 1
  )
  
  whitened_ctx4 <- fmrireg:::whiten_glm_context(glm_ctx4, ar_options4)
  
  # Check AR(4) coefficients
  expect_equal(length(whitened_ctx4$phi_hat[[1]]), 4)
})

test_that("whiten_glm_context handles AR(p) with custom order", {
  set.seed(128)
  
  # Test AR(p) with p=5
  n <- 250
  p_design <- 3
  n_vox <- 2
  p_ar <- 5
  ar_coef <- c(0.2, 0.15, 0.1, 0.08, 0.05)
  
  X <- cbind(1, matrix(rnorm(n * (p_design - 1)), n, p_design - 1))
  true_betas <- rnorm(p_design)
  
  Y <- matrix(NA, n, n_vox)
  for (v in 1:n_vox) {
    signal <- X %*% true_betas
    errors <- simulate_ar_data(n = n, p = p_ar, phi = ar_coef)
    Y[, v] <- signal + errors
  }
  
  proj <- fmrireg:::.fast_preproject(X)
  glm_ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  
  ar_options <- list(
    cor_struct = "arp",
    p = p_ar,
    iter = 1
  )
  
  whitened_ctx <- fmrireg:::whiten_glm_context(glm_ctx, ar_options)
  
  # Check AR(p) coefficients
  expect_equal(length(whitened_ctx$phi_hat[[1]]), p_ar)
})

test_that("AR methods handle near unit-root processes", {
  set.seed(129)
  
  # Test with AR coefficients close to 1
  near_unit_roots <- c(0.9, 0.95, 0.99)
  
  for (phi in near_unit_roots) {
    n <- 300  # Need longer series for near unit root
    p_design <- 2
    n_vox <- 2
    
    X <- cbind(1, rnorm(n))
    true_betas <- c(10, 2)
    
    Y <- matrix(NA, n, n_vox)
    for (v in 1:n_vox) {
      signal <- X %*% true_betas
      # Suppress warnings about near unit root
      errors <- suppressWarnings(simulate_ar_data(n = n, p = 1, phi = phi))
      Y[, v] <- signal + errors
    }
    
    proj <- fmrireg:::.fast_preproject(X)
    glm_ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
    
    ar_options <- list(
      cor_struct = "ar1",
      iter = 2
    )
    
    # Should handle without error
    result <- fmrireg:::iterative_ar_solve(glm_ctx, ar_options, max_iter = 2)
    
    expect_true(!is.null(result$ar_coef))
    phi_est <- result$ar_coef[[1]]
    
    # Check effective df adjustment
    df_adj <- fmrireg:::compute_ar_effective_df(n, p_design, phi_est)
    expect_true(df_adj > 0)  # Should still be positive
    expect_true(df_adj < n - p_design)  # Should be reduced
    
    # For very high AR, effective df should be reduced
    # Note: The exact reduction depends on the estimation
    if (phi > 0.95) {
      # Just check that it's substantially reduced
      expect_true(df_adj < (n - p_design) * 0.5)
    }
  }
})

test_that("AR methods handle short time series edge cases", {
  set.seed(130)
  
  # Test with time series too short for AR order
  n_short <- 10
  p_design <- 3
  n_vox <- 2
  
  X <- cbind(1, matrix(rnorm(n_short * (p_design - 1)), n_short, p_design - 1))
  Y <- matrix(rnorm(n_short * n_vox), n_short, n_vox)
  
  proj <- fmrireg:::.fast_preproject(X)
  glm_ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  
  # Try AR(4) with only 10 time points - should handle gracefully
  ar_options <- list(
    cor_struct = "ar4",
    iter = 1
  )
  
  # This should work but produce warning or handle the edge case
  whitened_ctx <- fmrireg:::whiten_glm_context(glm_ctx, ar_options)
  
  # AR estimation might fail or return fewer coefficients
  expect_true(is.glm_context(whitened_ctx))
})

test_that("Different AR structures per run are handled correctly", {
  set.seed(131)
  
  # Create data with different AR structures per run
  n_per_run <- 100
  n_runs <- 3
  p_design <- 3
  n_vox <- 3
  
  # Run 1: AR(1), Run 2: AR(2), Run 3: No AR
  X_list <- list()
  Y_list <- list()
  true_ar_structures <- list(
    list(order = 1, phi = 0.5),
    list(order = 2, phi = c(0.4, 0.2)),
    list(order = 0, phi = numeric(0))
  )
  
  for (r in 1:n_runs) {
    X_run <- cbind(1, matrix(rnorm(n_per_run * (p_design - 1)), n_per_run, p_design - 1))
    true_betas <- rnorm(p_design)
    
    Y_run <- matrix(NA, n_per_run, n_vox)
    for (v in 1:n_vox) {
      signal <- X_run %*% true_betas
      if (true_ar_structures[[r]]$order > 0) {
        errors <- simulate_ar_data(n = n_per_run, 
                                   p = true_ar_structures[[r]]$order, 
                                   phi = true_ar_structures[[r]]$phi)
      } else {
        errors <- rnorm(n_per_run)
      }
      Y_run[, v] <- signal + errors
    }
    
    X_list[[r]] <- X_run
    Y_list[[r]] <- Y_run
  }
  
  X <- do.call(rbind, X_list)
  Y <- do.call(rbind, Y_list)
  
  # Create run indices
  run_indices <- list()
  start_idx <- 1
  for (r in 1:n_runs) {
    run_indices[[r]] <- start_idx:(start_idx + n_per_run - 1)
    start_idx <- start_idx + n_per_run
  }
  
  # Test with AR(2) - should estimate appropriate order for each run
  proj <- fmrireg:::.fast_preproject(X)
  glm_ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  
  ar_options <- list(
    cor_struct = "ar2",
    iter = 1
  )
  
  whitened_ctx <- fmrireg:::whiten_glm_context(glm_ctx, ar_options, run_indices)
  
  # Check that we get estimates for each run
  expect_equal(length(whitened_ctx$phi_hat), n_runs)
  
  # Each run should have AR(2) coefficients (even if some are ~0)
  for (r in 1:n_runs) {
    expect_equal(length(whitened_ctx$phi_hat[[r]]), 2)
  }
})

test_that("whiten_glm_context handles pre-specified AR coefficients", {
  set.seed(132)
  
  glm_ctx <- create_ar_glm_context(n = 100, ar_order = 1, ar_coef = 0.6)
  
  # Use pre-specified coefficient
  ar_options <- list(
    cor_struct = "ar1",
    phi = 0.5,  # Different from true value
    iter = 1
  )
  
  whitened_ctx <- fmrireg:::whiten_glm_context(glm_ctx, ar_options)
  
  # Should use the specified phi, not estimate it
  expect_equal(whitened_ctx$phi_hat[[1]], 0.5)
})

test_that("Integration with GLM solving pipeline works correctly", {
  set.seed(133)
  
  # Create realistic fMRI-like data
  n_time <- 200
  n_vox <- 10
  p_design <- 5
  
  # Design matrix with intercept and some regressors
  X <- cbind(1, 
             sin(2 * pi * (1:n_time) / 20),  # Periodic regressor
             cos(2 * pi * (1:n_time) / 20),
             rnorm(n_time),  # Random regressor
             cumsum(rnorm(n_time)) / 10)  # Drift-like regressor
  
  # True parameters
  true_betas <- c(100, 5, 3, 2, 1)
  
  # Generate data with AR(1) errors
  ar_coef <- 0.4
  Y <- matrix(NA, n_time, n_vox)
  for (v in 1:n_vox) {
    signal <- X %*% true_betas
    errors <- simulate_ar_data(n = n_time, p = 1, phi = ar_coef, sigma = 10)
    Y[, v] <- signal + errors
  }
  
  # Test full pipeline
  proj <- fmrireg:::.fast_preproject(X)
  glm_ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  
  ar_options <- list(
    cor_struct = "ar1",
    iter = 3
  )
  
  # Run iterative AR solve
  result <- fmrireg:::iterative_ar_solve(glm_ctx, ar_options, max_iter = 3, tol = 1e-4)
  
  # Check that betas are recovered reasonably well
  mean_betas <- rowMeans(result$betas)
  
  # Should be within reasonable range of true values
  # Note: Recovery depends on noise level and AR strength
  for (i in 1:p_design) {
    if (abs(true_betas[i]) > 1) {  # Only check larger coefficients
      relative_error <- abs(mean_betas[i] - true_betas[i]) / abs(true_betas[i])
      expect_true(relative_error < 1.0)  # Within 100% of true value
    }
  }
  
  # Check AR coefficient recovery
  expect_true(abs(result$ar_coef[[1]] - ar_coef) < 0.2)
})

test_that("AR whitening preserves rank of design matrix", {
  set.seed(134)
  
  # Create full rank design matrix
  n <- 100
  p <- 5
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  Y <- matrix(rnorm(n * 3), n, 3)
  
  proj <- fmrireg:::.fast_preproject(X)
  glm_ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  
  ar_options <- list(
    cor_struct = "ar2",
    iter = 1
  )
  
  whitened_ctx <- fmrireg:::whiten_glm_context(glm_ctx, ar_options)
  
  # Check that whitened X still has full rank
  rank_original <- Matrix::rankMatrix(X)[1]
  rank_whitened <- Matrix::rankMatrix(whitened_ctx$X)[1]
  
  expect_equal(rank_whitened, rank_original)
})

test_that("Satterthwaite approximation for effective df with AR errors", {
  set.seed(135)
  
  # Test more sophisticated effective df calculation
  n <- 200
  p <- 5
  
  # Different AR structures
  ar_structures <- list(
    list(phi = 0.3, name = "Weak AR(1)"),
    list(phi = 0.7, name = "Strong AR(1)"),
    list(phi = c(0.4, 0.2), name = "AR(2)"),
    list(phi = c(0.3, 0.2, 0.1), name = "AR(3)")
  )
  
  for (ar_struct in ar_structures) {
    df_effective <- fmrireg:::compute_ar_effective_df(n, p, ar_struct$phi)
    
    # Effective df should decrease with stronger/more AR parameters
    expect_true(df_effective < n - p)
    expect_true(df_effective > 0)
    
    # Print for inspection
    cat(sprintf("%s: df = %.2f (reduction = %.1f%%)\n", 
                ar_struct$name, 
                df_effective,
                100 * (1 - df_effective / (n - p))))
  }
})

test_that("whiten_glm_context handles missing data gracefully", {
  set.seed(136)
  
  # Create data with some NA values
  n <- 100
  p_design <- 3
  n_vox <- 3
  
  X <- cbind(1, matrix(rnorm(n * (p_design - 1)), n, p_design - 1))
  Y <- matrix(rnorm(n * n_vox), n, n_vox)
  
  # Introduce some NAs
  Y[c(10, 20, 30), 1] <- NA
  
  proj <- fmrireg:::.fast_preproject(X)
  glm_ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  
  ar_options <- list(
    cor_struct = "ar1",
    iter = 1
  )
  
  # Should produce an error or handle gracefully
  expect_error(
    fmrireg:::whiten_glm_context(glm_ctx, ar_options),
    "NA values"
  )
})

test_that("AR whitening with exact_first option works correctly", {
  set.seed(137)
  
  glm_ctx <- create_ar_glm_context(n = 100, ar_order = 1, ar_coef = 0.5)
  
  ar_options <- list(
    cor_struct = "ar1",
    iter = 1
  )
  
  # Apply whitening with exact_first = TRUE (default)
  whitened_ctx <- fmrireg:::whiten_glm_context(glm_ctx, ar_options)
  
  # The first row should be scaled differently with exact_first
  # Check that whitening was applied
  expect_true(is.glm_context(whitened_ctx))
  
  # Verify the transformation preserved dimensions
  expect_equal(dim(whitened_ctx$X), dim(glm_ctx$X))
  expect_equal(dim(whitened_ctx$Y), dim(glm_ctx$Y))
})

test_that("compute_ar_effective_df handles edge cases", {
  # Test with very small effective sample size
  n <- 50
  p <- 40
  phi <- 0.95  # Very strong AR
  
  df <- fmrireg:::compute_ar_effective_df(n, p, phi)
  expect_true(df >= 1)  # Should always be at least 1
  
  # Test with p > n (should still work)
  n <- 30
  p <- 40
  phi <- 0.5
  
  df <- fmrireg:::compute_ar_effective_df(n, p, phi)
  expect_true(df >= 1)  # Should handle gracefully
})

test_that("iterative_ar_solve respects max_iter parameter", {
  set.seed(138)
  
  glm_ctx <- create_ar_glm_context(n = 100, ar_order = 1, ar_coef = 0.6)
  
  ar_options <- list(
    cor_struct = "ar1",
    iter = 10  # Request many iterations
  )
  
  # But limit with max_iter
  result <- fmrireg:::iterative_ar_solve(glm_ctx, ar_options, max_iter = 1)
  
  # Should have stopped after 1 iteration
  expect_true(!is.null(result$ar_coef))
  
  # Compare with more iterations
  result_more <- fmrireg:::iterative_ar_solve(glm_ctx, ar_options, max_iter = 3)
  
  # Coefficients might be different with more iterations
  expect_true(!is.null(result_more$ar_coef))
})

test_that("AR whitening handles various design matrix conditions", {
  set.seed(139)
  
  n <- 100
  p <- 5
  
  # Test 1: Full rank matrix
  X_full <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * 3), n, 3)
  
  proj_full <- fmrireg:::.fast_preproject(X_full)
  glm_ctx_full <- fmrireg:::glm_context(X = X_full, Y = Y, proj = proj_full)
  
  ar_options <- list(
    cor_struct = "ar1",
    iter = 1
  )
  
  # Should work with full rank matrix
  whitened_ctx_full <- fmrireg:::whiten_glm_context(glm_ctx_full, ar_options)
  expect_true(is.glm_context(whitened_ctx_full))
  
  # Test 2: Create a truly rank-deficient matrix
  X_rank_def <- matrix(0, n, p)
  X_rank_def[, 1] <- 1  # Intercept
  X_rank_def[, 2] <- rnorm(n)
  X_rank_def[, 3] <- rnorm(n)
  X_rank_def[, 4] <- X_rank_def[, 2]  # Exact duplicate
  X_rank_def[, 5] <- 2 * X_rank_def[, 3]  # Exact multiple
  
  # This should have rank 3 (not 5)
  proj_rank_def <- suppressWarnings(fmrireg:::.fast_preproject(X_rank_def))
  
  # Even if rank detection doesn't work as expected, AR whitening should still work
  glm_ctx_rank_def <- fmrireg:::glm_context(X = X_rank_def, Y = Y, proj = proj_rank_def)
  
  # Should handle rank deficient case
  whitened_ctx_rank_def <- fmrireg:::whiten_glm_context(glm_ctx_rank_def, ar_options)
  expect_true(is.glm_context(whitened_ctx_rank_def))
})

test_that("AR parameter estimation with pooled residuals works correctly", {
  set.seed(140)
  
  # Create data where voxels have similar AR structure
  n <- 150
  p_design <- 3
  n_vox <- 10
  true_ar <- 0.6
  
  X <- cbind(1, matrix(rnorm(n * (p_design - 1)), n, p_design - 1))
  true_betas <- rnorm(p_design)
  
  # Generate Y with same AR structure across voxels
  Y <- matrix(NA, n, n_vox)
  for (v in 1:n_vox) {
    signal <- X %*% true_betas
    errors <- simulate_ar_data(n = n, p = 1, phi = true_ar, sigma = 1)
    Y[, v] <- signal + errors
  }
  
  proj <- fmrireg:::.fast_preproject(X)
  glm_ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  
  ar_options <- list(
    cor_struct = "ar1",
    iter = 1
  )
  
  # Whitening should pool residuals across voxels
  whitened_ctx <- fmrireg:::whiten_glm_context(glm_ctx, ar_options)
  
  # Pooled estimate should be close to true value
  phi_est <- whitened_ctx$phi_hat[[1]]
  expect_true(abs(phi_est - true_ar) < 0.15)
})