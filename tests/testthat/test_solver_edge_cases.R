# Test edge cases for GLM solver

library(fmrireg)
library(testthat)

test_that("solve_glm_core handles rank-deficient matrices", {
  # Create rank-deficient design matrix
  n <- 50
  X <- cbind(
    1,
    rnorm(n),
    rnorm(n),
    1  # Duplicate of intercept
  )
  Y <- matrix(rnorm(n * 3), n, 3)
  
  # Check that matrix is rank deficient
  expect_lt(qr(X)$rank, ncol(X))
  
  # Should handle rank deficiency gracefully
  # The warning may or may not appear depending on implementation details
  proj <- suppressWarnings(fmrireg:::.fast_preproject(X))
  ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  
  result <- fmrireg:::solve_glm_core(ctx)
  
  # Should still return results
  expect_equal(dim(result$betas), c(ncol(X), ncol(Y)))
  # Should produce finite results even with rank deficiency
  expect_true(all(is.finite(result$betas)))
  expect_true(all(!is.na(result$betas)))
})

test_that("solve_glm_core handles near-singular matrices", {
  # Create nearly collinear predictors
  n <- 100
  x1 <- rnorm(n)
  x2 <- x1 + rnorm(n, sd = 0.01)  # Almost identical to x1
  X <- cbind(1, x1, x2)
  Y <- matrix(rnorm(n * 2), n, 2)
  
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  
  # Should complete without error
  result <- fmrireg:::solve_glm_core(ctx)
  
  expect_equal(dim(result$betas), c(3, 2))
  expect_true(all(!is.na(result$sigma2)))
  
  # Check condition number
  expect_gt(kappa(X), 100)  # High condition number indicates near-singularity
})

test_that("solve_glm_core handles extreme values", {
  n <- 50
  p <- 3
  
  # Test with very large values
  X_large <- cbind(1, matrix(rnorm(n * (p-1)) * 1e6, n, p-1))
  Y_large <- matrix(rnorm(n * 2) * 1e6, n, 2)
  
  proj_large <- fmrireg:::.fast_preproject(X_large)
  ctx_large <- fmrireg:::glm_context(X = X_large, Y = Y_large, proj = proj_large)
  result_large <- fmrireg:::solve_glm_core(ctx_large)
  
  expect_false(any(is.infinite(result_large$betas)))
  expect_false(any(is.nan(result_large$betas)))
  
  # Test with very small values
  X_small <- cbind(1, matrix(rnorm(n * (p-1)) * 1e-6, n, p-1))
  Y_small <- matrix(rnorm(n * 2) * 1e-6, n, 2)
  
  proj_small <- fmrireg:::.fast_preproject(X_small)
  ctx_small <- fmrireg:::glm_context(X = X_small, Y = Y_small, proj = proj_small)
  result_small <- fmrireg:::solve_glm_core(ctx_small)
  
  expect_false(any(is.infinite(result_small$betas)))
  expect_false(any(is.nan(result_small$betas)))
})

test_that("robust fitting handles extreme outliers", {
  n <- 100
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  
  # Create data with extreme outliers
  Y <- X %*% c(2, 1, -0.5) + rnorm(n, sd = 0.1)
  # Add extreme outliers
  outlier_idx <- sample(n, 10)
  Y[outlier_idx] <- Y[outlier_idx] + sample(c(-20, 20), 10, replace = TRUE)
  
  # First create initial OLS fit
  proj <- fmrireg:::.fast_preproject(X)
  initial_ctx <- fmrireg:::glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  initial_fit <- fmrireg:::solve_glm_core(initial_ctx)
  initial_ctx$betas <- initial_fit$betas
  
  cfg <- fmrireg:::fmri_lm_control(
    robust_options = list(
      type = "bisquare",
      c_tukey = 4.685,
      max_iter = 10
    )
  )
  
  # Robust fitting should downweight outliers
  result <- fmrireg:::robust_iterative_fitter(initial_ctx, cfg$robust, X)
  
  expect_true(all(result$robust_weights_final[outlier_idx] < 0.5))  # Outliers downweighted
  expect_true(mean(result$robust_weights_final[-outlier_idx]) > 0.8)  # Non-outliers have high weight
  
  # Estimates should be close to true values despite outliers
  expect_equal(as.vector(result$betas_robust), c(2, 1, -0.5), tolerance = 0.2)
})

test_that("robust fitting handles convergence failures", {
  # Create pathological case
  n <- 20
  X <- cbind(1, rnorm(n))
  Y <- matrix(c(rep(1, 10), rep(1000, 10)), ncol = 1)  # Bimodal response
  
  # First create initial OLS fit
  proj <- fmrireg:::.fast_preproject(X)
  initial_ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  initial_fit <- fmrireg:::solve_glm_core(initial_ctx)
  initial_ctx$betas <- initial_fit$betas
  
  cfg <- fmrireg:::fmri_lm_control(
    robust_options = list(
      type = "bisquare",
      max_iter = 2  # Force early stopping
    )
  )
  
  # Should complete even with limited iterations
  result <- fmrireg:::robust_iterative_fitter(initial_ctx, cfg$robust, X)
  
  expect_true(!is.null(result$betas_robust))
  expect_true(!is.null(result$robust_weights_final))
  expect_equal(length(result$robust_weights_final), n)
})

test_that("AR fitting handles short time series", {
  # Time series too short for AR order
  n <- 5  # Very short
  X <- cbind(1, rnorm(n))
  Y <- matrix(rnorm(n), ncol = 1)
  
  # Get residuals for AR estimation
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  fit <- fmrireg:::solve_glm_core(ctx)
  residuals <- Y - X %*% fit$betas
  
  # AR(2) with only 5 observations - should work but may not be reliable
  phi <- fmrireg:::.estimate_ar(as.vector(residuals), p_order = 2)
  
  # Should return coefficients
  expect_equal(length(phi), 2)
  # Check stationarity (all roots outside unit circle)
  if (length(phi) == 2) {
    # For AR(2), check characteristic polynomial roots
    # 1 - phi1*z - phi2*z^2 = 0
    # Stationarity requires |roots| > 1
    poly_coef <- c(1, -phi[1], -phi[2])
    roots <- polyroot(poly_coef)
    expect_true(all(Mod(roots) > 0.99))  # Allow slight numerical tolerance
  }
})

test_that("mixed_solve handles edge cases", {
  # Test with proper mixed model structure
  n <- 30
  # Fixed effects design matrix
  X <- cbind(1, rnorm(n))
  # Response
  y <- rnorm(n)
  # Random effects design matrix  
  Z <- diag(n)
  # Kinship/correlation matrix for random effects
  K <- diag(n)
  
  # Normal case should work
  result <- mixed_solve_cpp(y = y, Z = Z, K = K, X = X)
  expect_true(!is.null(result$beta))
  expect_equal(length(result$beta), ncol(X))
  expect_true(!is.null(result$Vu))
  expect_true(!is.null(result$Ve))
  
  # Edge case: Very small variance components
  y_small <- y * 1e-10
  result_small <- mixed_solve_cpp(y = y_small, Z = Z, K = K, X = X)
  expect_true(all(is.finite(result_small$beta)))
  
  # Edge case: Single observation per random effect level
  Z_single <- matrix(1, n, 1)
  K_single <- matrix(1, 1, 1)
  result_single <- mixed_solve_cpp(y = y, Z = Z_single, K = K_single, X = X)
  expect_equal(length(result_single$beta), ncol(X))
})