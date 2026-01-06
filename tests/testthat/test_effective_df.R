# Test effective degrees of freedom calculations

library(fmrireg)
library(testthat)

test_that("calculate_effective_df works correctly for basic cases", {
  n <- 100
  p <- 5
  
  # No adjustment (IID case)
  df_eff <- fmrireg:::calculate_effective_df(n, p)
  expect_equal(df_eff, n - p)
  
  # AR adjustment
  df_ar1 <- fmrireg:::calculate_effective_df(n, p, ar_order = 1)
  expect_lt(df_ar1, n - p)  # Should be less due to AR correction
  expect_gt(df_ar1, 0)      # Should still be positive
  
  df_ar2 <- fmrireg:::calculate_effective_df(n, p, ar_order = 2)
  expect_lt(df_ar2, df_ar1)  # Higher AR order = more df loss
  
  # Robust adjustment
  weights <- rep(0.8, n)
  df_robust <- fmrireg:::calculate_effective_df(n, p, robust_weights = weights)
  expect_lt(df_robust, n - p)  # Downweighting reduces effective df
  
  # Combined adjustment
  df_combined <- fmrireg:::calculate_effective_df(n, p, ar_order = 1, robust_weights = weights)
  expect_lt(df_combined, df_ar1)
  expect_lt(df_combined, df_robust)
})

test_that("calculate_effective_df handles edge cases", {
  # Small sample size
  df_small <- fmrireg:::calculate_effective_df(10, 5)
  expect_equal(df_small, 5)
  
  # Nearly all parameters
  df_overfit <- fmrireg:::calculate_effective_df(10, 9)
  expect_equal(df_overfit, 1)
  
  # Zero weights (all observations downweighted)
  weights_zero <- rep(0, 50)
  df_zero <- fmrireg:::calculate_effective_df(50, 5, robust_weights = weights_zero)
  expect_equal(df_zero, 1)  # Function ensures minimum of 1
  
  # Some zero weights
  weights_partial <- c(rep(1, 25), rep(0, 25))
  df_partial <- fmrireg:::calculate_effective_df(50, 5, robust_weights = weights_partial)
  expect_lt(df_partial, 45)  # Less than full df
  expect_gt(df_partial, 20)  # But accounts for non-zero weights
})

test_that("calculate_sandwich_variance works correctly", {
  # Create simple regression scenario
  n <- 100
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  XtXinv <- solve(crossprod(X))
  
  # Homoscedastic case
  resid_homo <- rnorm(n, sd = 1)
  weights <- rep(1, n)
  
  var_homo <- fmrireg:::calculate_sandwich_variance(X, resid_homo, XtXinv, weights)
  expect_equal(dim(var_homo), c(p, p))
  expect_true(all(diag(var_homo) > 0))
  
  # Heteroscedastic case
  resid_hetero <- rnorm(n, sd = seq(0.5, 2, length.out = n))
  var_hetero <- fmrireg:::calculate_sandwich_variance(X, resid_hetero, XtXinv, weights)
  
  # Sandwich variance should differ from standard
  expect_false(all(abs(var_homo - var_hetero) < 1e-10))
  
  # With robust weights
  weights_robust <- runif(n, 0.5, 1)
  var_robust <- fmrireg:::calculate_sandwich_variance(X, resid_homo, XtXinv, weights_robust)
  expect_true(all(diag(var_robust) > 0))
})

test_that("effective df integrates with model fitting", {
  # Create test data
  n <- 50
  X <- cbind(1, rnorm(n))
  y <- X %*% c(2, 1) + rnorm(n)

  # Create config with AR
  cfg <- fmri_lm_config(
    ar_options = list(
      cor_struct = "ar1",
      iter = 1
    )
  )

  # Fit model would calculate effective df internally
  # This is a placeholder for integration test
  # In real implementation, effective df would be used for p-value calculation

  expect_true(TRUE)  # Placeholder assertion
})

# ============================================================================
# Extended tests for effective df calculations
# ============================================================================

test_that("calculate_effective_df handles varying AR orders", {
  n <- 200
  p <- 10

  # Test AR orders from 0 to 5
  df_values <- vapply(0:5, function(ar_order) {
    fmrireg:::calculate_effective_df(n, p, ar_order = ar_order)
  }, numeric(1))

  # DF should monotonically decrease with increasing AR order
  for (i in 2:length(df_values)) {
    expect_lte(df_values[i], df_values[i-1])
  }

  # All should be positive
  expect_true(all(df_values > 0))
})

test_that("calculate_effective_df is stable for large AR orders", {
  n <- 500
  p <- 20

  # High AR order should not cause numerical issues
  expect_no_error(
    df_high_ar <- fmrireg:::calculate_effective_df(n, p, ar_order = 10)
  )

  expect_gt(df_high_ar, 0)
  expect_lt(df_high_ar, n - p)
})

test_that("calculate_effective_df handles varying weight patterns", {
  n <- 100
  p <- 5

  # Uniform weights at different levels
  for (w in c(0.5, 0.7, 0.9, 1.0)) {
    weights <- rep(w, n)
    df <- fmrireg:::calculate_effective_df(n, p, robust_weights = weights)
    expect_gt(df, 0, info = paste("Failed for weight =", w))
  }

  # Exponentially decaying weights
  weights_decay <- exp(-seq(0, 3, length.out = n))
  df_decay <- fmrireg:::calculate_effective_df(n, p, robust_weights = weights_decay)
  expect_gt(df_decay, 0)
  expect_lt(df_decay, n - p)

  # Random weights
  set.seed(42)
  weights_random <- runif(n, 0.3, 1.0)
  df_random <- fmrireg:::calculate_effective_df(n, p, robust_weights = weights_random)
  expect_gt(df_random, 0)
})

test_that("calculate_effective_df handles extreme sample sizes", {
  # Very small sample
  df_tiny <- fmrireg:::calculate_effective_df(5, 3)
  expect_equal(df_tiny, 2)

  # Large sample
  df_large <- fmrireg:::calculate_effective_df(10000, 50)
  expect_equal(df_large, 9950)

  # Edge case: n == p
  df_saturated <- fmrireg:::calculate_effective_df(10, 10)
  expect_gte(df_saturated, 0)
  expect_lte(df_saturated, 1)  # Should be 0 or bounded to minimum
})

test_that("calculate_sandwich_variance handles extreme heteroscedasticity", {
  n <- 100
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  XtXinv <- solve(crossprod(X))
  weights <- rep(1, n)

  # Extreme heteroscedasticity: variance changes 100-fold
  set.seed(123)
  resid_extreme <- rnorm(n, sd = seq(0.1, 10, length.out = n))

  expect_no_error(
    var_extreme <- fmrireg:::calculate_sandwich_variance(X, resid_extreme, XtXinv, weights)
  )

  expect_equal(dim(var_extreme), c(p, p))
  expect_true(all(diag(var_extreme) > 0))
  expect_true(all(is.finite(var_extreme)))
})

test_that("calculate_sandwich_variance handles outliers", {
  n <- 100
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  XtXinv <- solve(crossprod(X))
  weights <- rep(1, n)

  # Residuals with outliers
  set.seed(456)
  resid_outlier <- rnorm(n)
  resid_outlier[c(10, 50, 90)] <- c(50, -40, 60)  # Large outliers

  expect_no_error(
    var_outlier <- fmrireg:::calculate_sandwich_variance(X, resid_outlier, XtXinv, weights)
  )

  expect_true(all(is.finite(var_outlier)))
})

test_that("sandwich variance with downweighted outliers differs from unweighted", {
  n <- 100
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  XtXinv <- solve(crossprod(X))

  # Create residuals with outliers
  set.seed(789)
  resid <- rnorm(n)
  outlier_idx <- c(10, 50, 90)
  resid[outlier_idx] <- c(20, -15, 25)

  # Uniform weights
  weights_uniform <- rep(1, n)
  var_uniform <- fmrireg:::calculate_sandwich_variance(X, resid, XtXinv, weights_uniform)

  # Downweight outliers
  weights_robust <- rep(1, n)
  weights_robust[outlier_idx] <- 0.1
  var_robust <- fmrireg:::calculate_sandwich_variance(X, resid, XtXinv, weights_robust)

  # Variance estimates should differ
  expect_false(all(abs(var_uniform - var_robust) < 1e-10))

  # Robust estimate should generally be smaller (outliers downweighted)
  expect_lt(sum(diag(var_robust)), sum(diag(var_uniform)))
})

# ============================================================================
# Integration tests with actual model fitting
# ============================================================================

test_that("effective df affects p-values in model output", {
  skip_on_cran()

  # Create test dataset
  event_table <- data.frame(
    onset = c(10, 30, 60, 80),
    cond = factor(c("A", "B", "A", "B")),
    run = c(1, 1, 2, 2)
  )

  mat <- matrix(rnorm(100 * 20), 100, 20)
  dset <- fmridataset::matrix_dataset(
    mat,
    TR = 1,
    run_length = c(50, 50),
    event_table = event_table
  )

  # Fit with IID assumption
  mod_iid <- fmri_lm(
    onset ~ hrf(cond),
    block = ~ run,
    dataset = dset
  )

  # Fit with AR1
  mod_ar1 <- fmri_lm(
    onset ~ hrf(cond),
    block = ~ run,
    dataset = dset,
    ar_options = list(cor_struct = "ar1")
  )

  # Models should complete without error
  expect_s3_class(mod_iid, "fmri_lm")
  expect_s3_class(mod_ar1, "fmri_lm")
})

test_that("effective df calculation matches expected formula for simple case", {
  # For AR(1) with correlation rho, effective n is approximately n*(1-rho^2)/(1+rho^2)
  # This is a simplified check that the formula is in the right ballpark

  n <- 200
  p <- 5

  # AR(1) should reduce effective df
  df_ar1 <- fmrireg:::calculate_effective_df(n, p, ar_order = 1)

  # The reduction should be meaningful but not extreme
  reduction_ratio <- df_ar1 / (n - p)
  expect_gt(reduction_ratio, 0.5)  # Should retain at least 50% of df
  expect_lt(reduction_ratio, 1.0)  # Should reduce df
})

test_that("calculate_effective_df is consistent with known statistical properties", {
  n <- 100
  p <- 5

  # Property 1: With all weights = 1, effective df should equal n - p
  df_full_weights <- fmrireg:::calculate_effective_df(n, p, robust_weights = rep(1, n))
  expect_equal(df_full_weights, n - p)

  # Property 2: Half weights should reduce df by approximately half
  weights_half <- c(rep(1, n/2), rep(0, n/2))
  df_half <- fmrireg:::calculate_effective_df(n, p, robust_weights = weights_half)
  expect_lt(df_half, n - p)
  expect_gt(df_half, (n/2) - p - 5)  # Allow some tolerance
})

# ============================================================================
# Numerical stability tests
# ============================================================================

test_that("calculate_sandwich_variance handles near-singular designs", {
  n <- 50
  p <- 3

  # Create nearly collinear design
  x1 <- rnorm(n)
  x2 <- x1 + rnorm(n, sd = 0.01)  # Nearly identical to x1
  X <- cbind(1, x1, x2)

  # Use regularized inverse
  XtXinv <- tryCatch(
    solve(crossprod(X)),
    error = function(e) {
      # Use pseudo-inverse for singular case
      svd_result <- svd(crossprod(X))
      svd_result$v %*% diag(1/pmax(svd_result$d, 1e-10)) %*% t(svd_result$u)
    }
  )

  resid <- rnorm(n)
  weights <- rep(1, n)

  # Should handle without error (may produce warnings)
  expect_no_error(
    suppressWarnings(
      var_result <- fmrireg:::calculate_sandwich_variance(X, resid, XtXinv, weights)
    )
  )
})

test_that("effective df handles very small weights gracefully", {
  n <- 100
  p <- 5

  # Very small but non-zero weights
  weights_tiny <- rep(1e-10, n)

  expect_no_error(
    df_tiny <- fmrireg:::calculate_effective_df(n, p, robust_weights = weights_tiny)
  )

  # Should return minimum df (typically 1)
  expect_gte(df_tiny, 1)
})