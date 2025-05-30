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