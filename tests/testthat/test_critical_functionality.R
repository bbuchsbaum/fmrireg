# Test critical functionality that was identified as missing

library(fmrireg)
library(testthat)

test_that("AR whitening integration works", {
  # Simple AR(1) data
  n <- 100
  X <- cbind(1, rnorm(n))
  ar_coef <- 0.6
  e <- arima.sim(list(ar = ar_coef), n = n)
  Y <- X %*% c(1, 2) + as.vector(e)
  
  # Create config with AR
  config <- fmri_lm_control(noise = noise_spec("ar1", iter_gls = 2L))
  
  # Test integrated solver
  result <- fmrireg:::solve_integrated_glm(
    X = X,
    Y = matrix(Y, ncol = 1),
    config = config,
    run_indices = list(1:n)
  )
  
  # Should have AR coefficients
  expect_true(!is.null(result$ar_coef))
  expect_true(length(result$ar_coef) > 0)
  
  # AR coefficient should be reasonable
  ar_est <- unlist(result$ar_coef)[1]
  # AR estimation from residuals is inherently noisy, especially with small samples
  expect_equal(ar_est, ar_coef, tolerance = 0.5)
})

test_that("Robust fitting works", {
  # Data with outliers
  n <- 100
  X <- cbind(1, rnorm(n))
  Y <- 2 + 3*X[,2] + rnorm(n)
  
  # Add outliers - make them more extreme and deterministic
  outliers <- c(10, 20, 30)
  Y[outliers] <- Y[outliers] + c(8, -8, 8)  # More extreme, deterministic outliers
  
  # Robust config
  config <- fmri_lm_control(robust = robust_spec("huber"))
  
  # Fit
  result <- fmrireg:::solve_integrated_glm(
    X = X,
    Y = matrix(Y, ncol = 1),
    config = config
  )
  
  # Should have robust weights
  expect_true(!is.null(result$robust_weights))
  expect_equal(length(result$robust_weights), n)
  
  # Outliers should have lower weights
  expect_true(all(result$robust_weights[outliers] < 0.8))
})

test_that("Effective df calculations work", {
  n <- 100
  p <- 3
  
  # Test AR adjustment
  # Review 2026-10 (D4): ported from the deleted compute_effective_df_ar() /
  # compute_effective_df_robust() to the live implementations.
  df_ar <- fmrireg:::compute_ar_effective_df(n, p, phi = 0.6)
  expect_true(df_ar < n - p)  # Should be reduced
  expect_true(df_ar > 0)
  
  # Test robust adjustment with downweighting
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  # Create weights that actually downweight some observations
  weights <- rep(1, n)
  weights[1:20] <- 0.3  # Downweight 20% of observations
  
  # Weighted covariance matrix
  W <- diag(weights)
  XtWX <- t(X) %*% W %*% X
  XtWXinv <- solve(XtWX)
  
  df_robust <- fmrireg:::calculate_effective_df(n, p, robust_weights = weights)
  expect_equal(df_robust, sum(weights) - p)
  expect_true(df_robust < n - p)  # Should be reduced due to downweighting
  expect_true(df_robust > 0)
})


test_that("Sandwich variance estimator works", {
  # Heteroscedastic data
  n <- 100
  X <- cbind(1, rnorm(n))
  # Variance increases with X
  Y <- 2 + 3*X[,2] + rnorm(n) * (1 + abs(X[,2]))
  
  residuals <- Y - X %*% solve(crossprod(X), crossprod(X, Y))
  
  # Compute sandwich variance
  sandwich_vcov <- fmrireg:::compute_sandwich_variance(X, matrix(residuals, ncol = 1), type = "HC1")
  
  expect_equal(dim(sandwich_vcov), c(2, 2))
  expect_true(all(diag(sandwich_vcov) > 0))
  
  # Should be different from standard OLS variance
  sigma2 <- sum(residuals^2) / (n - 2)
  ols_vcov <- solve(crossprod(X)) * sigma2
  expect_true(any(abs(sandwich_vcov - ols_vcov) > 0.01))
})
