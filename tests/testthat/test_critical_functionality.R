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
  config <- fmrireg:::fmri_lm_config(
    robust = FALSE,
    ar_options = list(cor_struct = "ar1", iter = 2)
  )
  
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
  config <- fmrireg:::fmri_lm_config(
    robust = "huber",
    ar_options = list(cor_struct = "none")
  )
  
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

test_that("Contrast computation works", {
  # Simple design
  n <- 60
  X <- cbind(1, rep(c(0, 1), each = 30), rnorm(n))
  Y <- X %*% c(1, 2, 0.5) + rnorm(n)
  
  config <- fmrireg:::fmri_lm_config(robust = FALSE)
  
  # Fit model
  fit_result <- solve_integrated_glm(X, matrix(Y, ncol = 1), config)
  
  # Define contrast (group difference)
  contrast_mat <- matrix(c(0, 1, 0), nrow = 1)
  
  # Compute contrast
  con_result <- fmrireg:::compute_contrast(fit_result, contrast_mat)
  
  expect_true(!is.null(con_result$estimate))
  expect_true(!is.null(con_result$stderr))
  expect_true(!is.null(con_result$tstat))
  expect_true(!is.null(con_result$pvalue))
  
  # Estimate should be close to true difference (2)
  expect_equal(con_result$estimate[1,1], 2, tolerance = 0.5)
})

test_that("Effective df calculations work", {
  n <- 100
  p <- 3
  
  # Test AR adjustment
  df_ar <- fmrireg:::compute_effective_df_ar(n, p, ar_coef = 0.6)
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
  
  df_robust <- fmrireg:::compute_effective_df_robust(X, weights, XtWXinv)
  expect_true(df_robust < n - p)  # Should be reduced due to downweighting
  expect_true(df_robust > 0)
})

test_that("Bootstrap functionality works", {
  # Small example
  n <- 50
  X <- cbind(1, rnorm(n))
  Y <- X %*% c(1, 2) + rnorm(n)
  
  config <- fmrireg:::fmri_lm_config(robust = FALSE)
  fit_result <- solve_integrated_glm(X, matrix(Y, ncol = 1), config)
  
  # Run bootstrap
  boot_result <- fmrireg:::bootstrap_glm_inference(
    fit_result = fit_result,
    X = X,
    Y = matrix(Y, ncol = 1),
    config = config,
    nboot = 100,
    block_size = 10
  )
  
  expect_equal(dim(boot_result$boot_betas), c(100, 2, 1))
  expect_true(!is.null(boot_result$beta_ci))
  expect_equal(dim(boot_result$beta_ci), c(2, 2, 1))  # 2 percentiles, 2 params, 1 voxel
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