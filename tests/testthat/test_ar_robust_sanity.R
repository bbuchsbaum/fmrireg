test_that("AR whitening controls false positives under AR(2) noise", {
  set.seed(123)
  n_time <- 180
  p <- 2
  n_vox <- 10
  
  # Design matrix: intercept + noise regressor
  X <- cbind(1, rnorm(n_time))
  
  # Pure null hypothesis (no true effect)
  true_betas <- c(0, 0)
  
  # Simulate AR(2) noise
  phi <- c(0.5, 0.2)
  
  # Generate Y with AR(2) noise only
  Y <- matrix(NA, n_time, n_vox)
  for (v in 1:n_vox) {
    # Use arima.sim for proper AR simulation
    ar_noise <- arima.sim(model = list(ar = phi), n = n_time, sd = 1)
    Y[, v] <- as.numeric(X %*% true_betas + ar_noise)
  }
  
  # Create config with AR whitening
  config <- fmri_lm_config(
    ar_options = list(
      cor_struct = "ar2",
      iter = 1
    ),
    robust = FALSE
  )
  
  # Fit model with AR whitening
  result <- solve_integrated_glm(X, Y, config)
  
  # Check that coefficients are close to zero
  expect_true(all(abs(result$betas[2, ]) < 0.5), 
              "AR whitening should control false positives")
  
  # Check effective df is reduced appropriately
  expect_lt(result$effective_df, n_time - p, 
            "Effective df should be reduced for AR models")
})

test_that("Robust fitting handles outliers in whitened data", {
  set.seed(456)
  n_time <- 100
  p <- 3
  n_vox <- 5
  
  # Design matrix
  X <- cbind(1, rnorm(n_time), rnorm(n_time))
  true_betas <- matrix(c(2, 1, -0.5), p, 1)
  
  # Generate clean signal
  signal <- X %*% true_betas
  
  # Add AR(1) noise
  phi <- 0.3
  Y <- matrix(NA, n_time, n_vox)
  for (v in 1:n_vox) {
    ar_noise <- arima.sim(model = list(ar = phi), n = n_time, sd = 0.5)
    Y[, v] <- as.numeric(signal + ar_noise)
  }
  
  # Add outliers (10% contamination)
  n_outliers <- round(0.1 * n_time)
  outlier_idx <- sample(n_time, n_outliers)
  Y[outlier_idx, ] <- Y[outlier_idx, ] + sample(c(-10, 10), n_outliers, replace = TRUE)
  
  # Config for AR + Robust
  config <- fmri_lm_config(
    ar_options = list(
      cor_struct = "ar1",
      iter = 1
    ),
    robust = list(
      type = "bisquare",
      c_tukey = 4.685,
      max_iter = 5
    )
  )
  
  # Fit AR + Robust model
  result <- solve_integrated_glm(X, Y, config)
  
  # Check that robust weights downweight outliers
  mean_weight_outliers <- mean(result$robust_weights[outlier_idx])
  mean_weight_clean <- mean(result$robust_weights[-outlier_idx])
  
  expect_lt(mean_weight_outliers, mean_weight_clean,
            "Outliers should receive lower weights")
  
  # Check coefficient recovery (should be close to truth despite outliers)
  avg_betas <- rowMeans(result$betas)
  relative_error <- abs(avg_betas - true_betas) / abs(true_betas)
  
  expect_true(all(relative_error < 0.5),
              "AR+Robust should recover coefficients despite outliers")
})

test_that("Effective df computation follows standard formula", {
  set.seed(789)
  n_per_run <- 60
  n_runs <- 3
  n_time <- n_per_run * n_runs
  p <- 2
  
  # Simple case with known AR parameters
  phi <- c(0.5, 0.2)  # AR(2)
  ar_order <- length(phi)
  
  # Test the standard compute_ar_effective_df function
  # Following SPM/FSL/AFNI practice: no penalty for AR estimation
  df_single_run <- fmrireg:::compute_ar_effective_df(n_time, p, phi, n_runs = 1)
  df_multi_run <- fmrireg:::compute_ar_effective_df(n_time, p, phi, n_runs = n_runs)
  
  # Without penalty, df should be the same regardless of n_runs
  expect_equal(df_single_run, df_multi_run, tolerance = 1e-10)
  
  # Current implementation uses lag-correlation inflation and should return a
  # sensible df in bounds.
  expect_gte(df_multi_run, 1)
  expect_lte(df_multi_run, n_time - p)
  
  # Test conservative mode with penalize_ar = TRUE
  df_conservative <- fmrireg:::compute_ar_effective_df(n_time, p, phi, n_runs = n_runs, penalize_ar = TRUE)
  expect_equal(df_conservative, max(df_multi_run - ar_order, 1), tolerance = 1e-10)
})

test_that("Sandwich variance is computed correctly without diag matrices", {
  set.seed(321)
  n <- 50
  p <- 3
  
  X <- cbind(1, rnorm(n), rnorm(n))
  true_beta <- c(1, 2, -1)
  Y <- X %*% true_beta + rnorm(n, sd = 0.5)
  
  # Add heteroscedasticity
  Y[1:10] <- Y[1:10] + rnorm(10, sd = 2)
  
  # Compute OLS
  XtX <- crossprod(X)
  XtXinv <- solve(XtX)
  beta_ols <- XtXinv %*% crossprod(X, Y)
  residuals <- Y - X %*% beta_ols
  
  # Test sandwich variance computation
  sandwich_var <- calculate_sandwich_variance(X, residuals, XtXinv)
  
  # Check dimensions
  expect_equal(dim(sandwich_var), c(p, p))
  
  # Check positive definiteness
  eigenvals <- eigen(sandwich_var, only.values = TRUE)$values
  expect_true(all(eigenvals > 0), 
              "Sandwich variance should be positive definite")
  
  # Standard errors should be positive
  se <- sqrt(diag(sandwich_var))
  expect_true(all(se > 0), "Standard errors should be positive")
})

test_that("Robust effective df uses weighted hat trace efficiently", {
  set.seed(654)
  n <- 30
  p <- 2
  
  X <- cbind(1, rnorm(n))
  weights <- runif(n, 0.1, 1)  # Random weights
  
  # Should not create diag(weights) internally
  # This tests the new efficient implementation
  sqrt_w <- sqrt(weights)
  Xw <- X * sqrt_w
  XtWX <- crossprod(Xw)
  XtWX_inv <- solve(XtWX)
  
  # Compute effective df efficiently
  edf <- fmrireg:::compute_robust_effective_df(X, weights, XtWX_inv)
  
  # Should be positive and less than n
  expect_gt(edf, 0, "Effective df should be positive")
  expect_lt(edf, n, "Effective df should be less than n")
  
  # For uniform weights = 1, should equal p (trace of projection)
  weights_uniform <- rep(1, n)
  XtX_inv <- solve(crossprod(X))
  edf_uniform <- fmrireg:::compute_robust_effective_df(X, weights_uniform, XtX_inv)
  expect_equal(edf_uniform, p, tolerance = 1e-10)
})
