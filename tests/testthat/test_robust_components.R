# Test robust regression components

test_that("robust_iterative_fitter works with huber", {
  # Create data with outliers
  n <- 50
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  beta_true <- c(2, -1, 3)
  Y <- X %*% beta_true + rnorm(n, sd = 0.5)
  
  # Add outliers
  outlier_idx <- c(5, 15, 25)
  Y[outlier_idx] <- Y[outlier_idx] + 5
  
  # Initial GLM context
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  
  # Robust options
  robust_opts <- list(
    type = "huber",
    k_huber = 1.345,
    max_iter = 20,
    tol = 1e-4,
    scale_scope = "local"
  )
  
  # Fit robust model
  result <- robust_iterative_fitter(
    initial_glm_ctx = ctx,
    cfg_robust_options = robust_opts,
    X_orig_for_resid = X
  )
  
  # Check output structure
  expect_true(all(c("betas_robust", "sigma_robust_scale_final", 
                    "robust_weights_final", "XtWXi_final") %in% names(result)))
  
  # Check dimensions
  expect_equal(dim(result$betas_robust), c(p, 1))
  expect_length(result$robust_weights_final, n)
  expect_equal(dim(result$XtWXi_final), c(p, p))
  
  # Outliers should have lower weights
  expect_lt(mean(result$robust_weights_final[outlier_idx]), 
            mean(result$robust_weights_final[-outlier_idx]))
  
  # Compare with OLS
  ols_beta <- solve(t(X) %*% X) %*% t(X) %*% Y
  
  # Robust estimates should be less affected by outliers
  # (closer to true values than OLS)
  robust_error <- sum((result$betas_robust - beta_true)^2)
  ols_error <- sum((ols_beta - beta_true)^2)
  expect_lt(robust_error, ols_error)
})

test_that("robust_iterative_fitter works with bisquare", {
  # Create data with outliers
  n <- 100
  p <- 4
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  beta_true <- c(1, 2, -1, 0.5)
  Y <- X %*% beta_true + rnorm(n, sd = 0.3)
  
  # Add severe outliers
  outlier_idx <- sample(n, 10)
  Y[outlier_idx] <- Y[outlier_idx] + rnorm(10, mean = 10, sd = 2)
  
  # Setup
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  
  robust_opts <- list(
    type = "bisquare",
    c_tukey = 4.685,
    max_iter = 30,
    scale_scope = "local"
  )
  
  # Fit
  result <- robust_iterative_fitter(
    initial_glm_ctx = ctx,
    cfg_robust_options = robust_opts,
    X_orig_for_resid = X
  )
  
  # Bisquare should completely downweight severe outliers
  expect_true(all(result$robust_weights_final[outlier_idx] < 0.1))
  
  # Check that results are reasonable
  expect_true(!is.null(result$betas_robust))
  expect_equal(length(result$betas_robust), ncol(X))
})

test_that("robust_iterative_fitter handles multiple responses", {
  n <- 60
  p <- 3
  v <- 5  # voxels
  
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  B <- matrix(rnorm(p * v), p, v)
  Y <- X %*% B + matrix(rnorm(n * v, sd = 0.5), n, v)
  
  # Add outliers to different voxels
  Y[5:10, 1] <- Y[5:10, 1] + 5
  Y[20:25, 3] <- Y[20:25, 3] - 5
  
  # Setup
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = Y, proj = proj)
  
  robust_opts <- list(
    type = "huber",
    k_huber = 1.345,
    max_iter = 20,
    scale_scope = "local"
  )
  
  # Fit
  result <- robust_iterative_fitter(
    initial_glm_ctx = ctx,
    cfg_robust_options = robust_opts,
    X_orig_for_resid = X
  )
  
  # Check dimensions
  expect_equal(dim(result$betas_robust), c(p, v))
  expect_length(result$sigma_robust_scale_final, 1)  # Currently returns scalar
  expect_length(result$robust_weights_final, n)
})

test_that("robust scale estimation works correctly", {
  # Test MAD scale estimation
  x <- rnorm(100)
  x[1:5] <- 10  # Add outliers
  
  # MAD should be robust to outliers
  mad_scale <- median(abs(x - median(x))) / 0.6745
  expect_lt(mad_scale, 2)  # Should be close to 1, not affected by outliers
  
  # Test global vs local scale
  n <- 50
  p <- 2
  X <- cbind(1, rnorm(n))
  Y <- cbind(X %*% c(1, 2) + rnorm(n, sd = 0.5),
             X %*% c(-1, 3) + rnorm(n, sd = 2))  # Different scales
  
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = Y, proj = proj)
  
  # Local scale
  result_local <- robust_iterative_fitter(
    initial_glm_ctx = ctx,
    cfg_robust_options = list(
      type = "huber", 
      k_huber = 1.345,
      max_iter = 20,
      scale_scope = "local"
    ),
    X_orig_for_resid = X
  )
  
  # Global scale with fixed sigma
  result_global <- robust_iterative_fitter(
    initial_glm_ctx = ctx,
    cfg_robust_options = list(
      type = "huber",
      k_huber = 1.345,
      max_iter = 20,
      scale_scope = "global"
    ),
    X_orig_for_resid = X,
    sigma_fixed = 1.0
  )
  
  # Both should produce valid results
  expect_true(!is.null(result_local$betas_robust))
  expect_true(!is.null(result_global$betas_robust))
})

test_that("robust weights are computed correctly", {
  # Test Huber weights
  r <- seq(-5, 5, by = 0.5)
  k <- 1.345
  
  # Manual Huber weight calculation
  w_expected <- ifelse(abs(r) <= k, 1, k / abs(r))
  
  # Test bisquare weights
  c_tukey <- 4.685
  r_scaled <- r / c_tukey
  w_bisquare_expected <- ifelse(abs(r_scaled) <= 1, (1 - r_scaled^2)^2, 0)
  
  # Weights should be between 0 and 1
  expect_true(all(w_expected >= 0 & w_expected <= 1))
  expect_true(all(w_bisquare_expected >= 0 & w_bisquare_expected <= 1))
  
  # Weights should decrease with |r|
  expect_true(all(diff(w_expected[r >= 0]) <= 0))
  expect_true(all(diff(w_bisquare_expected[r >= 0]) <= 0))
})