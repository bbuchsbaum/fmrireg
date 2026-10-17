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
  config <- fmri_lm_control(noise = noise_spec("ar2", iter_gls = 1L))
  
  # Fit model with AR whitening
  result <- solve_integrated_glm(X, Y, config)
  
  # Check that coefficients are close to zero
  expect_true(all(abs(result$betas[2, ]) < 0.5), 
              "AR whitening should control false positives")

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
  config <- fmri_lm_control(
    noise = noise_spec("ar1", iter_gls = 1L),
    robust = robust_spec("bisquare", c_tukey = 4.685, max_iter = 5L)
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

test_that("Robust effective df is sum(w) - rank", {
  # Review 2026-10 (D4): compute_robust_effective_df() (unreachable; it
  # returned the weighted hat trace, which is always p) was deleted. The live
  # robust df is .robust_effective_df(): sum(w) - rank.
  set.seed(654)
  n <- 30
  p <- 2

  X <- cbind(1, rnorm(n))
  weights <- runif(n, 0.1, 1)  # Random weights

  edf <- fmrireg:::.robust_effective_df(weights, rank = qr(X)$rank)
  expect_gt(edf, 0, "Effective df should be positive")
  expect_lt(edf, n - p, "Effective df should be less than n - p")
  expect_equal(edf, sum(weights) - p)

  # The weighted hat trace is still the rank, for any weights.
  Xw <- X * sqrt(weights)
  expect_equal(sum((Xw %*% solve(crossprod(Xw))) * Xw), p, tolerance = 1e-10)

  # Uniform weights recover the OLS residual df.
  expect_equal(fmrireg:::.robust_effective_df(rep(1, n), rank = p), n - p)
})
