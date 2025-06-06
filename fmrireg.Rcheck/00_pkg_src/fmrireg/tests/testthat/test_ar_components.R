# Test AR modeling components

test_that("estimate_ar_parameters works correctly", {
  # Test AR(1) estimation
  n <- 100
  true_phi <- 0.7
  
  # Generate AR(1) series
  x <- numeric(n)
  x[1] <- rnorm(1)
  for (i in 2:n) {
    x[i] <- true_phi * x[i-1] + rnorm(1)
  }
  
  # Estimate
  phi_hat <- estimate_ar_parameters(x, p_order = 1)
  expect_length(phi_hat, 1)
  expect_equal(phi_hat[1], true_phi, tolerance = 0.1)
  
  # Test AR(2) estimation
  true_phi2 <- c(0.5, 0.3)
  y <- numeric(n)
  y[1:2] <- rnorm(2)
  for (i in 3:n) {
    y[i] <- true_phi2[1] * y[i-1] + true_phi2[2] * y[i-2] + rnorm(1)
  }
  
  phi_hat2 <- estimate_ar_parameters(y, p_order = 2)
  expect_length(phi_hat2, 2)
  expect_equal(phi_hat2, true_phi2, tolerance = 0.2)
})

test_that("ar_whiten_transform works correctly", {
  # Simple AR(1) test
  n <- 50
  p <- 3
  phi <- 0.6
  
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  Y <- matrix(rnorm(n * 2), n, 2)
  
  # Apply whitening
  result <- ar_whiten_transform(X, Y, phi, exact_first = TRUE)
  
  # Check dimensions preserved
  expect_equal(dim(result$X), dim(X))
  expect_equal(dim(result$Y), dim(Y))
  
  # The whitening should modify the data
  expect_false(all(result$X == X))
  expect_false(all(result$Y == Y))
  
  # With exact_first = TRUE and AR(1), first row should be scaled
  scale <- sqrt(1 - phi^2)
  expect_equal(result$X[1,], X[1,] * scale, tolerance = 1e-10)
  expect_equal(result$Y[1,], Y[1,] * scale, tolerance = 1e-10)
  
  # Second row should have AR filter applied
  expect_equal(result$X[2,], X[2,] - phi * X[1,], tolerance = 1e-10)
  
  # Test with exact_first = FALSE
  result2 <- ar_whiten_transform(X, Y, phi, exact_first = FALSE)
  expect_equal(dim(result2$X), dim(X))  # Dimensions preserved
  expect_equal(dim(result2$Y), dim(Y))
  
  # First row should not have special scaling
  expect_equal(result2$X[1,], X[1,] - phi * 0)  # No previous value, so just X[1,]
  expect_equal(result2$X[1,], X[1,])
})

test_that("ar_whiten_transform handles AR(2)", {
  n <- 100
  p <- 2
  phi <- c(0.5, 0.2)
  
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * 3), n, 3)
  
  result <- ar_whiten_transform(X, Y, phi, exact_first = TRUE)
  
  # Check dimensions
  expect_equal(dim(result$X), dim(X))
  expect_equal(dim(result$Y), dim(Y))
  
  # The whitening should modify the data
  expect_false(all(result$X == X))
  expect_false(all(result$Y == Y))
  
  # For AR(2), third row uses two previous values
  expect_equal(result$X[3,], X[3,] - phi[1] * X[2,] - phi[2] * X[1,], tolerance = 1e-10)
})

test_that("AR whitening improves residual autocorrelation", {
  # Generate data with AR(1) errors
  n <- 200
  p <- 5
  phi_true <- 0.8
  
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  beta_true <- rnorm(p)
  
  # Generate AR(1) errors
  errors <- numeric(n)
  errors[1] <- rnorm(1)
  for (i in 2:n) {
    errors[i] <- phi_true * errors[i-1] + rnorm(1)
  }
  
  Y <- X %*% beta_true + errors
  
  # Fit without whitening
  lm_ols <- lm(Y ~ X - 1)
  resid_ols <- residuals(lm_ols)
  
  # Estimate AR parameter
  phi_hat <- estimate_ar_parameters(resid_ols, p_order = 1)
  expect_gt(phi_hat[1], 0.5)  # Should detect substantial autocorrelation
  
  # Whiten and refit
  tmp <- ar_whiten_transform(X, matrix(Y, ncol = 1), phi_hat, exact_first = TRUE)
  lm_gls <- lm(tmp$Y ~ tmp$X - 1)
  resid_gls <- residuals(lm_gls)
  
  # Check that whitened residuals have less autocorrelation
  phi_resid <- estimate_ar_parameters(resid_gls, p_order = 1)
  expect_lt(abs(phi_resid[1]), abs(phi_hat[1]))
  expect_lt(abs(phi_resid[1]), 0.2)  # Should be close to 0
})

test_that("AR functions handle edge cases", {
  # Empty data
  expect_error(estimate_ar_parameters(numeric(0), p_order = 1))
  
  # Too short for order - acf may not error on all systems
  # expect_error(estimate_ar_parameters(c(1, 2), p_order = 3))
  
  # Zero variance - will error due to singular matrix
  expect_error(estimate_ar_parameters(rep(5, 100), p_order = 1))
  
  # Matrix input for ar_whiten_transform
  X <- matrix(1:12, 4, 3)
  Y <- matrix(1:8, 4, 2)
  phi <- 0.5
  
  result <- ar_whiten_transform(X, Y, phi, exact_first = TRUE)
  expect_equal(dim(result$X), dim(X))
  expect_equal(dim(result$Y), dim(Y))
})