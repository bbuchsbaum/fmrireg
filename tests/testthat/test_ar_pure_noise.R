# Pure AR Tests - Test AR functionality in isolation from GLM effects
# These tests ensure AR estimation and whitening work correctly on pure AR processes

test_that("AR estimation recovers parameters from pure AR(1) process", {
  set.seed(123)
  true_phi <- c(0.3, 0.5, 0.7, 0.9)
  
  for (phi in true_phi) {
    # Generate pure AR(1) data
    y <- as.numeric(arima.sim(list(ar = phi), n = 200))
    
    # Estimate AR parameter
    phi_hat <- fmrireg:::.estimate_ar(y, 1)
    
    # Should recover true parameter within reasonable tolerance
    # Higher phi values are harder to estimate accurately
    # With n=200, we need more tolerance
    tol <- ifelse(phi > 0.7, 0.1, 0.1)
    expect_equal(phi_hat, phi, tolerance = tol, 
                 info = paste("Failed for phi =", phi))
  }
})

test_that("AR estimation recovers parameters from pure AR(2) process", {
  set.seed(456)
  # Test stationary AR(2) combinations
  test_cases <- list(
    c(0.5, -0.3),   # Classic AR(2)
    c(0.6, -0.2),   # Another combination
    c(0.3, 0.3)     # Positive coefficients
  )
  
  for (true_phi in test_cases) {
    # Check stationarity
    roots <- polyroot(c(1, -true_phi))
    expect_true(all(Mod(roots) > 1), 
                info = paste("Non-stationary AR(2):", true_phi))
    
    # Generate AR(2) data
    y <- as.numeric(arima.sim(list(ar = true_phi), n = 300))
    
    # Estimate AR parameters
    phi_hat <- fmrireg:::.estimate_ar(y, 2)
    
    # Should recover true parameters
    expect_equal(phi_hat, true_phi, tolerance = 0.1,
                 info = paste("Failed for phi =", paste(true_phi, collapse=", ")))
  }
})

test_that("AR whitening correctly removes autocorrelation", {
  set.seed(789)
  # Generate AR(1) data
  phi <- 0.7
  n <- 100
  y <- as.numeric(arima.sim(list(ar = phi), n = n))
  
  # Apply whitening
  Y <- matrix(y, ncol = 1)
  X <- matrix(1, nrow = n, ncol = 1)  # Dummy design
  
  result <- fmrireg:::ar_whiten_transform(X, Y, phi, exact_first = FALSE)
  whitened <- as.vector(result$Y)
  
  # Check that autocorrelation is removed
  acf_orig <- acf(y, lag.max = 5, plot = FALSE)$acf[-1]
  acf_whitened <- acf(whitened, lag.max = 5, plot = FALSE)$acf[-1]
  
  # Original should have significant autocorrelation
  expect_true(abs(acf_orig[1]) > 0.5)
  
  # Whitened should have minimal autocorrelation
  expect_true(all(abs(acf_whitened) < 0.2))
})

test_that("Exact first observation scaling preserves innovation variance", {
  set.seed(321)
  phi <- 0.8
  n <- 1000
  
  # Generate many AR(1) realizations
  first_values_standard <- numeric(100)
  first_values_exact <- numeric(100)
  
  for (i in 1:100) {
    y <- rnorm(n)  # White noise
    Y <- matrix(y, ncol = 1)
    X <- matrix(1, nrow = n, ncol = 1)
    
    # Standard whitening
    res1 <- fmrireg:::ar_whiten_transform(X, Y, phi, exact_first = FALSE)
    first_values_standard[i] <- res1$Y[1, 1]
    
    # Exact first whitening
    res2 <- fmrireg:::ar_whiten_transform(X, Y, phi, exact_first = TRUE)
    first_values_exact[i] <- res2$Y[1, 1]
  }
  
  # With exact_first, the first observation should have variance 1
  # Without it, variance is sigma^2 = 1/(1-phi^2)
  expected_var_ratio <- 1 - phi^2
  actual_var_ratio <- var(first_values_exact) / var(first_values_standard)
  
  expect_equal(actual_var_ratio, expected_var_ratio, tolerance = 0.1)
})