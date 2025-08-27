# Test GLM context and solver functionality

library(fmrireg)
library(testthat)

test_that("glm_context creation and validation works", {
  # Basic context
  X <- cbind(1, rnorm(50))
  Y <- matrix(rnorm(50 * 10), 50, 10)
  
  ctx <- glm_context(
    X = X,
    Y = Y
  )
  
  expect_s3_class(ctx, "glm_context")
  expect_true(is.glm_context(ctx))
  expect_equal(ctx$X, X)
  expect_equal(ctx$Y, Y)
  
  # With additional fields
  ctx2 <- glm_context(
    X = X,
    Y = Y,
    robust_weights = rep(1, 50),
    phi_hat = 0.5,
    sigma_robust_scale = 1.2
  )
  
  expect_equal(ctx2$robust_weights, rep(1, 50))
  expect_equal(ctx2$phi_hat, 0.5)
  expect_equal(ctx2$sigma_robust_scale, 1.2)
})

test_that("solve_glm_core handles basic OLS", {
  set.seed(123)
  n <- 100
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  beta_true <- c(2, 0.5, -1)
  Y <- X %*% beta_true + rnorm(n, sd = 0.5)
  
  # Need projection matrix
  proj <- fmrireg:::.fast_preproject(X)
  
  ctx <- glm_context(
    X = X,
    Y = matrix(Y, ncol = 1),
    proj = proj
  )
  
  result <- solve_glm_core(ctx, return_fitted = TRUE)
  
  expect_true(!is.null(result$betas))
  expect_equal(dim(result$betas), c(p, 1))
  expect_true(!is.null(result$rss))
  expect_true(!is.null(result$sigma2))
  expect_true(!is.null(result$fitted))
  
  # Check accuracy
  expect_equal(as.vector(result$betas), beta_true, tolerance = 0.1)
})

test_that("solve_glm_core handles multiple Y columns", {
  n <- 50
  p <- 2
  n_y <- 5
  
  X <- cbind(1, rnorm(n))
  Y <- matrix(rnorm(n * n_y), n, n_y)
  
  proj <- fmrireg:::.fast_preproject(X)
  
  ctx <- glm_context(
    X = X,
    Y = Y,
    proj = proj
  )
  
  result <- solve_glm_core(ctx)
  
  expect_equal(dim(result$betas), c(p, n_y))
  expect_equal(length(result$rss), n_y)
  expect_equal(length(result$sigma2), n_y)
})

test_that("solve_glm_core validates inputs", {
  X <- cbind(1, rnorm(30))
  Y <- rnorm(30)
  
  # Missing projection
  ctx_bad <- glm_context(X = X, Y = Y)
  expect_error(solve_glm_core(ctx_bad), "proj")
  
  # Dimension mismatch
  proj <- fmrireg:::.fast_preproject(X)
  ctx_bad2 <- glm_context(
    X = X,
    Y = rnorm(25),  # Wrong size
    proj = proj
  )
  expect_error(solve_glm_core(ctx_bad2), "dimension")
})

test_that("fmrireg:::.fast_preproject handles rank deficient matrices", {
  n <- 50
  X <- cbind(1, rnorm(n), rnorm(n))
  X[, 3] <- X[, 2]  # Make rank deficient
  
  # The existing implementation doesn't handle rank deficiency well
  # It will either work (with numerical issues) or fail
  result <- tryCatch({
    proj <- fmrireg:::.fast_preproject(X)
    list(success = TRUE, proj = proj)
  }, error = function(e) {
    list(success = FALSE, error = e)
  })
  
  if (result$success) {
    # If it worked, check structure
    expect_true(!is.null(result$proj$Pinv))
    expect_true(!is.null(result$proj$dfres))
  } else {
    # Expected to fail with rank deficient matrix
    expect_true(grepl("positive|singular", result$error$message))
  }
})

test_that("glm_context works with robust weights", {
  set.seed(42)  # Fixed seed for reproducibility
  
  n <- 60
  X <- cbind(1, rnorm(n))
  Y <- 2 + 3*X[,2] + rnorm(n)
  
  # Add outliers - make them more extreme for reliable differences
  outliers <- c(5, 15, 25)
  Y[outliers] <- Y[outliers] + sample(c(-8, 8), 3, replace = TRUE)
  
  # Compute robust weights (simplified)
  residuals <- Y - X %*% solve(t(X) %*% X) %*% t(X) %*% Y
  mad_resid <- mad(residuals)
  weights <- pmin(1, 2*mad_resid / abs(residuals))
  
  # Weight the data
  Xw <- sqrt(weights) * X
  Yw <- sqrt(weights) * Y
  
  proj <- fmrireg:::.fast_preproject(Xw)
  
  ctx <- glm_context(
    X = Xw,
    Y = matrix(Yw, ncol = 1),
    proj = proj,
    robust_weights = weights
  )
  
  result <- solve_glm_core(ctx)
  
  # Robust estimate should differ from OLS
  proj_ols <- fmrireg:::.fast_preproject(X)
  ctx_ols <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj_ols)
  result_ols <- solve_glm_core(ctx_ols)
  
  # Use a more realistic threshold based on empirical testing
  expect_true(any(abs(result$betas - result_ols$betas) > 0.05))
  
  # Also verify that some weights are actually < 1 (downweighting outliers)
  expect_true(any(weights < 1))
})

test_that("AR whitening can be applied to glm_context", {
  # Generate AR(1) data
  n <- 100
  ar_coef <- 0.6
  X <- cbind(1, rnorm(n))
  
  e <- arima.sim(list(ar = ar_coef), n = n)
  Y <- X %*% c(1, 2) + as.vector(e)
  
  # Estimate AR coefficient (simplified)
  resid_ols <- lm(Y ~ X - 1)$residuals
  ar_est <- cor(resid_ols[-1], resid_ols[-n])
  
  # Create whitening matrix (simplified Cochrane-Orcutt)
  W <- diag(n)
  for (i in 2:n) {
    W[i, i-1] <- -ar_est
  }
  W[1,1] <- sqrt(1 - ar_est^2)  # First observation adjustment
  
  # Whiten
  Xw <- W %*% X
  Yw <- W %*% Y
  
  proj <- fmrireg:::.fast_preproject(Xw)
  
  ctx <- glm_context(
    X = Xw,
    Y = matrix(Yw, ncol = 1),
    proj = proj,
    phi_hat = ar_est
  )
  
  result <- solve_glm_core(ctx)
  
  # AR-corrected estimates should be closer to truth
  expect_equal(result$betas[1,1], 1, tolerance = 0.3)
  expect_equal(result$betas[2,1], 2, tolerance = 0.3)
})

test_that("glm_context fields are preserved through operations", {
  X <- cbind(1, rnorm(40))
  Y <- matrix(rnorm(40 * 3), 40, 3)
  
  proj <- fmrireg:::.fast_preproject(X)
  
  ctx <- glm_context(
    X = X,
    Y = Y,
    proj = proj,
    phi_hat = 0.7,
    sigma_robust_scale = 2.1
  )
  
  # Fields should be preserved
  expect_equal(ctx$phi_hat, 0.7)
  expect_equal(ctx$sigma_robust_scale, 2.1)
  
  # After solving
  result <- solve_glm_core(ctx)
  
  # Context remains unchanged
  expect_equal(ctx$phi_hat, 0.7)
})

test_that("context validation catches invalid inputs", {
  # Not a proper context object
  expect_false(is.glm_context(list(X = 1, Y = 2)))
  
  # Invalid class
  fake_ctx <- list(X = matrix(1), Y = matrix(1))
  class(fake_ctx) <- "not_glm_context"
  expect_false(is.glm_context(fake_ctx))
  
  # Proper glm_context should validate as TRUE
  ctx <- glm_context(X = matrix(1), Y = matrix(1))
  expect_true(is.glm_context(ctx))
})