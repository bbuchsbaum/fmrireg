# Test GLM context and solver components

library(fmrireg)

test_that("glm_context creates valid objects", {
  # Create test data
  n <- 100
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * 10), n, 10)
  
  # Create context
  ctx <- glm_context(X = X, Y = Y)
  
  expect_s3_class(ctx, "glm_context")
  expect_equal(dim(ctx$X), c(n, p))
  expect_equal(dim(ctx$Y), c(n, 10))
  expect_null(ctx$proj)  # Not pre-computed
  
  # Create context with projection
  proj <- .fast_preproject(X)
  ctx2 <- glm_context(X = X, Y = Y, proj = proj)
  
  expect_equal(ctx2$proj$dfres, n - qr(X)$rank)
  expect_equal(dim(ctx2$proj$XtXinv), c(p, p))
})

test_that("solve_glm_core produces correct results", {
  # Simple regression test
  n <- 50
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  true_beta <- c(2, -1, 0.5)
  Y <- X %*% true_beta + rnorm(n, sd = 0.1)
  
  # Create context and solve
  proj <- .fast_preproject(X)
  ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  result <- solve_glm_core(ctx)
  
  # Check dimensions
  expect_equal(dim(result$betas), c(p, 1))
  expect_length(result$sigma2, 1)
  expect_length(result$rss, 1)
  
  # Check accuracy (should be close to true values)
  expect_equal(as.vector(result$betas), true_beta, tolerance = 0.2)
  
  # Compare with lm()
  lm_fit <- lm(Y ~ X - 1)
  expect_equal(as.vector(result$betas), unname(coef(lm_fit)), tolerance = 1e-10)
  expect_equal(sqrt(result$sigma2), summary(lm_fit)$sigma, tolerance = 1e-10)
})

test_that("solve_glm_core handles multiple responses", {
  n <- 40
  p <- 4
  v <- 20  # number of voxels
  
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  B <- matrix(rnorm(p * v), p, v)  # Different betas for each voxel
  Y <- X %*% B + matrix(rnorm(n * v, sd = 0.5), n, v)
  
  # Solve
  proj <- .fast_preproject(X)
  ctx <- glm_context(X = X, Y = Y, proj = proj)
  result <- solve_glm_core(ctx, return_fitted = TRUE)
  
  # Check dimensions
  expect_equal(dim(result$betas), c(p, v))
  expect_length(result$sigma2, v)
  expect_length(result$rss, v)
  expect_equal(dim(result$fitted), c(n, v))
  
  # Check accuracy for first voxel
  lm1 <- lm(Y[,1] ~ X - 1)
  expect_equal(as.vector(result$betas[,1]), unname(coef(lm1)), tolerance = 1e-10)
})

test_that(".fast_preproject handles edge cases", {
  # Square matrix
  n <- 10
  X <- diag(n)
  proj <- .fast_preproject(X)
  
  expect_equal(proj$dfres, 0)
  expect_equal(proj$XtXinv, diag(n))
  
  # Rank deficient matrix
  n <- 20
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  X[, 5] <- X[, 1] + X[, 2]  # Make rank deficient
  
  # Should not error
  expect_error(proj <- .fast_preproject(X), NA)
  expect_true(qr(X)$rank < p)
})

test_that("solve_glm_core with weights works correctly", {
  # Test weighted least squares
  n <- 50
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  true_beta <- c(1, -2, 3)
  Y <- X %*% true_beta + rnorm(n)
  
  # Create weights
  weights <- runif(n, 0.1, 2)
  sqrtw <- sqrt(weights)
  
  # Weight the data
  Xw <- X * sqrtw
  Yw <- Y * sqrtw
  
  # Solve
  proj_w <- .fast_preproject(Xw)
  ctx_w <- glm_context(X = Xw, Y = matrix(Yw, ncol = 1), proj = proj_w)
  result_w <- solve_glm_core(ctx_w)
  
  # Compare with lm() with weights
  lm_w <- lm(Y ~ X - 1, weights = weights)
  expect_equal(as.vector(result_w$betas), unname(coef(lm_w)), tolerance = 1e-10)
})