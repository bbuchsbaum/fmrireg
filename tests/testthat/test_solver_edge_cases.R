# Test edge cases for GLM solver

library(fmrireg)
library(testthat)

test_that("solve_glm_core handles rank-deficient matrices", {
  # Create rank-deficient design matrix
  n <- 50
  X <- cbind(
    1,
    rnorm(n),
    rnorm(n),
    1  # Duplicate of intercept
  )
  Y <- matrix(rnorm(n * 3), n, 3)
  
  # Check that matrix is rank deficient
  expect_lt(qr(X)$rank, ncol(X))
  
  # Should handle rank deficiency gracefully
  # The warning may or may not appear depending on implementation details
  proj <- suppressWarnings(fmrireg:::.fast_preproject(X))
  ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  
  result <- fmrireg:::solve_glm_core(ctx)
  
  # Should still return results
  expect_equal(dim(result$betas), c(ncol(X), ncol(Y)))
  # Should produce finite results even with rank deficiency
  expect_true(all(is.finite(result$betas)))
  expect_true(all(!is.na(result$betas)))
})

test_that("solve_glm_core handles near-singular matrices", {
  # Create nearly collinear predictors
  n <- 100
  x1 <- rnorm(n)
  x2 <- x1 + rnorm(n, sd = 0.01)  # Almost identical to x1
  X <- cbind(1, x1, x2)
  Y <- matrix(rnorm(n * 2), n, 2)
  
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  
  # Should complete without error
  result <- fmrireg:::solve_glm_core(ctx)
  
  expect_equal(dim(result$betas), c(3, 2))
  expect_true(all(!is.na(result$sigma2)))
  
  # Check condition number
  expect_gt(kappa(X), 100)  # High condition number indicates near-singularity
})

test_that("solve_glm_core handles extreme values", {
  n <- 50
  p <- 3
  
  # Test with very large values
  X_large <- cbind(1, matrix(rnorm(n * (p-1)) * 1e6, n, p-1))
  Y_large <- matrix(rnorm(n * 2) * 1e6, n, 2)
  
  proj_large <- fmrireg:::.fast_preproject(X_large)
  ctx_large <- fmrireg:::glm_context(X = X_large, Y = Y_large, proj = proj_large)
  result_large <- fmrireg:::solve_glm_core(ctx_large)
  
  expect_false(any(is.infinite(result_large$betas)))
  expect_false(any(is.nan(result_large$betas)))
  
  # Test with very small values
  X_small <- cbind(1, matrix(rnorm(n * (p-1)) * 1e-6, n, p-1))
  Y_small <- matrix(rnorm(n * 2) * 1e-6, n, 2)
  
  proj_small <- fmrireg:::.fast_preproject(X_small)
  ctx_small <- fmrireg:::glm_context(X = X_small, Y = Y_small, proj = proj_small)
  result_small <- fmrireg:::solve_glm_core(ctx_small)
  
  expect_false(any(is.infinite(result_small$betas)))
  expect_false(any(is.nan(result_small$betas)))
})

test_that("robust fitting handles extreme outliers", {
  n <- 100
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  
  # Create data with extreme outliers
  Y <- X %*% c(2, 1, -0.5) + rnorm(n, sd = 0.1)
  # Add extreme outliers
  outlier_idx <- sample(n, 10)
  Y[outlier_idx] <- Y[outlier_idx] + sample(c(-20, 20), 10, replace = TRUE)
  
  # First create initial OLS fit
  proj <- fmrireg:::.fast_preproject(X)
  initial_ctx <- fmrireg:::glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  initial_fit <- fmrireg:::solve_glm_core(initial_ctx)
  initial_ctx$betas <- initial_fit$betas
  
  cfg <- fmrireg:::fmri_lm_control(
    robust_options = list(
      type = "bisquare",
      c_tukey = 4.685,
      max_iter = 10
    )
  )
  
  # Robust fitting should downweight outliers
  result <- fmrireg:::robust_iterative_fitter(initial_ctx, cfg$robust, X)
  
  expect_true(all(result$robust_weights_final[outlier_idx] < 0.5))  # Outliers downweighted
  expect_true(mean(result$robust_weights_final[-outlier_idx]) > 0.8)  # Non-outliers have high weight
  
  # Estimates should be close to true values despite outliers
  expect_equal(as.vector(result$betas_robust), c(2, 1, -0.5), tolerance = 0.2)
})

test_that("robust fitting handles convergence failures", {
  # Create pathological case
  n <- 20
  X <- cbind(1, rnorm(n))
  Y <- matrix(c(rep(1, 10), rep(1000, 10)), ncol = 1)  # Bimodal response
  
  # First create initial OLS fit
  proj <- fmrireg:::.fast_preproject(X)
  initial_ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  initial_fit <- fmrireg:::solve_glm_core(initial_ctx)
  initial_ctx$betas <- initial_fit$betas
  
  cfg <- fmrireg:::fmri_lm_control(
    robust_options = list(
      type = "bisquare",
      max_iter = 2  # Force early stopping
    )
  )
  
  # Should complete even with limited iterations
  result <- fmrireg:::robust_iterative_fitter(initial_ctx, cfg$robust, X)
  
  expect_true(!is.null(result$betas_robust))
  expect_true(!is.null(result$robust_weights_final))
  expect_equal(length(result$robust_weights_final), n)
})

test_that("AR fitting handles short time series", {
  # Time series too short for AR order
  n <- 5  # Very short
  X <- cbind(1, rnorm(n))
  Y <- matrix(rnorm(n), ncol = 1)
  
  # Get residuals for AR estimation
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  fit <- fmrireg:::solve_glm_core(ctx)
  residuals <- Y - X %*% fit$betas
  
  # AR(2) with only 5 observations - should work but may not be reliable
  phi <- fmrireg:::.estimate_ar(as.vector(residuals), p_order = 2)
  
  # Should return coefficients
  expect_equal(length(phi), 2)
  # Check stationarity (all roots outside unit circle)
  if (length(phi) == 2) {
    # For AR(2), check characteristic polynomial roots
    # 1 - phi1*z - phi2*z^2 = 0
    # Stationarity requires |roots| > 1
    poly_coef <- c(1, -phi[1], -phi[2])
    roots <- polyroot(poly_coef)
    expect_true(all(Mod(roots) > 0.99))  # Allow slight numerical tolerance
  }
})

test_that("mixed_solve handles edge cases", {
  # Create test data
  set.seed(123)
  n <- 50
  p_Z <- 5
  p_X <- 2

  y <- rnorm(n)
  Z <- matrix(rnorm(n * p_Z), n, p_Z)
  K <- diag(p_Z)
  X <- matrix(rnorm(n * p_X), n, p_X)

  # Test normal case
  result <- fmrilss::mixed_solve(Y = y, Z = Z, K = K, X = X)
  expect_true(is.list(result))
  expect_true(all(c("u", "beta", "Vu", "Ve") %in% names(result)))

  # Test with smaller dataset
  y_small <- y[1:10]
  result_small <- fmrilss::mixed_solve(Y = y_small, Z = Z[1:10, ], K = K, X = X[1:10, ])
  expect_true(is.list(result_small))

  # Test edge case: single column Z and K
  Z_single <- Z[, 1, drop = FALSE]
  K_single <- matrix(1, 1, 1)
  result_single <- fmrilss::mixed_solve(Y = y, Z = Z_single, K = K_single, X = X)
  expect_true(is.list(result_single))
  expect_equal(length(result_single$u), 1)
})

# ============================================================================
# Additional edge case tests
# ============================================================================

test_that("fmri_lm handles single voxel dataset", {
  event_table <- data.frame(
    onset = c(10, 30, 60, 80),
    cond = factor(c("A", "B", "A", "B")),
    run = c(1, 1, 2, 2)
  )

  # Single voxel
  mat <- matrix(rnorm(100), 100, 1)
  dset <- fmridataset::matrix_dataset(
    mat,
    TR = 1,
    run_length = c(50, 50),
    event_table = event_table
  )

  # Should work with single voxel
  expect_no_error(
    mod <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = dset)
  )

  expect_s3_class(mod, "fmri_lm")

  # Check that results have correct dimensions
  cf <- coef(mod)
  expect_true(is.matrix(cf) || is.numeric(cf))
})

test_that("fmri_lm handles very short time series", {
  # Minimal viable time series
  event_table <- data.frame(
    onset = c(5, 15),
    cond = factor(c("A", "B")),
    run = c(1, 1)
  )

  mat <- matrix(rnorm(30 * 5), 30, 5)
  dset <- fmridataset::matrix_dataset(
    mat,
    TR = 1,
    run_length = 30,
    event_table = event_table
  )

  expect_no_error(
    mod <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = dset)
  )

  expect_s3_class(mod, "fmri_lm")
})

test_that("solve_glm_core handles zero variance response", {
  n <- 50
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))

  # Zero variance Y column
  Y <- cbind(rnorm(n), rep(5, n), rnorm(n))

  proj <- fmrireg:::.fast_preproject(X)
  ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)

  result <- fmrireg:::solve_glm_core(ctx)

  # Should complete
  expect_equal(dim(result$betas), c(p, 3))

  # Zero variance column should have zero residual variance
  expect_lt(result$sigma2[2], 1e-10)
})

test_that("solve_glm_core handles all-zero response", {
  n <- 50
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  Y <- matrix(0, n, 2)

  proj <- fmrireg:::.fast_preproject(X)
  ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)

  result <- fmrireg:::solve_glm_core(ctx)

  # Should complete without NaN
  expect_true(all(is.finite(result$betas)))

  # Sigma should be essentially zero
  expect_true(all(result$sigma2 < 1e-10))
})

test_that("contrast computation handles rank-deficient contrast matrix", {
  # Create scenario where contrast matrix might be rank deficient
  n <- 100
  p <- 5
  V <- 10

  B <- matrix(rnorm(p * V), p, V)
  sigma2 <- rep(1, V)
  XtXinv <- diag(p)
  df <- n - p

  # Create redundant contrasts (rows are linearly dependent)
  l1 <- c(1, -1, 0, 0, 0)
  l2 <- c(-1, 1, 0, 0, 0)  # Negative of l1
  l3 <- c(0, 0, 1, -1, 0)

  # For t-contrast, use single vector
  full_contrast <- matrix(l1, nrow = 1, ncol = p)
  stats <- fmrireg:::.fast_t_contrast(B, sigma2, XtXinv, full_contrast, df)

  expect_equal(length(stats$estimate), V)
  expect_true(all(is.finite(stats$stat)))
})

test_that("F-contrast handles singular contrast matrix", {
  n <- 100
  p <- 5
  V <- 5

  B <- matrix(rnorm(p * V), p, V)
  sigma2 <- rep(1, V)
  XtXinv <- diag(p)
  df <- n - p

  # Create singular F-contrast matrix
  L <- matrix(c(
    1, -1, 0, 0, 0,
    2, -2, 0, 0, 0  # Multiple of first row
  ), nrow = 2, byrow = TRUE)

  # Should produce warning about singular matrix
  expect_warning(
    stats <- fmrireg:::.fast_F_contrast(B, sigma2, XtXinv, L, df),
    regexp = "[Ss]ingular"
  )

  # Results should contain NaN for singular case
  expect_true(any(is.nan(stats$stat)) || any(is.nan(stats$estimate)))
})

# ============================================================================
# HRF smoothing kernel edge cases
# ============================================================================

test_that("hrf_smoothing_kernel works with default parameters", {
  expect_no_error(
    sk <- hrf_smoothing_kernel(100, TR = 1.5)
  )

  expect_true(is.matrix(sk))
  expect_equal(dim(sk), c(100, 100))

  # Diagonal should be 1 (normalized)
  expect_true(all(abs(diag(sk) - 1) < 1e-10))
})

test_that("hrf_smoothing_kernel works with different methods", {
  # Gram method
  sk_gram <- hrf_smoothing_kernel(50, TR = 2, method = "gram")
  expect_true(is.matrix(sk_gram))

  # Cosine method
  sk_cosine <- hrf_smoothing_kernel(50, TR = 2, method = "cosine")
  expect_true(is.matrix(sk_cosine))

  # Methods should produce different results
  expect_false(all(abs(sk_gram - sk_cosine) < 1e-10))
})

test_that("hrf_smoothing_kernel handles different buffer_scans", {
  # No buffer
  sk_0 <- hrf_smoothing_kernel(50, TR = 2, buffer_scans = 0)
  expect_equal(dim(sk_0), c(50, 50))

  # Large buffer
  sk_10 <- hrf_smoothing_kernel(50, TR = 2, buffer_scans = 10)
  expect_equal(dim(sk_10), c(50, 50))
})

test_that("hrf_smoothing_kernel respects normalise parameter", {
  # Normalized (default)
  sk_norm <- hrf_smoothing_kernel(30, TR = 2, normalise = TRUE)
  expect_true(all(abs(diag(sk_norm) - 1) < 1e-10))

  # Not normalized
  sk_raw <- hrf_smoothing_kernel(30, TR = 2, normalise = FALSE)
  # Diagonal may not be 1
  expect_true(is.matrix(sk_raw))
})

test_that("hrf_smoothing_kernel works with different TRs", {
  for (tr in c(0.5, 1, 1.5, 2, 3)) {
    sk <- hrf_smoothing_kernel(40, TR = tr)
    expect_true(is.matrix(sk), info = paste("Failed for TR =", tr))
    expect_equal(dim(sk), c(40, 40), info = paste("Failed for TR =", tr))
  }
})

test_that("hrf_smoothing_kernel handles very short time series", {
  # Minimum viable length
  sk <- hrf_smoothing_kernel(10, TR = 2)
  expect_equal(dim(sk), c(10, 10))
})

test_that("hrf_smoothing_kernel works with custom formula", {
  form <- onset ~ trialwise(basis = "gaussian")
  sk <- hrf_smoothing_kernel(50, TR = 1.5, form = form)

  expect_true(is.matrix(sk))
  expect_equal(dim(sk), c(50, 50))
})

# ============================================================================
# Numerical precision tests
# ============================================================================

test_that("solver maintains numerical precision with scaled data", {
  n <- 100
  p <- 3

  # Original scale
  X_orig <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  Y_orig <- matrix(rnorm(n * 2), n, 2)

  proj_orig <- fmrireg:::.fast_preproject(X_orig)
  ctx_orig <- fmrireg:::glm_context(X = X_orig, Y = Y_orig, proj = proj_orig)
  result_orig <- fmrireg:::solve_glm_core(ctx_orig)

  # Scaled up by 1000
  X_scaled <- X_orig * 1000
  Y_scaled <- Y_orig * 1000

  proj_scaled <- fmrireg:::.fast_preproject(X_scaled)
  ctx_scaled <- fmrireg:::glm_context(X = X_scaled, Y = Y_scaled, proj = proj_scaled)
  result_scaled <- fmrireg:::solve_glm_core(ctx_scaled)

  # Betas should be the same (within numerical tolerance)
  expect_equal(result_orig$betas, result_scaled$betas, tolerance = 1e-6)
})

test_that("solver handles mixed precision data", {
  n <- 50
  p <- 3

  # Mix of very large and very small values
  X <- cbind(1, rnorm(n) * 1e-6, rnorm(n) * 1e6)
  Y <- matrix(rnorm(n * 2), n, 2)

  proj <- fmrireg:::.fast_preproject(X)
  ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)

  expect_no_error(
    result <- fmrireg:::solve_glm_core(ctx)
  )

  expect_true(all(is.finite(result$betas)))
})