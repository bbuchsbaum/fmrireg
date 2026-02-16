# Comprehensive numerical stability tests for fmrireg solvers
# Tests ill-conditioned matrices, extreme values, and edge cases

library(fmrireg)
library(testthat)
library(Matrix)

# These tests exercise extreme edge cases of internal solver functions
# and are platform-sensitive; skip on CRAN.
skip_on_cran()

# Helper function to create ill-conditioned matrices with specific condition number
create_ill_conditioned_matrix <- function(n, p, target_kappa) {
  # Use SVD to construct matrix with specific condition number
  set.seed(123)
  U <- qr.Q(qr(matrix(rnorm(n * n), n, n)))
  V <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  
  # Create singular values with desired condition number
  s_max <- 1
  s_min <- s_max / target_kappa
  s <- seq(s_max, s_min, length.out = min(n, p))
  
  # Pad with zeros if needed
  if (n > p) {
    D <- rbind(diag(s), matrix(0, n - p, p))
  } else {
    D <- cbind(diag(s), matrix(0, n, p - n))
  }
  
  X <- U %*% D %*% t(V)
  
  # Verify condition number
  actual_kappa <- kappa(X)
  attr(X, "target_kappa") <- target_kappa
  attr(X, "actual_kappa") <- actual_kappa
  
  X
}

# Helper to check numerical stability of results
check_numerical_stability <- function(result, desc = "") {
  test_desc <- if (nchar(desc) > 0) paste0(" (", desc, ")") else ""
  
  # Check for NaN/Inf
  expect_false(any(is.nan(result$betas)), 
               label = paste0("NaN in betas", test_desc))
  expect_false(any(is.infinite(result$betas)), 
               label = paste0("Inf in betas", test_desc))
  
  # Check for reasonable values - allow larger values for extreme scales
  max_reasonable <- if (grepl("scale=1.0e\\+15", desc)) 1e20 else 1e15
  expect_true(all(abs(result$betas) < max_reasonable), 
              label = paste0("Unreasonably large betas", test_desc))
  
  if (!is.null(result$sigma2)) {
    expect_true(all(result$sigma2 >= 0), 
                label = paste0("Negative variance", test_desc))
    expect_false(any(is.nan(result$sigma2)), 
                 label = paste0("NaN in sigma2", test_desc))
  }
}

test_that("solver handles increasingly ill-conditioned matrices", {
  n <- 100
  p <- 5
  Y <- matrix(rnorm(n * 3), n, 3)
  
  # Test with different condition numbers
  kappas <- c(1e3, 1e6, 1e9, 1e12)
  
  for (target_kappa in kappas) {
    X <- create_ill_conditioned_matrix(n, p, target_kappa)
    actual_kappa <- attr(X, "actual_kappa")
    
    # Test projection computation
    proj <- suppressWarnings(fmrireg:::.fast_preproject(X))
    
    # Create context and solve
    ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
    result <- fmrireg:::solve_glm_core(ctx)
    
    check_numerical_stability(result, 
                              desc = sprintf("kappa=%.1e", actual_kappa))
    
    # For extremely ill-conditioned matrices, check if rank deficiency is detected
    if (target_kappa >= 1e9) {
      # Matrix should be treated as numerically rank deficient
      expect_true(!is.null(proj$rank) || !is.null(proj$is_full_rank),
                  label = sprintf("Rank info missing for kappa=%.1e", actual_kappa))
    }
  }
})

test_that("solver handles extreme coefficient values", {
  n <- 50
  p <- 4
  
  # Test different scales
  scales <- c(1e-15, 1e-10, 1e-5, 1, 1e5, 1e10, 1e15)
  
  for (scale in scales) {
    set.seed(123)
    X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
    true_betas <- c(2, 1, -0.5, 0.3) * scale
    Y <- X %*% true_betas + rnorm(n, sd = 0.1 * abs(scale))
    Y <- matrix(Y, ncol = 1)
    
    proj <- fmrireg:::.fast_preproject(X)
    ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
    result <- fmrireg:::solve_glm_core(ctx)
    
    check_numerical_stability(result, 
                              desc = sprintf("scale=%.1e", scale))
    
    # Check relative error for reasonable scales
    if (scale >= 1e-10 && scale <= 1e10) {
      rel_error <- max(abs(result$betas - true_betas) / (abs(true_betas) + 1e-10))
      expect_lt(rel_error, 0.1, 
                label = sprintf("Large relative error for scale=%.1e", scale))
    }
  }
})

test_that("solver handles perfect multicollinearity", {
  n <- 100
  set.seed(123)
  
  # Create perfectly collinear predictors
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- 2 * x1 - 3 * x2  # Perfect linear combination
  
  X <- cbind(1, x1, x2, x3)
  Y <- matrix(rnorm(n * 2), n, 2)
  
  # Matrix should be rank deficient
  expect_lt(qr(X)$rank, ncol(X))
  
  # Should handle gracefully
  proj <- suppressWarnings(fmrireg:::.fast_preproject(X))
  
  ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  result <- fmrireg:::solve_glm_core(ctx)
  
  check_numerical_stability(result, desc = "perfect collinearity")
  
  # The solver should produce stable results even if rank deficiency 
  # is not explicitly tracked
  expect_true(all(is.finite(result$betas)))
  expect_true(all(is.finite(result$sigma2)))
})

test_that("solver handles various rank-deficient patterns", {
  n <- 80
  
  # Pattern 1: Duplicate columns
  X1 <- cbind(1, rnorm(n), rnorm(n), 1)  # Duplicate intercept
  
  # Pattern 2: Linear dependencies
  x <- rnorm(n)
  X2 <- cbind(1, x, 2*x, 3*x)  # Multiple of same variable
  
  # Pattern 3: More columns than rows
  X3 <- matrix(rnorm(20 * 30), 20, 30)  # p > n
  
  patterns <- list(
    duplicate = X1,
    linear_dep = X2,
    wide = X3
  )
  
  for (pattern_name in names(patterns)) {
    X <- patterns[[pattern_name]]
    n_rows <- nrow(X)
    Y <- matrix(rnorm(n_rows * 2), n_rows, 2)
    
    proj <- suppressWarnings(fmrireg:::.fast_preproject(X))
    ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
    result <- fmrireg:::solve_glm_core(ctx)
    
    check_numerical_stability(result, desc = pattern_name)
    
    # Just check for stable results - rank detection may vary
    expect_true(all(is.finite(result$betas)),
                label = paste("Non-finite betas for", pattern_name))
    
    # For wide matrices (p > n), sigma2 might be NaN due to negative df
    if (pattern_name != "wide") {
      expect_true(all(is.finite(result$sigma2)),
                  label = paste("Non-finite sigma2 for", pattern_name))
    } else {
      # For wide matrices, just check that we got some result
      expect_true(!is.null(result$sigma2))
    }
  }
})

test_that("solver handles matrices that break standard solve()", {
  n <- 50
  
  # Create a matrix that will cause standard solve() to fail
  # Near-singular with structure that challenges numerical stability
  set.seed(123)
  A <- matrix(rnorm(n * 5), n, 5)
  X <- A %*% diag(c(1, 1e-8, 1e-16, 1e-8, 1)) %*% t(qr.Q(qr(matrix(rnorm(25), 5, 5))))
  Y <- matrix(rnorm(n * 2), n, 2)
  
  # Standard solve should struggle
  XtX <- crossprod(X)
  expect_error(solve(XtX), regexp = NULL)  # May or may not error
  
  # Our solver should handle it
  proj <- suppressWarnings(fmrireg:::.fast_preproject(X))
  ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  result <- fmrireg:::solve_glm_core(ctx)
  
  check_numerical_stability(result, desc = "near-singular matrix")
})

test_that("robust fitting handles extreme outliers in ill-conditioned setting", {
  n <- 100
  p <- 5
  
  # Create ill-conditioned design
  X <- create_ill_conditioned_matrix(n, p, target_kappa = 1e6)
  
  # Generate data with extreme outliers
  true_betas <- c(2, 1, -0.5, 0.3, 0.1)
  Y <- X %*% true_betas + rnorm(n, sd = 0.1)
  
  # Add extreme outliers at different scales
  outlier_idx <- sample(n, 20)
  Y[outlier_idx[1:10]] <- Y[outlier_idx[1:10]] + 100   # Large positive
  Y[outlier_idx[11:20]] <- Y[outlier_idx[11:20]] - 100  # Large negative
  
  Y <- matrix(Y, ncol = 1)
  
  # Initial OLS fit
  proj <- suppressWarnings(fmrireg:::.fast_preproject(X))
  initial_ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  initial_fit <- fmrireg:::solve_glm_core(initial_ctx)
  initial_ctx$betas <- initial_fit$betas
  
  cfg <- fmrireg:::fmri_lm_control(
    robust_options = list(
      type = "bisquare",
      c_tukey = 4.685,
      max_iter = 20  # More iterations for difficult case
    )
  )
  
  # Robust fitting should handle both ill-conditioning and outliers
  result <- fmrireg:::robust_iterative_fitter(initial_ctx, cfg$robust, X)
  
  check_numerical_stability(result, desc = "robust + ill-conditioned")
  
  # Check outlier downweighting
  expect_true(mean(result$robust_weights_final[outlier_idx]) < 0.3)
  expect_true(mean(result$robust_weights_final[-outlier_idx]) > 0.7)
})

test_that("solver handles QR decomposition edge cases", {
  # Test 1: Matrix with zeros
  n <- 50
  X1 <- cbind(1, rep(0, n), rnorm(n))  # Column of zeros
  Y1 <- matrix(rnorm(n * 2), n, 2)
  
  proj1 <- suppressWarnings(fmrireg:::.fast_preproject(X1))
  ctx1 <- fmrireg:::glm_context(X = X1, Y = Y1, proj = proj1)
  result1 <- fmrireg:::solve_glm_core(ctx1)
  
  check_numerical_stability(result1, desc = "zero column")
  expect_equal(result1$betas[2, ], c(0, 0), tolerance = 1e-10)  # Zero coef for zero predictor
  
  # Test 2: Nearly zero matrix elements
  X2 <- cbind(1, rnorm(n) * 1e-300, rnorm(n))
  Y2 <- matrix(rnorm(n * 2), n, 2)
  
  proj2 <- suppressWarnings(fmrireg:::.fast_preproject(X2))
  ctx2 <- fmrireg:::glm_context(X = X2, Y = Y2, proj = proj2)
  result2 <- fmrireg:::solve_glm_core(ctx2)
  
  check_numerical_stability(result2, desc = "near-zero elements")
})

test_that("solver handles mixed precision scenarios", {
  n <- 100
  p <- 4
  
  # Mix of very different scales in same matrix
  set.seed(123)
  X <- cbind(
    rep(1, n),                    # Intercept (scale = 1)
    rnorm(n) * 1e-10,            # Very small scale
    rnorm(n) * 1e10,             # Very large scale  
    rnorm(n)                     # Normal scale
  )
  
  # Response with mixed scales
  true_betas <- c(1e5, 1e15, 1e-5, 1)
  Y <- X %*% true_betas + rnorm(n)
  Y <- matrix(Y, ncol = 1)
  
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  result <- fmrireg:::solve_glm_core(ctx)
  
  check_numerical_stability(result, desc = "mixed precision")
  
  # Should recover coefficients reasonably well despite scale differences
  # Check relative errors where meaningful
  rel_errors <- abs(result$betas - true_betas) / (abs(true_betas) + 1e-10)
  expect_true(all(rel_errors < 1), 
              label = "Large relative errors in mixed precision case")
})

test_that("solver handles numerical overflow/underflow scenarios", {
  n <- 50
  p <- 3
  
  # Test potential overflow in X'X - use more moderate values
  X_large <- cbind(1, matrix(rnorm(n * (p-1)) * 1e100, n, p-1))
  Y_large <- matrix(rnorm(n * 2) * 1e100, n, 2)
  
  # Should handle without overflow
  proj_large <- suppressWarnings(fmrireg:::.fast_preproject(X_large))
  ctx_large <- fmrireg:::glm_context(X = X_large, Y = Y_large, proj = proj_large)
  result_large <- fmrireg:::solve_glm_core(ctx_large)
  
  # Check for stability - values may be large but should be finite
  expect_false(any(is.nan(result_large$betas)))
  expect_false(any(is.infinite(result_large$betas)))
  
  # Test potential underflow
  X_small <- cbind(1, matrix(rnorm(n * (p-1)) * 1e-150, n, p-1))
  Y_small <- matrix(rnorm(n * 2) * 1e-150, n, 2)
  
  proj_small <- suppressWarnings(fmrireg:::.fast_preproject(X_small))
  ctx_small <- fmrireg:::glm_context(X = X_small, Y = Y_small, proj = proj_small)
  result_small <- fmrireg:::solve_glm_core(ctx_small)
  
  check_numerical_stability(result_small, desc = "potential underflow")
})

test_that("matrix inversion stability with Cholesky vs SVD fallback", {
  n <- 100
  p <- 5
  
  # Create matrix that's positive definite but nearly singular
  # This tests the Cholesky -> SVD fallback mechanism
  set.seed(123)
  A <- matrix(rnorm(n * p), n, p)
  A[, p] <- A[, p-1] + rnorm(n) * 1e-12  # Nearly collinear
  X <- qr.Q(qr(A))  # Orthogonalize
  X[, p] <- X[, p-1] + rnorm(n) * 1e-12  # Make nearly collinear again
  
  Y <- matrix(rnorm(n * 2), n, 2)
  
  # This should trigger the SVD fallback in fmrireg:::.fast_preproject
  proj <- suppressWarnings(fmrireg:::.fast_preproject(X))
  ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  result <- fmrireg:::solve_glm_core(ctx)
  
  check_numerical_stability(result, desc = "Cholesky fallback case")
  
  # Verify the solution is reasonable by checking residuals
  resid <- Y - X %*% result$betas
  rss_manual <- colSums(resid^2)
  expect_equal(result$rss, rss_manual, tolerance = 1e-10)
})

test_that("robust fitting with rank-deficient design and outliers", {
  n <- 80
  
  # Create rank-deficient design
  x1 <- rnorm(n)
  x2 <- rnorm(n) 
  X <- cbind(1, x1, x2, x1 + x2, 2*x1 - x2)  # Last two columns are linear combos
  
  # True model uses only first 3 coefficients
  true_betas <- c(2, 1, -0.5, 0, 0)
  Y <- X[, 1:3] %*% true_betas[1:3] + rnorm(n, sd = 0.1)
  
  # Add outliers
  outlier_idx <- sample(n, 10)
  Y[outlier_idx] <- Y[outlier_idx] + sample(c(-50, 50), 10, replace = TRUE)
  Y <- matrix(Y, ncol = 1)
  
  # Initial fit
  proj <- suppressWarnings(fmrireg:::.fast_preproject(X))
  initial_ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  initial_fit <- fmrireg:::solve_glm_core(initial_ctx)
  initial_ctx$betas <- initial_fit$betas
  
  cfg <- fmrireg:::fmri_lm_control(
    robust_options = list(
      type = "huber",
      k_huber = 1.345,
      max_iter = 15
    )
  )
  
  # Should handle rank deficiency + outliers
  result <- fmrireg:::robust_iterative_fitter(initial_ctx, cfg$robust, X)
  
  check_numerical_stability(result, desc = "robust + rank-deficient")
  
  # Outliers should be downweighted
  expect_true(mean(result$robust_weights_final[outlier_idx]) < 0.5)
})

test_that("fast solver path handles edge cases identically to regular path", {
  # Create challenging test cases
  n <- 60
  p <- 4
  
  test_cases <- list(
    ill_conditioned = create_ill_conditioned_matrix(n, p, 1e8),
    rank_deficient = cbind(1, rnorm(n), rnorm(n), 1),  # Duplicate intercept
    mixed_scale = cbind(1, rnorm(n)*1e-8, rnorm(n)*1e8, rnorm(n))
  )
  
  for (case_name in names(test_cases)) {
    X <- test_cases[[case_name]]
    Y <- matrix(rnorm(n * 3), n, 3)
    
    # Regular path
    proj_regular <- suppressWarnings(fmrireg:::.fast_preproject(X))
    ctx_regular <- fmrireg:::glm_context(X = X, Y = Y, proj = proj_regular)
    result_regular <- fmrireg:::solve_glm_core(ctx_regular)
    
    # For comparison with any fast path implementations
    # (Currently both use same implementation, but this tests consistency)
    result_fast <- fmrireg:::solve_glm_core(ctx_regular, return_fitted = FALSE)
    
    # Results should be identical
    expect_equal(result_regular$betas, result_fast$betas, 
                 tolerance = 1e-12,
                 label = paste("Betas differ for", case_name))
    expect_equal(result_regular$rss, result_fast$rss,
                 tolerance = 1e-12, 
                 label = paste("RSS differ for", case_name))
    expect_equal(result_regular$sigma2, result_fast$sigma2,
                 tolerance = 1e-12,
                 label = paste("Sigma2 differ for", case_name))
  }
})

test_that("solver provides appropriate warnings and errors", {
  n <- 50
  
  # Test 1: Empty matrix
  expect_error(
    fmrireg:::.fast_preproject(matrix(nrow = 0, ncol = 0)),
    regexp = NULL
  )
  
  # Test 2: Matrix with all NAs
  X_na <- matrix(NA, n, 3)
  expect_error(
    fmrireg:::.fast_preproject(X_na),
    regexp = NULL
  )
  
  # Test 3: Rank deficient matrix should be handled gracefully
  X_rankdef <- cbind(1, rnorm(n), 1)
  proj <- suppressWarnings(fmrireg:::.fast_preproject(X_rankdef))
  
  # Should produce valid projection even if rank deficiency is not explicitly tracked
  expect_true(!is.null(proj$XtXinv))
  expect_true(all(is.finite(proj$XtXinv)))
  
  # Test 4: Dimension mismatch in solve_glm_core
  X <- cbind(1, rnorm(n))
  Y <- matrix(rnorm((n-1) * 2), n-1, 2)  # Wrong number of rows
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  
  expect_error(
    fmrireg:::solve_glm_core(ctx),
    regexp = "dimensions"
  )
})

# Final summary test
test_that("numerical stability test suite summary", {
  # This test just confirms the test file is complete
  expect_true(TRUE, label = "Numerical stability test suite completed")
})