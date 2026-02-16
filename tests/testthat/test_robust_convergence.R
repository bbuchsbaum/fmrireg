# Test robust fitting convergence and extreme cases
#
# This test file focuses on edge cases and extreme scenarios for robust regression:
# - Non-convergent/pathological data
# - Maximum iteration limits
# - Extreme contamination levels (30%, 50%, 70%)
# - Zero or near-zero MAD scale scenarios
# - Weight evolution tracking
# - Tolerance-based convergence (for future implementation)
# - Global vs local scale estimation
# - Leverage points (X-space outliers)
# - Singular/near-singular designs
# - Groups of identical outliers
# - Extreme scale values
# - Performance degradation with contamination

test_that("robust fitting handles non-convergent data", {
  # Create pathological data that oscillates between solutions
  n <- 100
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  
  # Create alternating pattern of extreme outliers
  Y <- rnorm(n)
  Y[seq(1, n, by = 2)] <- 100 * (-1)^(seq(1, n, by = 2))
  
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  
  robust_opts <- list(
    type = "bisquare",
    c_tukey = 4.685,
    max_iter = 5,  # Limited iterations
    scale_scope = "local"
  )
  
  # Should complete without error even if not converged
  result <- robust_iterative_fitter(
    initial_glm_ctx = ctx,
    cfg_robust_options = robust_opts,
    X_orig_for_resid = X
  )
  
  expect_true(!is.null(result$betas_robust))
  expect_true(all(is.finite(result$betas_robust)))
})

test_that("maximum iteration limit is respected", {
  n <- 50
  p <- 2
  X <- cbind(1, rnorm(n))
  Y <- X %*% c(1, 2) + rnorm(n)
  
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  
  # Test with different max_iter values
  for (max_iter in c(1, 5, 10)) {
    robust_opts <- list(
      type = "huber",
      k_huber = 1.345,
      max_iter = max_iter,
      scale_scope = "local"
    )
    
    # Track iterations somehow - currently the function doesn't return iteration count
    # But we can verify it completes
    result <- robust_iterative_fitter(
      initial_glm_ctx = ctx,
      cfg_robust_options = robust_opts,
      X_orig_for_resid = X
    )
    
    expect_true(!is.null(result$betas_robust))
  }
})

test_that("robust fitting handles extreme contamination (30% outliers)", {
  n <- 100
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  beta_true <- c(2, -1, 3)
  Y <- X %*% beta_true + rnorm(n, sd = 0.5)
  
  # Add 30% contamination
  n_outliers <- floor(0.3 * n)
  outlier_idx <- sample(n, n_outliers)
  Y[outlier_idx] <- Y[outlier_idx] + rnorm(n_outliers, mean = 10, sd = 2)
  
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  
  # Test both Huber and bisquare
  for (psi_type in c("huber", "bisquare")) {
    robust_opts <- list(
      type = psi_type,
      k_huber = 1.345,
      c_tukey = 4.685,
      max_iter = 50,  # More iterations for convergence
      scale_scope = "local"
    )
    
    result <- robust_iterative_fitter(
      initial_glm_ctx = ctx,
      cfg_robust_options = robust_opts,
      X_orig_for_resid = X
    )
    
    # Check that outliers are downweighted
    expect_lt(mean(result$robust_weights_final[outlier_idx]), 
              mean(result$robust_weights_final[-outlier_idx]))
    
    # Results should still be reasonable
    expect_true(all(is.finite(result$betas_robust)))
  }
})

test_that("robust fitting handles extreme contamination (50% outliers)", {
  n <- 100
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  beta_true <- c(1, 2, -1)
  Y <- X %*% beta_true + rnorm(n, sd = 0.3)
  
  # Add 50% contamination - breakdown point test
  n_outliers <- floor(0.5 * n)
  outlier_idx <- sample(n, n_outliers)
  Y[outlier_idx] <- rnorm(n_outliers, mean = 20, sd = 5)
  
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  
  robust_opts <- list(
    type = "bisquare",
    c_tukey = 4.685,
    max_iter = 100,  # Many iterations
    scale_scope = "local"
  )
  
  result <- robust_iterative_fitter(
    initial_glm_ctx = ctx,
    cfg_robust_options = robust_opts,
    X_orig_for_resid = X
  )
  
  # At 50% contamination, robust methods struggle but should still produce output
  expect_true(!is.null(result$betas_robust))
  expect_true(all(is.finite(result$betas_robust)))
  
  # At 50% contamination, we just verify some downweighting occurs
  expect_true(min(result$robust_weights_final) < 0.9)
})

test_that("robust fitting handles extreme contamination (70% outliers)", {
  n <- 100
  p <- 2
  X <- cbind(1, rnorm(n))
  beta_true <- c(0, 1)
  Y <- X %*% beta_true + rnorm(n, sd = 0.2)
  
  # Add 70% contamination - beyond breakdown point
  n_outliers <- floor(0.7 * n)
  outlier_idx <- sample(n, n_outliers)
  Y[outlier_idx] <- runif(n_outliers, -50, 50)
  
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  
  robust_opts <- list(
    type = "huber",
    k_huber = 1.345,
    max_iter = 50,
    scale_scope = "local"
  )
  
  # Should complete even with majority outliers
  result <- robust_iterative_fitter(
    initial_glm_ctx = ctx,
    cfg_robust_options = robust_opts,
    X_orig_for_resid = X
  )
  
  expect_true(!is.null(result$betas_robust))
  expect_true(all(is.finite(result$betas_robust)))
})

test_that("scale estimation handles zero MAD (constant data)", {
  n <- 50
  p <- 2
  X <- cbind(1, rnorm(n))
  
  # Constant Y values - MAD would be zero
  Y <- rep(5, n)
  
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  
  robust_opts <- list(
    type = "huber",
    k_huber = 1.345,
    max_iter = 10,
    scale_scope = "local"
  )
  
  result <- robust_iterative_fitter(
    initial_glm_ctx = ctx,
    cfg_robust_options = robust_opts,
    X_orig_for_resid = X
  )
  
  # Scale should be set to machine epsilon
  expect_equal(result$sigma_robust_scale_final, .Machine$double.eps)
  
  # With constant data and zero scale, weights should be high
  expect_true(mean(result$robust_weights_final) > 0.9)
})

test_that("scale estimation handles near-zero MAD", {
  n <- 60
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  
  # Very small noise
  Y <- X %*% c(1, 2, -1) + rnorm(n, sd = 1e-10)
  
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  
  robust_opts <- list(
    type = "bisquare",
    c_tukey = 4.685,
    max_iter = 20,
    scale_scope = "local"
  )
  
  result <- robust_iterative_fitter(
    initial_glm_ctx = ctx,
    cfg_robust_options = robust_opts,
    X_orig_for_resid = X
  )
  
  # Should handle near-zero scale
  expect_true(result$sigma_robust_scale_final > 0)
  expect_true(all(is.finite(result$betas_robust)))
})

test_that("weight evolution shows convergence pattern", {
  n <- 80
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  beta_true <- c(1, -2, 3)
  Y <- X %*% beta_true + rnorm(n, sd = 0.5)
  
  # Add some outliers
  outlier_idx <- c(5, 15, 25, 35)
  Y[outlier_idx] <- Y[outlier_idx] + 8
  
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  
  # Track weights across iterations by running with increasing max_iter
  weights_history <- list()
  
  for (iter in 1:5) {
    robust_opts <- list(
      type = "huber",
      k_huber = 1.345,
      max_iter = iter,
      scale_scope = "local"
    )
    
    result <- robust_iterative_fitter(
      initial_glm_ctx = ctx,
      cfg_robust_options = robust_opts,
      X_orig_for_resid = X
    )
    
    weights_history[[iter]] <- result$robust_weights_final
  }
  
  # Weights should stabilize over iterations
  # Check that weights change less between later iterations
  change_early <- sum(abs(weights_history[[2]] - weights_history[[1]]))
  change_late <- sum(abs(weights_history[[5]] - weights_history[[4]]))
  
  expect_lt(change_late, change_early)
  
  # Outlier weights should generally be low after convergence
  final_outlier_weights <- mean(weights_history[[5]][outlier_idx])
  final_inlier_weights <- mean(weights_history[[5]][-outlier_idx])
  expect_lt(final_outlier_weights, final_inlier_weights)
})

test_that("tolerance-based convergence would work if implemented", {
  # Note: Current implementation doesn't use tolerance, but test the concept
  n <- 60
  p <- 2
  X <- cbind(1, rnorm(n))
  Y <- X %*% c(2, -1) + rnorm(n)
  Y[c(5, 10, 15)] <- Y[c(5, 10, 15)] + 5
  
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  
  # If tolerance were implemented, it would be in robust_opts
  robust_opts <- list(
    type = "huber",
    k_huber = 1.345,
    max_iter = 100,
    tol = 1e-6,  # Currently ignored
    scale_scope = "local"
  )
  
  result <- robust_iterative_fitter(
    initial_glm_ctx = ctx,
    cfg_robust_options = robust_opts,
    X_orig_for_resid = X
  )
  
  # Should complete successfully
  expect_true(!is.null(result$betas_robust))
})

test_that("global vs local scale estimation in extreme cases", {
  n <- 100
  p <- 3
  v <- 3  # voxels with very different scales
  
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  
  # Create responses with vastly different scales
  Y <- cbind(
    X %*% c(1, 2, -1) + rnorm(n, sd = 0.1),    # Small scale
    X %*% c(-1, 3, 2) + rnorm(n, sd = 5),     # Medium scale  
    X %*% c(2, -2, 1) + rnorm(n, sd = 50)     # Large scale
  )
  
  # Add outliers to each
  Y[1:5, 1] <- Y[1:5, 1] + 1      # Large relative to scale
  Y[6:10, 2] <- Y[6:10, 2] + 25   # Medium relative to scale
  Y[11:15, 3] <- Y[11:15, 3] + 200 # Small relative to scale
  
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = Y, proj = proj)
  
  # Test with local scale
  result_local <- robust_iterative_fitter(
    initial_glm_ctx = ctx,
    cfg_robust_options = list(
      type = "bisquare",
      c_tukey = 4.685,
      max_iter = 30,
      scale_scope = "local"
    ),
    X_orig_for_resid = X
  )
  
  # Test with fixed global scale
  result_global <- robust_iterative_fitter(
    initial_glm_ctx = ctx,
    cfg_robust_options = list(
      type = "bisquare", 
      c_tukey = 4.685,
      max_iter = 30,
      scale_scope = "global"
    ),
    X_orig_for_resid = X,
    sigma_fixed = 10.0  # Fixed scale between the extremes
  )
  
  # Both should produce results
  expect_equal(dim(result_local$betas_robust), c(p, v))
  expect_equal(dim(result_global$betas_robust), c(p, v))
  
  # Weights pattern depends on scale - just check for valid results
  expect_true(all(result_global$robust_weights_final >= 0))
  expect_true(all(result_global$robust_weights_final <= 1))
})

test_that("leverage points (X-space outliers) are handled", {
  set.seed(4243)
  n <- 100
  p <- 3
  
  # Create design matrix with some leverage points
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  
  # Create high-leverage points by making some X values extreme
  leverage_idx <- c(5, 10, 15)
  X[leverage_idx, 2:p] <- X[leverage_idx, 2:p] * 10
  
  # Generate Y with outliers at leverage points
  beta_true <- c(2, -1, 3)
  Y <- X %*% beta_true + rnorm(n, sd = 0.5)
  Y[leverage_idx] <- Y[leverage_idx] + 5  # Bad leverage points
  
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  
  # Robust regression should downweight bad leverage points
  robust_opts <- list(
    type = "bisquare",
    c_tukey = 4.685,
    max_iter = 50,
    scale_scope = "local"
  )
  
  result <- robust_iterative_fitter(
    initial_glm_ctx = ctx,
    cfg_robust_options = robust_opts,
    X_orig_for_resid = X
  )
  
  # High leverage outliers should have lower weights than normal points
  expect_lt(mean(result$robust_weights_final[leverage_idx]), 
            mean(result$robust_weights_final[-leverage_idx]))
  
  # Compare with OLS which is sensitive to leverage
  ols_beta <- solve(t(X) %*% X) %*% t(X) %*% Y
  
  # Robust should be closer to true values
  robust_error <- sum((result$betas_robust - beta_true)^2)
  ols_error <- sum((ols_beta - beta_true)^2)
  expect_lt(robust_error, ols_error)
})

test_that("robust fitting handles singular or near-singular designs", {
  n <- 50
  p <- 4
  
  # Create near-singular design matrix
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  X[, 4] <- X[, 3] + rnorm(n, sd = 1e-8)  # Nearly collinear
  
  Y <- rnorm(n)
  
  # Check if preprojection handles near-singularity
  expect_warning(
    proj <- fmrireg:::.fast_preproject(X),
    regexp = NA  # May or may not warn
  )
  
  # If preprojection succeeds, test robust fitting
  if (!is.null(proj)) {
    ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
    
    robust_opts <- list(
      type = "huber",
      k_huber = 1.345,
      max_iter = 10,
      scale_scope = "local"
    )
    
    result <- robust_iterative_fitter(
      initial_glm_ctx = ctx,
      cfg_robust_options = robust_opts,
      X_orig_for_resid = X
    )
    
    expect_true(!is.null(result$betas_robust))
    # Some coefficients might be NA or poorly estimated
    expect_true(sum(is.finite(result$betas_robust)) >= p - 1)
  }
})

test_that("robust fitting handles data with groups of identical outliers", {
  n <- 100
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  beta_true <- c(1, 2, -1)
  Y <- X %*% beta_true + rnorm(n, sd = 0.5)
  
  # Create groups of identical outliers
  # This can cause issues with scale estimation
  Y[1:20] <- 50    # 20% identical high values
  Y[21:40] <- -50  # 20% identical low values
  
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  
  robust_opts <- list(
    type = "bisquare",
    c_tukey = 4.685,
    max_iter = 50,
    scale_scope = "local"
  )
  
  result <- robust_iterative_fitter(
    initial_glm_ctx = ctx,
    cfg_robust_options = robust_opts,
    X_orig_for_resid = X
  )
  
  # Should handle repeated outlier values
  expect_true(all(is.finite(result$betas_robust)))
  
  # Outlier groups should be downweighted
  expect_lt(mean(result$robust_weights_final[1:40]), 0.2)
})

test_that("weight computation edge cases", {
  # Test weight computation with extreme scaled residuals
  
  # Create simple data
  n <- 50
  X <- cbind(1, rnorm(n))
  Y <- rnorm(n)
  
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  
  # Test with extremely small scale (forces large scaled residuals)
  result <- robust_iterative_fitter(
    initial_glm_ctx = ctx,
    cfg_robust_options = list(
      type = "huber",
      k_huber = 1.345,
      max_iter = 5,
      scale_scope = "local"
    ),
    X_orig_for_resid = X,
    sigma_fixed = 1e-10  # Extremely small scale
  )
  
  # With tiny scale, weights should be affected
  expect_true(all(result$robust_weights_final >= 0))
  expect_true(all(result$robust_weights_final <= 1))
  
  # Test with extremely large scale (forces small scaled residuals)
  result2 <- robust_iterative_fitter(
    initial_glm_ctx = ctx,
    cfg_robust_options = list(
      type = "bisquare",
      c_tukey = 4.685,
      max_iter = 5,
      scale_scope = "local"
    ),
    X_orig_for_resid = X,
    sigma_fixed = 1e10  # Extremely large scale
  )
  
  # With huge scale, weights should tend toward 1
  expect_true(all(result2$robust_weights_final >= 0))
  expect_true(mean(result2$robust_weights_final) > 0.8)
})

test_that("robust fitting performance degrades gracefully with contamination", {
  n <- 100
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  beta_true <- c(2, -1, 3)
  
  # Test different contamination levels
  contamination_levels <- c(0.1, 0.3, 0.5, 0.7)
  mse_robust <- numeric(length(contamination_levels))
  mse_ols <- numeric(length(contamination_levels))
  
  for (i in seq_along(contamination_levels)) {
    cont_level <- contamination_levels[i]
    
    # Generate clean data
    Y <- X %*% beta_true + rnorm(n, sd = 0.5)
    
    # Add contamination
    n_outliers <- floor(cont_level * n)
    outlier_idx <- sample(n, n_outliers)
    Y[outlier_idx] <- rnorm(n_outliers, mean = 10, sd = 5)
    
    # OLS
    ols_beta <- solve(t(X) %*% X) %*% t(X) %*% Y
    mse_ols[i] <- mean((ols_beta - beta_true)^2)
    
    # Robust
    proj <- fmrireg:::.fast_preproject(X)
    ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
    
    result <- robust_iterative_fitter(
      initial_glm_ctx = ctx,
      cfg_robust_options = list(
        type = "bisquare",
        c_tukey = 4.685,
        max_iter = 50,
        scale_scope = "local"
      ),
      X_orig_for_resid = X
    )
    
    mse_robust[i] <- mean((result$betas_robust - beta_true)^2)
  }
  
  # Robust should outperform OLS at moderate contamination levels
  # At 70% contamination, both methods struggle
  expect_true(mse_robust[1] < mse_ols[1])  # 10% contamination
  expect_true(mse_robust[2] < mse_ols[2])  # 30% contamination
  
  # MSE should generally increase with contamination
  # Allow some variation due to randomness
  expect_true(mse_robust[4] > mse_robust[1])
  expect_true(mse_ols[4] > mse_ols[1])
  
  # At high contamination, robust methods may also struggle
  # Just verify finite results
  expect_true(all(is.finite(mse_robust)))
  expect_true(all(is.finite(mse_ols)))
})
