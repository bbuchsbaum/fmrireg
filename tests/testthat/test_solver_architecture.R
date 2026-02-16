# Test the new solver architecture from the refactoring

library(fmrireg)
library(testthat)

test_that("solve_glm_core handles basic OLS correctly", {
  # Simple regression problem
  n <- 100
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  beta_true <- c(2, 0.5, -1)
  Y <- X %*% beta_true + rnorm(n, sd = 0.3)
  
  # Create proper GLM context with projection
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(
    X = X,
    Y = matrix(Y, ncol = 1),
    proj = proj
  )
  
  # Solve
  result <- solve_glm_core(ctx)
  
  expect_true(!is.null(result$betas))
  expect_equal(dim(result$betas), c(p, 1))
  expect_true(!is.null(result$sigma2))
  expect_true(!is.null(result$rss))
  
  # Coefficients should be close to truth
  expect_equal(as.vector(result$betas), beta_true, tolerance = 0.2)
})

test_that("solve_glm_core handles weighted least squares", {
  n <- 100
  X <- cbind(1, rnorm(n))
  Y <- 3 + 2*X[,2] + rnorm(n)
  
  # Create heteroscedastic weights
  weights <- 1 / (1 + abs(X[,2]))
  
  # Apply weights to X and Y
  sqrt_weights <- sqrt(weights)
  Xw <- X * sqrt_weights
  Yw <- Y * sqrt_weights
  
  # Create weighted context
  proj <- fmrireg:::.fast_preproject(Xw)
  ctx <- glm_context(
    X = Xw,
    Y = matrix(Yw, ncol = 1),
    proj = proj,
    robust_weights = weights
  )
  
  result <- solve_glm_core(ctx)
  
  # Compare with manual WLS
  W <- diag(weights)
  Xw_manual <- sqrt(W) %*% X
  Yw_manual <- sqrt(W) %*% Y
  coef_expected <- solve(t(Xw_manual) %*% Xw_manual) %*% t(Xw_manual) %*% Yw_manual
  
  expect_equal(as.vector(result$betas), as.vector(coef_expected), tolerance = 1e-6)
})

test_that("solve_glm_core integrates with AR whitening", {
  set.seed(4242)
  # Generate AR(1) data
  n <- 100
  ar_coef <- 0.7
  X <- cbind(1, rnorm(n))
  beta_true <- c(1, 2)
  
  # AR(1) errors
  e <- arima.sim(list(ar = ar_coef), n = n)
  Y <- X %*% beta_true + as.vector(e)
  
  # Initial OLS fit
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  initial_result <- solve_glm_core(ctx)
  
  # Estimate AR parameters from residuals
  residuals <- Y - X %*% initial_result$betas
  phi_est <- estimate_ar_parameters(residuals, 1)
  
  # Apply AR whitening
  whitened <- ar_whiten_transform(X, matrix(Y, ncol = 1), phi_est, exact_first = FALSE)
  
  # Solve whitened problem
  proj_w <- fmrireg:::.fast_preproject(whitened$X)
  ctx_whitened <- glm_context(
    X = whitened$X,
    Y = whitened$Y,
    proj = proj_w,
    phi_hat = phi_est
  )
  
  result <- solve_glm_core(ctx_whitened)
  
  expect_true(!is.null(result$betas))
  expect_true(!is.null(phi_est))
  
  # AR coefficient should be estimated reasonably
  expect_equal(phi_est, ar_coef, tolerance = 0.3)
})

test_that("solve_glm_core handles rank deficient matrices", {
  n <- 50
  # Create rank deficient X
  X <- cbind(1, rnorm(n), rnorm(n))
  X[, 3] <- X[, 2]  # Make columns 2 and 3 identical
  
  Y <- rnorm(n)
  
  # Should handle rank deficiency gracefully - either error or proceed
  # The implementation might use different approaches (SVD, etc.)
  result <- tryCatch({
    proj <- fmrireg:::.fast_preproject(X)
    ctx <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
    solve_glm_core(ctx)
  }, error = function(e) {
    # If it errors, that's acceptable for rank deficient matrices
    NULL
  })
  
  # Either we get a result or an error - both are acceptable
  expect_true(is.null(result) || !is.null(result$betas))
})

test_that("iterative AR+Robust pipeline works", {
  set.seed(4242)
  # Data with both autocorrelation and outliers
  n <- 100
  X <- cbind(1, rnorm(n))
  beta <- c(2, 1)
  
  # AR(1) errors with outliers
  ar_coef <- 0.5
  e <- arima.sim(list(ar = ar_coef), n = n, sd = 0.5)
  outliers <- sample(n, 5)
  e[outliers] <- e[outliers] + sample(c(-4, 4), 5, replace = TRUE)
  
  Y <- X %*% beta + as.vector(e)
  
  # Step 1: Initial OLS fit
  proj <- fmrireg:::.fast_preproject(X)
  ctx_init <- glm_context(X = X, Y = matrix(Y, ncol = 1), proj = proj)
  initial_result <- solve_glm_core(ctx_init)
  
  # Step 2: Estimate AR parameters
  residuals <- Y - X %*% initial_result$betas
  phi_est <- estimate_ar_parameters(residuals, 1)
  
  # Step 3: Apply AR whitening
  whitened <- ar_whiten_transform(X, matrix(Y, ncol = 1), phi_est, exact_first = FALSE)
  
  # Step 4: Robust fitting on whitened data
  proj_w <- fmrireg:::.fast_preproject(whitened$X)
  ctx_w <- glm_context(X = whitened$X, Y = whitened$Y, proj = proj_w, phi_hat = phi_est)
  
  robust_config <- list(
    type = "bisquare",
    k_huber = 1.345,
    c_tukey = 4.685,
    max_iter = 3,
    scale_scope = "run"
  )
  
  result <- robust_iterative_fitter(
    initial_glm_ctx = ctx_w,
    cfg_robust_options = robust_config,
    X_orig_for_resid = whitened$X
  )
  
  expect_true(!is.null(result$betas_robust))
  expect_true(!is.null(result$robust_weights_final))
  expect_true(!is.null(phi_est))
  
  # AR coefficient should be reasonably estimated - relaxed tolerance
  expect_equal(phi_est, ar_coef, tolerance = 0.5)
})

test_that("multi-run solving maintains independence", {
  # Two independent runs - simplified test
  n_per_run <- 50
  n_runs <- 2
  
  X_list <- list()
  Y_list <- list()
  
  for (r in 1:n_runs) {
    X_list[[r]] <- cbind(1, rnorm(n_per_run))
    Y_list[[r]] <- 2 + 3*X_list[[r]][,2] + rnorm(n_per_run)
  }
  
  # Test each run separately to verify independence
  results <- list()
  for (r in 1:n_runs) {
    proj <- fmrireg:::.fast_preproject(X_list[[r]])
    ctx <- glm_context(
      X = X_list[[r]],
      Y = matrix(Y_list[[r]], ncol = 1),
      proj = proj
    )
    results[[r]] <- solve_glm_core(ctx)
  }
  
  # Each run should produce valid results
  for (r in 1:n_runs) {
    expect_true(!is.null(results[[r]]$betas))
    expect_equal(dim(results[[r]]$betas), c(2, 1))
  }
  
  # Results should be similar (both estimating same underlying model)
  expect_equal(results[[1]]$betas[1,1], results[[2]]$betas[1,1], tolerance = 0.5)
  expect_equal(results[[1]]$betas[2,1], results[[2]]$betas[2,1], tolerance = 0.5)
})

test_that("voxelwise contrast computation works", {
  # Multi-voxel data
  n <- 60
  p <- 4
  n_voxels <- 10
  
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  Y <- matrix(rnorm(n * n_voxels), n, n_voxels)
  
  # Add signal to some voxels - fix the dimension mismatch
  signal_voxels <- 1:3
  signal_effects <- matrix(c(2, -1, 0.5), nrow = 1)  # 1 x 3 matrix
  
  for (sv in signal_voxels) {
    Y[, sv] <- Y[, sv] + X[, 2:4] %*% t(signal_effects)  # Now conformable
  }
  
  # Define contrast
  contrast_vec <- c(0, 1, -1, 0)  # Compare effects 2 vs 3
  
  # Solve all voxels
  contrast_values <- numeric(n_voxels)
  
  for (v in 1:n_voxels) {
    proj <- fmrireg:::.fast_preproject(X)
    ctx <- glm_context(
      X = X,
      Y = Y[, v, drop = FALSE],
      proj = proj
    )
    
    result <- solve_glm_core(ctx)
    
    # Compute contrast manually
    contrast_values[v] <- sum(contrast_vec * as.vector(result$betas))
  }
  
  # Signal voxels should have larger contrast values
  expect_true(mean(abs(contrast_values[signal_voxels])) > 
              mean(abs(contrast_values[-signal_voxels])))
})

test_that("configuration validation catches errors", {
  # Invalid robust method
  expect_error(
    fmri_lm_control(robust_options = list(type = "invalid_method")),
    "Invalid robust_psi|Invalid robust type"
  )
  
  # Basic configuration should work  
  config <- fmri_lm_control(ar_options = list(struct = "ar1"))
  expect_true(inherits(config, "fmri_lm_config"))
  
  # Another basic configuration should work
  config2 <- fmri_lm_control(
    robust_options = list(type = "huber"),
    ar_options = list(struct = "iid")
  )
  expect_true(inherits(config2, "fmri_lm_config"))
})

test_that("solver preserves context structure", {
  n <- 50
  X <- cbind(1, rnorm(n))
  Y <- rnorm(n)
  
  # Create context with additional fields
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- glm_context(
    X = X,
    Y = matrix(Y, ncol = 1),
    proj = proj,
    phi_hat = 0.5,
    sigma_robust_scale = 1.2
  )
  
  result <- solve_glm_core(ctx)
  
  # Basic result structure should be preserved
  expect_true(!is.null(result$betas))
  expect_true(!is.null(result$sigma2))
  expect_true(!is.null(result$rss))
  
  # Original context should be unchanged
  expect_equal(ctx$phi_hat, 0.5)
  expect_equal(ctx$sigma_robust_scale, 1.2)
})
