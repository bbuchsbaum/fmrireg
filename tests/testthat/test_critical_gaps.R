# Test critical functionality gaps

library(fmrireg)
library(testthat)

test_that("GLM solver handles rank deficiency", {
  # Create rank-deficient matrix
  n <- 50
  X <- cbind(1, rnorm(n), rnorm(n))
  X <- cbind(X, X[,2] + X[,3])  # Linear combination
  Y <- matrix(rnorm(n * 2), n, 2)
  
  # Check that matrix is rank deficient
  expect_lt(qr(X)$rank, ncol(X))
  
  # Should handle gracefully
  expect_silent({
    proj <- fmrireg:::.fast_preproject(X)
    ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
    result <- fmrireg:::solve_glm_core(ctx)
  })
  
  expect_equal(dim(result$betas), c(ncol(X), ncol(Y)))
})

test_that("Robust fitting downweights outliers", {
  n <- 100
  X <- cbind(1, rnorm(n))
  y <- 2 + 3 * X[,2] + rnorm(n, sd = 0.5)
  
  # Add outliers
  outlier_idx <- 1:10
  y[outlier_idx] <- y[outlier_idx] + 10
  
  # Create proper context with projection
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- fmrireg:::glm_context(X = X, Y = matrix(y, ncol = 1), proj = proj)
  cfg <- list(
    type = "huber",
    k_huber = 1.345,
    max_iter = 5,
    scale_scope = "global"
  )
  
  result <- fmrireg:::robust_iterative_fitter(ctx, cfg, X_orig_for_resid = X)
  
  # Outliers should have lower weights
  expect_true(mean(result$robust_weights_final[outlier_idx]) < mean(result$robust_weights_final[-outlier_idx]))
  expect_true(all(result$robust_weights_final >= 0))
  expect_true(all(result$robust_weights_final <= 1))
})

test_that("AR whitening transforms data correctly", {
  n <- 100
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  Y <- matrix(rnorm(n * 2), n, 2)
  
  # Create AR(1) structure
  phi <- 0.6
  
  # Apply whitening - note function signature is (Y, X, phi) and returns list
  result <- fmrireg:::ar_whiten_inplace(Y, X, phi)
  X_white <- result$X
  Y_white <- result$Y
  
  # Check dimensions preserved
  expect_equal(dim(X_white), dim(X))
  expect_equal(dim(Y_white), dim(Y))
  
  # Check that the result has the expected structure
  expect_true(is.matrix(X_white))
  expect_true(is.matrix(Y_white))
  
  # Basic sanity check - first element should be unchanged (no previous value)
  expect_equal(X_white[1,], X[1,])
  expect_equal(Y_white[1,], Y[1,])
})

test_that("Mixed solve handles weighted regression", {
  n <- 50
  X <- cbind(1, rnorm(n))
  y <- X %*% c(2, 3) + rnorm(n)
  
  # For mixed model: need Z (random effects design), K (kinship), and proper response vector
  Z <- diag(n)  # Identity matrix for random effects
  K <- diag(n)  # Identity kinship matrix
  
  result <- mixed_solve_cpp(y = y, Z = Z, K = K, X = X)
  
  expect_equal(length(result$beta), ncol(X))
  expect_true(!is.null(result$Vu))  # Variance component for random effects
  expect_true(!is.null(result$Ve))  # Variance component for residuals
  
  # Coefficients should be reasonable
  expect_true(all(abs(result$beta - c(2, 3)) < 1))
})

test_that("Data chunking works correctly", {
  # Create test dataset
  n_time <- 200
  n_voxels <- 10
  data_mat <- matrix(rnorm(n_time * n_voxels), n_time, n_voxels)
  
  dset <- matrix_dataset(
    data_mat,
    TR = 2,
    run_length = c(100, 100)
  )
  
  # Test chunking - may return fewer chunks than requested
  chunks <- data_chunks(dset, nchunks = 4, runwise = FALSE)
  
  expect_true(length(chunks) > 0)
  expect_true(is.list(chunks))
  
  # Check that chunks are data_chunk objects
  expect_true(all(sapply(chunks, function(ch) inherits(ch, "data_chunk"))))
})

test_that("Event model integration works", {
  # Create simple event data
  event_data <- data.frame(
    onsets = c(5, 15, 25, 35),
    condition = factor(c("A", "B", "A", "B")),
    run = rep(1, 2)
  )
  
  sframe <- fmrihrf::sampling_frame(blocklens = 50, TR = 2)
  
  # Create event model
  em <- event_model(
    onsets ~ hrf(condition),
    data = event_data,
    block = ~ run,
    sampling_frame = sframe
  )
  
  expect_s3_class(em, "event_model")
  expect_true(!is.null(em$terms))
  
  # Get design matrix
  dm <- design_matrix(em)
  expect_equal(nrow(dm), 50)  # Number of time points
  expect_true(ncol(dm) >= 2)  # At least 2 conditions
})

test_that("Baseline model creates proper drift terms", {
  sframe <- fmrihrf::sampling_frame(blocklens = c(100, 100), TR = 2)
  
  bmodel <- baseline_model(
    basis = "poly",
    degree = 2,
    sframe = sframe
  )
  
  expect_s3_class(bmodel, "baseline_model")
  
  dm <- design_matrix(bmodel)
  expect_equal(nrow(dm), 200)  # Total time points
  expect_true(ncol(dm) >= 4)   # At least intercept + 2 poly terms per run
})

test_that("fmri_lm basic functionality works", {
  # Create minimal test case
  n <- 100
  event_data <- data.frame(
    onsets = seq(10, 90, by = 20),
    condition = factor(rep(c("A", "B"), length.out = 5)),
    run = rep(1, 5)
  )
  
  # Simple data
  Y <- matrix(rnorm(n * 5), n, 5)
  
  dset <- matrix_dataset(
    Y,
    TR = 1,
    run_length = n,
    event_table = event_data
  )
  
  # Fit model
  result <- fmri_lm(
    onsets ~ hrf(condition),
    block = ~ run,
    dataset = dset
  )
  
  expect_s3_class(result, "fmri_lm")
  expect_true(!is.null(result$result))
  
  # Extract coefficients
  betas <- coef(result, type = "betas")
  expect_true(!is.null(betas))
})
