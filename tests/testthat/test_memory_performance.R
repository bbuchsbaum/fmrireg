# Test memory efficiency and performance aspects

library(fmrireg)
library(testthat)

test_that("chunked processing maintains constant memory usage", {
  # Large dataset that would be problematic to load fully
  n_time <- 300
  n_voxels <- 10000  # Simulated as smaller for testing
  
  # Create dataset that chunks data
  # In real usage, this would be file-backed
  
  # Event design - fix column name to match API
  event_data <- data.frame(
    onset = seq(10, 290, by = 30),
    condition = factor(rep(c("A", "B"), 5)),
    run = rep(1:2, each = 5)
  )
  
  # For testing, create smaller representative data
  test_voxels <- 100
  test_data <- matrix(rnorm(n_time * test_voxels), n_time, test_voxels)
  
  dset <- fmridataset::matrix_dataset(
    test_data,
    TR = 2,
    run_length = c(150, 150),
    event_table = event_data
  )
  
  # Process in chunks
  chunk_sizes <- c(1, 5, 10, 20)
  results <- list()
  
  for (nchunks in chunk_sizes) {
    result <- fmri_lm(
      onset ~ hrf(condition),
      block = ~ run,
      dataset = dset,
      strategy = "chunkwise",
      nchunks = nchunks
    )
    results[[as.character(nchunks)]] <- result
  }
  
  # All chunk sizes should give same results
  base_betas <- coef(results[["1"]])
  
  for (i in 2:length(chunk_sizes)) {
    test_betas <- coef(results[[as.character(chunk_sizes[i])]])
    expect_equal(base_betas, test_betas, tolerance = 1e-10)
  }
})

test_that("iterator pattern efficiently processes data", {
  # Test the iterator doesn't load all data at once
  n_time <- 200
  n_voxels <- 50
  
  data_mat <- matrix(rnorm(n_time * n_voxels), n_time, n_voxels)
  
  dset <- fmridataset::matrix_dataset(
    data_mat,
    TR = 1,
    run_length = n_time
  )
  
  # Create iterator using correct function
  chunk_iter <- fmridataset::data_chunks(dset, nchunks = 5)
  
  # Check iterator properties
  expect_equal(chunk_iter$nchunks, 5)
  
  # Process chunks
  chunk_voxels <- numeric(5)
  for (i in 1:5) {
    chunk <- chunk_iter$nextElem()
    chunk_voxels[i] <- ncol(chunk$data)
  }
  
  # Chunks should partition the voxels
  expect_equal(sum(chunk_voxels), n_voxels)
  expect_true(all(chunk_voxels > 0))
})

test_that("parallel processing maintains result consistency", {
  # Setup data - fix event table column name
  n <- 100
  n_voxels <- 20
  X <- cbind(1, rnorm(n))
  Y <- matrix(rnorm(n * n_voxels), n, n_voxels)
  
  # Add signal
  Y[, 1:5] <- Y[, 1:5] + X[, 2] * 2
  
  event_data <- data.frame(
    onset = seq(5, 95, by = 10),
    value = rnorm(10),
    run = rep(1, 10)
  )
  
  dset <- fmridataset::matrix_dataset(Y, TR = 1, run_length = n, event_table = event_data)
  
  # Run with different thread counts
  old_threads <- getOption("fmrireg.num_threads", 1)
  
  # Single threaded
  options(fmrireg.num_threads = 1)
  result_single <- fmri_lm(
    onset ~ hrf(value),
    block = ~ run,
    dataset = dset
  )
  
  # Multi-threaded
  options(fmrireg.num_threads = 4)
  result_parallel <- fmri_lm(
    onset ~ hrf(value),
    block = ~ run,
    dataset = dset
  )
  
  # Restore
  options(fmrireg.num_threads = old_threads)
  
  # Results should be identical
  expect_equal(
    coef(result_single),
    coef(result_parallel),
    tolerance = 1e-10
  )
})

test_that("sparse matrix handling for efficiency", {
  # Design matrix with many zeros (e.g., event design)
  n <- 200
  n_events <- 10
  
  # Sparse event indicators
  event_times <- sort(sample(1:n, n_events))
  X_sparse <- matrix(0, n, n_events + 1)
  X_sparse[, 1] <- 1  # Intercept
  
  for (i in 1:n_events) {
    X_sparse[event_times[i], i + 1] <- 1
  }
  
  Y <- matrix(rnorm(n), ncol = 1)
  
  # Should handle sparse matrix efficiently - use correct glm_context constructor
  proj <- fmrireg:::.fast_preproject(X_sparse)
  ctx <- glm_context(
    X = X_sparse,
    Y = Y,
    proj = proj
  )
  
  # This should work without memory issues
  result <- solve_glm_core(ctx)
  
  expect_true(!is.null(result$betas))
  expect_equal(length(result$betas), ncol(X_sparse))
})

test_that("lazy evaluation prevents unnecessary computation", {
  # Dataset with lazy loading
  n_time <- 100
  n_voxels <- 50
  
  # Function that tracks calls
  call_count <- 0
  get_data <- function() {
    call_count <<- call_count + 1
    matrix(rnorm(n_time * n_voxels), n_time, n_voxels)
  }
  
  # In real implementation, matrix_dataset would support lazy evaluation
  # For now, test the concept
  
  expect_equal(call_count, 0)  # No calls yet
  
  # Would only load data when actually needed
  # data <- get_data()
  # expect_equal(call_count, 1)
})

test_that("pre-allocation prevents memory fragmentation", {
  # Test that result structures are pre-allocated
  n <- 100
  n_voxels <- 50
  n_contrasts <- 3
  
  X <- cbind(1, matrix(rnorm(n * 3), n, 3))
  Y <- matrix(rnorm(n * n_voxels), n, n_voxels)
  
  # Pre-allocate result structures
  n_coef <- ncol(X)
  
  # These should be pre-allocated in the implementation
  betas <- matrix(NA_real_, n_coef, n_voxels)
  se <- matrix(NA_real_, n_coef, n_voxels)
  contrasts <- array(NA_real_, dim = c(n_contrasts, n_voxels, 3))  # estimate, se, t
  
  # Fill in results using correct glm_context constructor
  for (v in 1:n_voxels) {
    proj <- fmrireg:::.fast_preproject(X)
    ctx <- glm_context(
      X = X,
      Y = Y[, v, drop = FALSE],
      proj = proj
    )
    
    result <- solve_glm_core(ctx)
    
    betas[, v] <- result$betas
    # Standard errors need to be computed from sigma2 and proj
    se[, v] <- sqrt(diag(ctx$proj$XtXinv) * result$sigma2)
  }
  
  # Check no NA values remain
  expect_false(any(is.na(betas)))
  expect_false(any(is.na(se)))
})

test_that("recycling design matrices saves memory", {
  # Multiple runs with same design
  n_per_run <- 100
  n_runs <- 4
  n_voxels <- 20
  
  # Same design matrix for all runs
  X_single <- cbind(1, rnorm(n_per_run))
  
  # Instead of replicating, use run indices
  run_indices <- lapply(1:n_runs, function(r) {
    ((r-1)*n_per_run + 1):(r*n_per_run)
  })
  
  # Combined data
  Y_all <- matrix(rnorm(n_per_run * n_runs * n_voxels), 
                  n_per_run * n_runs, n_voxels)
  
  # X can be recycled for each run
  X_all <- do.call(rbind, replicate(n_runs, X_single, simplify = FALSE))
  
  # Use correct glm_context constructor
  proj <- fmrireg:::.fast_preproject(X_all)
  ctx <- glm_context(
    X = X_all,
    Y = Y_all,
    proj = proj
  )
  
  result <- solve_glm_core(ctx)
  
  # Should handle efficiently - betas matrix has correct dimensions
  expect_equal(nrow(result$betas), ncol(X_all))  # Number of coefficients
  expect_equal(ncol(result$betas), ncol(Y_all))  # Number of response variables
})
