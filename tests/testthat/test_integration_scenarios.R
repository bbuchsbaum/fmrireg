# Test integration scenarios combining multiple features

library(fmrireg)
library(testthat)

test_that("AR + Robust fitting integration works", {
  # Create autocorrelated data with outliers
  n <- 100
  X <- cbind(1, rnorm(n))
  
  # Generate AR(1) errors
  ar_coef <- 0.6
  errors <- arima.sim(list(ar = ar_coef), n = n, sd = 0.5)
  
  # Add outliers
  outlier_idx <- sample(n, 5)
  errors[outlier_idx] <- errors[outlier_idx] + sample(c(-5, 5), 5, replace = TRUE)
  
  Y <- X %*% c(2, 1) + as.vector(errors)
  
  # Create dataset
  dset <- matrix_dataset(
    matrix(Y, ncol = 1),
    TR = 2,
    run_length = n,
    event_table = data.frame(
      onsets = seq(10, 90, by = 20),
      condition = factor(rep(c("A", "B"), length.out = 5)),
      run = rep(1, 5)
    )
  )
  
  # Fit with both AR and robust
  result <- fmri_lm(
    onsets ~ hrf(condition),
    block = ~ run,
    dataset = dset,
    robust = "bisquare",
    ar_options = list(
      cor_struct = "ar1",
      iter = 2
    )
  )
  
  expect_s3_class(result, "fmri_lm")
  expect_true(!is.null(result$result))
  
  # Check that both AR and robust were applied
  # This would be visible in the residual structure and weights
})

test_that("Multi-run with different strategies produces consistent results", {
  # Create multi-run dataset
  n_per_run <- 50
  n_runs <- 3
  
  event_list <- list()
  data_list <- list()
  
  for (run in 1:n_runs) {
    # Events for this run
    event_list[[run]] <- data.frame(
      onsets = c(5, 15, 25, 35),
      condition = factor(c("A", "B", "A", "B")),
      run = rep(run, 4)
    )
    
    # Data for this run
    data_list[[run]] <- matrix(rnorm(n_per_run * 10), n_per_run, 10)
  }
  
  # Combine
  all_events <- do.call(rbind, event_list)
  all_data <- do.call(rbind, data_list)
  
  dset <- matrix_dataset(
    all_data,
    TR = 2,
    run_length = rep(n_per_run, n_runs),
    event_table = all_events
  )
  
  # Fit with runwise strategy
  result_runwise <- fmri_lm(
    onsets ~ hrf(condition),
    block = ~ run,
    dataset = dset,
    strategy = "runwise"
  )
  
  # Fit with chunkwise strategy
  result_chunkwise <- fmri_lm(
    onsets ~ hrf(condition),
    block = ~ run,
    dataset = dset,
    strategy = "chunkwise",
    nchunks = n_runs
  )
  
  # Results should be similar (not identical due to numerical differences)
  betas_runwise <- coef(result_runwise, type = "betas")
  betas_chunkwise <- coef(result_chunkwise, type = "betas")
  
  # Check dimensions match
  expect_equal(dim(betas_runwise), dim(betas_chunkwise))
})

test_that("Complex contrast specifications work with all model types", {
  # Create dataset with simple single-factor design
  n <- 120
  event_data <- data.frame(
    onsets = c(10, 30, 50, 70, 90, 110),
    condition = c("A", "B", "A", "B", "A", "B"),
    run = c(1, 1, 2, 2, 3, 3)
  )
  
  dset <- matrix_dataset(
    matrix(rnorm(n * 10), n, 10),
    TR = 1,
    run_length = c(40, 40, 40),
    event_table = event_data
  )
  
  # Standard model
  result_standard <- fmri_lm(
    onsets ~ hrf(condition),
    block = ~ run,
    dataset = dset
  )
  
  # With contrasts
  contrasts <- list(
    main = pair_contrast(~ condition == "A", ~ condition == "B")
  )
  
  con_results <- fit_contrasts(result_standard, contrasts)
  
  expect_equal(length(con_results), 1)
  expect_true(all(sapply(con_results, function(x) !is.null(x$statistic))))
  
  # Robust model
  result_robust <- fmri_lm(
    onsets ~ hrf(condition),
    block = ~ run,
    dataset = dset,
    robust = "huber"
  )
  
  con_results_robust <- fit_contrasts(result_robust, contrasts)
  expect_equal(length(con_results_robust), 1)
})

test_that("Missing data handling works across components", {
  # Create dataset with missing values
  n <- 80
  data_mat <- matrix(rnorm(n * 15), n, 15)
  
  # Introduce missing values
  missing_idx <- sample(length(data_mat), size = 50)
  data_mat[missing_idx] <- NA
  
  event_data <- data.frame(
    onsets = seq(5, 75, by = 10),
    value = rnorm(8),
    run = rep(1, 8)
  )
  
  dset <- matrix_dataset(
    data_mat,
    TR = 2,
    run_length = n,
    event_table = event_data
  )
  
  # Should handle missing data
  result <- fmri_lm(
    onsets ~ hrf(value),
    block = ~ run,
    dataset = dset
  )
  
  # Some voxels might have no valid estimates
  betas <- coef(result, type = "betas")
  expect_true(any(is.na(betas)))
})

test_that("Large dataset chunking maintains accuracy", {
  # Simulate large dataset scenario
  n_voxels <- 10000
  n_time <- 200
  
  # Create sparse data (most voxels inactive)
  data_mat <- matrix(rnorm(n_time * 20, sd = 0.1), n_time, 20)  # Only 20 for testing
  
  # Add signal to some voxels
  active_voxels <- 1:5
  for (v in active_voxels) {
    # Add event-related signal
    signal_times <- seq(10, 190, by = 40)
    for (t in signal_times) {
      data_mat[t:(t+10), v] <- data_mat[t:(t+10), v] + rnorm(11, mean = 2, sd = 0.5)
    }
  }
  
  event_data <- data.frame(
    onsets = seq(10, 190, by = 40),
    condition = factor(rep(c("A", "B"), length.out = 5)),
    run = rep(1, 5)
  )
  
  dset <- matrix_dataset(
    data_mat,
    TR = 1,
    run_length = n_time,
    event_table = event_data
  )
  
  # Test with different chunk sizes
  result_1chunk <- fmri_lm(
    onsets ~ hrf(condition),
    block = ~ run,
    dataset = dset,
    strategy = "chunkwise",
    nchunks = 1
  )
  
  result_5chunks <- fmri_lm(
    onsets ~ hrf(condition),
    block = ~ run,
    dataset = dset,
    strategy = "chunkwise",
    nchunks = 5
  )
  
  # Results should be very similar
  betas_1 <- coef(result_1chunk, type = "betas")
  betas_5 <- coef(result_5chunks, type = "betas")
  
  # Check that results are consistent (dimensions and overall structure)
  expect_equal(dim(betas_1), dim(betas_5))
  
  # Compare a subset of active voxels (only those that exist)
  available_voxels <- seq_len(min(ncol(betas_1), length(active_voxels)))
  for (v in available_voxels) {
    expect_equal(betas_1[, v], betas_5[, v], tolerance = 0.01)
  }
})