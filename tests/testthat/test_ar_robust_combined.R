# Test combined AR + Robust regression

test_that("AR + Robust combination works in runwise", {
  skip_if_not_installed("neuroim2")
  
  # Generate data with AR errors and outliers
  n_time <- 100
  n_vox <- 5
  n_runs <- 2
  true_ar <- 0.6
  
  # Create design
  onsets <- c(10, 30, 50, 70, 90)
  blockvar <- rep(1:n_runs, each = n_time/n_runs)
  
  # Generate AR(1) errors with outliers
  Y <- matrix(0, n_time, n_vox)
  for (run in 1:n_runs) {
    run_idx <- which(blockvar == run)
    for (v in 1:n_vox) {
      # AR(1) noise
      noise <- numeric(length(run_idx))
      noise[1] <- rnorm(1)
      for (t in 2:length(run_idx)) {
        noise[t] <- true_ar * noise[t-1] + rnorm(1)
      }
      
      # Add outliers
      if (v <= 2) {  # Add outliers to first 2 voxels
        outlier_times <- sample(length(run_idx), 3)
        noise[outlier_times] <- noise[outlier_times] + rnorm(3, mean = 5, sd = 1)
      }
      
      Y[run_idx, v] <- noise
    }
  }
  
  # Add signal
  ev_df <- data.frame(
    onset = rep(onsets, n_runs),
    run = rep(1:n_runs, each = length(onsets))
  )
  
  dset <- matrix_dataset(Y, TR = 1, run_length = rep(n_time/n_runs, n_runs), event_table = ev_df)
  
  # Fit with AR + Robust  
  fit <- fmri_lm(
    onset ~ hrf(onset, hrf_spmg1()),
    block = ~ run,
    dataset = dset,
    strategy = "runwise",
    robust_options = list(type = "bisquare", tuning = 4.685),
    ar_options = list(struct = "ar1", iter_gls = 1),
    use_fast_path = TRUE
  )
  
  expect_s3_class(fit, "fmri_lm")
  expect_true(!is.null(fit$result$betas))
  
  # Check that both AR and robust were applied
  cfg <- attr(fit, "config")
  expect_equal(cfg$ar$struct, "ar1")
  expect_equal(cfg$robust$type, "bisquare")
})

test_that("AR + Robust with re-estimation works", {
  skip_if_not_installed("neuroim2")
  
  # Simple single-run test
  n_time <- 80
  n_vox <- 3
  
  # Generate strong AR(1) signal
  true_ar <- 0.8
  Y <- matrix(0, n_time, n_vox)
  
  for (v in 1:n_vox) {
    noise <- numeric(n_time)
    noise[1] <- rnorm(1)
    for (t in 2:n_time) {
      noise[t] <- true_ar * noise[t-1] + rnorm(1, sd = 0.5)
    }
    Y[, v] <- noise
  }
  
  # Add outliers that might affect AR estimation
  Y[c(10, 20, 30), 1] <- Y[c(10, 20, 30), 1] + 10
  
  dset <- matrix_dataset(Y, TR = 1, run_length = n_time, event_table = data.frame(onset = c(15, 35, 55), run = 1))
  
  # Fit with re-estimation
  fit_reest <- fmri_lm(
    onset ~ hrf(onset, hrf_spmg1()),
    block = ~ run,
    dataset = dset,
    robust_options = list(
      type = "huber",
      k_huber = 1.345,
      reestimate_phi = TRUE
    ),
    ar_options = list(struct = "ar1"),
    use_fast_path = TRUE
  )
  
  # Fit without re-estimation
  fit_no_reest <- fmri_lm(
    onset ~ hrf(onset, hrf_spmg1()),
    block = ~ run,
    dataset = dset,
    robust_options = list(
      type = "huber",
      k_huber = 1.345,
      reestimate_phi = FALSE
    ),
    ar_options = list(struct = "ar1"),
    use_fast_path = TRUE
  )
  
  # Both should work
  expect_s3_class(fit_reest, "fmri_lm")
  expect_s3_class(fit_no_reest, "fmri_lm")
})

test_that("process_run_ar_robust handles edge cases", {
  # Small run with few time points
  n_time <- 20
  n_vox <- 2
  
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  dset <- matrix_dataset(Y, TR = 2, run_length = n_time, event_table = data.frame(onset = c(5, 15), run = 1))
  
  # Create minimal model
  sframe <- sampling_frame(n_time, TR = 2)
  ev <- event_model(onset ~ hrf(onset, hrf_spmg1()),
                    data = data.frame(onset = c(5, 15), run = 1),
                    block = ~ run)
  bmodel <- baseline_model(sframe, degree = 1)
  fmodel <- fmri_model(ev, bmodel)
  
  # Create config
  cfg <- fmri_lm_control(
    robust_options = list(type = "bisquare"),
    ar_options = list(struct = "ar1")
  )
  
  # Get run chunk
  chunks <- exec_strategy("runwise")(dset)
  
  # Process with AR + Robust
  result <- process_run_ar_robust(
    run_chunk = chunks[[1]],
    model = fmodel,
    cfg = cfg
  )
  
  # Check output structure
  expect_true(all(c("betas", "sigma2", "XtXinv", "robust_weights", 
                    "phi_hat", "ar_order") %in% names(result)))
  expect_equal(result$ar_order, 1)
  expect_length(result$phi_hat, 1)
  expect_length(result$robust_weights, n_time)
})

test_that("Chunkwise AR + Robust works", {
  skip_if_not_installed("neuroim2")
  
  # Multi-run data for chunkwise
  n_time <- 60
  n_vox <- 20
  n_runs <- 2
  
  # Generate AR data
  Y <- matrix(0, n_time, n_vox)
  blockvar <- rep(1:n_runs, each = n_time/n_runs)
  
  for (run in 1:n_runs) {
    run_idx <- which(blockvar == run)
    ar_coef <- 0.5 + 0.2 * (run - 1)  # Different AR per run
    
    for (v in 1:n_vox) {
      noise <- numeric(length(run_idx))
      noise[1] <- rnorm(1)
      for (t in 2:length(run_idx)) {
        noise[t] <- ar_coef * noise[t-1] + rnorm(1)
      }
      
      # Add outliers to some voxels
      if (v %% 5 == 0) {
        outliers <- sample(length(run_idx), 2)
        noise[outliers] <- noise[outliers] + rnorm(2, mean = 4)
      }
      
      Y[run_idx, v] <- noise
    }
  }
  
  dset <- matrix_dataset(Y, TR = 1, run_length = rep(n_time/n_runs, n_runs), 
                         event_table = data.frame(
                           onset = rep(c(10, 20), n_runs),
                           run = rep(1:n_runs, each = 2)
                         ))
  
  # Test chunkwise with AR + Robust
  fit_chunk <- fmri_lm(
    onset ~ hrf(onset, hrf_spmg1()),
    block = ~ run,
    dataset = dset,
    strategy = "chunkwise",
    nchunks = 2,
    robust_options = list(type = "huber"),
    ar_options = list(struct = "ar1"),
    use_fast_path = TRUE
  )
  
  expect_s3_class(fit_chunk, "fmri_lm")
  expect_equal(attr(fit_chunk, "strategy"), "chunkwise")
  
  # Check results structure
  expect_true(!is.null(fit_chunk$result$betas))
  expect_true(!is.null(fit_chunk$result$contrasts))
})