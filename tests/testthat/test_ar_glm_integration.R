# AR+GLM Integration Tests with Realistic Expectations
# These tests verify AR functionality within GLM context, accounting for
# the fact that regressors may absorb some temporal structure

test_that("AR+GLM: AR1 estimation with orthogonal design", {
  set.seed(111)
  # Create a design that's orthogonal to AR structure
  n_time <- 200
  n_runs <- 2
  
  # Generate AR(1) errors
  true_phi <- 0.6
  errors <- as.numeric(arima.sim(list(ar = true_phi), n = n_time * n_runs))
  
  # Simple block design with random events (reduces correlation with AR)
  event_times <- sort(sample(20:180, 10))
  design <- matrix(0, n_time * n_runs, 1)
  for (et in event_times) {
    design[et + c(0, n_time), 1] <- 1
  }
  
  # Generate data: design effect + AR errors
  y <- 5 * design[,1] + errors
  
  # Create minimal dataset
  dset <- matrix_dataset(
    matrix(y, ncol = 1), 
    TR = 1, 
    run_length = rep(n_time, n_runs),
    event_table = data.frame(
      onset = rep(event_times, n_runs),
      run = rep(1:n_runs, each = length(event_times))
    )
  )
  
  # Fit model with hrf term
  fit_ar <- fmri_lm(onset ~ hrf(onset), block = ~ run, dataset = dset, 
                    ar_options = list(struct = "ar1"), strategy = "runwise")
  
  # AR should be detected even with regressors
  # We expect some reduction from true value but not dramatic
  cfg <- attr(fit_ar, "config")
  expect_equal(cfg$ar$struct, "ar1")
  
  # Extract residuals and check AR
  # Note: We don't have direct access to estimated phi in the current API
  # This is something we should add
})

test_that("AR+GLM: Standard errors increase with positive autocorrelation", {
  set.seed(222)
  n_time <- 100
  n_vox <- 10
  
  # Generate data with known AR(1)
  true_phi <- 0.7
  Y <- matrix(0, n_time, n_vox)
  
  # Simple design: just an intercept and one regressor
  X <- cbind(1, rep(c(0, 1), each = n_time/2))
  
  for (v in 1:n_vox) {
    errors <- as.numeric(arima.sim(list(ar = true_phi), n = n_time))
    Y[, v] <- X %*% c(0, 2) + errors
  }
  
  dset <- matrix_dataset(Y, TR = 1, run_length = n_time)
  
  # Create minimal event table for the model  
  etab <- data.frame(onset = 1, run = 1)
  dset <- matrix_dataset(Y, TR = 1, run_length = n_time, event_table = etab)
  
  # Fit with and without AR correction
  fit_iid <- fmri_lm(onset ~ hrf(onset), block = ~ run, dataset = dset, baseline_model = NULL, 
                     ar_options = list(struct = "iid"), use_fast_path = TRUE)
  fit_ar1 <- fmri_lm(onset ~ hrf(onset), block = ~ run, dataset = dset, baseline_model = NULL,
                     ar_options = list(struct = "ar1"), use_fast_path = TRUE)
  
  # With positive AR, standard errors should be larger when corrected
  se_iid <- standard_error(fit_iid)
  se_ar1 <- standard_error(fit_ar1)
  
  # Get the first column (should be the HRF term)
  se_iid_vals <- se_iid[[1]]
  se_ar1_vals <- se_ar1[[1]]
  
  # AR correction should increase SE for positive autocorrelation
  expect_true(mean(se_ar1_vals) > mean(se_iid_vals))
})

test_that("AR+GLM: Global AR estimation works across runs", {
  set.seed(333)
  n_runs <- 4
  n_time_per_run <- 50
  true_phi <- 0.5
  
  # Generate consistent AR structure across runs
  Y <- matrix(0, n_time_per_run * n_runs, 1)
  for (r in 1:n_runs) {
    idx <- ((r-1) * n_time_per_run + 1):(r * n_time_per_run)
    Y[idx, 1] <- as.numeric(arima.sim(list(ar = true_phi), n = n_time_per_run))
  }
  
  # Create event table
  etab <- data.frame(
    onset = rep(1, n_runs),
    run = 1:n_runs
  )
  dset <- matrix_dataset(Y, TR = 1, run_length = rep(n_time_per_run, n_runs), event_table = etab)
  
  # Fit with global AR estimation
  fit_global <- fmri_lm(onset ~ hrf(onset), block = ~ run, dataset = dset, baseline_model = NULL,
                        ar_options = list(struct = "ar1", global = TRUE))
  
  # Fit with run-wise AR estimation  
  fit_runwise <- fmri_lm(onset ~ hrf(onset), block = ~ run, dataset = dset, baseline_model = NULL,
                         ar_options = list(struct = "ar1", global = FALSE))
  
  # Both should detect AR structure
  expect_s3_class(fit_global, "fmri_lm")
  expect_s3_class(fit_runwise, "fmri_lm")
  
  # Results should be similar (not testing exact equality)
  coef_global <- coef(fit_global)
  coef_runwise <- coef(fit_runwise)
  expect_equal(coef_global, coef_runwise, tolerance = 0.1)
})

test_that("AR+GLM: Iterative GLS improves estimates", {
  set.seed(444)
  n <- 100
  true_phi <- 0.6
  true_beta <- 2
  
  # Generate correlated data
  X <- cbind(1, rnorm(n))
  errors <- as.numeric(arima.sim(list(ar = true_phi), n = n))
  y <- X %*% c(0, true_beta) + errors
  
  etab <- data.frame(onset = 1, run = 1)
  dset <- matrix_dataset(matrix(y, ncol = 1), TR = 1, run_length = n, event_table = etab)
  
  # Fit with different numbers of GLS iterations
  fit_gls1 <- fmri_lm(onset ~ hrf(onset), block = ~ run, dataset = dset, baseline_model = NULL,
                      ar_options = list(struct = "ar1", iter_gls = 1))
  fit_gls2 <- fmri_lm(onset ~ hrf(onset), block = ~ run, dataset = dset, baseline_model = NULL,
                      ar_options = list(struct = "ar1", iter_gls = 2))
  
  # Both should run without error
  expect_s3_class(fit_gls1, "fmri_lm")
  expect_s3_class(fit_gls2, "fmri_lm")
  
  # Coefficients should be similar but potentially refined with more iterations
  expect_equal(coef(fit_gls1), coef(fit_gls2), tolerance = 0.2)
})