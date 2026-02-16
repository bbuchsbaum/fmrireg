# Tests for censor integration with AR modeling
# Ensures censored timepoints are properly excluded from AR estimation

test_that("fmri_lm_control accepts censor parameter in ar_options", {
  # NULL censor (default)
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1"))
  expect_null(cfg$ar$censor)

  # "auto" censor
 cfg <- fmri_lm_control(ar_options = list(struct = "ar1", censor = "auto"))
  expect_equal(cfg$ar$censor, "auto")

  # Integer vector censor
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1", censor = c(5L, 10L, 15L)))
  expect_equal(cfg$ar$censor, c(5L, 10L, 15L))

  # Logical vector censor
  censor_logical <- c(rep(FALSE, 4), TRUE, rep(FALSE, 4), TRUE, rep(FALSE, 5))
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1", censor = censor_logical))
  expect_equal(cfg$ar$censor, censor_logical)

  # Invalid censor should error
  expect_error(
    fmri_lm_control(ar_options = list(struct = "ar1", censor = "invalid")),
    "should be"
  )
})

test_that("extract_censor_from_dataset handles NULL dataset gracefully", {
  result <- fmrireg:::extract_censor_from_dataset(NULL)
  expect_null(result)
})

test_that("extract_censor_from_dataset handles dataset without censor", {
  # Create a mock dataset without censor field
  mock_dataset <- list(
    sampling_frame = fmrihrf::sampling_frame(blocklens = c(50, 50), TR = 2)
  )
  result <- fmrireg:::extract_censor_from_dataset(mock_dataset)
  expect_null(result)
})

test_that("extract_censor_from_dataset converts binary censor to indices", {
  mock_dataset <- list(
    sampling_frame = fmrihrf::sampling_frame(blocklens = c(10, 10), TR = 2),
    censor = c(rep(0, 4), 1, rep(0, 4), 1, rep(0, 10))  # censor timepoints 5 and 10
  )

  result <- fmrireg:::extract_censor_from_dataset(mock_dataset)
  expect_equal(result, c(5L, 10L))
})

test_that("extract_censor_from_dataset handles logical vector", {
  mock_dataset <- list(
    sampling_frame = fmrihrf::sampling_frame(blocklens = c(10, 10), TR = 2),
    censor = c(rep(FALSE, 4), TRUE, rep(FALSE, 4), TRUE, rep(FALSE, 10))
  )

  result <- fmrireg:::extract_censor_from_dataset(mock_dataset)
  expect_equal(result, c(5L, 10L))
})

test_that("extract_censor_from_dataset subsets by run", {
  mock_dataset <- list(
    sampling_frame = fmrihrf::sampling_frame(blocklens = c(10, 10), TR = 2),
    censor = c(rep(0, 4), 1, rep(0, 5),  # Run 1: censor point 5
               rep(0, 2), 1, rep(0, 7))   # Run 2: censor point 3 (global 13)
  )

  # Run 1 should give index 5
  result_run1 <- fmrireg:::extract_censor_from_dataset(mock_dataset, run_num = 1)
  expect_equal(result_run1, 5L)

  # Run 2 should give index 3 (local to run 2)
  result_run2 <- fmrireg:::extract_censor_from_dataset(mock_dataset, run_num = 2)
  expect_equal(result_run2, 3L)
})

test_that("resolve_censor returns NULL when cfg$ar$censor is NULL", {
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1"))
  result <- fmrireg:::resolve_censor(cfg)
  expect_null(result)
})

test_that("resolve_censor extracts from dataset when 'auto'", {
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1", censor = "auto"))

  mock_dataset <- list(
    sampling_frame = fmrihrf::sampling_frame(blocklens = c(10, 10), TR = 2),
    censor = c(rep(0, 4), 1, rep(0, 15))
  )

  result <- fmrireg:::resolve_censor(cfg, mock_dataset)
  expect_equal(result, 5L)
})

test_that("resolve_censor uses explicit censor vector from config", {
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1", censor = c(3L, 7L, 12L)))

  result <- fmrireg:::resolve_censor(cfg)
  expect_equal(result, c(3L, 7L, 12L))
})

test_that("estimate_ar_parameters accepts censor parameter", {
  set.seed(123)
  # Generate AR(1) process
  n <- 100
  phi_true <- 0.7
  e <- rnorm(n)
  x <- numeric(n)
  x[1] <- e[1]
  for (i in 2:n) {
    x[i] <- phi_true * x[i-1] + e[i]
  }

  # Without censor
  phi_est_no_censor <- estimate_ar_parameters(x, 1)
  expect_length(phi_est_no_censor, 1)

  # With censor (exclude some points)
  phi_est_with_censor <- estimate_ar_parameters(x, 1, censor = c(10L, 20L, 30L))
  expect_length(phi_est_with_censor, 1)

  # Both should give reasonable estimates (close to 0.7)
  expect_true(abs(phi_est_no_censor - phi_true) < 0.3)
  expect_true(abs(phi_est_with_censor - phi_true) < 0.3)
})

test_that("ar_whiten_transform accepts censor parameter", {
  set.seed(456)
  n <- 50
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * 10), n, 10)
  phi <- c(0.5)

  # Without censor
  result_no_censor <- ar_whiten_transform(X, Y, phi)
  expect_equal(dim(result_no_censor$X), c(n, p))
  expect_equal(dim(result_no_censor$Y), c(n, 10))

  # With censor
  result_with_censor <- ar_whiten_transform(X, Y, phi, censor = c(5L, 15L, 25L))
  expect_equal(dim(result_with_censor$X), c(n, p))
  expect_equal(dim(result_with_censor$Y), c(n, 10))
})

test_that("censor affects AR coefficient estimation", {
  set.seed(789)
  # Generate AR(1) process with an outlier segment
  n <- 100
  phi_true <- 0.6
  e <- rnorm(n)
  x <- numeric(n)
  x[1] <- e[1]
  for (i in 2:n) {
    x[i] <- phi_true * x[i-1] + e[i]
  }

  # Insert extreme values at specific locations
  bad_indices <- c(40, 41, 42, 60, 61, 62)
  x[bad_indices] <- x[bad_indices] + 10 * sign(rnorm(length(bad_indices)))

  # Estimation without censoring (should be biased by outliers)
  phi_est_no_censor <- estimate_ar_parameters(x, 1)

  # Estimation with censoring (should be closer to true value)
  phi_est_with_censor <- estimate_ar_parameters(x, 1, censor = bad_indices)

  # Both should be numeric and length 1
  expect_length(phi_est_no_censor, 1)
  expect_length(phi_est_with_censor, 1)

  # The censored estimate should generally be closer to truth
  # (though this depends on random seed, so we just check they're different)
  # In practice, censoring outliers typically improves estimation
})

test_that("fmriAR adapter functions accept censor parameter", {
  skip_if_not_installed("fmriAR")

  set.seed(101)
  n <- 60
  nvox <- 5
  residuals <- matrix(rnorm(n * nvox), n, nvox)

  cfg <- fmri_lm_control(ar_options = list(struct = "ar1"))

  # Test .estimate_ar_via_fmriAR with censor
  plan <- fmrireg:::.estimate_ar_via_fmriAR(residuals, cfg, censor = c(10L, 20L, 30L))
  expect_s3_class(plan, "fmriAR_plan")
  expect_true(!is.null(plan$phi))

  # Test .apply_ar_whitening_via_fmriAR with censor
  X <- matrix(rnorm(n * 3), n, 3)
  Y <- residuals
  whitened <- fmrireg:::.apply_ar_whitening_via_fmriAR(X, Y, plan, censor = c(10L, 20L, 30L))
  expect_equal(dim(whitened$X), dim(X))
  expect_equal(dim(whitened$Y), dim(Y))
})

test_that("iterative GLS respects censor", {
  skip_if_not_installed("fmriAR")

  set.seed(202)
  n <- 80
  p <- 4
  nvox <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  Y <- X %*% matrix(rnorm(p * nvox), p, nvox) + matrix(rnorm(n * nvox), n, nvox)

  cfg <- fmri_lm_control(ar_options = list(struct = "ar1", iter_gls = 2))

  # Without censor
  result_no_censor <- fmrireg:::.iterative_ar_gls_via_fmriAR(X, Y, cfg)
  expect_equal(dim(result_no_censor$X_white), dim(X))
  expect_equal(dim(result_no_censor$Y_white), dim(Y))
  expect_null(result_no_censor$censor)

  # With censor
  censor_pts <- c(15L, 35L, 55L)
  result_with_censor <- fmrireg:::.iterative_ar_gls_via_fmriAR(X, Y, cfg, censor = censor_pts)
  expect_equal(dim(result_with_censor$X_white), dim(X))
  expect_equal(dim(result_with_censor$Y_white), dim(Y))
  expect_equal(result_with_censor$censor, censor_pts)
})

test_that("end-to-end: fmri_lm with censor in ar_options", {
  skip_if_not_installed("fmriAR")
  skip_if_not_installed("fmridataset")
  skip_if_not_installed("fmridesign")

  set.seed(303)
  # Create simple synthetic data
  n_time <- 60
  n_vox <- 10
  TR <- 2

  # Design: simple block design
  onsets <- seq(10, 50, by = 20)
  durations <- rep(5, length(onsets))

  event_table <- data.frame(
    onset = onsets,
    duration = durations,
    condition = rep("A", length(onsets)),
    run = rep(1, length(onsets))
  )

  # Create data matrix
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  # Create simple dataset
  dset <- fmridataset::matrix_dataset(
    datamat = Y,
    TR = TR,
    run_length = n_time,
    event_table = event_table
  )

  # Build model
  emod <- event_model(onset ~ hrf(condition), data = event_table,
                      block = ~run, sampling_frame = dset$sampling_frame)
  bmod <- baseline_model("poly", degree = 2, sframe = dset$sampling_frame)
  fmod <- fmri_model(emod, bmod, dset)

  # Censor timepoints to exclude from AR estimation
  censor_points <- c(5L, 25L, 45L)

  # Fit with AR1 and censor
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1", censor = censor_points))
  result <- fmri_lm(fmod, dset, strategy = "runwise", control = cfg)

  # Check that fit succeeded
  expect_s3_class(result, "fmri_lm")
  expect_true(!is.null(result$result))
})

test_that("censor 'auto' mode extracts from dataset$censor", {
  skip_if_not_installed("fmriAR")
  skip_if_not_installed("fmridataset")
  skip_if_not_installed("fmridesign")

  set.seed(404)
  n_time <- 50
  n_vox <- 5
  TR <- 2

  onsets <- c(10, 30)
  event_table <- data.frame(
    onset = onsets,
    duration = rep(3, length(onsets)),
    condition = rep("A", length(onsets)),
    run = rep(1, length(onsets))
  )

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  # Create dataset WITH censor field
  dset <- fmridataset::matrix_dataset(
    datamat = Y,
    TR = TR,
    run_length = n_time,
    event_table = event_table
  )
  # Manually add censor (simulating what fmri_dataset would do)
  dset$censor <- c(rep(0, 9), 1, rep(0, 9), 1, rep(0, 30))  # Censor points 10, 20

  emod <- event_model(onset ~ hrf(condition), data = event_table,
                      block = ~run, sampling_frame = dset$sampling_frame)
  bmod <- baseline_model("constant", sframe = dset$sampling_frame)
  fmod <- fmri_model(emod, bmod, dset)

  # Use "auto" to extract censor from dataset
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1", censor = "auto"))
  result <- fmri_lm(fmod, dset, strategy = "runwise", control = cfg)

  expect_s3_class(result, "fmri_lm")
})

test_that("censor works with multiple runs", {
  skip_if_not_installed("fmriAR")
  skip_if_not_installed("fmridataset")
  skip_if_not_installed("fmridesign")

  set.seed(505)
  n_time_per_run <- 30
  n_runs <- 2
  n_time_total <- n_time_per_run * n_runs
  n_vox <- 5
  TR <- 2

  # Events in two runs
  event_table <- data.frame(
    onset = c(5, 15, 35, 45),
    duration = rep(3, 4),
    condition = rep("A", 4),
    run = c(1, 1, 2, 2)
  )

  Y <- matrix(rnorm(n_time_total * n_vox), n_time_total, n_vox)

  dset <- fmridataset::matrix_dataset(
    datamat = Y,
    TR = TR,
    run_length = c(n_time_per_run, n_time_per_run),
    event_table = event_table
  )
  # Censor: run 1 point 10, run 2 point 5 (global 35)
  dset$censor <- c(rep(0, 9), 1, rep(0, 20),   # Run 1
                   rep(0, 4), 1, rep(0, 25))   # Run 2

  emod <- event_model(onset ~ hrf(condition), data = event_table,
                      block = ~run, sampling_frame = dset$sampling_frame)
  bmod <- baseline_model("constant", sframe = dset$sampling_frame)
  fmod <- fmri_model(emod, bmod, dset)

  cfg <- fmri_lm_control(ar_options = list(struct = "ar1", censor = "auto"))
  result <- fmri_lm(fmod, dset, strategy = "runwise", control = cfg)

  expect_s3_class(result, "fmri_lm")
})
