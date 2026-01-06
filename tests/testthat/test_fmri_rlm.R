# Test fmri_rlm wrapper function

library(testthat)
library(fmrireg)

# ============================================================================
# Setup: Create test datasets
# ============================================================================

create_test_dataset <- function(n_timepoints = 100, n_voxels = 50, n_runs = 2) {
  run_length <- n_timepoints / n_runs

  # Create event table
  event_table <- data.frame(
    onset = c(10, 30, 60, 80),
    fac = factor(c("A", "B", "A", "B")),
    run = c(1, 1, 2, 2)
  )

  # Create random data matrix
  mat <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  # Create dataset
  fmridataset::matrix_dataset(
    mat,
    TR = 1,
    run_length = rep(run_length, n_runs),
    event_table = event_table
  )
}

# ============================================================================
# Basic functionality tests
# ============================================================================

test_that("fmri_rlm runs without error", {
  dset <- create_test_dataset()

  expect_no_error(
    lm_robust <- fmri_rlm(onset ~ hrf(fac), block = ~ run, dataset = dset)
  )

  expect_s3_class(lm_robust, "fmri_rlm")
  expect_s3_class(lm_robust, "fmri_lm")  # Should inherit from fmri_lm
})

test_that("fmri_rlm produces equivalent results to fmri_lm with robust=TRUE", {
  dset <- create_test_dataset(n_voxels = 20)

  # Fit using fmri_rlm wrapper
  lm_rlm <- fmri_rlm(
    onset ~ hrf(fac),
    block = ~ run,
    dataset = dset,
    robust_psi = "huber"
  )

  # Fit using fmri_lm with robust options
  lm_robust <- fmri_lm(
    onset ~ hrf(fac),
    block = ~ run,
    dataset = dset,
    robust = TRUE,
    robust_options = list(type = "huber")
  )

  # Results should be very similar
  expect_equal(coef(lm_rlm), coef(lm_robust), tolerance = 1e-6)
})

test_that("fmri_rlm has correct class structure", {
  dset <- create_test_dataset()
  lm_robust <- fmri_rlm(onset ~ hrf(fac), block = ~ run, dataset = dset)

  classes <- class(lm_robust)
  expect_true("fmri_rlm" %in% classes)
  expect_true("fmri_lm" %in% classes)

  # fmri_rlm should be first (for method dispatch)
  expect_equal(classes[1], "fmri_rlm")
})

# ============================================================================
# Test robust options
# ============================================================================

test_that("fmri_rlm works with Huber psi function", {
  dset <- create_test_dataset()

  expect_no_error(
    lm_huber <- fmri_rlm(
      onset ~ hrf(fac),
      block = ~ run,
      dataset = dset,
      robust_psi = "huber",
      robust_k_huber = 1.345
    )
  )

  expect_s3_class(lm_huber, "fmri_rlm")
})

test_that("fmri_rlm works with bisquare psi function", {
  dset <- create_test_dataset()

  expect_no_error(
    lm_bisquare <- fmri_rlm(
      onset ~ hrf(fac),
      block = ~ run,
      dataset = dset,
      robust_psi = "bisquare",
      robust_c_tukey = 4.685
    )
  )

  expect_s3_class(lm_bisquare, "fmri_rlm")
})

test_that("fmri_rlm respects robust_max_iter parameter", {
  dset <- create_test_dataset()

  # Single iteration
  lm_1iter <- fmri_rlm(
    onset ~ hrf(fac),
    block = ~ run,
    dataset = dset,
    robust_max_iter = 1
  )

  # Multiple iterations
  lm_5iter <- fmri_rlm(
    onset ~ hrf(fac),
    block = ~ run,
    dataset = dset,
    robust_max_iter = 5
  )

  # Both should complete
  expect_s3_class(lm_1iter, "fmri_rlm")
  expect_s3_class(lm_5iter, "fmri_rlm")
})

test_that("fmri_rlm handles different robust_scale_scope", {
  dset <- create_test_dataset()

  # Run-level scale estimation
  lm_run <- fmri_rlm(
    onset ~ hrf(fac),
    block = ~ run,
    dataset = dset,
    robust_scale_scope = "run"
  )

  # Global scale estimation
  lm_global <- fmri_rlm(
    onset ~ hrf(fac),
    block = ~ run,
    dataset = dset,
    robust_scale_scope = "global"
  )

  expect_s3_class(lm_run, "fmri_rlm")
  expect_s3_class(lm_global, "fmri_rlm")
})

# ============================================================================
# Test AR correlation structure options
# ============================================================================

test_that("fmri_rlm works with AR1 correlation structure", {
  dset <- create_test_dataset()

  expect_no_error(
    lm_ar1 <- fmri_rlm(
      onset ~ hrf(fac),
      block = ~ run,
      dataset = dset,
      cor_struct = "ar1"
    )
  )

  expect_s3_class(lm_ar1, "fmri_rlm")
})

test_that("fmri_rlm works with AR2 correlation structure", {
  dset <- create_test_dataset()

  expect_no_error(
    lm_ar2 <- fmri_rlm(
      onset ~ hrf(fac),
      block = ~ run,
      dataset = dset,
      cor_struct = "ar2"
    )
  )

  expect_s3_class(lm_ar2, "fmri_rlm")
})

test_that("fmri_rlm respects cor_iter parameter", {
  dset <- create_test_dataset()

  lm_iter1 <- fmri_rlm(
    onset ~ hrf(fac),
    block = ~ run,
    dataset = dset,
    cor_struct = "ar1",
    cor_iter = 1
  )

  lm_iter3 <- fmri_rlm(
    onset ~ hrf(fac),
    block = ~ run,
    dataset = dset,
    cor_struct = "ar1",
    cor_iter = 3
  )

  expect_s3_class(lm_iter1, "fmri_rlm")
  expect_s3_class(lm_iter3, "fmri_rlm")
})

test_that("fmri_rlm handles cor_global option", {
  dset <- create_test_dataset()

  # Local AR estimation
  lm_local <- fmri_rlm(
    onset ~ hrf(fac),
    block = ~ run,
    dataset = dset,
    cor_struct = "ar1",
    cor_global = FALSE
  )

  # Global AR estimation
  lm_global <- fmri_rlm(
    onset ~ hrf(fac),
    block = ~ run,
    dataset = dset,
    cor_struct = "ar1",
    cor_global = TRUE
  )

  expect_s3_class(lm_local, "fmri_rlm")
  expect_s3_class(lm_global, "fmri_rlm")
})

# ============================================================================
# Test strategy options
# ============================================================================

test_that("fmri_rlm works with runwise strategy", {
  dset <- create_test_dataset()

  expect_no_error(
    lm_runwise <- fmri_rlm(
      onset ~ hrf(fac),
      block = ~ run,
      dataset = dset,
      strategy = "runwise"
    )
  )

  expect_s3_class(lm_runwise, "fmri_rlm")
})

test_that("fmri_rlm works with chunkwise strategy", {
  dset <- create_test_dataset()

  expect_no_error(
    lm_chunkwise <- fmri_rlm(
      onset ~ hrf(fac),
      block = ~ run,
      dataset = dset,
      strategy = "chunkwise",
      nchunks = 5
    )
  )

  expect_s3_class(lm_chunkwise, "fmri_rlm")
})

# ============================================================================
# Test print method
# ============================================================================

test_that("print.fmri_rlm works correctly", {
  dset <- create_test_dataset()
  lm_robust <- fmri_rlm(onset ~ hrf(fac), block = ~ run, dataset = dset)

  # Capture print output
  output <- capture.output(print(lm_robust))

  # Should contain "Robust" indicator
  expect_true(any(grepl("[Rr]obust", output)))
})

# ============================================================================
# Test with single run dataset
# ============================================================================

test_that("fmri_rlm works with single run dataset", {
  event_table <- data.frame(
    onset = c(10, 30, 50, 70),
    fac = factor(c("A", "B", "A", "B")),
    run = rep(1, 4)
  )

  mat <- matrix(rnorm(100 * 30), 100, 30)
  dset <- fmridataset::matrix_dataset(
    mat,
    TR = 1,
    run_length = 100,
    event_table = event_table
  )

  expect_no_error(
    lm_single <- fmri_rlm(onset ~ hrf(fac), block = ~ run, dataset = dset)
  )

  expect_s3_class(lm_single, "fmri_rlm")
})

# ============================================================================
# Test combined robust + AR options
# ============================================================================

test_that("fmri_rlm handles combined robust and AR options", {
  dset <- create_test_dataset()

  expect_no_error(
    lm_combined <- fmri_rlm(
      onset ~ hrf(fac),
      block = ~ run,
      dataset = dset,
      robust_psi = "huber",
      robust_k_huber = 1.5,
      robust_max_iter = 3,
      cor_struct = "ar1",
      cor_iter = 2,
      cor_global = TRUE
    )
  )

  expect_s3_class(lm_combined, "fmri_rlm")
})

# ============================================================================
# Test accessor methods work on fmri_rlm objects
# ============================================================================

test_that("accessor methods work on fmri_rlm objects",
{
  dset <- create_test_dataset(n_voxels = 20)
  lm_robust <- fmri_rlm(onset ~ hrf(fac), block = ~ run, dataset = dset)

  # Test that standard accessors work
  expect_no_error(cf <- coef(lm_robust))
  expect_true(is.matrix(cf) || is.numeric(cf))

  expect_no_error(se <- standard_error(lm_robust))

  expect_no_error(st <- stats(lm_robust))
})

# ============================================================================
# Edge cases
# ============================================================================

test_that("fmri_rlm handles outliers appropriately", {
  # Create dataset with outliers
  event_table <- data.frame(
    onset = c(10, 30, 60, 80),
    fac = factor(c("A", "B", "A", "B")),
    run = c(1, 1, 2, 2)
  )

  mat <- matrix(rnorm(100 * 20), 100, 20)
  # Add some outliers
  mat[c(5, 25, 45), ] <- mat[c(5, 25, 45), ] * 10

  dset <- fmridataset::matrix_dataset(
    mat,
    TR = 1,
    run_length = c(50, 50),
    event_table = event_table
  )

  # Robust fitting should handle outliers
  expect_no_error(
    lm_robust <- fmri_rlm(onset ~ hrf(fac), block = ~ run, dataset = dset)
  )

  expect_s3_class(lm_robust, "fmri_rlm")
})
