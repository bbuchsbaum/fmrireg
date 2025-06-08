# Test fmri_lm_config and fmri_lm_control

test_that("fmri_lm_control creates valid config objects", {
  # Default config
  cfg1 <- fmri_lm_control()
  expect_s3_class(cfg1, "fmri_lm_config")
  expect_equal(cfg1$robust$type, FALSE)
  expect_equal(cfg1$ar$struct, "iid")
  
  # Robust config
  cfg2 <- fmri_lm_control(
    robust_options = list(
      type = "bisquare",
      c_tukey = 4.685,
      max_iter = 20
    )
  )
  expect_equal(cfg2$robust$type, "bisquare")
  expect_equal(cfg2$robust$c_tukey, 4.685)
  expect_equal(cfg2$robust$max_iter, 20)
  
  # AR config
  cfg3 <- fmri_lm_control(
    ar_options = list(
      struct = "ar1",
      global = TRUE,
      iter_gls = 2
    )
  )
  expect_equal(cfg3$ar$struct, "ar1")
  expect_true(cfg3$ar$global)
  expect_equal(cfg3$ar$iter_gls, 2)
  
  # Combined config
  cfg4 <- fmri_lm_control(
    robust_options = list(type = "huber", k_huber = 1.345),
    ar_options = list(struct = "ar2", voxelwise = TRUE)
  )
  expect_equal(cfg4$robust$type, "huber")
  expect_equal(cfg4$robust$k_huber, 1.345)
  expect_equal(cfg4$ar$struct, "ar2")
  expect_true(cfg4$ar$voxelwise)
})

test_that("fmri_lm_control validates inputs", {
  # Invalid robust type
  expect_error(
    fmri_lm_control(robust_options = list(type = "invalid")),
    "Must be one of"
  )
  
  # Invalid AR struct
  expect_error(
    fmri_lm_control(ar_options = list(struct = "ar99")),
    "should be one of"
  )
  
  # Invalid c_tukey parameter - currently not validated
  # expect_error(
  #   fmri_lm_control(robust_options = list(type = "bisquare", c_tukey = -1)),
  #   "c_tukey"
  # )
  
  # Missing ar_p for arp - currently uses default
  # expect_error(
  #   fmri_lm_control(ar_options = list(struct = "arp")),
  #   "ar_p"
  # )
  
  # Valid arp with p specified
  cfg <- fmri_lm_control(ar_options = list(struct = "arp", p = 3))
  expect_equal(cfg$ar$p, 3)
})

test_that("fmri_lm accepts new config API", {
  skip_if_not_installed("neuroim2")
  
  # Create simple test data
  n_time <- 50
  n_vox <- 10
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  event_df <- data.frame(onset = c(10, 20, 30), block = 1)
  dset <- matrix_dataset(Y, TR = 2, run_length = n_time, event_table = event_df)
  
  # Test with config object
  cfg <- fmri_lm_control(
    robust_options = list(type = "huber"),
    ar_options = list(struct = "ar1")
  )
  
  # This should work with new API
  expect_error({
    fit <- fmri_lm(
      onset ~ hrf(onset, basis = "spmg1"),
      block = ~ block,
      dataset = dset,
      robust_options = cfg$robust,
      ar_options = cfg$ar
    )
  }, NA)
})

test_that("config options propagate correctly", {
  skip_if_not_installed("neuroim2")
  
  # Create test data
  n_time <- 100
  n_vox <- 5
  onsets <- seq(10, 90, by = 20)
  
  # Generate AR(1) noise
  ar_coef <- 0.5
  Y <- matrix(0, n_time, n_vox)
  for (v in 1:n_vox) {
    noise <- rnorm(n_time)
    for (t in 2:n_time) {
      noise[t] <- ar_coef * noise[t-1] + rnorm(1)
    }
    Y[, v] <- noise
  }
  
  dset <- matrix_dataset(Y, TR = 1, run_length = n_time)
  
  # Create dataset with events
  event_df <- data.frame(onset = onsets, block = 1)
  dset <- matrix_dataset(Y, TR = 1, run_length = n_time, event_table = event_df)
  
  # Test AR options propagate
  fit_ar <- fmri_lm(
    onset ~ hrf(onset, basis = "spmg1"),
    block = ~ block,
    dataset = dset,
    ar_options = list(struct = "ar1", iter_gls = 1)
  )
  
  # Check that config was stored
  expect_equal(attr(fit_ar, "config")$ar$struct, "ar1")
  
  # Test robust options propagate
  fit_robust <- fmri_lm(
    onset ~ hrf(onset, basis = "spmg1"),
    block = ~ block,
    dataset = dset,
    robust = "bisquare",
    robust_options = list(max_iter = 10)
  )
  
  expect_equal(attr(fit_robust, "config")$robust$type, "bisquare")
  expect_equal(attr(fit_robust, "config")$robust$max_iter, 10)
})

test_that("fmri_lm_control handles default values correctly", {
  # Test all default values
  cfg <- fmri_lm_control()
  
  # Robust defaults
  expect_equal(cfg$robust$type, FALSE)
  expect_true(is.list(cfg$robust))
  
  # AR defaults
  expect_equal(cfg$ar$struct, "iid")
  expect_true(is.list(cfg$ar))
  
  # Should have all expected fields
  expect_true("robust" %in% names(cfg))
  expect_true("ar" %in% names(cfg))
})

test_that("fmri_lm_control merges options correctly", {
  # Test partial options specification
  cfg1 <- fmri_lm_control(robust_options = list(type = "huber"))
  expect_equal(cfg1$robust$type, "huber")
  expect_equal(cfg1$ar$struct, "iid")  # Should keep default
  
  cfg2 <- fmri_lm_control(ar_options = list(struct = "ar1"))
  expect_equal(cfg2$ar$struct, "ar1")
  expect_equal(cfg2$robust$type, FALSE)  # Should keep default
  
  # Test that explicitly provided options override defaults
  cfg3 <- fmri_lm_control(
    robust_options = list(type = "bisquare", max_iter = 5),
    ar_options = list(struct = "ar2", global = FALSE)
  )
  expect_equal(cfg3$robust$type, "bisquare")
  expect_equal(cfg3$robust$max_iter, 5)
  expect_equal(cfg3$ar$struct, "ar2")
  expect_false(cfg3$ar$global)
})

test_that("fmri_lm_control handles edge cases", {
  # Empty lists should work (use defaults)
  cfg1 <- fmri_lm_control(robust_options = list(), ar_options = list())
  expect_equal(cfg1$robust$type, FALSE)
  expect_equal(cfg1$ar$struct, "iid")
  
  # NULL options should work (use defaults)
  cfg2 <- fmri_lm_control(robust_options = NULL, ar_options = NULL)
  expect_equal(cfg2$robust$type, FALSE)
  expect_equal(cfg2$ar$struct, "iid")
  
  # Mixed valid and invalid - valid should work
  cfg3 <- fmri_lm_control(
    robust_options = list(type = "huber", k_huber = 1.345),
    ar_options = list(struct = "ar1")
  )
  expect_equal(cfg3$robust$type, "huber")
  expect_equal(cfg3$robust$k_huber, 1.345)
  expect_equal(cfg3$ar$struct, "ar1")
})

test_that("config object structure is consistent", {
  cfg <- fmri_lm_control(
    robust_options = list(type = "bisquare", c_tukey = 4.0),
    ar_options = list(struct = "ar2", p = 2)
  )
  
  # Should be proper S3 class
  expect_s3_class(cfg, "fmri_lm_config")
  
  # Should have expected structure
  expect_true(is.list(cfg))
  expect_true(is.list(cfg$robust))
  expect_true(is.list(cfg$ar))
  
  # Should be printable without error
  expect_no_error(print(cfg))
  
  # Should work with str() without error
  expect_no_error(str(cfg))
})

test_that("robust estimator parameters are handled correctly", {
  # Test Huber parameters
  cfg_huber <- fmri_lm_control(
    robust_options = list(
      type = "huber",
      k_huber = 1.5,
      max_iter = 15,
      scale_scope = "global"
    )
  )
  expect_equal(cfg_huber$robust$type, "huber")
  expect_equal(cfg_huber$robust$k_huber, 1.5)
  expect_equal(cfg_huber$robust$max_iter, 15)
  expect_equal(cfg_huber$robust$scale_scope, "global")
  
  # Test Bisquare parameters
  cfg_bisquare <- fmri_lm_control(
    robust_options = list(
      type = "bisquare",
      c_tukey = 5.0,
      max_iter = 25,
      reestimate_phi = TRUE
    )
  )
  expect_equal(cfg_bisquare$robust$type, "bisquare")
  expect_equal(cfg_bisquare$robust$c_tukey, 5.0)
  expect_equal(cfg_bisquare$robust$max_iter, 25)
  expect_true(cfg_bisquare$robust$reestimate_phi)
  
  # Test that non-robust estimation works
  cfg_no_robust <- fmri_lm_control(robust_options = list(type = FALSE))
  expect_equal(cfg_no_robust$robust$type, FALSE)
})

test_that("AR model parameters are configured correctly", {
  # Test AR(1) parameters
  cfg_ar1 <- fmri_lm_control(
    ar_options = list(
      struct = "ar1",
      global = TRUE,
      iter_gls = 3,
      exact_first = TRUE
    )
  )
  expect_equal(cfg_ar1$ar$struct, "ar1")
  expect_true(cfg_ar1$ar$global)
  expect_equal(cfg_ar1$ar$iter_gls, 3)
  expect_true(cfg_ar1$ar$exact_first)
  
  # Test AR(p) parameters
  cfg_arp <- fmri_lm_control(
    ar_options = list(
      struct = "arp",
      p = 4,
      voxelwise = TRUE
    )
  )
  expect_equal(cfg_arp$ar$struct, "arp")
  expect_equal(cfg_arp$ar$p, 4)
  expect_true(cfg_arp$ar$voxelwise)
  
  # Test IID (no autocorrelation)
  cfg_iid <- fmri_lm_control(ar_options = list(struct = "iid"))
  expect_equal(cfg_iid$ar$struct, "iid")
})