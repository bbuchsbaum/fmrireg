test_that("calibration harness reruns the complete inference pipeline", {
  source(validation_script("fmri_lm_inference_calibration.R"), local = TRUE)
  out <- run_fmri_lm_inference_calibration(
    nsim = 4L, n = 70L, phi = 0.3,
    variance_method = "model", seed = 8110L,
    absolute_tolerance = 1
  )

  expect_s3_class(out, "fmri_lm_calibration")
  expect_identical(out$settings$nsim, 4L)
  expect_equal(nrow(out$draws), 4L)
  expect_named(out$metrics,
               c("metric", "target", "estimate", "lower_99", "upper_99", "pass"))
  expect_true(all(is.finite(out$draws)))
})
