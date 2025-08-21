context("afni hrf aliases")

test_that("afni_hrf accepts sin/gam aliases", {
  # afni_hrf needs a variable name
  spec_sin <- afni_hrf(cond, basis = "sin")
  spec_gam <- afni_hrf(cond, basis = "gam")

  # AFNI_HRF objects store the name as the base object value itself
  expect_equal(as.vector(spec_sin$hrf), "SIN")
  expect_equal(as.vector(spec_gam$hrf), "GAM")
})
