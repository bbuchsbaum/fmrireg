context("afni hrf aliases")

test_that("afni_hrf accepts sin/gam aliases", {
  spec_sin <- afni_hrf(x, basis = "sin")
  spec_gam <- afni_hrf(x, basis = "gam")

  expect_equal(attr(spec_sin$hrf, "name"), "SIN")
  expect_equal(attr(spec_gam$hrf, "name"), "GAM")
})
