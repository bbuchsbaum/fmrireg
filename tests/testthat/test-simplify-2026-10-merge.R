# Simplification 2026-10: checks added while merging the parallel streams.

test_that("a failed AR estimate warns instead of silently dropping whitening", {
  local_mocked_bindings(.estimate_ar_via_fmriAR = function(...) stop("boom"))
  set.seed(1)
  r <- matrix(rnorm(200), 50, 4)
  expect_warning(
    phi <- .estimate_ar_parameters_routed(r, 1L),
    class = "fmrireg_ar_estimation_failed"
  )
  expect_identical(phi, 0)
})
