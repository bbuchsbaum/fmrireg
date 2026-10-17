# Coverage for the fmrireg spatial-FDR posthoc and plugin resolve_basis paths.

test_that("fmrigds spatial posthoc helper derives q from p or z", {
  set.seed(3)
  P <- 20
  group <- rep(1:4, each = 5)
  p_arr <- array(runif(P * 1 * 2, 0.001, 0.2), dim = c(P, 1, 2))
  q1 <- fmrireg:::.fmrigds_spatial_posthoc_compute(list(p = p_arr), group, alpha = 0.1)
  expect_equal(dim(q1), dim(p_arr))
  expect_true(all(is.finite(q1) | is.na(q1)))

  z_arr <- array(rnorm(P * 1 * 1), dim = c(P, 1, 1))
  q2 <- fmrireg:::.fmrigds_spatial_posthoc_compute(list(z = z_arr), group, alpha = 0.1)
  expect_equal(dim(q2), dim(z_arr))

  expect_error(
    fmrireg:::.fmrigds_spatial_posthoc_compute(list(), group, 0.05),
    "requires p"
  )
  expect_error(
    fmrireg:::.fmrigds_spatial_posthoc(list(p = p_arr), opts = list()),
    "requires options\\$group"
  )
})

test_that("resolve_basis and register_basis cover plugin registry errors", {
  expect_error(resolve_basis("definitely_not_registered_xyz"), "Unknown basis")

  register_basis("cov_test_basis", function(...) {
    fmrihrf::HRF_SPMG1
  })
  resolved <- resolve_basis("cov_test_basis")
  expect_true(!is.null(resolved))

  expect_error(register_basis("", function(...) NULL))
  expect_error(register_basis("cov_test_basis2", "not-a-function"), "must be a function")
})
