# Coverage for the plugin resolve_basis paths.

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
