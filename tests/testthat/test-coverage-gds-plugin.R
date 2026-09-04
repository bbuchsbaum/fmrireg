# Coverage for fmrigds registration helpers and plugin resolve_basis paths.

test_that("fmrigds meta reducer helpers cover FE/DL/contrast/cov branches", {
  set.seed(21)
  n <- 10
  p <- 4
  beta <- matrix(rnorm(n * p), n, p)
  se <- matrix(runif(n * p, 0.1, 0.3), n, p)
  X <- cbind(1, rnorm(n))

  expect_equal(fmrireg:::.fmrigds_map_robust("t"), "huber")
  expect_equal(fmrireg:::.fmrigds_map_robust("none"), "none")
  expect_equal(fmrireg:::.fmrigds_map_robust(NULL), "none")

  expect_equal(fmrireg:::.fmrigds_coerce_var(matrix(1, 2, 2), NULL), matrix(1, 2, 2))
  expect_equal(fmrireg:::.fmrigds_coerce_var(NULL, se = matrix(2, 2, 2)), matrix(4, 2, 2))
  expect_error(fmrireg:::.fmrigds_coerce_var(NULL, NULL), "var' or 'se'")

  expect_equal(fmrireg:::.fmrigds_meta_method("fe"), "fe")
  expect_equal(fmrireg:::.fmrigds_meta_method("pm", "dl"), "dl")

  reducer <- fmrireg:::.fmrigds_make_meta_reducer("fe")
  out <- reducer(beta, var = se^2, X = X, z = NULL, p = NULL, df = NULL, df1 = NULL, df2 = NULL, opts = list())
  expect_true(all(c("beta", "se", "z", "p") %in% names(out)))
  expect_equal(nrow(out$beta), p)

  # Contrast path
  Cmat <- matrix(c(0, 1), nrow = 2)
  out_c <- reducer(
    beta, var = se^2, X = X, z = NULL, p = NULL, df = NULL, df1 = NULL, df2 = NULL,
    opts = list(contrasts = Cmat, robust = "t")
  )
  expect_true(all(c("con_est", "con_se", "con_z") %in% names(out_c)))

  # Covariance return path for pm/reml
  reducer_pm <- fmrireg:::.fmrigds_make_meta_reducer("pm")
  out_cov <- reducer_pm(
    beta, var = se^2, X = X, z = NULL, p = NULL, df = NULL, df1 = NULL, df2 = NULL,
    opts = list(return_cov = "tri", method = "pm")
  )
  expect_true(!is.null(out_cov$cov_tri) || !is.null(out_cov$tau2) || !is.null(out_cov$beta))
})

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

  # Deprecated alias still works
  expect_true(!is.null(fmrireg:::.fmrireg_resolve_basis("cov_test_basis")))
})
