# fmriAR_adapter.R: fixed-order early returns and zero fallback.

test_that(".estimate_phi_fixed_order covers short blocks and zero fallback", {
  skip_if_not_installed("fmriAR")

  # n_eff <= 1
  tiny <- matrix(1, 1, 3)
  expect_equal(
    fmrireg:::.estimate_phi_fixed_order(tiny, order = 2L, pooling = "global")[[1]],
    c(0, 0)
  )

  # Run pooling with a 1-row slice
  resid <- matrix(rnorm(20 * 2), 20, 2)
  phi_run <- fmrireg:::.estimate_phi_fixed_order(
    resid, order = 1L, pooling = "run",
    run_indices = list(1L, 2:20)
  )
  expect_equal(length(phi_run), 2L)
  expect_equal(phi_run[[1]], 0)

  # Force fmriAR failure: the fixed-order estimator falls back to zeros.
  suppressWarnings(
    trace(
      fmriAR::fit_noise,
      tracer = quote(stop("forced failure for coverage")),
      print = FALSE
    )
  )
  on.exit(try(untrace(fmriAR::fit_noise), silent = TRUE), add = TRUE)

  zeros <- fmrireg:::.estimate_phi_fixed_order(
    matrix(rnorm(40 * 3), 40, 3), order = 2L, pooling = "global"
  )
  expect_equal(zeros[[1]], c(0, 0))

  short <- matrix(rnorm(6), 3, 2)
  short_phi <- fmrireg:::.estimate_phi_fixed_order(
    short, order = 3L, pooling = "global"
  )
  expect_equal(short_phi[[1]], c(0, 0, 0))

  try(untrace(fmriAR::fit_noise), silent = TRUE)
})

test_that(".iterative_ar_gls_via_fmriAR runs multi-iter AR1 whitening", {
  skip_if_not_installed("fmriAR")
  set.seed(15)
  n <- 40L
  X <- cbind(1, rnorm(n))
  Y <- matrix(rnorm(n * 3), n, 3)

  gls <- fmrireg:::.iterative_ar_gls_via_fmriAR(
    X, Y,
    cfg = list(struct = "ar1", iter_gls = 2L),
    max_iter = 2L
  )
  expect_equal(dim(gls$X_white), dim(X))
  expect_equal(dim(gls$Y_white), dim(Y))
  expect_true(gls$iterations >= 1L)
  expect_true(!is.null(gls$ar_coef) || !is.null(gls$plan) || is.null(gls$plan))
})
