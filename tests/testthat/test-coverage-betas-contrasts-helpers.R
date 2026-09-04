# Twelfth wave: ols_betas, gen_beta_design, construct_contrast_matrix, sandwich HC2/HC3.

test_that("ols_betas and gen_beta_design cover fixed/null-fixed layouts", {
  set.seed(211)
  n <- 60L
  V <- 3L
  X <- cbind(1, rnorm(n), rnorm(n))
  Y <- matrix(as.numeric(X %*% c(0.2, 0.8, -0.4)), n, V) +
    matrix(rnorm(n * V, sd = 0.2), n, V)
  B <- fmrireg:::ols_betas(X, Y)
  expect_equal(dim(B), c(3L, V))
  expect_equal(unname(B[, 1]), unname(qr.coef(qr(X), Y[, 1])), tolerance = 1e-8)

  etab <- data.frame(
    onset = seq(5, n - 10, length.out = 8),
    condition = factor(rep(c("A", "B"), 4)),
    nuisance = rnorm(8),
    run = 1L
  )
  dset <- matrix_dataset(
    matrix(rnorm(n * 2), n, 2),
    TR = 1, run_length = n, event_table = etab
  )
  bmod <- baseline_model(basis = "constant", sframe = dset$sampling_frame)

  with_fixed <- fmrireg:::gen_beta_design(
    fixed = onset ~ hrf(nuisance),
    ran = onset ~ hrf(condition),
    block = ~ run,
    bmod = bmod,
    dset = dset
  )
  expect_true(!is.null(with_fixed$emod_fixed))
  expect_true(is.matrix(as.matrix(with_fixed$dmat_fixed)) ||
                inherits(with_fixed$dmat_fixed, "tbl_df"))
  expect_true(length(with_fixed$fixed_ind) >= 1L)
  expect_true(length(with_fixed$ran_ind) >= 1L)
  expect_true(length(with_fixed$base_ind) >= 1L)

  no_fixed <- fmrireg:::gen_beta_design(
    fixed = NULL,
    ran = onset ~ hrf(condition),
    block = ~ run,
    bmod = bmod,
    dset = dset
  )
  expect_null(no_fixed$emod_fixed)
  expect_null(no_fixed$fixed_ind)
  expect_equal(no_fixed$ran_ind, seq_len(ncol(no_fixed$dmat_ran)))
})

test_that("construct_contrast_matrix and compute_contrast/f-stat edges", {
  expect_equal(
    fmrireg:::construct_contrast_matrix(diag(2), list(), c("a", "b")),
    diag(2)
  )
  expect_equal(
    fmrireg:::construct_contrast_matrix(c(1, -1, 0), list(), letters[1:3]),
    matrix(c(1, -1, 0), nrow = 1)
  )
  expect_error(
    fmrireg:::construct_contrast_matrix(list(a = 1), list(), "a"),
    "not yet implemented"
  )

  set.seed(212)
  n <- 40L
  p <- 3L
  V <- 5L
  X <- cbind(1, rnorm(n), rnorm(n))
  Betas <- matrix(rnorm(p * V), p, V)
  sigma2 <- runif(V, 0.2, 0.8)
  fit_result <- list(
    betas = Betas,
    XtXinv = chol2inv(chol(crossprod(X))),
    sigma2 = sigma2,
    dfres = n - p
  )

  # Vector contrast expands to 1-row matrix
  tcon <- fmrireg:::compute_contrast(fit_result, c(0, 1, -1), contrast_names = "AvB")
  expect_equal(rownames(tcon$estimate), "AvB")
  expect_equal(dim(tcon$estimate), c(1L, V))
  expect_true(all(is.finite(tcon$tstat)))

  # Single sigma2 scalar branch
  fit_scalar <- fit_result
  fit_scalar$sigma2 <- 0.5
  t2 <- fmrireg:::compute_contrast(fit_scalar, matrix(c(0, 1, 0), 1, 3))
  expect_equal(ncol(t2$stderr), V)

  expect_error(
    fmrireg:::compute_contrast(fit_result, matrix(1, 1, 2)),
    "must match number of parameters"
  )
  expect_error(
    fmrireg:::compute_contrast(list(betas = Betas, sigma2 = 1), c(1, 0, 0)),
    "variance-covariance"
  )

  fstat <- fmrireg:::compute_f_statistic(fit_result, rbind(c(0, 1, 0), c(0, 0, 1)))
  expect_equal(dim(fstat$fstat), c(1L, V))
  expect_equal(fstat$df1, 2)

  # Parallel message branch
  old <- getOption("fmrireg.num_threads")
  on.exit(options(fmrireg.num_threads = old), add = TRUE)
  options(fmrireg.num_threads = 2)
  expect_message(
    fmrireg:::compute_voxelwise_contrasts(
      fit_result,
      list(matrix(c(0, 1, -1), 1, 3), rbind(c(0, 1, 0), c(0, 0, 1))),
      parallel = TRUE
    ),
    "Parallel contrast"
  )
})

test_that("compute_sandwich_variance covers HC0/HC2/HC3 multi-voxel branches", {
  set.seed(213)
  n <- 30L
  p <- 3L
  X <- cbind(1, rnorm(n), rnorm(n))
  resid1 <- matrix(rnorm(n), n, 1)
  hc0 <- fmrireg:::compute_sandwich_variance(X, resid1, type = "HC0")
  expect_equal(dim(hc0), c(p, p))
  expect_true(all(is.finite(hc0)))

  resid_multi <- matrix(rnorm(n * 4), n, 4)
  expect_warning(
    hc2 <- fmrireg:::compute_sandwich_variance(X, resid_multi, type = "HC2"),
    "HC2 not implemented"
  )
  expect_equal(dim(hc2), c(p, p))
  expect_warning(
    hc3 <- fmrireg:::compute_sandwich_variance(X, resid_multi, type = "HC3"),
    "HC3 not implemented"
  )
  expect_equal(dim(hc3), c(p, p))
})
