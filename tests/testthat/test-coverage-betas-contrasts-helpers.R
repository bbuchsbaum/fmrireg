# Twelfth wave: ols_betas, gen_beta_design, sandwich HC2/HC3.

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
  dset <- matrix_frame(
    matrix(rnorm(n * 2), n, 2),
    TR = 1, run_length = n, event_table = etab
  )
  bmod <- baseline_model(basis = "constant", sframe = frame_sframe(dset))

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

