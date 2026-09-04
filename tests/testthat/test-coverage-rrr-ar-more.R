# Additional rrr bootstrap SE helper and fmriAR adapter branches.

test_that(".rrr_bootstrap_task_se returns finite SEs with and without blocks/contrasts", {
  set.seed(21)
  n <- 30L
  k <- 3L
  V <- 5L
  X0 <- cbind(1, rnorm(n), rnorm(n))
  B <- matrix(c(0.2, 0.5, -0.3), k, 1)
  Y0 <- X0 %*% matrix(rep(B, V), k, V) + matrix(rnorm(n * V, sd = 0.4), n, V)
  B_hat <- solve(crossprod(X0), crossprod(X0, Y0))
  rank_control <- fmrireg:::.rrr_normalize_args(list(rank = 2L, rank_mode = "fixed"))

  # Contrast remapped onto task columns
  w <- c(0, 1, -1)
  attr(w, "colind") <- 1:3
  t_cons <- list(A_vs_B = w)

  se_plain <- fmrireg:::.rrr_bootstrap_task_se(
    X0, Y0, B_hat, rank_control,
    t_contrasts_remapped = t_cons,
    nboot = 12L, block_size = 1L, seed = 5L
  )
  expect_equal(dim(se_plain$se_beta), dim(B_hat))
  expect_true(all(is.finite(se_plain$se_beta)))
  expect_true(all(se_plain$se_beta >= 0))
  expect_equal(names(se_plain$se_t), "A_vs_B")
  expect_equal(length(se_plain$se_t$A_vs_B), V)

  se_block <- fmrireg:::.rrr_bootstrap_task_se(
    X0, Y0, B_hat, rank_control,
    t_contrasts_remapped = list(),
    nboot = 8L, block_size = 5L, seed = 6L
  )
  expect_equal(dim(se_block$se_beta), dim(B_hat))
  expect_equal(length(se_block$se_t), 0L)
})

test_that("fmriAR adapter order/df/lag/run-label/whiten branches", {
  expect_equal(fmrireg:::.get_ar_order(plan = list(order = c(p = 2L, q = 0L))), 2L)
  expect_equal(fmrireg:::.get_ar_order(cfg = list(struct = "ar3")), 3L)
  expect_equal(fmrireg:::.get_ar_order(cfg = list(struct = "arp", p = 5L)), 5L)
  expect_equal(fmrireg:::.get_ar_order(), 0L)

  # Effective df: no AR -> n - p
  expect_equal(
    fmrireg:::.compute_ar_effective_df_compat(n = 50, p = 4, plan = NULL),
    46
  )
  df_ar1 <- fmrireg:::.compute_ar_effective_df_compat(
    n = 50, p = 4, ar_coef = 0.4
  )
  expect_true(df_ar1 < 46 && df_ar1 >= 1)

  df_multi <- fmrireg:::.compute_ar_effective_df_compat(
    n = 40, p = 3,
    plan = list(phi = list(c(0.3, 0.1), c(0.25, 0.05)), theta = list(0.1, 0.05))
  )
  expect_true(is.finite(df_multi) && df_multi >= 1)

  df_list_ar <- fmrireg:::.compute_ar_effective_df_compat(
    n = 40, p = 3, ar_coef = list(0.2, 0.3), ma_coef = list(0.05, 0.02)
  )
  expect_true(is.finite(df_list_ar))

  X <- cbind(1, rnorm(40))
  budget <- fmrireg:::.ar_correction_lag_budget(X, target_order = 2L)
  expect_true(is.integer(budget) || is.numeric(budget))
  expect_true(budget >= 1L)

  budget_runs <- fmrireg:::.ar_correction_lag_budget(
    X, target_order = 1L,
    run_indices = list(1:20, 21:40),
    censor = c(1L, 21L),
    max_lag = 10L
  )
  expect_true(budget_runs >= 1L && budget_runs <= 10L)

  # Invalid / edge run labels
  expect_null(fmrireg:::.build_run_labels(5L, list()))
  expect_error(
    fmrireg:::.build_run_labels(4L, list(1:2, 2:4)),
    "partition"
  )
  expect_error(
    fmrireg:::.build_run_labels(4L, list(1:2, c(4L, 3L))),
    "partition|contiguous"
  )

  expect_null(fmrireg:::.shared_ar_correction_design(matrix(1, 3, 2), list(), NULL))
  expect_equal(
    fmrireg:::.shared_ar_correction_design(matrix(1, 3, 2), list(struct = "ar1"), X[1:3, ]),
    X[1:3, ]
  )

  expect_error(
    fmrireg:::.block_diagonal_design(list()),
    "non-empty"
  )

  # phi fixed-order: order 0 and short blocks
  resid <- matrix(rnorm(20 * 3), 20, 3)
  expect_equal(
    fmrireg:::.estimate_phi_fixed_order(resid, order = 0L, pooling = "global"),
    list(numeric(0))
  )
  phi0_run <- fmrireg:::.estimate_phi_fixed_order(
    resid, order = 0L, pooling = "run", run_indices = list(1:10, 11:20)
  )
  expect_equal(length(phi0_run), 2L)

  skip_if_not_installed("fmriAR")

  plan <- fmrireg:::.estimate_ar_via_fmriAR(
    resid, list(struct = "ar1"), run_indices = NULL
  )
  expect_true(!is.null(plan$phi) || !is.null(plan$order))

  whitened <- fmrireg:::.apply_ar_whitening_via_fmriAR(
    X = cbind(1, rnorm(20)),
    Y = resid,
    plan = plan
  )
  expect_equal(dim(whitened$X), c(20L, 2L))
  expect_equal(dim(whitened$Y), dim(resid))

  # Missing theta filled in whitener
  plan_no_theta <- plan
  plan_no_theta$theta <- NULL
  whitened2 <- fmrireg:::.apply_ar_whitening_via_fmriAR(
    X = cbind(1, rnorm(20)), Y = resid, plan = plan_no_theta
  )
  expect_equal(dim(whitened2$Y), dim(resid))

  gls <- fmrireg:::.iterative_ar_gls_via_fmriAR(
    X = cbind(1, rnorm(20)),
    Y = resid,
    cfg = list(struct = "ar1", iter_gls = 1L),
    max_iter = 1L
  )
  expect_equal(dim(gls$X_white), c(20L, 2L))
  expect_equal(dim(gls$Y_white), dim(resid))

  # IID / n_iter=0 passthrough
  passthrough <- fmrireg:::.iterative_ar_gls_via_fmriAR(
    X = cbind(1, rnorm(20)),
    Y = resid,
    cfg = list(struct = "iid"),
    max_iter = 0L
  )
  expect_null(passthrough$plan)
  expect_equal(passthrough$iterations, 0L)
})
