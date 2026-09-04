# RRR engine: bootstrap+contrast success, filter policies, restore, preflight errors.

tiny_rrr <- function(n = 64L, V = 6L, seed = 51L, with_contrast = TRUE) {
  set.seed(seed)
  n_ev <- 8L
  onsets <- seq(5, by = max(4, floor((n - 10) / n_ev)), length.out = n_ev)
  onsets <- pmin(onsets, n - 3)
  etab <- data.frame(
    onset = onsets,
    condition = factor(rep(c("A", "B"), length.out = n_ev)),
    run = 1L
  )
  Y <- matrix(rnorm(n * V), n, V)
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = etab)
  form <- if (with_contrast) {
    con <- pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B")
    onset ~ hrf(condition, contrasts = list(con))
  } else {
    onset ~ hrf(condition)
  }
  emod <- event_model(form, data = etab, block = ~ run, sampling_frame = dset$sampling_frame)
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = dset$sampling_frame)
  list(model = fmri_model(emod, bmod, dset), dataset = dset)
}

test_that(".fit_rrr_gls_engine bootstrap SE overwrites contrast se/stat/prob", {
  fx <- tiny_rrr(with_contrast = TRUE)
  fit <- fmrireg:::.fit_rrr_gls_engine(
    fx$model, fx$dataset,
    args = list(
      rank = 2L,
      se_mode = "bootstrap",
      bootstrap_n = 8L,
      bootstrap_block_size = 4L,
      bootstrap_seed = 9L,
      contrast_policy = "warn_drop"
    ),
    cfg = fmri_lm_control()
  )
  expect_s3_class(fit, "fmri_lm")
  expect_equal(attr(fit, "engine"), "rrr_gls")
  expect_true(nrow(fit$result$contrasts) >= 1L)
  expect_equal(fit$rrr$se_mode, "bootstrap")
  expect_true(!is.null(fit$rrr$bootstrap$nboot))
  expect_equal(fit$result$variance_model$method, "bootstrap")
})

test_that(".fit_rrr_gls_engine rss_budget / energy rank + AR whitening", {
  skip_if_not_installed("fmriAR")
  fx <- tiny_rrr(n = 70L, V = 5L, seed = 52L, with_contrast = FALSE)
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1"))

  fit_rss <- fmrireg:::.fit_rrr_gls_engine(
    fx$model, fx$dataset,
    args = list(rank_mode = "rss_budget", rss_budget = 0.25, se_mode = "conditional"),
    cfg = cfg
  )
  expect_s3_class(fit_rss, "fmri_lm")
  expect_true(fit_rss$rrr$rank_used >= 0L)
  expect_true(!is.null(fit_rss$ar_coef) || !is.null(fit_rss$result$ar_coef))

  fit_en <- fmrireg:::.fit_rrr_gls_engine(
    fx$model, fx$dataset,
    args = list(rank_mode = "energy", energy_keep = 0.9, se_mode = "conditional"),
    cfg = fmri_lm_control()
  )
  expect_s3_class(fit_en, "fmri_lm")
  expect_true(isTRUE(fit_en$result$variance_model$metadata$adaptive_rank))
})

test_that(".rrr_filter_task_contrasts policies and restore metadata", {
  event_idx <- 1:2
  ok <- list(A_vs_B = structure(c(1, -1), colind = 1:2))
  bad <- list(bad_base = structure(c(1, 0, 1), colind = c(1L, 2L, 5L)))
  fbad <- list(Fbad = structure(diag(2), colind = c(1L, 4L)))

  keep <- fmrireg:::.rrr_filter_task_contrasts(ok, list(), event_idx, policy = "warn_drop")
  expect_equal(names(keep$simple), "A_vs_B")
  expect_equal(attr(keep$simple$A_vs_B, "colind"), 1:2)

  expect_warning(
    dropped <- fmrireg:::.rrr_filter_task_contrasts(bad, fbad, event_idx, policy = "warn_drop"),
    "dropped contrasts"
  )
  expect_true("bad_base" %in% dropped$dropped)
  expect_true("Fbad" %in% dropped$dropped)
  expect_equal(length(dropped$simple), 0L)

  silent <- fmrireg:::.rrr_filter_task_contrasts(bad, list(), event_idx, policy = "drop")
  expect_equal(silent$dropped, "bad_base")

  expect_error(
    fmrireg:::.rrr_filter_task_contrasts(bad, list(), event_idx, policy = "error"),
    "supports contrasts on event/task"
  )

  # Weights without colind default to sequential indices
  no_ci <- list(plain = c(1, -1))
  remapped <- fmrireg:::.rrr_filter_task_contrasts(no_ci, list(), event_idx, policy = "drop")
  expect_equal(names(remapped$simple), "plain")

  # Restore original conmat/colind on packaged contrast results
  packed <- list(
    dplyr::tibble(
      type = "contrast",
      name = "A_vs_B",
      conmat = list(structure(c(1, -1), colind = 1:2)),
      colind = list(1:2),
      data = list(dplyr::tibble(estimate = 1, se = 0.1, stat = 10, prob = 0.01))
    )
  )
  orig <- list(A_vs_B = structure(c(1, -1), colind = c(3L, 4L)))
  restored <- fmrireg:::.rrr_restore_contrast_metadata(packed, orig, list())
  expect_equal(attr(restored[[1]]$conmat[[1]], "colind"), c(3L, 4L))
  expect_equal(restored[[1]]$colind[[1]], c(3L, 4L))
  expect_equal(fmrireg:::.rrr_restore_contrast_metadata(list(), orig, list()), list())
})

test_that(".preflight_rrr_gls_engine and extract/whiten helpers", {
  fx <- tiny_rrr(n = 40L, V = 4L, seed = 53L, with_contrast = FALSE)
  expect_true(fmrireg:::.preflight_rrr_gls_engine(
    fx$model, fx$dataset, list(rank = 1L), fmri_lm_control()
  ))

  Y <- fmrireg:::.rrr_extract_response_matrix(fx$dataset)
  expect_equal(dim(Y), c(40L, 4L))
  expect_error(
    fmrireg:::.rrr_extract_response_matrix(list()),
    "matrix|dataset|data"
  )

  runs <- fmrireg:::.rrr_derive_run_indices(fx$model, 40L)
  expect_true(is.list(runs) || is.null(runs) || is.integer(unlist(runs)))

  expect_equal(fmrireg:::.rrr_ar_order_from_cfg(fmri_lm_control()), 0L)
  cfg_ar <- fmri_lm_control(ar_options = list(struct = "ar2"))
  expect_equal(fmrireg:::.rrr_ar_order_from_cfg(cfg_ar), 2L)

  # Residualize with and without nuisance
  X0 <- cbind(1, rnorm(30))
  Y0 <- matrix(rnorm(30 * 3), 30, 3)
  Zw <- cbind(1, rnorm(30))
  out <- fmrireg:::.rrr_residualize_against_nuisance(X0, Y0, Zw)
  expect_equal(dim(out$X0), dim(X0))
  expect_equal(dim(out$Y0), dim(Y0))
  out2 <- fmrireg:::.rrr_residualize_against_nuisance(X0, Y0, NULL)
  expect_equal(out2$rank_z, 0L)
})

test_that(".rrr_build_beta_stats and bootstrap_task_se edge seeds", {
  B <- matrix(rnorm(6), 3, 2)
  se <- matrix(abs(rnorm(4)), 2, 2)
  stats <- fmrireg:::.rrr_build_beta_stats(
    B_full = B,
    event_indices = 1:2,
    se_event = se,
    dfres = 20,
    df_inference = 18,
    varnames = paste0("b", 1:3),
    sigma = c(0.5, 0.6)
  )
  expect_equal(stats$type, "beta")
  expect_equal(ncol(stats$data[[1]]$estimate[[1]]), 3L)

  set.seed(3)
  X0 <- cbind(1, rnorm(24), rnorm(24))
  Y0 <- X0 %*% matrix(c(0.2, 0.5, -0.2), 3, 1) %*% matrix(1, 1, 3) +
    matrix(rnorm(24 * 3, sd = 0.3), 24, 3)
  Bhat <- solve(crossprod(X0), crossprod(X0, Y0))
  rc <- list(rank_mode = "fixed", rank = 1L)
  boot <- fmrireg:::.rrr_bootstrap_task_se(
    X0, Y0, Bhat, rc,
    t_contrasts_remapped = list(),
    nboot = 5L, block_size = 1L, seed = NULL
  )
  expect_equal(dim(boot$se_beta), dim(Bhat))
})
