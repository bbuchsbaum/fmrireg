# Coverage for rrr/lowrank engine helpers and fmriAR adapter normalization.

tiny_design_fixture <- function(n_time = 60L, n_vox = 8L, n_events = 6L) {
  set.seed(17)
  onsets <- seq(5, n_time - 8, length.out = n_events)
  etab <- data.frame(
    onset = onsets,
    condition = factor(rep(c("A", "B"), length.out = n_events)),
    run = rep(1:2, length.out = n_events)
  )
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  dset <- matrix_dataset(
    Y, TR = 1,
    run_length = c(n_time / 2, n_time / 2),
    event_table = etab
  )
  emod <- event_model(
    onset ~ hrf(condition),
    data = etab,
    block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = dset$sampling_frame)
  list(
    model = fmri_model(emod, bmod, dset),
    dataset = dset,
    Y = Y
  )
}

test_that(".rrr_normalize_args validates rank modes and aliases", {
  dflt <- fmrireg:::.rrr_normalize_args(NULL)
  expect_equal(dflt$rank_mode, "fixed")
  expect_equal(dflt$se_mode, "conditional")
  expect_equal(dflt$contrast_policy, "warn_drop")

  aliased <- fmrireg:::.rrr_normalize_args(list(energy = 0.9, nboot = 25L, rank_mode = "energy"))
  expect_equal(aliased$energy_keep, 0.9)
  expect_equal(aliased$bootstrap_n, 25L)

  expect_error(fmrireg:::.rrr_normalize_args(list(rank = 0L)), "positive integer")
  expect_error(
    fmrireg:::.rrr_normalize_args(list(rank_mode = "energy", energy_keep = 1.5)),
    "energy_keep"
  )
  expect_error(
    fmrireg:::.rrr_normalize_args(list(rank_mode = "rss_budget")),
    "rss_budget"
  )
  expect_error(
    fmrireg:::.rrr_normalize_args(list(se_mode = "bootstrap", bootstrap_n = 1L)),
    "bootstrap_n"
  )
  expect_error(
    fmrireg:::.rrr_normalize_args(list(bootstrap_block_size = 0L)),
    "bootstrap_block_size"
  )

  ok <- fmrireg:::.rrr_normalize_args(list(
    rank = 2L, se_mode = "bootstrap", bootstrap_n = 10L,
    bootstrap_block_size = 5L, bootstrap_seed = 1L,
    rank_mode = "rss_budget", rss_budget = 0.5
  ))
  expect_equal(ok$rank, 2L)
  expect_equal(ok$rss_budget, 0.5)
})

test_that(".rrr_choose_rank / residualize / contrast filter helpers", {
  d <- c(5, 3, 1, 0.2)
  expect_equal(
    fmrireg:::.rrr_choose_rank(d, list(rank_mode = "fixed", rank = 2L)),
    2L
  )
  expect_equal(
    fmrireg:::.rrr_choose_rank(d, list(rank_mode = "fixed", rank = NULL)),
    4L
  )
  r_energy <- fmrireg:::.rrr_choose_rank(
    d, list(rank_mode = "energy", energy_keep = 0.95)
  )
  expect_true(r_energy >= 1L && r_energy <= 4L)
  r_rss <- fmrireg:::.rrr_choose_rank(
    d, list(rank_mode = "rss_budget", rss_budget = 1e6)
  )
  expect_true(r_rss >= 1L)
  expect_equal(fmrireg:::.rrr_choose_rank(c(0, NA), list(rank_mode = "fixed")), 0L)

  X <- cbind(1, rnorm(40), rnorm(40))
  Y <- matrix(rnorm(40 * 5), 40, 5)
  Z <- cbind(rnorm(40), rnorm(40))
  res0 <- fmrireg:::.rrr_residualize_against_nuisance(X[, 1:2], Y, NULL)
  expect_equal(res0$rank_z, 0L)
  expect_equal(dim(res0$X0), c(40L, 2L))
  resZ <- fmrireg:::.rrr_residualize_against_nuisance(X[, 1:2], Y, Z)
  expect_true(resZ$rank_z >= 1L)
  expect_equal(dim(resZ$Y0), dim(Y))

  w_ok <- structure(c(1, -1), colind = c(1L, 2L))
  w_bad <- structure(1, colind = 5L)
  filtered <- fmrireg:::.rrr_filter_task_contrasts(
    list(good = w_ok, bad = w_bad),
    list(),
    event_indices = 1:3,
    policy = "drop"
  )
  expect_true("good" %in% names(filtered$simple))
  expect_true("bad" %in% filtered$dropped)

  expect_warning(
    fmrireg:::.rrr_filter_task_contrasts(
      list(bad = w_bad), list(), event_indices = 1:2, policy = "warn_drop"
    ),
    "dropped contrasts"
  )
  expect_error(
    fmrireg:::.rrr_filter_task_contrasts(
      list(bad = w_bad), list(), event_indices = 1:2, policy = "error"
    ),
    "event/task"
  )

  # Restore metadata
  fake_cres <- list(tibble::tibble(
    name = "good",
    conmat = list(c(1, -1)),
    colind = list(1:2)
  ))
  restored <- fmrireg:::.rrr_restore_contrast_metadata(
    fake_cres,
    original_simple = list(good = w_ok),
    original_f = list()
  )
  expect_equal(restored[[1]]$colind[[1]], c(1L, 2L))
})

test_that(".rrr extract/ar-order/whiten/build-beta-stats cover", {
  fx <- tiny_design_fixture()
  Y <- fmrireg:::.rrr_extract_response_matrix(fx$dataset)
  expect_equal(dim(Y), dim(fx$Y))

  runs <- fmrireg:::.rrr_derive_run_indices(fx$model, nrow(Y))
  expect_equal(length(runs), 2L)
  expect_equal(sum(lengths(runs)), nrow(Y))

  cfg <- fmri_lm_control()
  expect_equal(fmrireg:::.rrr_ar_order_from_cfg(cfg), 0L)
  cfg_ar <- fmri_lm_control(ar_options = list(struct = "ar1"))
  expect_equal(fmrireg:::.rrr_ar_order_from_cfg(cfg_ar), 1L)
  cfg_arp <- list(ar = list(struct = "arp", p = 3L))
  expect_equal(fmrireg:::.rrr_ar_order_from_cfg(cfg_arp), 3L)

  X_task <- as.matrix(design_matrix(fx$model))[, 1:2, drop = FALSE]
  Z_nuis <- as.matrix(design_matrix(fx$model))[, -(1:2), drop = FALSE]
  whitened <- fmrireg:::.rrr_whiten_data(X_task, Z_nuis, Y, cfg, run_indices = runs)
  expect_equal(whitened$ar_order, 0L)
  expect_equal(dim(whitened$Yw), dim(Y))

  if (requireNamespace("fmriAR", quietly = TRUE)) {
    whitened_ar <- fmrireg:::.rrr_whiten_data(
      X_task, Z_nuis, Y, cfg_ar, run_indices = runs
    )
    expect_equal(whitened_ar$ar_order, 1L)
    expect_equal(dim(whitened_ar$Xw_task), dim(X_task))
  }

  # Fit subspace + beta stats
  X0 <- matrix(rnorm(50 * 3), 50, 3)
  Y0 <- matrix(rnorm(50 * 6), 50, 6)
  subspace <- fmrireg:::.rrr_fit_task_subspace(
    X0, Y0, list(rank_mode = "fixed", rank = 1L)
  )
  expect_equal(nrow(subspace$B_task), 3L)
  expect_equal(ncol(subspace$B_task), 6L)
  expect_equal(subspace$rank_used, 1L)

  B_full <- rbind(subspace$B_task, matrix(0, 2, 6))
  se_event <- matrix(0.2, 3, 6)
  bstats <- fmrireg:::.rrr_build_beta_stats(
    B_full, event_indices = 1:3, se_event = se_event,
    dfres = 40, varnames = paste0("b", 1:5), sigma = rep(1, 6)
  )
  expect_s3_class(bstats, "tbl_df")
  expect_equal(bstats$type, "beta")
})

test_that("lowrank group AR helper and whiten initial OLS", {
  set.seed(9)
  n <- 40L
  X <- cbind(1, rnorm(n), rnorm(n))
  Z <- matrix(rnorm(n * 4), n, 4)
  ar_opts <- list(struct = "ar1", exact_first = TRUE)

  if (requireNamespace("fmriAR", quietly = TRUE)) {
    whitened <- fmrireg:::.lowrank_whiten_initial_ols(
      X, Z, ar_order = 1L, ar_opts = ar_opts
    )
    expect_true(!is.null(whitened$phi))
    expect_equal(dim(whitened$X), dim(X))
    expect_equal(dim(whitened$Y), dim(Z))

    resid_list <- list(
      g1 = as.numeric(whitened$residuals[, 1]),
      g2 = as.numeric(whitened$residuals[, 2])
    )
    grouped <- fmrireg:::.lowrank_group_ar_estimates(
      residuals = resid_list,
      sizes = list(g1 = 10, g2 = 5),
      ar_order = 1L,
      ar_opts = ar_opts,
      shrink_c0 = 5,
      design = X
    )
    expect_equal(length(grouped$phi_groups), 2L)
    expect_true(!is.null(grouped$global_phi))
  }

  expect_error(
    fmrireg:::.lowrank_group_ar_estimates(
      residuals = list(), sizes = list(), ar_order = 1L,
      ar_opts = list(struct = "ar1"), shrink_c0 = 1, design = X
    ),
    "non-empty list"
  )
})

test_that("fmriAR adapter config normalization and run helpers", {
  expect_equal(
    fmrireg:::.normalize_ar_options(NULL)$struct,
    "iid"
  )
  from_cor <- fmrireg:::.normalize_ar_options(list(cor_struct = "ar2"))
  expect_equal(from_cor$struct, "ar2")
  expect_equal(from_cor$cor_struct, "ar2")

  none <- fmrireg:::.normalize_ar_options(list(struct = "none"))
  expect_equal(none$struct, "iid")

  ar5 <- fmrireg:::.normalize_ar_options(list(struct = "ar5"))
  expect_equal(ar5$struct, "arp")
  expect_equal(ar5$p, 5L)

  ar0 <- fmrireg:::.normalize_ar_options(list(struct = "ar0"))
  expect_equal(ar0$struct, "iid")

  expect_error(
    fmrireg:::.normalize_ar_options(list(struct = "arp")),
    "p must be specified"
  )
  arp_from_cor <- fmrireg:::.normalize_ar_options(
    list(struct = "arp", cor_struct = "ar3")
  )
  expect_equal(arp_from_cor$p, 3L)

  iter_map <- fmrireg:::.normalize_ar_options(list(struct = "ar1", iter = 3L))
  expect_equal(iter_map$iter_gls, 3L)

  expect_error(
    fmrireg:::.normalize_ar_options(list(struct = "ar1", shared_estimator = "bad")),
    "shared_estimator"
  )

  expect_false(fmrireg:::.temporal_noise_enabled(list(struct = "iid")))
  expect_true(fmrireg:::.temporal_noise_enabled(list(struct = "ar1")))
  expect_true(fmrireg:::.temporal_noise_enabled(list(struct = "iid", q = 1L)))

  expect_equal(fmrireg:::.target_ar_order(list(struct = "ar2")), 2L)
  expect_equal(fmrireg:::.target_ar_order(list(struct = "arp", p = 4L)), 4L)
  expect_equal(fmrireg:::.target_ar_order(list(struct = "iid")), 0L)

  expect_null(fmrireg:::.build_run_labels(10L, NULL))
  labs <- fmrireg:::.build_run_labels(6L, list(1:3, 4:6))
  expect_equal(labs, c(1L, 1L, 1L, 2L, 2L, 2L))

  expect_equal(
    fmrireg:::.split_run_slices(6L, "global", list(1:3, 4:6)),
    list(1:6)
  )
  expect_equal(
    fmrireg:::.split_run_slices(6L, "run", list(1:3, 4:6)),
    list(1:3, 4:6)
  )

  block <- fmrireg:::.block_diagonal_design(list(matrix(1:4, 2), matrix(5:8, 2)))
  expect_equal(dim(block), c(4L, 4L))
  expect_equal(block[1:2, 1:2], matrix(1:4, 2))
  expect_equal(block[3:4, 3:4], matrix(5:8, 2))

  cfg <- fmrireg:::.fmrireg_to_fmriAR_config(list(struct = "ar1", global = TRUE))
  expect_equal(cfg$method, "ar")
  expect_equal(cfg$p, 1L)
  expect_equal(cfg$pooling, "global")

  cfg_arma <- fmrireg:::.fmrireg_to_fmriAR_config(
    list(struct = "ar1", q = 1L), n_runs = 2L
  )
  expect_equal(cfg_arma$method, "arma")
  expect_equal(cfg_arma$q, 1L)

  cfg_run <- fmrireg:::.fmrireg_to_fmriAR_config(
    list(struct = "ar2", global = FALSE), n_runs = 3L
  )
  expect_equal(cfg_run$pooling, "run")
  expect_equal(cfg_run$p, 2L)

  expect_error(
    fmrireg:::.fmrireg_to_fmriAR_config(list(struct = "ar1", q = -1L)),
    "nonnegative"
  )

  expect_message(
    fmrireg:::.fmrireg_to_fmriAR_config(list(struct = "ar1", voxelwise = TRUE)),
    "voxelwise"
  )

  # mean_series residual input collapses columns
  resid <- matrix(rnorm(20 * 4), 20, 4)
  collapsed <- fmrireg:::.shared_ar_residual_input(
    resid, list(struct = "ar1", shared_estimator = "mean_series")
  )
  expect_equal(ncol(collapsed), 1L)
  expect_equal(as.numeric(collapsed), rowMeans(resid), tolerance = 1e-12)

  pooled <- fmrireg:::.shared_ar_residual_input(
    resid, list(struct = "ar1", shared_estimator = "pooled_acvf")
  )
  expect_equal(dim(pooled), dim(resid))
})
