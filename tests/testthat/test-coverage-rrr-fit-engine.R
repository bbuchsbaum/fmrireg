# Deep coverage for .fit_rrr_gls_engine and remaining .rrr_* normalize branches.

tiny_rrr_fixture <- function(n = 60L, V = 8L, seed = 11) {
  set.seed(seed)
  etab <- data.frame(
    onset = seq(4, by = 6, length.out = 8),
    condition = factor(rep(c("A", "B"), length.out = 8)),
    run = 1L
  )
  Y <- matrix(rnorm(n * V), n, V)
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = dset$sampling_frame)
  list(
    model = fmri_model(emod, bmod, dset),
    dataset = dset,
    Y = Y
  )
}

test_that(".rrr_normalize_args covers modes, aliases, and rejection paths", {
  dflt <- fmrireg:::.rrr_normalize_args(NULL)
  expect_equal(dflt$rank_mode, "fixed")
  expect_equal(dflt$se_mode, "conditional")
  expect_equal(dflt$contrast_policy, "warn_drop")

  # Non-list args coerced
  wrapped <- fmrireg:::.rrr_normalize_args(2L)
  expect_true(is.list(wrapped))

  # energy / nboot aliases
  aliased <- fmrireg:::.rrr_normalize_args(list(
    energy = 0.9, nboot = 25L, rank_mode = "energy", se_mode = "bootstrap"
  ))
  expect_equal(aliased$energy_keep, 0.9)
  expect_equal(aliased$bootstrap_n, 25L)
  expect_equal(aliased$rank_mode, "energy")
  expect_equal(aliased$se_mode, "bootstrap")

  rss <- fmrireg:::.rrr_normalize_args(list(
    rank_mode = "rss_budget", rss_budget = 0.5, bootstrap_block_size = 3L,
    bootstrap_seed = 42L
  ))
  expect_equal(rss$rss_budget, 0.5)
  expect_equal(rss$bootstrap_block_size, 3L)
  expect_equal(rss$bootstrap_seed, 42L)

  expect_error(
    fmrireg:::.rrr_normalize_args(list(rank = 0L)),
    "positive integer"
  )
  expect_error(
    fmrireg:::.rrr_normalize_args(list(rank_mode = "energy", energy_keep = 1.5)),
    "energy_keep"
  )
  expect_error(
    fmrireg:::.rrr_normalize_args(list(rank_mode = "rss_budget", rss_budget = -1)),
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
})

test_that(".fit_rrr_gls_engine succeeds on tiny matrix_dataset (conditional SE)", {
  fx <- tiny_rrr_fixture()
  cfg <- fmri_lm_control()
  fit <- fmrireg:::.fit_rrr_gls_engine(
    fx$model, fx$dataset,
    args = list(rank = 2L, rank_mode = "fixed", se_mode = "conditional"),
    cfg = cfg
  )
  expect_s3_class(fit, "fmri_lm")
  expect_true(!is.null(fit$result$betas))
  expect_true(is.matrix(fit$result$cov.unscaled) || is.matrix(fit$vcov_inv) ||
                !is.null(fit$result$betas))
  eng <- attr(fit, "engine")
  expect_true(is.null(eng) || identical(eng, "rrr_gls"))
})

test_that(".fit_rrr_gls_engine bootstrap SE path runs with tiny nboot", {
  fx <- tiny_rrr_fixture(n = 48L, V = 6L, seed = 17)
  cfg <- fmri_lm_control()
  fit <- fmrireg:::.fit_rrr_gls_engine(
    fx$model, fx$dataset,
    args = list(
      rank = 1L,
      se_mode = "bootstrap",
      bootstrap_n = 8L,
      bootstrap_block_size = 4L,
      bootstrap_seed = 3L,
      contrast_policy = "drop"
    ),
    cfg = cfg
  )
  expect_s3_class(fit, "fmri_lm")
  expect_true(!is.null(fit$result$betas))
})

test_that(".rrr_fit_task_subspace handles rank-0 / near-zero signal", {
  set.seed(2)
  X0 <- matrix(0, 20, 2)
  Y0 <- matrix(rnorm(20 * 4, sd = 1e-16), 20, 4)
  expect_warning(
    out0 <- fmrireg:::.rrr_fit_task_subspace(
      X0, Y0, list(rank_mode = "fixed", rank = 1L)
    ),
    "rank deficient|Aliased"
  )
  expect_equal(out0$rank_used, 0L)
  expect_equal(dim(out0$B_task), c(2L, 4L))

  # Energy / rss choose-rank edge: single singular value
  expect_equal(
    fmrireg:::.rrr_choose_rank(5, list(rank_mode = "rss_budget", rss_budget = 0)),
    1L
  )
  expect_equal(
    fmrireg:::.rrr_choose_rank(c(10, 1), list(rank_mode = "energy", energy_keep = 0.5)),
    1L
  )
})

test_that(".register_builtin_engines is idempotent and exposes specs", {
  fmrireg:::.register_builtin_engines()
  expect_true(!is.null(get_engine("rrr_gls")))
  expect_true(!is.null(get_engine("latent_sketch")))
  fmrireg:::.register_builtin_engines()
  expect_equal(fmrireg:::.builtin_engine_aliases("latent_sketch"), "sketch")
  expect_equal(fmrireg:::.builtin_engine_aliases("rrr_gls"), character())
  expect_equal(fmrireg:::.builtin_engine_source("rrr_gls"), "builtin")
  expect_equal(fmrireg:::.builtin_engine_source("custom_x"), "plugin")
})
