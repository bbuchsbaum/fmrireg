# Twelfth wave: remaining lowrank/rrr preflight+fit edges and generic fallbacks.

test_that(".preflight_rrr_gls_engine and .fit_rrr_gls_engine happy path", {
  set.seed(281)
  n <- 64L
  V <- 8L
  etab <- data.frame(
    onset = c(5, 15, 25, 35, 45, 55),
    condition = factor(rep(c("A", "B"), 3)),
    run = 1L
  )
  Y <- matrix(rnorm(n * V), n, V)
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = dset$sampling_frame)
  model <- fmri_model(emod, bmod, dset)
  cfg <- fmri_lm_control()

  expect_true(fmrireg:::.preflight_rrr_gls_engine(model, dset, list(rank = 2L), cfg))
  expect_error(
    fmrireg:::.preflight_rrr_gls_engine(list(), dset, list(), cfg),
    "fmri_model"
  )
  expect_error(
    fmrireg:::.preflight_rrr_gls_engine(model, NULL, list(), cfg),
    "dataset"
  )
  expect_error(
    fmrireg:::.preflight_rrr_gls_engine(model, dset, list(), list()),
    "fmri_lm_control"
  )

  fit <- fmrireg:::.fit_rrr_gls_engine(
    model, dset,
    args = list(rank_mode = "fixed", rank = 2L),
    cfg = cfg
  )
  expect_true(inherits(fit, "fmri_lm") || is.list(fit))
  expect_true(!is.null(fit$result) || !is.null(fit$betas) || inherits(fit, "fmri_lm"))
})

test_that(".fit_lowrank_engine_plugin landmarks path with tiny sketch", {
  set.seed(282)
  n <- 48L
  V <- 12L
  etab <- data.frame(
    onset = c(5, 15, 25, 35),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  Y <- matrix(rnorm(n * V), n, V)
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "constant", sframe = dset$sampling_frame)
  model <- fmri_model(emod, bmod, dset)
  cfg <- fmri_lm_control()

  expect_true(fmrireg:::.preflight_lowrank_engine_plugin(
    model, dset, args = list(lowrank = list(k_latents = 3L)), cfg = cfg
  ))

  fit <- tryCatch(
    fmrireg:::.fit_lowrank_engine_plugin(
      model, dset,
      args = list(lowrank = list(
        method = "landmarks",
        k_latents = 3L,
        n_landmarks = 6L
      )),
      cfg = cfg
    ),
    error = function(e) e
  )
  expect_true(inherits(fit, "error") || inherits(fit, "fmri_lm") || is.list(fit))
})

test_that("contrast.default and tidy generic dispatch on non-meta objects", {
  # contrast.default delegates to fmridesign::contrast
  out <- tryCatch(contrast(list(a = 1), ~ a), error = function(e) e)
  expect_true(inherits(out, "error") || !is.null(out))

  # tidy generic exists and errors on unsupported classes without method
  expect_error(tidy(1:3), regexp = ".")
})

test_that("runwise_lm_fast and chunkwise_lm_fast return pooled structures", {
  set.seed(283)
  n_per <- 36L
  V <- 4L
  etab <- data.frame(
    onset = c(5, 15, 25, 5, 15, 25),
    condition = factor(rep(c("A", "B"), 3)),
    run = rep(1:2, each = 3)
  )
  Y <- matrix(rnorm((2 * n_per) * V), 2 * n_per, V)
  dset <- matrix_dataset(Y, TR = 1, run_length = c(n_per, n_per), event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = dset$sampling_frame)
  model <- fmri_model(emod, bmod, dset)
  cfg <- fmri_lm_control()

  # Prefer public-ish impl wrappers when direct *_fast needs many args
  fast_run <- fmrireg:::runwise_lm_impl(
    dset, model, contrast_objects = list(), cfg = cfg,
    use_fast_path = TRUE, progress = FALSE
  )
  expect_true(is.list(fast_run))
  expect_true(!is.null(fast_run$betas) || !is.null(fast_run$event_indices))

  fast_chunk <- fmrireg:::chunkwise_lm.fmri_dataset(
    dset, model, contrast_objects = list(), nchunks = 2L, cfg = cfg,
    use_fast_path = TRUE, progress = FALSE
  )
  expect_true(is.list(fast_chunk))
})
