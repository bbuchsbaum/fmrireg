# chunkwise/runwise strategy internals: slow path, progress, voxelwise AR edges.

make_two_run_fx <- function(n_per = 40L, V = 4L, seed = 17L) {
  set.seed(seed)
  etab <- data.frame(
    onset = c(5, 15, 25, 35, 5, 15, 25, 35),
    condition = factor(rep(c("A", "B"), 4)),
    run = rep(1:2, each = 4)
  )
  Y <- matrix(rnorm((2 * n_per) * V), 2 * n_per, V)
  dset <- matrix_dataset(Y, TR = 1, run_length = c(n_per, n_per), event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = dset$sampling_frame)
  list(model = fmri_model(emod, bmod, dset), dataset = dset)
}

test_that("chunkwise_lm.fmri_dataset fast and slow paths with progress", {
  fx <- make_two_run_fx()
  cfg <- fmri_lm_control()
  cons <- list()

  fast <- fmrireg:::chunkwise_lm.fmri_dataset(
    fx$dataset, fx$model, cons, nchunks = 2L, cfg = cfg,
    use_fast_path = TRUE, progress = FALSE, verbose = FALSE
  )
  expect_true(is.list(fast))
  expect_true(!is.null(fast$betas) || !is.null(fast$event_indices) || !is.null(fast$rss))

  slow <- fmrireg:::chunkwise_lm.fmri_dataset(
    fx$dataset, fx$model, cons, nchunks = 2L, cfg = cfg,
    use_fast_path = FALSE, progress = FALSE, verbose = TRUE
  )
  expect_true(is.list(slow))

  expect_error(
    fmrireg:::chunkwise_lm.fmri_dataset(
      fx$dataset, fx$model, cons, nchunks = 2L, cfg = list(),
      use_fast_path = TRUE
    ),
    "fmri_lm_control"
  )
  expect_error(
    fmrireg:::chunkwise_lm.fmri_dataset(
      fx$dataset, fx$model, cons, nchunks = 2L, cfg = cfg,
      parallel_chunks = NA
    ),
    "parallel_chunks"
  )
})

test_that("runwise_lm_impl fast/slow and progress bar branches", {
  fx <- make_two_run_fx(n_per = 36L, V = 3L, seed = 18L)
  cfg <- fmri_lm_control()
  cons <- list()

  fast <- fmrireg:::runwise_lm_impl(
    fx$dataset, fx$model, cons, cfg = cfg,
    use_fast_path = TRUE, progress = FALSE, verbose = TRUE
  )
  expect_true(is.list(fast))
  expect_true(!is.null(fast$betas) || !is.null(fast$event_indices))

  slow <- fmrireg:::runwise_lm_impl(
    fx$dataset, fx$model, cons, cfg = cfg,
    use_fast_path = FALSE, progress = FALSE, verbose = FALSE
  )
  expect_true(is.list(slow))

  # AR runwise (global) if fmriAR available
  skip_if_not_installed("fmriAR")
  cfg_ar <- fmri_lm_control(ar_options = list(struct = "ar1"))
  ar_fit <- fmrireg:::runwise_lm_impl(
    fx$dataset, fx$model, cons, cfg = cfg_ar,
    use_fast_path = TRUE, progress = FALSE
  )
  expect_true(is.list(ar_fit))
})

test_that("prepare_fmri_lm_contrasts processes model contrasts and missing metadata", {
  fx <- make_two_run_fx(seed = 19L)
  # No contrasts on demo model -> empty processed list
  out <- fmrireg:::prepare_fmri_lm_contrasts(fx$model)
  expect_true(is.list(out))

  # Inject contrast_weights on event model and strip col_indices to hit error
  emod <- fx$model$event_model
  dm <- emod$design_matrix
  attr(dm, "col_indices") <- NULL
  emod$design_matrix <- dm
  # Fake contrast_weights method via attribute if present
  cw <- tryCatch(contrast_weights(fx$model$event_model), error = function(e) list())
  if (length(cw) == 0L) {
    # Manually attach a contrast_weights-like list through a stub model
    fmod <- fx$model
    fmod$event_model$contrast_weights <- list(
      A_vs_B = list(name = "A_vs_B", term = "condition", weights = c(1, -1))
    )
    # prepare uses contrast_weights() S3; skip if we can't inject
    expect_true(is.list(fmrireg:::prepare_fmri_lm_contrasts(fx$model)))
  } else {
    fmod <- fx$model
    fmod$event_model$design_matrix <- dm
    expect_error(fmrireg:::prepare_fmri_lm_contrasts(fmod), "column-index metadata")
  }
})

test_that("coef_names.fmri_lm and tidy-ish methods cover remaining branches", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  cn <- coef_names(fit)
  expect_true(is.character(cn))
  expect_true(length(cn) >= 1L)

  # include_baseline path via coef
  cf_all <- coef(fit, type = "betas", include_baseline = TRUE)
  expect_true(is.matrix(cf_all) || inherits(cf_all, "tbl_df") || is.data.frame(cf_all))

  st <- stats(fit)
  expect_true(!is.null(st))
  se <- standard_error(fit)
  expect_true(!is.null(se))
  pv <- p_values(fit)
  expect_true(!is.null(pv))
})
