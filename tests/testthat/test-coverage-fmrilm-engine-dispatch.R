# Twelfth wave: fmrilm engine-ignore warn + dispatch + attach config metadata.

test_that(".fmri_lm_warn_engine_ignores flags strategy/nchunks/progress args", {
  call_obj <- quote(fmri_lm(
    onset ~ hrf(condition),
    block = ~ run,
    dataset = dset,
    engine = "rrr_gls",
    strategy = "runwise",
    nchunks = 4L,
    progress = TRUE,
    use_fast_path = FALSE
  ))
  expect_warning(
    ignored <- fmrireg:::.fmri_lm_warn_engine_ignores(call_obj, "rrr_gls"),
    "ignores"
  )
  expect_true(all(c("strategy", "nchunks", "progress", "use_fast_path") %in% ignored))

  call_clean <- quote(fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset, engine = "rrr_gls"))
  expect_equal(
    fmrireg:::.fmri_lm_warn_engine_ignores(call_clean, "rrr_gls"),
    character(0)
  )
})

test_that(".fmri_lm_dispatch_engine runs rrr_gls on matrix dataset", {
  set.seed(291)
  n <- 56L
  V <- 6L
  etab <- data.frame(
    onset = c(5, 15, 25, 35, 45),
    condition = factor(c("A", "B", "A", "B", "A")),
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

  fit <- fmrireg:::.fmri_lm_dispatch_engine(
    model = model,
    dataset = dset,
    engine = "rrr_gls",
    cfg = cfg,
    engine_args = list(rank_mode = "fixed", rank = 2L)
  )
  expect_s3_class(fit, "fmri_lm")
  expect_true(!is.null(fit$result) || !is.null(coef(fit)))
})

test_that(".fmri_lm_attach_config_metadata stores cfg on fit objects", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1"))
  attached <- fmrireg:::.fmri_lm_attach_config_metadata(fit, cfg)
  expect_s3_class(attr(attached, "requested_config"), "fmri_lm_control")
  expect_s3_class(attr(attached, "executed_config"), "fmri_lm_control")
  expect_identical(attr(attached, "config"), cfg)
  expect_equal(attr(attached, "config")$ar$struct, "ar1")
})
