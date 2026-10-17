# fmrilm.R early helpers: extract engine context.

test_that(".fmri_lm_extract_engine_context peels engine args and rejects extras", {
  ctx <- fmrireg:::.fmri_lm_extract_engine_context(list(
    engine = "rrr_gls",
    lowrank = list(m = 10),
    cfg = fmri_lm_control(),
    engine_args = list(rank = 2L),
    rrr_gls = list(se_mode = "bootstrap")
  ))
  expect_equal(ctx$engine, "rrr_gls")
  expect_equal(ctx$lowrank$m, 10)
  expect_equal(ctx$engine_args$rank, 2L)
  expect_equal(ctx$engine_args$se_mode, "bootstrap")
  expect_true(inherits(ctx$engine_cfg, "fmri_lm_control"))

  # Non-list engine_args coerced
  ctx2 <- fmrireg:::.fmri_lm_extract_engine_context(list(
    engine = "latent_sketch",
    engine_args = 5L,
    latent_sketch = 7L
  ))
  expect_true(is.list(ctx2$engine_args))

  expect_error(
    fmrireg:::.fmri_lm_extract_engine_context(list(engine = "ols", bogus = 1)),
    "Unexpected arguments"
  )
})

test_that(".fmri_lm_normalize_engine_name and resolve_engine_spec aliases", {
  expect_equal(fmrireg:::.fmri_lm_normalize_engine_name("sketch"), "latent_sketch")
  expect_equal(fmrireg:::.fmri_lm_normalize_engine_name("latent_sketch"), "latent_sketch")
  expect_equal(fmrireg:::.fmri_lm_normalize_engine_name("rrr_gls"), "rrr_gls")

  spec <- fmrireg:::.fmri_lm_resolve_engine_spec("sketch")
  expect_true(is.list(spec) || inherits(spec, "fmrireg_engine_spec") || !is.null(spec))

  # Warn helper should warn when ignored args are supplied
  expect_warning(
    fmrireg:::.fmri_lm_warn_engine_ignores(
      call("fmri_lm", engine = "rrr_gls", strategy = "runwise", nchunks = 2),
      engine = "rrr_gls"
    ),
    "ignores"
  )
  expect_silent(
    fmrireg:::.fmri_lm_warn_engine_ignores(
      call("fmri_lm", engine = "rrr_gls"),
      engine = "rrr_gls"
    )
  )
})

test_that("create_fmri_model builds from formula and validates inputs", {
  etab <- data.frame(
    onset = c(5, 20, 35),
    condition = factor(c("A", "B", "A")),
    run = 1L
  )
  Y <- matrix(rnorm(40 * 3), 40, 3)
  dset <- matrix_frame(Y, TR = 1, run_length = 40, event_table = etab)
  fm <- fmrireg:::create_fmri_model(
    onset ~ hrf(condition), block = ~ run, dataset = dset
  )
  expect_s3_class(fm, "fmri_model")
  expect_error(
    fmrireg:::create_fmri_model("not a formula", block = ~ run, dataset = dset),
    "formula"
  )
  expect_error(
    fmrireg:::create_fmri_model(onset ~ hrf(condition), block = "run", dataset = dset),
    "formula"
  )
})
