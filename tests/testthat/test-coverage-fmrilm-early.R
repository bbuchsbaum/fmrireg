# fmrilm.R early helpers: extract engine context, fast_rlm_run branches.

test_that(".fmri_lm_extract_engine_context peels engine args and rejects extras", {
  ctx <- fmrireg:::.fmri_lm_extract_engine_context(list(
    engine = "rrr_gls",
    lowrank = list(m = 10),
    ar_options = list(struct = "ar1"),
    robust_options = list(type = "huber"),
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

test_that("fast_rlm_run covers huber/bisquare and sigma_fixed paths", {
  set.seed(12)
  n <- 40L
  p <- 3L
  V <- 4L
  X <- cbind(1, rnorm(n), rnorm(n))
  Y <- X %*% matrix(c(0.5, 1, -0.5), p, V) + matrix(rnorm(n * V, sd = 0.4), n, V)
  # Inject an outlier timepoint
  Y[1, ] <- Y[1, ] + 8

  proj <- fmrireg:::.fast_preproject(X)
  fit_h <- fmrireg:::fast_rlm_run(X, Y, proj, psi = "huber", max_it = 3L)
  expect_equal(dim(fit_h$betas), c(p, V))
  expect_true(length(fit_h$sigma) >= 1L)
  expect_equal(length(fit_h$weights), n)
  expect_true(all(is.finite(fit_h$betas)))
  expect_true(all(fit_h$sigma > 0))
  expect_equal(length(fit_h$se), p)

  fit_b <- fmrireg:::fast_rlm_run(X, Y, proj = NULL, psi = "bisquare", max_it = 2L)
  expect_equal(dim(fit_b$betas), c(p, V))
  expect_true(all(is.finite(fit_b$se)))

  # Global scale via sigma_fixed
  fit_g <- fmrireg:::fast_rlm_run(
    X, Y, proj, psi = "huber", max_it = 2L,
    sigma_fixed = as.numeric(fit_h$sigma)[1]
  )
  expect_equal(dim(fit_g$betas), c(p, V))
  expect_equal(length(fit_g$sigma), 1L)

  expect_error(
    fmrireg:::fast_rlm_run(X, Y + NA, proj),
    "NA/Inf"
  )
  Xbad <- X
  Xbad[1, 1] <- Inf
  expect_error(
    fmrireg:::fast_rlm_run(Xbad, Y, proj),
    "NA/Inf"
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
  dset <- matrix_dataset(Y, TR = 1, run_length = 40, event_table = etab)
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
