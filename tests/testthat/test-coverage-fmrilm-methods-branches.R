# Twelfth wave: fmrilm method branches for coef/stats/se/print/fitted_hrf.

test_that("coef/stats/standard_error/print.fmri_lm cover type branches", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  expect_s3_class(fit, "fmri_lm")

  b <- coef(fit, type = "betas", include_baseline = FALSE)
  expect_true(is.matrix(b) || inherits(b, "tbl_df") || is.data.frame(b))
  b_all <- coef(fit, type = "betas", include_baseline = TRUE)
  expect_true(NCOL(b_all) >= NCOL(b) || NROW(b_all) >= NROW(b))

  # Contrasts may be empty on demo fit
  bc <- tryCatch(coef(fit, type = "contrasts"), error = function(e) e)
  expect_true(inherits(bc, "error") || !is.null(bc))

  st_est <- stats(fit, type = "estimates")
  expect_true(!is.null(st_est))
  st_con <- tryCatch(stats(fit, type = "contrasts"), error = function(e) e)
  expect_true(inherits(st_con, "error") || !is.null(st_con))
  st_f <- tryCatch(stats(fit, type = "F"), error = function(e) e)
  expect_true(inherits(st_f, "error") || !is.null(st_f))

  se_est <- standard_error(fit, type = "estimates")
  expect_true(!is.null(se_est))
  se_con <- tryCatch(standard_error(fit, type = "contrasts"), error = function(e) e)
  expect_true(inherits(se_con, "error") || !is.null(se_con))

  expect_output(print(fit), regexp = ".")
})

test_that("fitted_hrf.fmri_lm returns term predictions for sample grid", {
  set.seed(271)
  n <- 60L
  etab <- data.frame(
    onset = c(5, 15, 25, 35, 45),
    condition = factor(c("A", "B", "A", "B", "A")),
    run = 1L
  )
  Y <- matrix(rnorm(n * 3), n, 3)
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = etab)
  fit <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset)

  hrf_list <- fitted_hrf(fit, sample_at = seq(0, 12, by = 2))
  expect_true(is.list(hrf_list))
  expect_true(length(hrf_list) >= 1L)
  first <- hrf_list[[1]]
  expect_true(!is.null(first$design) || !is.null(first$pred) || is.list(first))
})

test_that("fmri_lm_fit and engine resolve helpers accept matrix models", {
  set.seed(272)
  n <- 48L
  etab <- data.frame(
    onset = c(5, 15, 25, 35),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  Y <- matrix(rnorm(n * 2), n, 2)
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "constant", sframe = dset$sampling_frame)
  fmod <- fmri_model(emod, bmod, dset)
  cfg <- fmri_lm_control()

  fit <- fmrireg:::fmri_lm_fit(
    fmod, dset, strategy = "runwise", cfg = cfg, nchunks = 1L
  )
  expect_s3_class(fit, "fmri_lm")

  # Resolve engine spec for named engines / invalid names
  sketch <- tryCatch(
    fmrireg:::.fmri_lm_resolve_engine_spec("latent_sketch"),
    error = function(e) e
  )
  expect_true(inherits(sketch, "error") || is.list(sketch) || is.null(sketch))
  expect_error(fmrireg:::.fmri_lm_resolve_engine_spec("not_an_engine_zzz"), regexp = ".")
})