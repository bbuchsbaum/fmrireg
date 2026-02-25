test_that("rrr_gls engine fits and exposes diagnostics", {
  dset <- .demo_matrix_dataset()

  fit <- fmri_lm(
    onsets ~ hrf(condition),
    block = ~run,
    dataset = dset,
    engine = "rrr_gls",
    engine_args = list(rank = 1L)
  )

  expect_s3_class(fit, "fmri_lm")
  expect_equal(attr(fit, "engine"), "rrr_gls")
  expect_equal(attr(fit, "strategy"), "engine")
  expect_true(is.list(fit$rrr))
  expect_equal(fit$rrr$rank_used, 1L)

  se <- standard_error(fit, "estimates")
  expect_true(all(is.finite(as.matrix(se))))
})


test_that("rrr_gls supports bootstrap standard errors", {
  dset <- .demo_matrix_dataset()

  fit <- fmri_lm(
    onsets ~ hrf(condition),
    block = ~run,
    dataset = dset,
    engine = "rrr_gls",
    engine_args = list(
      rank = 1L,
      se_mode = "bootstrap",
      bootstrap_n = 10L,
      bootstrap_seed = 1L
    )
  )

  expect_s3_class(fit, "fmri_lm")
  expect_equal(fit$rrr$se_mode, "bootstrap")

  se <- standard_error(fit, "estimates")
  expect_true(all(is.finite(as.matrix(se))))
})


test_that("rrr_gls supports shared AR whitening", {
  set.seed(123)
  Tlen <- 80
  V <- 20

  ev <- data.frame(
    onsets = seq(4, 52, by = 8),
    condition = factor(rep(c("A", "B"), length.out = 7)),
    run = 1L
  )

  Y <- matrix(rnorm(Tlen * V), nrow = Tlen, ncol = V)
  dset <- fmridataset::matrix_dataset(
    Y,
    TR = 2,
    run_length = Tlen,
    event_table = ev
  )

  fit <- fmri_lm(
    onsets ~ hrf(condition),
    block = ~run,
    dataset = dset,
    ar_options = list(struct = "ar1"),
    engine = "rrr_gls",
    engine_args = list(rank = 2L)
  )

  expect_s3_class(fit, "fmri_lm")
  expect_equal(attr(fit, "engine"), "rrr_gls")
  expect_true(is.list(ar_parameters(fit, scope = "raw")))
})


test_that("rrr_gls enforces event-only contrast scope for post-hoc contrasts", {
  dset <- .demo_matrix_dataset()

  fit <- fmri_lm(
    onsets ~ hrf(condition),
    block = ~run,
    dataset = dset,
    engine = "rrr_gls",
    engine_args = list(rank = 1L)
  )

  expect_gt(length(fit$result$event_indices), 0L)
  expect_gt(length(fit$result$baseline_indices), 0L)

  good <- structure(1, colind = fit$result$event_indices[1])
  good_res <- fit_contrasts(fit, list(task = good))
  expect_true("task" %in% names(good_res))

  bad <- structure(1, colind = fit$result$baseline_indices[1])
  expect_error(
    fit_contrasts(fit, list(base = bad)),
    "unsupported coefficient indices"
  )
})
