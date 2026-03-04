library(testthat)
library(fmrireg)

test_that("fitted_hrf returns expected shape for demo fit", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  sample_at <- seq(0, 6, by = 2)

  out <- fitted_hrf(fit, sample_at = sample_at)

  expect_named(out, "condition")
  expect_true(is.matrix(out$condition$pred))
  expect_s3_class(out$condition$design, "tbl_df")

  eterm <- terms(fit$model$event_model)[[1]]
  ncond <- nrow(cells(eterm, exclude_basis = TRUE))
  nvox <- nrow(fit$result$betas$data[[1]]$estimate[[1]])

  expect_equal(dim(out$condition$pred), c(length(sample_at) * ncond, nvox))
  expect_equal(nrow(out$condition$design), length(sample_at) * ncond)
  expect_equal(out$condition$design$time, rep(sample_at, ncond))
})

test_that("fitted_hrf handles multi-basis event terms", {
  dset <- fmrireg:::.demo_matrix_dataset()
  sframe <- fmrireg:::.demo_sampling_frame()
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = sframe)

  fit <- suppressWarnings(
    fmri_lm(
      onsets ~ hrf(condition, basis = "spmg3"),
      block = ~run,
      dataset = dset,
      baseline_model = bmod,
      progress = FALSE
    )
  )

  sample_at <- seq(0, 8, by = 2)
  out <- suppressWarnings(fitted_hrf(fit, sample_at = sample_at))

  eterm <- terms(fit$model$event_model)[[1]]
  ncond <- nrow(cells(eterm, exclude_basis = TRUE))
  nvox <- nrow(fit$result$betas$data[[1]]$estimate[[1]])

  expect_equal(dim(out$condition$pred), c(length(sample_at) * ncond, nvox))
  expect_equal(nrow(out$condition$design), length(sample_at) * ncond)
})

test_that("fitted_hrf validates sample_at values", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())

  expect_error(fitted_hrf(fit, sample_at = c(0, 1, NA_real_)), "finite")
  expect_error(fitted_hrf(fit, sample_at = c(0, 1, Inf)), "finite")

  out <- fitted_hrf(fit, sample_at = numeric(0))
  expect_equal(nrow(out$condition$pred), 0L)
  expect_equal(nrow(out$condition$design), 0L)
})
