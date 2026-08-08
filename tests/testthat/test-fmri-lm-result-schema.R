test_that("built-in fits expose the versioned result and variance contract", {
  dset <- .demo_matrix_dataset()
  control <- fmri_lm_control(estimation = estimation_spec("joint"))
  fit <- suppressWarnings(fmri_lm(
    onsets ~ hrf(condition), block = ~run, dataset = dset,
    control = control, compute = compute_spec(voxel_chunks = 1L)
  ))

  expect_identical(fit$result$schema_version, 2L)
  expect_s3_class(variance_model(fit), "fmri_lm_variance_model")
  expect_identical(variance_model(fit)$method, "model")
  expect_identical(variance_model(fit)$metadata$estimation_scope, "joint")
  expect_identical(fit$result$fit_state$voxel_status,
                   fit$result$voxel_status)
  expect_length(fit$result$df$inference, ncol(coef(fit)))
  expect_s3_class(attr(fit, "compute"), "fmri_lm_compute_spec")
  expect_identical(attr(fit, "requested_control"), control)
  expect_s3_class(attr(fit, "executed_control"), "fmri_lm_control")
})

test_that("model-based residual-df fits do not retain inference residuals", {
  dset <- .demo_matrix_dataset()
  control <- fmri_lm_control(
    estimation = estimation_spec("joint"),
    noise = noise_spec("ar1", iter_gls = 2L),
    robust = robust_spec("huber", max_iter = 2L),
    variance = variance_spec("model", df = "residual")
  )
  fit <- suppressWarnings(fmri_lm(
    onsets ~ hrf(condition), block = ~run, dataset = dset,
    control = control, compute = compute_spec(voxel_chunks = 2L)
  ))

  expect_null(fit$result$inference_context)
  expect_null(fit$result$residuals)
  expect_true(length(fit$result$fit_state$robust_weights) > 0L)
  expect_identical(variance_model(fit)$covariance_scope, "shared")
})

test_that("joint results are invariant to voxel partitioning and chunk parallelism", {
  skip_if_not_installed("future.apply")
  dset <- .demo_matrix_dataset()
  control <- fmri_lm_control(estimation = estimation_spec("joint"))
  fit <- function(chunks, parallel) suppressWarnings(fmri_lm(
    onsets ~ hrf(condition), block = ~run, dataset = dset,
    control = control,
    compute = compute_spec(voxel_chunks = chunks, parallel = parallel)
  ))

  one <- fit(1L, "none")
  two <- fit(2L, "none")
  future_chunks <- suppressWarnings(fit(2L, "chunks"))

  expect_equal(coef(two), coef(one), tolerance = 1e-10)
  expect_equal(stats(two, "estimates"), stats(one, "estimates"), tolerance = 1e-10)
  expect_equal(coef(future_chunks), coef(two), tolerance = 1e-10)
  expect_equal(stats(future_chunks, "estimates"),
               stats(two, "estimates"), tolerance = 1e-10)
})

test_that("runwise backend is compute-only for supported model fits", {
  dset <- .demo_matrix_dataset()
  control <- fmri_lm_control(estimation = estimation_spec("runwise_meta"))

  matrix_fit <- suppressWarnings(fmri_lm(
    onsets ~ hrf(condition), block = ~run, dataset = dset,
    control = control, compute = compute_spec(backend = "matrix")
  ))
  reference_fit <- suppressWarnings(fmri_lm(
    onsets ~ hrf(condition), block = ~run, dataset = dset,
    control = control, compute = compute_spec(backend = "reference")
  ))

  expect_equal(coef(reference_fit), coef(matrix_fit), tolerance = 1e-8)
  expect_equal(stats(reference_fit, "estimates"),
               stats(matrix_fit, "estimates"), tolerance = 1e-8)
})
