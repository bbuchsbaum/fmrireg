test_that("runwise_lm wrapper matches modular implementation", {
  dset <- .demo_matrix_dataset()
  model <- create_fmri_model(
    formula = onsets ~ hrf(condition),
    block = ~run,
    dataset = dset
  )
  contrast_objects <- prepare_fmri_lm_contrasts(model)$standard
  cfg <- fmri_lm_control()

  res_wrapper <- suppressWarnings(
    fmrireg:::runwise_lm(
      dset = dset,
      model = model,
      contrast_objects = contrast_objects,
      cfg = cfg,
      use_fast_path = TRUE,
      progress = FALSE
    )
  )
  res_impl <- suppressWarnings(
    fmrireg:::runwise_lm_impl(
      dset = dset,
      model = model,
      contrast_objects = contrast_objects,
      cfg = cfg,
      use_fast_path = TRUE,
      progress = FALSE
    )
  )

  expect_equal(
    res_wrapper$betas$data[[1]]$estimate[[1]],
    res_impl$betas$data[[1]]$estimate[[1]],
    tolerance = 1e-10
  )
  expect_equal(res_wrapper$cov.unscaled, res_impl$cov.unscaled, tolerance = 1e-10)
  expect_equal(res_wrapper$event_indices, res_impl$event_indices)
  expect_equal(res_wrapper$baseline_indices, res_impl$baseline_indices)
})

test_that("chunkwise_lm generic dispatches fmri_frame inputs to the frame method", {
  dset <- .demo_matrix_dataset()
  model <- create_fmri_model(
    formula = onsets ~ hrf(condition),
    block = ~run,
    dataset = dset
  )
  contrast_objects <- prepare_fmri_lm_contrasts(model)$standard
  cfg <- fmri_lm_control()

  res_method <- suppressWarnings(
    fmrireg:::chunkwise_lm.fmri_frame(
      x = dset,
      model = model,
      contrast_objects = contrast_objects,
      nchunks = 1,
      cfg = cfg,
      use_fast_path = TRUE,
      progress = FALSE
    )
  )
  res_generic <- suppressWarnings(
    fmrireg:::chunkwise_lm(
      dset,
      model = model,
      contrast_objects = contrast_objects,
      nchunks = 1,
      cfg = cfg,
      use_fast_path = TRUE,
      progress = FALSE
    )
  )

  expect_equal(
    res_generic$betas$data[[1]]$estimate[[1]],
    res_method$betas$data[[1]]$estimate[[1]],
    tolerance = 1e-10
  )
  expect_equal(res_generic$cov.unscaled, res_method$cov.unscaled, tolerance = 1e-10)
})
