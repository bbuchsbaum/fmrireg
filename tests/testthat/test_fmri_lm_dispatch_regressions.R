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

test_that("legacy chunkwise alias dispatches to active implementation", {
  dset <- .demo_matrix_dataset()
  model <- create_fmri_model(
    formula = onsets ~ hrf(condition),
    block = ~run,
    dataset = dset
  )
  contrast_objects <- prepare_fmri_lm_contrasts(model)$standard
  cfg <- fmri_lm_control()

  res_active <- suppressWarnings(
    fmrireg:::chunkwise_lm.fmri_dataset(
      x = dset,
      model = model,
      contrast_objects = contrast_objects,
      nchunks = 1,
      cfg = cfg,
      use_fast_path = TRUE,
      progress = FALSE
    )
  )
  res_legacy <- suppressWarnings(
    fmrireg:::chunkwise_lm.fmri_dataset_old(
      x = dset,
      model = model,
      contrast_objects = contrast_objects,
      nchunks = 1,
      cfg = cfg,
      use_fast_path = TRUE,
      progress = FALSE
    )
  )

  expect_equal(
    res_legacy$betas$data[[1]]$estimate[[1]],
    res_active$betas$data[[1]]$estimate[[1]],
    tolerance = 1e-10
  )
  expect_equal(res_legacy$cov.unscaled, res_active$cov.unscaled, tolerance = 1e-10)
  expect_equal(res_legacy$event_indices, res_active$event_indices)
  expect_equal(res_legacy$baseline_indices, res_active$baseline_indices)
})
