test_that("latent_sketch rejects unsupported robust and preprocessing modes", {
  dset <- .demo_matrix_dataset()
  Tlen <- nrow(fmridataset::get_data_matrix(dset))
  low <- lowrank_control(time_sketch = list(method = "gaussian", m = min(16L, Tlen)))
  
  expect_error(
    fmri_lm(
      onsets ~ hrf(condition),
      block = ~run,
      dataset = dset,
      engine = "latent_sketch",
      lowrank = low,
      robust = TRUE
    ),
    "does not support robust fitting"
  )
  
  expect_error(
    fmri_lm(
      onsets ~ hrf(condition),
      block = ~run,
      dataset = dset,
      engine = "latent_sketch",
      lowrank = low,
      volume_weights = TRUE
    ),
    "does not support volume_weights or soft_subspace preprocessing"
  )
})

test_that("latent_sketch validates by_cluster requirements", {
  dset <- .demo_matrix_dataset()
  Tlen <- nrow(fmridataset::get_data_matrix(dset))
  low <- lowrank_control(time_sketch = list(method = "gaussian", m = min(16L, Tlen)))
  
  expect_error(
    fmri_lm(
      onsets ~ hrf(condition),
      block = ~run,
      dataset = dset,
      engine = "latent_sketch",
      lowrank = low,
      ar_options = list(by_cluster = TRUE, order = 1L)
    ),
    "requires `lowrank\\$parcels`"
  )

  skip_if_not_installed("fmristore")
  skip_if_not_installed("neuroim2")

  latent_fixture <- local({
    set.seed(11)
    n_comp <- 3L
    n_voxels <- 12L
    basis <- matrix(rnorm(Tlen * n_comp), nrow = Tlen, ncol = n_comp)
    loadings <- matrix(rnorm(n_voxels * n_comp), nrow = n_voxels, ncol = n_comp)
    lvec <- fmristore::LatentNeuroVec(
      basis = basis,
      loadings = loadings,
      space = neuroim2::NeuroSpace(c(3, 2, 2, Tlen)),
      mask = rep(TRUE, n_voxels),
      offset = rep(0, n_voxels)
    )
    ds <- fmridataset::latent_dataset(
      source = list(lvec),
      TR = 1,
      run_length = Tlen
    )
    ds$event_table <- data.frame(onsets = 1, condition = factor("A"), run = 1L)
    ds
  })

  expect_error(
    fmri_lm(
      onsets ~ hrf(condition),
      block = ~run,
      dataset = latent_fixture,
      engine = "latent_sketch",
      lowrank = lowrank_control(
        parcels = seq_len(ncol(fmridataset::get_data_matrix(dset))),
        time_sketch = list(method = "gaussian", m = min(16L, Tlen))
      ),
      ar_options = list(by_cluster = TRUE, order = 1L)
    ),
    "does not support by_cluster AR whitening for latent_dataset inputs"
  )
})

test_that("latent_sketch preserves model-defined contrasts and covariance", {
  dset <- .demo_matrix_dataset()
  Tlen <- nrow(fmridataset::get_data_matrix(dset))
  con <- contrast_set(pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B"))
  
  fit <- fmri_lm(
    onsets ~ hrf(condition, contrasts = con),
    block = ~run,
    dataset = dset,
    engine = "latent_sketch",
    lowrank = lowrank_control(time_sketch = list(method = "gaussian", m = min(16L, Tlen)))
  )
  
  expect_s3_class(fit, "fmri_lm")
  expect_gt(nrow(fit$result$contrasts), 0L)
  expect_true(any(grepl("A_vs_B", fit$result$contrasts$name, fixed = TRUE)))
  expect_gt(length(fit$bcons), 0L)
  expect_true(is.matrix(fit$result$cov.unscaled))
  
  posthoc <- fit_contrasts(fit, list(task = structure(1, colind = fit$result$event_indices[1])))
  expect_true("task" %in% names(posthoc))
})

test_that("rrr_gls rejects unsupported preprocessing modes", {
  dset <- .demo_matrix_dataset()
  
  expect_error(
    fmri_lm(
      onsets ~ hrf(condition),
      block = ~run,
      dataset = dset,
      engine = "rrr_gls",
      engine_args = list(rank = 1L),
      volume_weights = TRUE
    ),
    "does not support volume_weights or soft_subspace preprocessing"
  )
})
