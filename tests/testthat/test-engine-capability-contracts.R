make_engine_call <- function(entry, dset, engine, ..., lowrank = NULL) {
  if (identical(entry, "formula")) {
    return(fmri_lm(
      onsets ~ hrf(condition),
      block = ~run,
      dataset = dset,
      engine = engine,
      lowrank = lowrank,
      ...
    ))
  }

  model <- create_fmri_model(
    formula = onsets ~ hrf(condition),
    block = ~run,
    dataset = dset
  )
  fmri_lm(
    model,
    dataset = dset,
    engine = engine,
    lowrank = lowrank,
    ...
  )
}

test_that("engine capability gate is consistent across formula and fmri_model entry points", {
  dset <- .demo_matrix_dataset()
  entries <- c("formula", "model")

  for (entry in entries) {
    expect_error(
      make_engine_call(
        entry = entry,
        dset = dset,
        engine = "sketch",
        volume_weights = TRUE
      ),
      "latent_sketch does not support volume_weights or soft_subspace preprocessing"
    )

    expect_error(
      make_engine_call(
        entry = entry,
        dset = dset,
        engine = "rrr_gls",
        engine_args = list(rank = 1L),
        ar_voxelwise = TRUE
      ),
      "rrr_gls supports only shared \\(non-voxelwise\\) temporal covariance"
    )

    expect_error(
      make_engine_call(
        entry = entry,
        dset = dset,
        engine = "latent_sketch",
        ar_options = list(by_cluster = TRUE, order = 1L)
      ),
      "latent_sketch requires `lowrank\\$parcels` when ar_options\\$by_cluster = TRUE"
    )
  }
})

test_that("engines reject variance and df methods outside their contract", {
  dset <- .demo_matrix_dataset()
  model <- .demo_fmri_model()

  expect_error(
    fmri_lm(
      model, dataset = dset,
      control = fmri_lm_control(variance = variance_spec("hac", max_lag = 2L)),
      engine = "rrr_gls", engine_args = list(rank = 1L)
    ),
    "does not support variance method 'hac'"
  )
  expect_error(
    fmri_lm(
      model, dataset = dset,
      control = fmri_lm_control(
        variance = variance_spec("model", df = "satterthwaite")
      ),
      engine = "latent_sketch"
    ),
    "does not support df method 'satterthwaite'"
  )
})

test_that("derived engine controls keep canonical specs and kernel aliases aligned", {
  cfg <- fmri_lm_control(
    noise = noise_spec("ar1", pooling = "parcel", parcels = c(1L, 1L),
                       voxelwise = TRUE),
    robust = robust_spec("huber"),
    weights = weights_spec("tukey"),
    projection = projection_spec("soft_subspace", nuisance_matrix = matrix(1, 4, 1))
  )

  executed <- fmrireg:::.derive_engine_execution_config(
    cfg,
    list(robust = FALSE, preprocessing = FALSE,
         ar_voxelwise = FALSE, ar_by_cluster = FALSE)
  )

  expect_identical(executed$robust$type, "none")
  expect_identical(executed$noise$voxelwise, FALSE)
  expect_identical(executed$noise$pooling, "run")
  expect_identical(executed$noise, executed$ar)
  expect_identical(executed$weights, executed$volume_weights)
  expect_identical(executed$projection, executed$soft_subspace)
  expect_identical(executed$weights$method, "none")
  expect_identical(executed$projection$method, "none")
})
