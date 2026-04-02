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
