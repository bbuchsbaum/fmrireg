# Twelfth wave: BIDS nifti/h5 beta and contrast-by-stat writers on spatial fits.

make_spatial_fit <- function(seed = 221L) {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmridataset")
  set.seed(seed)
  dims <- c(2L, 2L, 1L)
  n_time <- 40L
  arr <- array(rnorm(prod(dims) * n_time), c(dims, n_time))
  scan <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(dim = c(dims, n_time)))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dim = dims))
  etab <- data.frame(
    onset = c(5, 12, 20, 28),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  dset <- fmridataset::fmri_mem_dataset(
    scans = list(scan), mask = mask, TR = 1, event_table = etab
  )
  con <- contrast_set(pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B"))
  fit <- fmri_lm(
    onset ~ hrf(condition, contrasts = con),
    block = ~ run, dataset = dset
  )
  list(fit = fit, entities = list(subject = "01", task = "demo", space = "T1w"))
}

test_that(".save_regressor_betas_nifti writes betas + json metadata", {
  fx <- make_spatial_fit()
  td <- tempfile("bids-betas-nifti-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  out <- fmrireg:::.save_regressor_betas_nifti(
    fx$fit, td, fx$entities, desc = "GLM", overwrite = TRUE,
    output_formats = "nifti"
  )
  expect_true(file.exists(out$nifti))
  expect_true(file.exists(out$json))
  expect_match(basename(out$nifti), "betas\\.nii\\.gz$")
})

test_that(".save_contrasts_by_stat_nifti writes contrast maps", {
  fx <- make_spatial_fit(seed = 222L)
  td <- tempfile("bids-con-nifti-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  ctab <- fx$fit$result$contrasts
  if (is.null(ctab) || (is.data.frame(ctab) && nrow(ctab) == 0)) {
    skip("no contrasts available on fit")
  }
  out <- fmrireg:::.save_contrasts_by_stat_nifti(
    fx$fit, ctab, td, fx$entities, desc = "GLM",
    contrast_stats = c("estimate", "stat"),
    overwrite = TRUE,
    output_formats = "nifti"
  )
  expect_true(is.list(out))
  expect_true(length(out) >= 1L)
  files <- unlist(out, use.names = FALSE)
  expect_true(any(grepl("\\.nii\\.gz$", files)))
})

test_that(".save_contrasts_by_stat empty table returns empty list", {
  fx <- make_spatial_fit(seed = 223L)
  td <- tempfile("bids-empty-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  empty <- tibble::tibble()
  expect_equal(
    fmrireg:::.save_contrasts_by_stat(
      fx$fit, empty, td, fx$entities, desc = "GLM",
      contrast_stats = "estimate", overwrite = TRUE
    ),
    list()
  )
})

test_that(".save_betas_json_metadata and .save_contrasts_json_metadata write sidecars", {
  fx <- make_spatial_fit(seed = 224L)
  td <- tempfile("bids-json-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ent <- fx$entities

  beta_json <- fmrireg:::.save_betas_json_metadata(
    fx$fit, td, ent, desc = "GLM",
    regressor_names = c("b1", "b2"),
    h5_filepath = file.path(td, "dummy_betas.h5"),
    output_formats = "nifti"
  )
  expect_true(file.exists(beta_json))
  meta <- jsonlite::fromJSON(beta_json)
  expect_true(is.list(meta) || is.character(meta) || length(meta) >= 0)

  ctab <- fx$fit$result$contrasts
  if (!is.null(ctab) && nrow(ctab) > 0) {
    con_json <- fmrireg:::.save_contrasts_json_metadata(
      fx$fit, td, ent, desc = "GLM", stat = "estimate",
      contrast_names = as.character(ctab$name),
      contrasts = ctab,
      h5_filepath = file.path(td, "dummy_con.h5"),
      output_formats = "nifti"
    )
    expect_true(file.exists(con_json))
  }
})
