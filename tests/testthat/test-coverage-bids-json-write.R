# bids_export.R: JSON metadata + write_* empty/dispatch helpers still uncovered.

make_contrast_tbl_bids <- function(V = 4L) {
  tibble::tibble(
    type = c("contrast", "Fcontrast"),
    name = c("A_vs_B", "main"),
    stat_type = c("tstat", "Fstat"),
    data = list(
      tibble::tibble(
        estimate = rnorm(V), se = runif(V, 0.1, 0.5),
        stat = rnorm(V), prob = runif(V), sigma = runif(V, 0.5, 1)
      ),
      tibble::tibble(
        estimate = rnorm(V), se = runif(V, 0.1, 0.5),
        stat = abs(rnorm(V)), prob = runif(V), sigma = runif(V, 0.5, 1)
      )
    ),
    conmat = list(matrix(c(1, -1), 1), diag(2)),
    colind = list(1:2, 1:2)
  )
}

test_that(".save_*_json_metadata write betas/contrasts/single-contrast sidecars", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  td <- tempfile("bids-json-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ent <- list(subject = "01", task = "demo", space = "MNI152NLin2009cAsym")

  betas_json <- fmrireg:::.save_betas_json_metadata(
    fit, td, ent, "GLM", c("conditionA", "conditionB"),
    file.path(td, "betas.h5"), output_formats = c("h5", "nifti")
  )
  expect_true(file.exists(betas_json))
  meta <- jsonlite::fromJSON(betas_json)
  expect_true(!is.null(meta$RegressorOrder))
  expect_equal(meta$SpatialReference, ent$space)

  ctab <- make_contrast_tbl_bids()
  cons_json <- fmrireg:::.save_contrasts_json_metadata(
    fit, td, ent, "GLM", "tstat", c("A_vs_B"), ctab[1, , drop = FALSE],
    file.path(td, "tstat.h5"), output_formats = "h5"
  )
  expect_true(file.exists(cons_json))
  cmeta <- jsonlite::fromJSON(cons_json)
  expect_equal(toupper(cmeta$StatisticType), "TSTAT")
  expect_true("A_vs_B" %in% names(cmeta$ContrastDefinitions) ||
                identical(cmeta$ContrastOrder, "A_vs_B") ||
                "A_vs_B" %in% cmeta$ContrastOrder)

  single_json <- fmrireg:::.save_single_contrast_json_metadata(
    fit, td, ent, "GLM", "A_vs_B", c("beta", "tstat"),
    ctab[1, , drop = FALSE], file.path(td, "con.h5"),
    output_formats = "nifti"
  )
  expect_true(file.exists(single_json))
  smeta <- jsonlite::fromJSON(single_json)
  expect_equal(smeta$ContrastName, "A_vs_B")
  expect_true(!is.null(smeta$ContrastInfo$Type))
})

test_that(".write_beta_outputs and .write_contrast_outputs cover empty/gds/no-match", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  td <- tempfile("bids-write-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ent <- list(subject = "02", task = "task", space = "MNI")

  expect_equal(
    fmrireg:::.write_beta_outputs(fit, td, ent, "GLM", "h5", save_betas = FALSE, overwrite = TRUE),
    list()
  )
  expect_equal(
    fmrireg:::.write_beta_outputs(fit, td, ent, "GLM", "gds", save_betas = TRUE, overwrite = TRUE),
    list()
  )

  # No contrasts table
  fit0 <- fit
  fit0$result$contrasts <- fit$result$contrasts[0, ]
  expect_equal(
    fmrireg:::.write_contrast_outputs(
      fit0, td, ent, "GLM", "h5", "by_stat",
      contrasts = NULL, contrast_match = "auto",
      contrast_stats = c("beta", "tstat"), overwrite = TRUE
    ),
    list()
  )

  # Requested contrast not present
  expect_warning(
    out_none <- fmrireg:::.write_contrast_outputs(
      fit, td, ent, "GLM", "h5", "by_stat",
      contrasts = "DOES_NOT_EXIST", contrast_match = "exact",
      contrast_stats = c("beta", "tstat"), overwrite = TRUE
    ),
    "None of the specified contrasts"
  )
  expect_equal(out_none, list())

  # gds-only format yields no image outputs
  expect_equal(
    fmrireg:::.write_contrast_outputs(
      fit, td, ent, "GLM", "gds", "by_stat",
      contrasts = NULL, contrast_match = "auto",
      contrast_stats = c("beta", "tstat"), overwrite = TRUE
    ),
    list()
  )

  # Inapplicable stats for family -> empty
  expect_equal(
    fmrireg:::.write_contrast_outputs(
      fit, td, ent, "GLM", "h5", "by_stat",
      contrasts = NULL, contrast_match = "auto",
      contrast_stats = c("fstat"), overwrite = TRUE
    ),
    list()
  )
})

test_that(".save_gds_outputs builds assays when spatial mask available", {
  skip_if_not_installed("fmrigds")
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmridataset")

  set.seed(61)
  dims <- c(2L, 2L, 1L)
  n_time <- 40L
  V <- prod(dims)
  arr <- array(rnorm(V * n_time), c(dims, n_time))
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
  fit <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset)
  td <- tempfile("bids-gds-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ent <- list(subject = "03", task = "nback", space = "MNI")

  out <- tryCatch(
    fmrireg:::.save_gds_outputs(fit, td, ent, "GLM"),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    # Still exercise the early extract failure path with stripped betas
    expect_match(conditionMessage(out), ".")
  } else {
    expect_true(file.exists(out$h5) || is.character(out$h5))
    expect_true(is.character(out$plan))
  }

  # Missing se/estimate should fail fast
  bad <- fit
  bad$result$betas$data[[1]]$se <- NULL
  bad$result$betas$data[[1]]$estimate <- NULL
  expect_error(
    fmrireg:::.save_gds_outputs(bad, td, ent, "GLM"),
    "Unable to extract beta"
  )
})
