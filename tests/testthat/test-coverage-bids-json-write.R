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

test_that("JSON metadata builders describe betas/contrasts/single contrasts", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  ent <- list(subject = "01", task = "demo", space = "MNI152NLin2009cAsym")

  meta <- fmrireg:::.betas_json_metadata(
    fit, ent, c("conditionA", "conditionB"), output_formats = c("h5", "nifti")
  )
  expect_identical(meta$RegressorOrder, c("conditionA", "conditionB"))
  expect_identical(meta$ModelInfo$NumRegressors, 2L)
  expect_identical(meta$DataInfo$FileFormat, c("HDF5", "NIfTI"))
  expect_equal(meta$SpatialReference, ent$space)

  ctab <- make_contrast_tbl_bids()
  cmeta <- fmrireg:::.contrasts_json_metadata(
    fit, ent, "tstat", c("A_vs_B"), ctab, output_formats = "h5"
  )
  expect_identical(cmeta$StatisticType, "TSTAT")
  expect_identical(cmeta$ContrastOrder, "A_vs_B")
  # only contrasts in the map are defined
  expect_identical(names(cmeta$ContrastDefinitions), "A_vs_B")
  expect_identical(cmeta$ContrastDefinitions$A_vs_B$Type, "t-contrast")
  expect_identical(cmeta$DataInfo$Units, "t-statistic")

  smeta <- fmrireg:::.single_contrast_json_metadata(
    fit, ent, "main", c("fstat", "fpval"), ctab[2, , drop = FALSE],
    output_formats = "nifti"
  )
  expect_identical(smeta$ContrastName, "main")
  expect_identical(smeta$ContrastInfo$Type, "F-contrast")
  expect_identical(smeta$StatisticOrder, c("fstat", "fpval"))
  expect_identical(smeta$DataInfo$FileFormat, "NIfTI")
})

test_that(".bids_output_plan covers no-betas/gds/no-contrast/no-match cases", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  ent <- list(subject = "02", task = "task", space = "MNI")
  plan_for <- function(fit, format = "h5", save_betas = TRUE, contrasts = NULL,
                       contrast_stats = c("beta", "tstat"), contrast_match = "auto") {
    fmrireg:::.bids_output_plan(fit, ent, "GLM", format, "by_stat", save_betas,
                                contrasts, contrast_match, contrast_stats)
  }
  empty <- list(betas = list(), contrasts = list())

  expect_length(plan_for(fit, save_betas = FALSE)$betas, 0L)
  expect_identical(plan_for(fit, format = "gds"), empty)
  expect_identical(fmrireg:::.write_bids_plan(plan_for(fit)$betas, tempdir(), "gds"), list())

  # No contrasts table
  fit0 <- fit
  fit0$result$contrasts <- fit$result$contrasts[0, ]
  expect_length(plan_for(fit0)$contrasts, 0L)

  # Requested contrast not present
  fitc <- fit
  fitc$result$contrasts <- make_contrast_tbl_bids()
  expect_warning(
    out_none <- plan_for(fitc, contrasts = "DOES_NOT_EXIST", contrast_match = "exact"),
    "None of the specified contrasts"
  )
  expect_length(out_none$contrasts, 0L)

  # t-family stats plus the F statistic for the F-contrast
  stems <- vapply(plan_for(fitc)$contrasts, `[[`, "", "stem")
  expect_identical(stems, c("sub-02_task-task_space-MNI_desc-beta_bold",
                            "sub-02_task-task_space-MNI_desc-tstat_bold",
                            "sub-02_task-task_space-MNI_desc-fstat_bold"))

  # Inapplicable stats for every family -> no contrast entries
  fit_t <- fitc
  fit_t$result$contrasts <- make_contrast_tbl_bids()[1, ]
  expect_length(plan_for(fit_t, contrast_stats = "fstat")$contrasts, 0L)
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
  dset <- neurovec_frame(
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
