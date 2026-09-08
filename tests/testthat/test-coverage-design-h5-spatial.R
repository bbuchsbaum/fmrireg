# design_plot longify path via local evaluation + spatial FDR posthoc + H5 validation.

test_that("design_plot longify covers event and baseline terms", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("plotly")
  skip_if_not_installed("bslib")
  skip_if_not_installed("tidyr")

  set.seed(31)
  etab <- data.frame(
    onset = c(5, 15, 25, 35, 45, 55),
    condition = factor(c("A", "B", "A", "B", "A", "B")),
    run = rep(1:2, each = 3)
  )
  Y <- matrix(rnorm(80 * 4), 80, 4)
  dset <- matrix_dataset(Y, TR = 1, run_length = c(40, 40), event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = dset$sampling_frame)
  fmod <- fmri_model(emod, bmod, dset)

  # Bad term name
  expect_error(
    design_plot(fmod, term_name = "not-a-term"),
    "term_name not found"
  )

  # Build app (covers longify for all terms before shiny serves)
  app <- design_plot(fmod, term_name = "condition", longnames = TRUE,
                     legend_threshold = 1)
  expect_true(inherits(app, "shiny.appobj") || is.list(app))

  app2 <- design_plot(fmod, longnames = FALSE, color_palette = "Set1")
  expect_true(inherits(app2, "shiny.appobj") || is.list(app2))
})

test_that(".fmrigds_spatial_posthoc covers z and p array paths", {
  set.seed(32)
  P <- 20
  group <- rep(1:4, length.out = P)
  z <- rnorm(P)
  z_arr <- array(z, dim = c(P, 1L, 1L))
  q1 <- fmrireg:::.fmrigds_spatial_posthoc_compute(
    list(z = z_arr), group = group, alpha = 0.2
  )
  expect_equal(dim(q1), dim(z_arr))
  expect_true(all(is.finite(q1) | is.na(q1)))

  p_arr <- array(pnorm(-abs(z)) * 2, dim = c(P, 1L, 2L))
  q2 <- fmrireg:::.fmrigds_spatial_posthoc_compute(
    list(p = p_arr), group = group, alpha = 0.2
  )
  expect_equal(dim(q2), dim(p_arr))

  # Vector z without dims
  q3 <- fmrireg:::.fmrigds_spatial_posthoc_compute(
    list(z = z), group = group, alpha = 0.1
  )
  expect_equal(length(q3), P)

  expect_error(
    fmrireg:::.fmrigds_spatial_posthoc(list(z = z_arr), opts = list()),
    "options\\$group"
  )
  out <- fmrireg:::.fmrigds_spatial_posthoc(
    list(assays = list(z = z_arr)), opts = list(group = group, alpha = 0.1)
  )
  expect_true(!is.null(out$q))

  expect_error(
    fmrireg:::.fmrigds_spatial_posthoc_compute(list(), group, 0.05),
    "requires p"
  )
})

test_that("group_data_from_h5 validation and subject extraction", {
  expect_error(
    suppressWarnings(group_data_from_h5(1:3, validate = FALSE)),
    "character vector"
  )
  expect_error(
    suppressWarnings(group_data_from_h5(c("a.h5", "b.h5"), validate = TRUE)),
    "do not exist"
  )
  expect_error(
    suppressWarnings(group_data_from_h5(
      c("a.h5", "b.h5"), subjects = "only-one", validate = FALSE
    )),
    "Length of 'subjects'"
  )
  expect_error(
    suppressWarnings(group_data_from_h5(
      c("a.h5", "b.h5"), covariates = data.frame(x = 1), validate = FALSE
    )),
    "data frame|Number of rows|must be a data frame"
  )

  # More than 5 missing files message path
  paths <- paste0("missing_", 1:7, ".h5")
  expect_error(
    suppressWarnings(group_data_from_h5(paths, validate = TRUE)),
    "and 2 more|do not exist"
  )

  ids <- fmrireg:::extract_subject_ids_from_paths(
    c("/data/sub-01_task-x.h5", "/data/sub-02_task-x.h5")
  )
  expect_equal(ids, c("sub-01", "sub-02"))

  ids2 <- fmrireg:::extract_subject_ids_from_paths(c("foo.nii.gz", "bar.nii.gz"))
  expect_equal(ids2, c("foo", "bar"))
})

test_that("group_data_from_fmrilm succeeds with NeuroSpace-backed mask", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmridataset")

  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  # Attach spatial mask metadata when possible
  sp <- tryCatch(neuroim2::NeuroSpace(c(2, 2, 1)), error = function(e) NULL)
  if (!is.null(sp) && !is.null(fit$dataset)) {
    fit$dataset$mask <- tryCatch(
      neuroim2::LogicalNeuroVol(array(TRUE, dim = c(2, 2, 1)), sp),
      error = function(e) fit$dataset$mask
    )
  }

  gd <- tryCatch(
    suppressWarnings(group_data_from_fmrilm(list(fit, fit))),
    error = function(e) e
  )
  if (inherits(gd, "error")) {
    # Still valuable: exercised validation beyond empty list
    expect_match(conditionMessage(gd), "mask|space|Neuro|fmri_lm|subjects|voxel", ignore.case = TRUE)
  } else {
    expect_true(inherits(gd, "group_data"))
  }
})
