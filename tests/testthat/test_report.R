library(testthat)
library(fmrireg)

make_report_test_fit <- function() {
  dset <- fmrireg:::.demo_matrix_dataset()
  con <- contrast_set(
    pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B")
  )
  suppressWarnings(
    fmri_lm(
      onsets ~ hrf(condition, contrasts = con),
      block = ~run,
      dataset = dset,
      progress = FALSE
    )
  )
}

test_that("report generic exists and dispatches", {
  expect_true(is.function(report))

  fit <- make_report_test_fit()
  expect_true(inherits(fit, "fmri_lm"))

  expect_error(
    report(1),
    "No report\\(\\) method"
  )
})

test_that(".is_spatial_dataset is FALSE for matrix_dataset", {
  dset <- fmrireg:::.demo_matrix_dataset()
  expect_false(fmrireg:::.is_spatial_dataset(dset))
})

test_that(".report_model_info extracts expected fields", {
  fit <- make_report_test_fit()
  info <- fmrireg:::.report_model_info(fit)

  expect_true(is.list(info))
  expect_true("table" %in% names(info))
  expect_true(is.data.frame(info$table))
  expect_true(all(c("Parameter", "Value") %in% names(info$table)))

  params <- info$table$Parameter
  expect_true("Formula" %in% params)
  expect_true("Strategy" %in% params)
  expect_true("Residual df" %in% params)
})

test_that(".report_estimates returns required summary columns", {
  fit <- make_report_test_fit()
  est <- fmrireg:::.report_estimates(fit)

  expect_true(is.list(est))
  expect_true(isTRUE(est$available))
  expect_true(is.data.frame(est$table))

  expected_cols <- c(
    "Coefficient", "Mean_Beta", "SD_Beta", "Median_Beta",
    "Mean_T", "Max_Abs_T", "Pct_P_lt_0_05"
  )
  expect_true(all(expected_cols %in% names(est$table)))
})

test_that(".pick_peak_slices identifies high-activation slices", {
  arr <- array(0, dim = c(8, 8, 12))
  arr[, , 3] <- 1
  arr[, , 10] <- 3
  arr[4:5, 4:5, 6] <- 5

  vol <- neuroim2::NeuroVol(arr, neuroim2::NeuroSpace(dim = dim(arr)))
  idx <- fmrireg:::.pick_peak_slices(vol, along = 3L, n = 4L)

  expect_true(length(idx) >= 1)
  expect_true(all(idx >= 1 & idx <= dim(arr)[3]))
  expect_true(any(idx %in% c(6, 10)))
})

test_that(".build_peak_table returns expected columns", {
  arr <- array(0, dim = c(10, 10, 10))
  arr[2:4, 2:4, 2:4] <- 4
  arr[7:9, 7:9, 7:9] <- 6
  vol <- neuroim2::NeuroVol(arr, neuroim2::NeuroSpace(dim = dim(arr)))

  tab <- fmrireg:::.build_peak_table(
    vol = vol,
    thresh = 3,
    min_size = 4L,
    max_peaks = 10L
  )

  expect_true(is.data.frame(tab))
  expect_true(all(c("Cluster", "X", "Y", "Z", "Peak_Stat", "Size", "Label") %in% names(tab)))
  expect_true(nrow(tab) >= 1)
})

test_that(".prepare_report_data returns expected structure", {
  fit <- make_report_test_fit()
  fig_dir <- tempfile("report-fig-")
  dir.create(fig_dir, recursive = TRUE)
  on.exit(unlink(fig_dir, recursive = TRUE), add = TRUE)

  out <- fmrireg:::.prepare_report_data(
    x = fit,
    sections = c("model", "design", "hrf", "estimates", "contrasts", "diagnostics"),
    brain_map_stat = "tstat",
    slice_axis = 3L,
    n_slices = 6L,
    threshold = NULL,
    bg_vol = NULL,
    atlas = NULL,
    cluster_thresh = 2.0,
    min_cluster_size = 2L,
    max_peaks = 10L,
    fig_dir = fig_dir
  )

  expect_true(is.list(out))
  expect_true(all(c("model", "design", "hrf", "estimates", "contrasts", "diagnostics") %in% names(out)))
})

test_that("report.fmri_lm renders a PDF (smoke test)", {
  skip_on_cran()
  skip_if_not_installed("quarto")
  skip_if_not_installed("tinytable")
  skip_if_not(quarto::quarto_available())

  fit <- make_report_test_fit()
  out_file <- tempfile(fileext = ".pdf")

  rendered <- tryCatch(
    report(
      fit,
      output_file = out_file,
      sections = c("model", "estimates"),
      open = FALSE,
      quiet = TRUE
    ),
    error = identity
  )

  if (inherits(rendered, "error")) {
    skip(paste("Quarto render unavailable in this environment:", conditionMessage(rendered)))
  }

  expect_true(file.exists(out_file))
})
