test_that("tidy_fitted_hrf returns plottable long table", {
  skip_if_not_installed("fmridataset")
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  set.seed(101)
  n <- 50L
  Y <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  ev <- data.frame(
    onsets = c(6, 16, 26, 36),
    run = 1,
    condition = factor(c("A", "B", "A", "B"))
  )

  dset <- fmridataset::matrix_dataset(
    datamat = Y,
    TR = 1,
    run_length = n,
    event_table = ev
  )

  fit <- fmri_lm(
    onsets ~ hrf(condition),
    block = ~run,
    dataset = dset,
    strategy = "runwise",
    nchunks = 1
  )

  tbl <- tidy_fitted_hrf(
    fit,
    sample_at = 0:10,
    term = "condition",
    term_match = "contains",
    voxel = 1L
  )

  expect_s3_class(tbl, "tbl_df")
  expect_gt(nrow(tbl), 0L)
  expect_true(all(c("term", "time", "condition", "estimate", "voxel") %in% names(tbl)))
  expect_true(all(is.finite(tbl$estimate)))
})
