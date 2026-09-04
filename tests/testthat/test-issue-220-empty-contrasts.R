# Regression tests for issue #220.
#
# With estimation_spec("runwise_meta"), more than one run, and no contrasts,
# meta_contrasts() used to synthesise a placeholder row named "meta_contrast"
# whose estimate was NULL. It reached fit$result$contrasts, leaked into the
# public accessors, and made write_results() abort at the write step -- after
# the model had already been estimated.

library(testthat)
library(fmrireg)

make_runwise_fit <- function(with_contrast = FALSE, n_vox = 8L) {
  set.seed(1)
  nvols <- 120L
  tr <- 2
  onsets <- seq(6, by = 12, length.out = 9)
  ev <- data.frame(
    run = rep(1:2, each = 9),
    onset = rep(onsets, 2),
    cond = factor(rep(c("a", "b", "c"), 6))
  )
  sframe <- fmrihrf::sampling_frame(blocklens = c(nvols, nvols), TR = tr)

  rhs <- if (with_contrast) {
    con <- contrast_set(pair_contrast(~ cond == "a", ~ cond == "b", name = "a_gt_b"))
    fmridesign::event_model(
      onset ~ hrf(cond, basis = "spmg1", contrasts = con),
      data = ev, block = ~ run, sampling_frame = sframe
    )
  } else {
    fmridesign::event_model(
      onset ~ hrf(cond, basis = "spmg1"),
      data = ev, block = ~ run, sampling_frame = sframe
    )
  }

  Y <- matrix(rnorm(2 * nvols * n_vox), 2 * nvols, n_vox)
  dset <- fmridataset::matrix_dataset(
    Y, TR = tr, run_length = c(nvols, nvols), event_table = ev
  )
  fmod <- fmri_model(
    rhs,
    fmridesign::baseline_model(basis = "bs", degree = 3, sframe = sframe),
    dset
  )

  fmri_lm(
    fmod,
    control = fmri_lm_control(estimation = estimation_spec("runwise_meta")),
    compute = compute_spec(voxel_chunks = 1, parallel = "none", progress = FALSE)
  )
}

test_that("runwise_meta with no contrasts leaves the contrast table empty", {
  fit <- make_runwise_fit(with_contrast = FALSE)
  ct <- fit$result$contrasts

  expect_equal(nrow(ct), 0L)
  expect_false("meta_contrast" %in% ct$name)
})

test_that("the empty contrast table keeps the columns real contrasts have", {
  fit <- make_runwise_fit(with_contrast = FALSE)
  ct <- fit$result$contrasts
  populated <- make_runwise_fit(with_contrast = TRUE)$result$contrasts

  # A zero-column tibble would make dplyr::filter(ct, type == "contrast") fail
  # with "object 'type' not found", so the columns must be present even when
  # there are no rows.
  expect_true(all(c("type", "name", "stat_type", "conmat", "colind", "data") %in% names(ct)))
  expect_equal(nrow(dplyr::filter(ct, type == "contrast")), 0L)
  expect_true(all(names(ct) %in% names(populated)))
})

test_that("no phantom contrast leaks into the public accessors", {
  fit <- make_runwise_fit(with_contrast = FALSE)

  expect_equal(coef_names(fit, "contrasts"), character(0))
  expect_equal(coef_names(fit, "F"), character(0))
  expect_equal(nrow(tidy(fit, type = "contrasts")), 0L)
})

test_that("single-run and multi-run runwise agree on the empty contrast table", {
  set.seed(2)
  nvols <- 120L
  tr <- 2
  ev <- data.frame(
    run = 1L,
    onset = seq(6, by = 12, length.out = 9),
    cond = factor(rep(c("a", "b", "c"), 3))
  )
  sframe <- fmrihrf::sampling_frame(blocklens = nvols, TR = tr)
  emod <- fmridesign::event_model(
    onset ~ hrf(cond, basis = "spmg1"),
    data = ev, block = ~ run, sampling_frame = sframe
  )
  Y <- matrix(rnorm(nvols * 8L), nvols, 8L)
  dset <- fmridataset::matrix_dataset(Y, TR = tr, run_length = nvols, event_table = ev)
  fit1 <- fmri_lm(
    fmri_model(emod, fmridesign::baseline_model(basis = "bs", degree = 3, sframe = sframe), dset),
    control = fmri_lm_control(estimation = estimation_spec("runwise_meta")),
    compute = compute_spec(voxel_chunks = 1, parallel = "none", progress = FALSE)
  )

  fit2 <- make_runwise_fit(with_contrast = FALSE)

  expect_equal(nrow(fit1$result$contrasts), 0L)
  expect_equal(names(fit1$result$contrasts), names(fit2$result$contrasts))
})

test_that("requested contrasts are still pooled across runs", {
  fit <- make_runwise_fit(with_contrast = TRUE, n_vox = 8L)
  ct <- fit$result$contrasts

  expect_equal(nrow(ct), 1L)
  expect_equal(ct$name, "a_gt_b")
  expect_equal(length(unlist(ct$data[[1]]$estimate)), 8L)
  expect_true(all(is.finite(unlist(ct$data[[1]]$estimate))))
})

test_that("meta_contrasts() has a single definition that honours weighting=", {
  # Two definitions existed (metafuns.R and fmri_lm_internal.R); the shadowed
  # one silently dropped the weighting argument.
  expect_true("weighting" %in% names(formals(fmrireg:::meta_contrasts)))
  expect_true("weighting" %in% names(formals(fmrireg:::meta_betas)))
  expect_equal(nrow(fmrireg:::meta_contrasts(list(), weighting = "equal")), 0L)
  expect_equal(nrow(fmrireg:::meta_contrasts(list())), 0L)
})

# ---------------------------------------------------------------------------
# write_results(): a contrast that cannot be written must not discard the
# outputs that were written successfully. Before this change the staging
# directory was deleted wholesale on any failure, so a bad contrast row threw
# away completed beta maps -- the expensive part of a first-level fit.
# ---------------------------------------------------------------------------

make_image_fit <- function() {
  set.seed(123)
  dims <- c(3, 3, 2)
  n_time <- 50
  scans <- lapply(1:2, function(run) {
    arr <- array(rnorm(prod(dims) * n_time), c(dims, n_time))
    neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(dim = c(dims, n_time)))
  })
  mask <- neuroim2::LogicalNeuroVol(
    array(TRUE, dims), neuroim2::NeuroSpace(dim = dims)
  )
  event_table <- data.frame(
    onset = rep(c(5, 15, 25, 35, 45), 2),
    condition = factor(rep(c("A", "B", "A", "B", "A"), 2)),
    run = rep(1:2, each = 5)
  )
  dset <- fmridataset::fmri_mem_dataset(
    scans = scans, mask = mask, TR = 1.5, event_table = event_table
  )
  con <- contrast_set(pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B"))
  fmri_lm(
    onset ~ hrf(condition, contrasts = con),
    block = ~ run, dataset = dset, durations = 0,
    strategy = "chunkwise", nchunks = 2, progress = FALSE
  )
}

test_that("a failing contrast write keeps the betas and reports what survived", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")

  mod <- make_image_fit()

  # Corrupt the contrast payload the way the #220 sentinel did: a statistic
  # whose length does not match the mask. by_contrast hard-errors on this.
  mod$result$contrasts$data[[1]]$estimate <- list(NULL)

  out_dir <- tempfile()
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  expect_error(
    write_results(
      mod, path = out_dir, subject = "01", task = "test",
      space = "MNI152NLin2009cAsym",
      strategy = "by_contrast", save_betas = TRUE
    ),
    "Failed to write BIDS contrast results"
  )

  # The betas were writable and must still be on disk.
  written <- list.files(out_dir, recursive = TRUE)
  expect_true(any(grepl("betas\\.h5$", written)))
  expect_false(any(grepl("contrast-", written)))

  # No staging directory is left behind.
  expect_false(any(grepl("^\\.", list.files(out_dir, all.files = TRUE, no.. = TRUE))))
})

test_that("a contrast-free fit writes betas without error", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")

  mod <- make_image_fit()
  mod$result$contrasts <- fmrireg:::empty_contrast_table()

  out_dir <- tempfile()
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  res <- write_results(
    mod, path = out_dir, subject = "01", task = "test",
    space = "MNI152NLin2009cAsym",
    strategy = "by_contrast", save_betas = TRUE
  )

  expect_true(file.exists(res$betas$h5))
  expect_false(any(grepl("contrast-", list.files(out_dir, recursive = TRUE))))
})
