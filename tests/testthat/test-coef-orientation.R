# Regression: coef.fmri_lm orientation is voxels × terms in every mode
# (event-only, baseline-inclusive, contrasts), including square-matrix cases
# where dim(coef) alone cannot detect a transpose.

.make_coef_orientation_fit <- function(n_vox, n_cond, seed = 227L) {
  set.seed(seed)
  n_time <- 100L
  onsets <- seq(10, 190, by = 20)
  n_evt <- length(onsets)
  stopifnot(n_cond >= 2L, n_evt >= n_cond)
  cond <- factor(rep(LETTERS[seq_len(n_cond)], length.out = n_evt))
  ev <- data.frame(onset = onsets, cond = cond, run = 1L)
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  dset <- matrix_frame(Y, TR = 2, run_length = n_time, event_table = ev)

  # Named pairwise contrast when exactly two conditions; otherwise no formula contrast
  fml <- if (n_cond == 2L) {
    onset ~ hrf(
      cond,
      contrasts = pair_contrast(~ cond == "A", ~ cond == "B", name = "AvB")
    )
  } else {
    onset ~ hrf(cond)
  }

  fmri_lm(
    fml,
    block = ~ run,
    dataset = dset,
    control = fmri_lm_control(estimation = estimation_spec("joint")),
    progress = FALSE
  )
}

.expect_voxels_by_terms <- function(mat, n_vox, term_names = NULL) {
  m <- as.matrix(mat)
  expect_equal(nrow(m), n_vox)
  expect_true(!is.null(colnames(m)))
  expect_null(rownames(m))
  if (!is.null(term_names)) {
    expect_equal(colnames(m), term_names)
  }
  invisible(m)
}

test_that("coef.fmri_lm is voxels × terms for event, baseline, and contrast modes", {
  n_vox <- 3L
  fit <- .make_coef_orientation_fit(n_vox = n_vox, n_cond = 2L)

  b <- .expect_voxels_by_terms(coef(fit), n_vox)
  expect_equal(ncol(b), 2L)
  expect_true(all(grepl("cond", colnames(b))))

  bb <- .expect_voxels_by_terms(coef(fit, include_baseline = TRUE), n_vox)
  dm_names <- colnames(design_matrix(fit$model))
  expect_equal(colnames(bb), dm_names)
  expect_gt(ncol(bb), ncol(b))
  # Event columns match the event-only accessor (same values, same margin)
  expect_equal(bb[, colnames(b), drop = FALSE], b, ignore_attr = FALSE)

  cc <- .expect_voxels_by_terms(coef(fit, type = "contrasts"), n_vox, "AvB")
  expect_equal(ncol(cc), 1L)
})

test_that("coef.fmri_lm square case cannot hide a transpose via dimensions alone", {
  # Equal number of voxels and event coefficients: dim() is identical either way.
  n <- 3L
  fit <- .make_coef_orientation_fit(n_vox = n, n_cond = n, seed = 228L)

  b <- as.matrix(coef(fit))
  expect_equal(dim(b), c(n, n))
  # Term names must be colnames (voxels × terms), not rownames (terms × voxels)
  expect_true(!is.null(colnames(b)))
  expect_null(rownames(b))
  expect_equal(ncol(b), n)
  expect_true(all(grepl("cond", colnames(b))))

  bb <- as.matrix(coef(fit, include_baseline = TRUE))
  expect_equal(nrow(bb), n)
  expect_equal(colnames(bb), colnames(design_matrix(fit$model)))
  expect_equal(bb[, colnames(b), drop = FALSE], b)

  # Sanity: event-only agrees with the stored beta matrix columns (no silent t())
  stored <- fit$result$betas$data[[1]]$estimate[[1]]
  evt <- fit$result$event_indices
  expect_equal(unname(b), unname(stored[, evt, drop = FALSE]))
})

test_that("coef.fmri_lm event columns align with design_matrix names", {
  fit <- .make_coef_orientation_fit(n_vox = 4L, n_cond = 2L, seed = 229L)
  dm <- design_matrix(fit$model)
  evt <- fit$result$event_indices
  expect_equal(colnames(coef(fit)), colnames(dm)[evt])
  expect_equal(colnames(coef(fit, include_baseline = TRUE)), colnames(dm))
})
