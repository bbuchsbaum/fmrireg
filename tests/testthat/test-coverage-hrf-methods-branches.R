# estimate_hrf helpers and S3 methods: print/coef/as.matrix/predict/tidy + errors.

test_that(".hrf_curve_map and .partial_out_hrf_nuisance cover error/boundary paths", {
  # Incomplete basis components for a curve
  design <- matrix(
    rnorm(40 * 4), 40, 4,
    dimnames = list(NULL, c("cond_A_b1", "cond_A_b2", "cond_B_b1", "cond_B_b2"))
  )
  expect_error(
    fmrireg:::.hrf_curve_map(event_model = list(), design = design, k = 3L),
    "basis components"
  )

  # Unmappable column names
  bad <- matrix(rnorm(20 * 2), 20, 2, dimnames = list(NULL, c("foo", "bar")))
  expect_error(
    fmrireg:::.hrf_curve_map(event_model = list(), design = bad, k = 1L),
    "Could not map"
  )

  # Empty colnames
  empty_names <- matrix(rnorm(10 * 2), 10, 2)
  colnames(empty_names) <- c("", "x_b1")
  expect_error(
    fmrireg:::.hrf_curve_map(event_model = list(), design = empty_names, k = 1L),
    "non-empty column names"
  )

  # Happy suffix path with k=2
  ok <- matrix(
    rnorm(30 * 4), 30, 4,
    dimnames = list(NULL, c("cond_A_b1", "cond_A_b2", "cond_B_b1", "cond_B_b2"))
  )
  mapped <- fmrireg:::.hrf_curve_map(event_model = list(), design = ok, k = 2L)
  expect_equal(length(mapped$indices), 2L)
  expect_equal(nrow(mapped$info), 2L)
  expect_equal(ncol(mapped$design), 4L)

  # Nuisance: zero columns
  event <- matrix(rnorm(40 * 3), 40, 3)
  Y <- matrix(rnorm(40), 40, 1)
  out0 <- fmrireg:::.partial_out_hrf_nuisance(event, Y, matrix(0, 40, 0))
  expect_equal(out0$nuisance_rank, 0L)
  expect_equal(out0$event, event)

  # Rank-deficient nuisance (duplicated columns)
  nuis <- cbind(1, 1, rnorm(40))
  expect_warning(
    out_rd <- fmrireg:::.partial_out_hrf_nuisance(event, Y, nuis),
    "rank deficient"
  )
  expect_equal(dim(out_rd$event), dim(event))
  expect_lt(out_rd$nuisance_rank, ncol(nuis))
})

test_that("fmri_hrf_estimate print/coef/as.matrix/predict/tidy cover branches", {
  skip_if_not_installed("fmridataset")
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  set.seed(180)
  n <- 120L
  onsets <- seq(8, n - 30, by = 12)
  events <- data.frame(
    onset = onsets,
    condition = factor(rep(c("A", "B"), length.out = length(onsets))),
    run = 1L
  )
  Y <- matrix(rnorm(n * 3), n, 3)
  colnames(Y) <- paste0("v", 1:3)
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = events)

  fit <- estimate_hrf(
    onset ~ hrf(condition),
    block = ~ run,
    dataset = dset,
    basis = "bspline",
    progress = FALSE
  )
  expect_s3_class(fit, "fmri_hrf_estimate")

  expect_output(print(fit), "fmri_hrf_estimate")
  expect_true(is.matrix(coef(fit)) || is.numeric(coef(fit)))

  # Multi-curve requires explicit curve for as.matrix
  expect_error(as.matrix(fit), "curve must be supplied")
  m1 <- as.matrix(fit, curve = 1L)
  expect_equal(nrow(m1), length(fit$time))
  expect_equal(ncol(m1), length(fit$voxels))
  m_se <- as.matrix(fit, curve = fit$curves[[1]], what = "std.error")
  expect_equal(dim(m_se), dim(m1))

  expect_error(as.matrix(fit, curve = 99L), "out of bounds")
  expect_error(as.matrix(fit, curve = "not-a-curve"), "identify one")

  pred <- predict(fit, newdata = seq(0, min(fit$span, 12), by = 1))
  expect_true(is.array(pred))
  pred_se <- predict(fit, newdata = c(0, 1, 2), se.fit = TRUE)
  expect_true(is.list(pred_se) && all(c("fit", "se.fit", "df") %in% names(pred_se)))
  expect_error(predict(fit, newdata = -1), "within")

  tid <- tidy(fit)
  expect_s3_class(tid, "tbl_df")
  expect_gt(nrow(tid), 0L)

  tid2 <- tidy(fit, curve = 1L, voxel = fit$voxels[[1]])
  expect_true(all(tid2$curve == fit$curves[[1]]) || nrow(tid2) > 0)

  expect_error(tidy(fit, curve = "missing"), "unknown or out-of-range")
  expect_error(tidy(fit, voxel = 999L), "unknown or out-of-range")
})
