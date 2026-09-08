# estimate_hrf.R: S3 methods and reorder/mapping internals still lightly covered.

make_hrf_fit <- function() {
  set.seed(44)
  n <- 120L
  onsets <- seq(8, n - 28, by = 14)
  events <- data.frame(
    onset = onsets,
    condition = factor(rep(c("A", "B"), length.out = length(onsets))),
    run = 1L
  )
  Y <- matrix(rnorm(n * 2), n, 2)
  colnames(Y) <- c("v1", "v2")
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = events)
  estimate_hrf(
    onset ~ hrf(condition),
    block = ~ run,
    dataset = dset,
    basis = "tent",
    k = 5L,
    lambda = 0.2,
    progress = FALSE
  )
}

test_that("fmri_hrf_estimate print/coef/as.matrix/predict methods", {
  fit <- make_hrf_fit()
  expect_s3_class(fit, "fmri_hrf_estimate")

  expect_output(print(fit), "fmri_hrf_estimate")
  expect_output(print(fit), "basis")

  cf <- coef(fit)
  expect_true(is.matrix(cf) || is.array(cf) || is.numeric(cf))

  # Multi-curve requires curve argument
  expect_error(as.matrix(fit), "curve must be supplied")
  m1 <- as.matrix(fit, curve = 1L)
  expect_equal(nrow(m1), length(fit$time))
  expect_equal(ncol(m1), length(fit$voxels))
  m2 <- as.matrix(fit, curve = fit$curves[[2]], what = "std.error")
  expect_equal(dim(m2), dim(m1))
  expect_true(all(is.finite(m2) | is.na(m2)))

  expect_error(as.matrix(fit, curve = 99L), "out of bounds")
  expect_error(as.matrix(fit, curve = "nope"), "identify one")

  pred <- predict(fit, newdata = fit$time[seq_len(min(5L, length(fit$time)))])
  expect_true(is.array(pred) || is.list(pred))

  pred_se <- predict(
    fit,
    newdata = fit$time[1:3],
    se.fit = TRUE
  )
  expect_true(is.list(pred_se) || is.array(pred_se))

  expect_error(predict(fit, newdata = c(-1, 0)), "within")
  expect_error(predict(fit, newdata = c(0, fit$span + 5)), "within")
  expect_error(predict(fit, newdata = numeric()), "within")
})

test_that(".hrf_curve_map rejects incomplete basis maps and reorders suffixes", {
  # Construct a design whose colnames do not all match _bNN suffix
  design <- matrix(rnorm(40), 10, 4)
  colnames(design) <- c("condA_b01", "condA_b02", "condB_b01", "badname")
  emod <- structure(list(), class = "event_model")
  expect_error(
    fmrireg:::.hrf_curve_map(emod, design, k = 2L),
    "Could not map all event-design columns"
  )

  # Incomplete basis set for a curve
  design2 <- matrix(rnorm(30), 10, 3)
  colnames(design2) <- c("condA_b01", "condA_b02", "condA_b03")
  expect_error(
    fmrireg:::.hrf_curve_map(emod, design2, k = 4L),
    "does not contain exactly basis components"
  )

  # Happy path with suffix names (no semantic colmap available)
  design3 <- matrix(rnorm(40), 10, 4)
  colnames(design3) <- c("condA_b02", "condA_b01", "condB_b01", "condB_b02")
  ordered <- fmrireg:::.hrf_curve_map(emod, design3, k = 2L)
  expect_equal(ncol(ordered$design), 4L)
  expect_equal(length(ordered$indices), 2L)
  expect_true(all(c("curve", "term", "condition") %in% names(ordered$info)))
})
