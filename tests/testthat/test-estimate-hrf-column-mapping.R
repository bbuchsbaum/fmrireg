# estimate_hrf() groups design columns into curves through
# fmridesign::design_colmap(), not through column names. fmridesign 0.6.0 fills
# multi-basis columns condition-major but names them basis-major, so a
# name-based mapping would assemble each curve from other conditions' columns.
# These tests pin the mapping to what the columns actually contain.

colmap_fixture <- function() {
  withr::local_seed(4)
  nscan <- 200L
  sf <- fmrihrf::sampling_frame(blocklens = c(nscan, nscan), TR = 1)
  ev <- do.call(rbind, lapply(1:2, function(r) {
    ons <- sort(sample(seq(5, nscan - 30, by = 4), 30))
    data.frame(run = r, onset = ons, condition = factor(rep(c("fast", "slow"), 15)))
  }))
  list(sf = sf, ev = ev, nscan = nscan, times = fmrihrf::samples(sf, global = TRUE))
}

event_regressor <- function(fx, cond, hrf, delay = 0) {
  e <- fx$ev[fx$ev$condition == cond, ]
  as.matrix(fmrihrf::evaluate(
    fmrihrf::regressor(fmrihrf::global_onsets(fx$sf, e$onset, e$run) + delay, hrf),
    fx$times
  ))
}

test_that("design_colmap() describes what multi-basis columns contain", {
  fx <- colmap_fixture()
  em <- fmridesign::event_model(onset ~ hrf(condition, basis = "spmg3"), data = fx$ev,
                                block = ~run, sampling_frame = fx$sf)
  DM <- as.matrix(fmridesign::design_matrix(em))
  cm <- as.data.frame(fmridesign::design_colmap(em))
  cm <- cm[order(as.integer(cm$col)), ]

  truth <- do.call(cbind, lapply(c("fast", "slow"), function(cond) {
    R <- event_regressor(fx, cond, fmrihrf::HRF_SPMG3)
    colnames(R) <- paste0(cond, "_", seq_len(ncol(R)))
    R
  }))
  best <- apply(abs(stats::cor(DM, truth)), 1, function(r) colnames(truth)[which.max(r)])
  mapped <- paste0(sub("^.*\\.", "", cm$condition), "_", cm$basis_ix)
  expect_equal(unname(mapped), unname(best))
})

test_that("estimate_hrf() recovers each condition's own response shape", {
  fx <- colmap_fixture()
  # "fast" peaks near 5 s; "slow" is the same response delayed by 4 s
  withr::local_seed(5)
  signal <- 2 * event_regressor(fx, "fast", fmrihrf::HRF_SPMG1) +
    2 * event_regressor(fx, "slow", fmrihrf::HRF_SPMG1, delay = 4)
  Y <- sapply(1:3, function(v) signal[, 1] + stats::rnorm(nrow(signal), sd = 0.3)) + 100
  ds <- matrix_frame(Y, TR = 1, run_length = c(fx$nscan, fx$nscan), event_table = fx$ev)

  est <- estimate_hrf(onset ~ hrf(condition), block = ~run, dataset = ds,
                      rsam = seq(0, 20, by = 1), k = 8, lambda = 1)
  tb <- tidy(est)
  mean_curve <- stats::aggregate(estimate ~ curve + time, data = tb, FUN = mean)
  curves <- split(mean_curve, mean_curve$curve)
  names(curves) <- sub("^.*\\.", "", names(curves))
  peak <- vapply(curves, function(d) d$time[which.max(d$estimate)], numeric(1))

  # absolute peak times, and each curve's shape against its true response
  expect_lte(abs(peak[["fast"]] - 5), 1)
  expect_lte(abs(peak[["slow"]] - 9), 1)
  t <- curves$fast$time
  expect_gt(stats::cor(curves$fast$estimate, fmrihrf::HRF_SPMG1(t)), 0.95)
  expect_gt(stats::cor(curves$slow$estimate, fmrihrf::HRF_SPMG1(pmax(t - 4, 0))), 0.95)
})
