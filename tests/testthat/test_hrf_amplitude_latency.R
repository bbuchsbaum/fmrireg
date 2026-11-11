skip_if_not(
  exists("amplitude_latency", asNamespace("fmrihrf"), inherits = FALSE),
  "fmrihrf::amplitude_latency is required"
)

.al_fixture <- function(shift_b = 1.5, basis = "spmg2", noise = 0.01) {
  set.seed(42)
  # Jittered, well-separated events, sampled where the model samples them
  events <- data.frame(
    onset = cumsum(runif(24, 11, 16)),
    condition = factor(sample(rep(c("A", "B"), 12))),
    run = 1
  )
  time <- fmrihrf::samples(fmrihrf::sampling_frame(200, TR = 2))
  resp <- function(on) as.numeric(fmrihrf::evaluate(fmrihrf::regressor(on), time))
  a <- resp(events$onset[events$condition == "A"])
  b <- resp(events$onset[events$condition == "B"] + shift_b)
  Y <- cbind(100 + a + b, 100 + 2 * a + 0.5 * b) +
    matrix(rnorm(400, sd = noise), 200, 2)
  dset <- matrix_frame(Y, TR = 2, run_length = 200, event_table = events)
  form <- if (identical(basis, "spmg2")) {
    onset ~ hrf(condition, basis = "spmg2")
  } else {
    onset ~ hrf(condition, basis = fmrihrf::HRF_SPMGQ)
  }
  fmri_lm(form, block = ~ run, dataset = dset)
}

test_that("hrf_amplitude_latency recovers amplitude and latency per condition", {
  fit <- .al_fixture()
  al <- hrf_amplitude_latency(fit)
  expect_s3_class(al, "hrf_amplitude_latency")
  expect_equal(dim(al$amplitude), c(2L, 2L))
  expect_equal(colnames(al$amplitude), c("condition_condition.A", "condition_condition.B"))
  expect_equal(al$conditions$basis, c("spmg2", "spmg2"))
  expect_equal(unname(al$amplitude), rbind(c(1, 1), c(2, 0.5)), tolerance = 0.05)
  expect_equal(unname(al$latency[, 1]), c(0, 0), tolerance = 0.2)
  expect_equal(unname(al$latency[, 2]), c(1.5, 1.5), tolerance = 0.2)
  expect_true(all(al$fit > 0.99))

  # The canonical coefficient alone underestimates the delayed condition
  b <- coef(fit)
  expect_lt(b[1, "condition_condition.B_b01"], 0.92)

  boost <- hrf_amplitude_latency(fit, method = "boost")
  expect_true(all(is.na(boost$latency)))
  expect_gt(boost$amplitude[1, 2], b[1, "condition_condition.B_b01"])
  expect_output(print(al), "2 voxels x 2 conditions")
})

test_that("hrf_amplitude_latency works with the quadrature basis", {
  skip_if_not(exists("HRF_SPMGQ", asNamespace("fmrihrf"), inherits = FALSE))
  fit <- .al_fixture(shift_b = 2.5, basis = "spmgq")
  al <- hrf_amplitude_latency(fit)
  expect_equal(unname(al$amplitude), rbind(c(1, 1), c(2, 0.5)), tolerance = 0.05)
  expect_equal(unname(al$latency[, 2]), c(2.5, 2.5), tolerance = 0.2)
})

test_that("hrf_amplitude_latency fails clearly on unsuitable input", {
  fit <- .al_fixture()
  expect_error(hrf_amplitude_latency(fit, term = "nope"), "No multi-basis event term")
  expect_error(hrf_amplitude_latency(list()), "fmri_lm")
  events <- data.frame(onset = seq(4, 180, by = 12),
                       condition = factor(rep(c("A", "B"), length.out = 15)), run = 1)
  dset <- matrix_frame(matrix(rnorm(200), 100, 2), TR = 2, run_length = 100,
                       event_table = events)
  fit1 <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset)
  expect_error(hrf_amplitude_latency(fit1), "no event term with a multi-basis HRF")
})
