# End-to-end tests that a censor column carried by an fmri_frame actually
# reaches the fit. The pre-existing test-censor-ar-integration.R exercises
# extract_censor_from_dataset() and resolve_censor() directly, which is why a
# frame-carried censor could be silently discarded by every fmri_lm() call
# while those unit tests stayed green.

# A 2-run frame whose censored volumes hold gross artefacts. Censoring has to
# change the AR fit; if it does not, the censor column was ignored.
make_censor_fixture <- function(bad = c(10L, 11L, 12L, 55L, 56L)) {
  set.seed(42)
  runs <- c(40L, 40L)
  TR <- 2
  nvox <- 4L
  total <- sum(runs)
  datamat <- matrix(rnorm(total * nvox), total, nvox)
  datamat[bad, ] <- datamat[bad, ] + 50
  ev <- do.call(rbind, lapply(seq_along(runs), function(r) {
    on <- seq(4, runs[r] * TR - 12, by = 8)
    data.frame(
      onset = on,
      condition = factor(rep(c("A", "B"), length.out = length(on))),
      run = r
    )
  }))
  cen <- rep(FALSE, total)
  cen[bad] <- TRUE
  list(
    censored = matrix_frame(datamat, TR = TR, run_length = runs,
                            event_table = ev, censor = cen),
    plain = matrix_frame(datamat, TR = TR, run_length = runs, event_table = ev),
    bad = bad
  )
}

fit_with <- function(ds, ar_opts) {
  fmri_lm(
    onset ~ hrf(condition), block = ~ run, dataset = ds,
    control = fmri_lm_control(
      estimation = estimation_spec("runwise"),
      ar_options = ar_opts
    ),
    compute = compute_spec()
  )
}

coefs <- function(fit) as.numeric(as.matrix(coef(fit)))

test_that("matrix_frame(censor=) reaches fmri_lm without opting in", {
  fx <- make_censor_fixture()

  # The whole point: no censor= in ar_options, yet the frame's column is used.
  with_censor <- coefs(fit_with(fx$censored, list(struct = "ar1")))
  without <- coefs(fit_with(fx$plain, list(struct = "ar1")))

  expect_false(isTRUE(all.equal(with_censor, without)))
})

test_that("an unset censor matches an explicit censor='auto'", {
  fx <- make_censor_fixture()

  implicit <- coefs(fit_with(fx$censored, list(struct = "ar1")))
  explicit <- coefs(fit_with(fx$censored, list(struct = "ar1", censor = "auto")))

  expect_equal(implicit, explicit)
})

test_that("censor='none' opts out of a frame-carried censor", {
  fx <- make_censor_fixture()

  opted_out <- coefs(fit_with(fx$censored, list(struct = "ar1", censor = "none")))
  no_column <- coefs(fit_with(fx$plain, list(struct = "ar1")))

  expect_equal(opted_out, no_column)
})

test_that("an explicit censor vector still overrides the frame's column", {
  fx <- make_censor_fixture()

  # Censoring different volumes than the frame records must change the answer,
  # proving the config spec wins rather than being merged with the column.
  explicit <- coefs(fit_with(fx$censored,
                             list(struct = "ar1", censor = c(20L, 21L, 70L))))
  from_frame <- coefs(fit_with(fx$censored, list(struct = "ar1")))

  expect_false(isTRUE(all.equal(explicit, from_frame)))
})

test_that("resolve_censor falls back to the dataset when cfg censor is unset", {
  fx <- make_censor_fixture()
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1"))

  expect_equal(resolve_censor(cfg, fx$censored), fx$bad)
  # No dataset and no spec is still NULL, not an error.
  expect_null(resolve_censor(cfg))
  # A frame with no censor column stays NULL.
  expect_null(resolve_censor(cfg, fx$plain))
})

test_that("resolve_censor honours the 'none' opt-out", {
  fx <- make_censor_fixture()
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1", censor = "none"))

  expect_null(resolve_censor(cfg, fx$censored))
})

test_that("noise_spec accepts 'none' and rejects other strings", {
  expect_equal(fmri_lm_control(ar_options = list(struct = "ar1", censor = "none"))$ar$censor,
               "none")
  expect_error(
    fmri_lm_control(ar_options = list(struct = "ar1", censor = "sometimes")),
    "must be NULL"
  )
})

test_that("a censor column under an iid model warns that it is inert", {
  fx <- make_censor_fixture()

  # Censoring only feeds AR estimation, so under iid it cannot do anything.
  # Silence there is the trap this warning exists to close.
  expect_warning(
    fit_with(fx$censored, list(struct = "iid")),
    "no effect under an iid noise model"
  )
  # No column, no warning.
  expect_no_warning(fit_with(fx$plain, list(struct = "iid")))
  # Opting out silences it.
  expect_no_warning(fit_with(fx$censored, list(struct = "iid", censor = "none")))
})
