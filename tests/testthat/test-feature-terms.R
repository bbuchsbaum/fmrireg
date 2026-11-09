skip_if_not(
  exists("feature", asNamespace("fmridesign"), inherits = FALSE),
  "fmridesign::feature is required"
)
skip_if_not(
  exists("feature_regressor", asNamespace("fmrihrf"), inherits = FALSE),
  "fmrihrf::feature_regressor is required"
)

.local_feature_dataset <- function() {
  ds <- make_test_matrix_dataset()
  tr <- frame_sframe(ds)$TR
  bl <- fmrihrf::blocklens(frame_sframe(ds))
  dt <- 0.5
  per_run <- lapply(bl, function(n) {
    tvec <- seq(0, by = dt, length.out = max(1L, floor(n * tr / dt)))
    abs(sin(2 * pi * tvec / 8))
  })
  list(ds = ds, dt = dt, rms = per_run)
}

test_that("mixed hrf + feature formula fits through fmri_lm", {
  fx <- .local_feature_dataset()
  rms <- fx$rms
  dt <- fx$dt

  fit <- fmri_lm(
    onset ~ hrf(condition) + feature(rms, dt = dt, id = "rms",
                                     center = FALSE, scale = "none"),
    block = ~ run,
    dataset = fx$ds,
    progress = FALSE
  )
  expect_s3_class(fit, "fmri_lm")

  dm <- design_matrix(fit$model)
  expect_true(any(grepl("rms", colnames(dm), fixed = TRUE)))
  expect_true(any(grepl("condition", colnames(dm), fixed = TRUE)))

  tm <- term_matrices(fit$model)
  expect_gt(length(attr(tm, "event_term_indices")), 0L)
  expect_equal(nrow(as.matrix(coef(fit))), ncol(frame_data(fx$ds)))
})

test_that("fitted_hrf skips feature terms instead of erroring", {
  fx <- .local_feature_dataset()
  rms <- fx$rms
  dt <- fx$dt

  fit <- fmri_lm(
    onset ~ hrf(condition) + feature(rms, dt = dt, id = "rms",
                                     center = FALSE, scale = "none"),
    block = ~ run,
    dataset = fx$ds,
    progress = FALSE
  )

  hrfs <- fitted_hrf(fit, sample_at = c(0, 4, 8))
  expect_type(hrfs, "list")
  expect_false(any(grepl("rms", names(hrfs), ignore.case = TRUE)))
  expect_true(length(hrfs) >= 1L)
  expect_true(all(c("pred", "design") %in% names(hrfs[[1]])))
})

test_that("shortnames and cells accept mixed feature models", {
  fx <- .local_feature_dataset()
  rms <- fx$rms
  dt <- fx$dt

  fit <- fmri_lm(
    onset ~ hrf(condition) + feature(rms, dt = dt, id = "rms",
                                     center = FALSE, scale = "none"),
    block = ~ run,
    dataset = fx$ds,
    progress = FALSE
  )

  expect_silent(sn <- shortnames(fit$model$event_model))
  expect_true(is.character(sn))
  expect_true(length(sn) >= 1L)
  expect_silent(cells(fit$model))
})

test_that("preflight does not require feature series as event columns", {
  fx <- .local_feature_dataset()
  rms <- fx$rms
  dt <- fx$dt
  ds <- fx$ds

  tmpl <- fmri_template(
    onset ~ hrf(condition) + feature(rms, dt = dt, id = "rms"),
    ~ run
  )
  job <- instantiate(tmpl, list(
    id = "sub-01",
    scans = frame_data(ds),
    TR = 2,
    run_length = c(40L, 40L),
    events = frame_events(ds)
  ))
  expect_silent(rep <- preflight(job, on_issue = "collect"))
  expect_true(rep$ok)
})

test_that("inject_basis rewrites feature() basis arguments", {
  rms <- 1:5
  form <- ~ feature(rms, dt = 0.5, id = "rms")
  newform <- fmrireg:::inject_basis(form, fmrihrf::HRF_SPMG2)
  expect_identical(rlang::f_rhs(newform)$basis, fmrihrf::HRF_SPMG2)
})

# Bands, sub-bands and contrasts on feature terms need an fmridesign that
# declares them (design_bands() arrived with that support).
.feature_bands_available <- function() {
  exists("design_bands", asNamespace("fmridesign"), inherits = FALSE) &&
    "durations" %in% names(formals(fmrihrf::feature_regressor))
}

.local_banded_feature_fit <- function() {
  set.seed(4172)
  ds <- make_test_matrix_dataset()
  sf <- frame_sframe(ds)
  bl <- fmrihrf::blocklens(sf)
  tr <- sf$TR[[1L]]
  on <- lapply(bl, function(n) sort(runif(30, 0, n * tr - 10)))
  du <- lapply(bl, function(n) runif(30, 0.05, 0.3))
  phon <- lapply(bl, function(n) {
    matrix(rnorm(30 * 4), 30, 4,
           dimnames = list(NULL, c("stop", "fric", "lab", "voi")))
  })
  fmri_lm(
    onset ~ hrf(condition) +
      feature(phon, times = on, durations = du, id = "phon",
              band = "phonological",
              subband = c("manner", "manner", "place", "voicing"),
              contrasts = column_contrast("stop", "fric", name = "stop_vs_fric")),
    block = ~ run, dataset = ds, progress = FALSE
  )
}

test_that("a contrast declared on a feature term is fitted with the model", {
  skip_if_not(.feature_bands_available(), "fmridesign feature bands are required")
  fit <- .local_banded_feature_fit()
  fitted <- stats(fit, "contrasts")
  expect_true("stop_vs_fric" %in% names(fitted))
  expect_true(all(is.finite(as.numeric(fitted[["stop_vs_fric"]]))))
})

test_that("post-hoc contrast specs resolve against feature terms", {
  skip_if_not(.feature_bands_available(), "fmridesign feature bands are required")
  fit <- .local_banded_feature_fit()
  spec <- fit_contrasts(
    fit, list(lab_vs_voi = column_contrast("lab", "voi", name = "lab_vs_voi"))
  )
  manual <- fit_contrasts(fit, list(m = c(phon_lab = 1, phon_voi = -1)))
  expect_equal(spec$lab_vs_voi$estimate, manual$m$estimate, tolerance = 1e-10)
  expect_equal(spec$lab_vs_voi$stat, manual$m$stat, tolerance = 1e-10)
  # a spec naming an unknown variable resolves against no term
  expect_error(
    fit_contrasts(fit, list(x = pair_contrast(~ zzz == "a", ~ zzz == "b", name = "x"))),
    "does not resolve"
  )
})

test_that("band and sub-band F-tests come from Fcontrasts() on the event model", {
  skip_if_not(.feature_bands_available(), "fmridesign feature bands are required")
  fit <- .local_banded_feature_fit()
  Fc <- fmridesign::Fcontrasts(fit$model$event_model)
  expect_true(all(c("phon#phon", "phon#manner", "phon#place", "phon#voicing") %in%
                    names(Fc)))
  tests <- fit_contrasts(fit, Fc)
  expect_identical(names(tests), names(Fc))
  expect_identical(tests[["phon#manner"]]$stat_type, "Fstat")

  # the manner F is the joint test of its two columns
  L <- rbind(c(phon_stop = 1, phon_fric = 0), c(phon_stop = 0, phon_fric = 1))
  manual <- fit_contrasts(fit, list(m = L))
  expect_equal(tests[["phon#manner"]]$stat, manual$m$stat, tolerance = 1e-8)

  # a one-column sub-band test is that column's t-test
  t_voi <- fit_contrasts(fit, list(m = c(phon_voi = 1)))$m$stat
  expect_equal(tests[["phon#voicing"]]$stat^2, t_voi^2, tolerance = 1e-8)
})
