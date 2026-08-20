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
  tr <- ds$sampling_frame$TR
  bl <- fmrihrf::blocklens(ds$sampling_frame)
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
  expect_equal(ncol(as.matrix(coef(fit))), ncol(fx$ds$datamat))
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
    scans = ds$datamat,
    TR = 2,
    run_length = c(40L, 40L),
    events = ds$event_table
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
