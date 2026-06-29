test_that("baseline_spec constructs and validates inputs", {
  bs <- baseline_spec(degree = 3, basis = "bs")
  expect_s3_class(bs, "baseline_spec")
  expect_equal(bs$degree, 3)
  expect_error(baseline_spec(degree = -1), "non-negative")
  expect_error(baseline_spec(confounds = 1:3), "confounds")
  expect_silent(baseline_spec(confounds = c("trans_x", "rot_y")))
  # bs/ns require degree >= 3 (fail fast, not on a worker); poly is fine at 2.
  expect_error(baseline_spec(basis = "bs", degree = 2), "degree >= 3")
  expect_silent(baseline_spec(basis = "poly", degree = 2))
})

test_that("fmri_template constructs and validates", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run,
                        baseline = baseline_spec(degree = 3))
  expect_s3_class(tmpl, "fmri_template")
  expect_s3_class(tmpl$baseline, "baseline_spec")
  expect_s3_class(tmpl$control, "fmri_lm_config")
  expect_true(validate_template(tmpl))
  expect_error(fmri_template("not a formula", ~ run), "formula")
  expect_error(fmri_template(onset ~ hrf(condition), ~ run, baseline = list()),
               "baseline_spec")
})

test_that("dataset_spec and fmri_job construct and validate", {
  ds <- dataset_spec("fmri_dataset",
                     args = list(scans = c("a.nii.gz", "b.nii.gz"),
                                 TR = 2, run_length = c(200, 200)),
                     source = "file")
  expect_s3_class(ds, "dataset_spec")
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
  job <- fmri_job("sub-01", tmpl, ds, meta = list(subject = "01", task = "x"))
  expect_s3_class(job, "fmri_job")
  expect_true(validate_job(job))
  expect_error(fmri_job("", tmpl, ds), "id")
  expect_error(fmri_job("sub-01", list(), ds), "template")
})

test_that("file-backed job serializes small and round-trips identically", {
  tmpl <- fmri_template(onset ~ hrf(condition) + hrf(modulation), ~ run,
                        baseline = baseline_spec(degree = 3,
                                                 confounds = c("trans_x", "rot_y")))
  ds <- dataset_spec("fmri_dataset",
                     args = list(scans = sprintf("sub-01/run-%d_bold.nii.gz", 1:3),
                                 TR = 2, run_length = c(200, 200, 200),
                                 base_path = "/study/derivatives"),
                     source = "file")
  job <- fmri_job("sub-01", tmpl, ds, meta = list(subject = "01", task = "stroop"))

  f <- withr::local_tempfile(fileext = ".rds")
  saveRDS(job, f)
  # A file-backed recipe carries paths only -> tiny on disk (no voxel data).
  expect_lt(file.info(f)$size, 50 * 1024)

  back <- readRDS(f)
  expect_s3_class(back, "fmri_job")
  # Round-trip fidelity of the recipe value. (Not expect_identical: formulas
  # carry an environment reference that saveRDS/readRDS necessarily renews.)
  expect_equal(job, back)
  expect_equal(deparse(back$template$formula), deparse(job$template$formula))
  expect_equal(back$dataset_spec$args$run_length, c(200, 200, 200))
})

test_that("template strips formula/block env so local scope does not bloat jobs", {
  build_in_scope <- function() {
    big_local <- numeric(300000)  # ~2.4 MB that must NOT ride along
    fmri_template(onset ~ hrf(condition), ~ run)
  }
  tmpl <- build_in_scope()
  expect_identical(environment(tmpl$formula), globalenv())
  expect_identical(environment(tmpl$block), globalenv())

  ds <- dataset_spec("fmri_dataset",
                     args = list(scans = "x.nii.gz", TR = 2, run_length = 100),
                     source = "file")
  job <- fmri_job("s1", tmpl, ds)
  f <- withr::local_tempfile(fileext = ".rds")
  saveRDS(job, f)
  expect_lt(file.info(f)$size, 50 * 1024)
})

test_that("reducer capturing local state warns; top-level reducer does not", {
  make_capturing_reducer <- function() {
    big_state <- seq_len(1000)
    function(fit, job) big_state
  }
  r <- make_capturing_reducer()
  expect_warning(
    fmri_template(onset ~ hrf(condition), ~ run, reducer = r),
    "serial"
  )
  # A top-level / base function captures nothing -> no warning.
  expect_silent(fmri_template(onset ~ hrf(condition), ~ run, reducer = identity))
})
