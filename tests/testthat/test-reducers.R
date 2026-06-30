test_that("builtin reducers are tagged and exempt from the capture warning", {
  for (r in list(reduce_identity(), reduce_contrasts(), reduce_betas(),
                 reduce_write_results())) {
    expect_s3_class(r, "fmri_reducer")
    expect_true(is.function(r))
  }
  # A builtin reducer in a template must NOT trigger the closure-capture warning.
  expect_silent(
    fmri_template(onset ~ hrf(condition), ~ run, reducer = reduce_contrasts())
  )
})

test_that("apply_reducer treats NULL as identity", {
  fit <- make_test_fit()
  job <- fmri_job("sub-01", fmri_template(onset ~ hrf(condition), ~ run),
                  dataset_spec("matrix_dataset", source = "inline"))
  expect_identical(fmrireg:::apply_reducer(NULL, fit, job), fit)
})

test_that("reduce_identity returns the fitted object", {
  fit <- make_test_fit()
  job <- fmri_job("sub-01", fmri_template(onset ~ hrf(condition), ~ run),
                  dataset_spec("matrix_dataset", source = "inline"))
  r <- reduce_identity()
  expect_identical(r(fit, job), fit)
})

test_that("reduce_contrasts / reduce_betas return tidy frames with a job id", {
  fit <- make_test_fit()
  job <- fmri_job("sub-07", fmri_template(onset ~ hrf(condition), ~ run),
                  dataset_spec("matrix_dataset", source = "inline"))

  b <- reduce_betas()(fit, job)
  expect_s3_class(b, "data.frame")
  expect_true("job_id" %in% names(b))
  expect_true(all(b$job_id == "sub-07"))
  expect_gt(nrow(b), 0)

  # No id column when add_id = FALSE
  b2 <- reduce_betas(add_id = FALSE)(fit, job)
  expect_false("job_id" %in% names(b2))
})

test_that("reduce_betas(include_baseline = TRUE) is orientation-robust (no crash)", {
  fit <- make_test_fit()
  job <- fmri_job("sub-bl", fmri_template(onset ~ hrf(condition), ~ run),
                  dataset_spec("matrix_dataset", source = "inline"))
  b <- reduce_betas(include_baseline = TRUE)(fit, job)
  expect_s3_class(b, "data.frame")
  expect_true(all(c("job_id", "term", "voxel", "estimate") %in% names(b)))
  # baseline terms are included and labelled
  expect_true(any(grepl("base", b$term)))
  # condition terms still present
  expect_true(any(grepl("condition", b$term)))
})
