# Regression tests for issues raised by the adversarial (council) review.

test_that("template-level contrasts actually drive reduce_contrasts output (A)", {
  cs <- contrast_set(pair_contrast(~ condition == "A", ~ condition == "B",
                                   name = "AvB"))
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run,
                        contrasts = cs, reducer = reduce_contrasts())
  ds <- make_test_matrix_dataset()
  job <- instantiate(tmpl, list(id = "sub-c", scans = ds$datamat, TR = 2,
                                run_length = c(40L, 40L), events = ds$event_table))
  out <- run_job(job)
  expect_s3_class(out, "data.frame")
  expect_gt(nrow(out), 0)                     # was 0 before the fix (dead field)
  expect_true(all(c("job_id", "term", "voxel", "estimate") %in% names(out)))
  expect_true(all(out$job_id == "sub-c"))
  expect_equal(nrow(out), ncol(ds$datamat))   # one row per voxel for the single contrast
})

test_that("baseline_spec(confounds=) selects only the named columns (B)", {
  ds <- make_test_matrix_dataset()
  conf <- matrix(rnorm(80 * 3), 80, 3,
                 dimnames = list(NULL, c("trans_x", "trans_y", "rot_z")))

  tmpl <- fmri_template(onset ~ hrf(condition), ~ run,
                        baseline = baseline_spec(degree = 3, confounds = "trans_x"))
  job <- instantiate(tmpl, list(id = "s1", scans = ds$datamat, TR = 2,
                                run_length = c(40L, 40L), events = ds$event_table,
                                confounds = conf))
  sel <- fmrireg:::.select_confounds(job$nuisance, "trans_x")
  expect_equal(colnames(sel), "trans_x")

  bd <- fmridesign::design_matrix(build_model(job)$baseline_model)
  expect_equal(sum(grepl("nuis", colnames(bd))), 2L)  # 1 selected col x 2 runs

  # missing selector column is an error (not silently dropped)
  tmpl2 <- fmri_template(onset ~ hrf(condition), ~ run,
                         baseline = baseline_spec(degree = 3, confounds = "nope"))
  job2 <- instantiate(tmpl2, list(id = "s2", scans = ds$datamat, TR = 2,
                                  run_length = c(40L, 40L), events = ds$event_table,
                                  confounds = conf))
  expect_error(build_model(job2), "not found")
})

test_that("job ids must be path-safe (C)", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
  ds <- dataset_spec("matrix_dataset", source = "inline")
  expect_error(fmri_job("a/b", tmpl, ds), "path separator")
  expect_error(fmri_job("..", tmpl, ds), "'\\.'")
  expect_silent(fmri_job("sub-01", tmpl, ds))
})

test_that("preflight validates per-run (list) nuisance (D)", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
  ds <- make_test_matrix_dataset()

  # too few run-matrices
  job1 <- instantiate(tmpl, list(id = "s1", scans = ds$datamat, TR = 2,
                                 run_length = c(40L, 40L), events = ds$event_table,
                                 confounds = list(matrix(rnorm(40 * 2), 40, 2))))
  rep1 <- preflight(job1, on_issue = "collect")
  expect_false(rep1$ok)
  expect_true(any(grepl("nuisance list length", rep1$issues$message)))

  # right count, wrong rows in run 2
  job2 <- instantiate(tmpl, list(id = "s2", scans = ds$datamat, TR = 2,
                                 run_length = c(40L, 40L), events = ds$event_table,
                                 confounds = list(matrix(0, 40, 2), matrix(0, 30, 2))))
  rep2 <- preflight(job2, on_issue = "collect")
  expect_false(rep2$ok)
  expect_true(any(grepl("nuisance run 2 has 30 rows", rep2$issues$message)))
})

test_that("realize_dataset only allows known constructors (F)", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
  ds <- make_test_matrix_dataset()
  job <- instantiate(tmpl, list(id = "s1", scans = ds$datamat, TR = 2,
                                run_length = c(40L, 40L), events = ds$event_table))
  job$dataset_spec$constructor <- "Sys.setenv"  # tamper with a non-dataset fn
  expect_error(realize_dataset(job), "not allowed")
})
