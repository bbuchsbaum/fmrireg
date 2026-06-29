test_that("preflight passes a well-formed job", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
  ds <- make_test_matrix_dataset()
  job <- instantiate(tmpl, list(id = "sub-01", scans = ds$datamat, TR = 2,
                                run_length = c(40L, 40L), events = ds$event_table))
  rep <- preflight(job)
  expect_s3_class(rep, "fmri_preflight")
  expect_true(rep$ok)
  expect_equal(nrow(rep$issues), 0)
})

test_that("preflight flags missing design columns", {
  tmpl <- fmri_template(onset ~ hrf(condition) + hrf(modulation), ~ run)
  ds <- make_test_matrix_dataset()  # has condition + run, but no 'modulation'
  job <- instantiate(tmpl, list(id = "sub-02", scans = ds$datamat, TR = 2,
                                run_length = c(40L, 40L), events = ds$event_table))
  expect_warning(rep <- preflight(job), "modulation")
  expect_false(rep$ok)
  expect_true(any(grepl("modulation", rep$issues$message)))
})

test_that("preflight flags run_length / data mismatch and errors on demand", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
  ds <- make_test_matrix_dataset()
  job <- instantiate(tmpl, list(id = "sub-03", scans = ds$datamat, TR = 2,
                                run_length = c(40L, 50L), events = ds$event_table))
  expect_error(preflight(job, on_issue = "error"), "run_length")
})

test_that("preflight checks file existence for file-backed jobs when asked", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
  job <- instantiate(tmpl, list(
    id = "sub-04", scans = c("/nope/run-1_bold.nii.gz", "/nope/run-2_bold.nii.gz"),
    TR = 2, run_length = c(100L, 100L),
    events = data.frame(onset = c(5, 105), condition = factor(c("A", "B")),
                        run = c(1, 2))))
  # Without check_files the missing files are not an issue.
  expect_true(preflight(job, on_issue = "collect")$ok)
  rep <- preflight(job, check_files = TRUE, on_issue = "collect")
  expect_false(rep$ok)
  expect_true(any(grepl("not found", rep$issues$message)))
})

test_that("preflight aggregates issues across a list of jobs", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
  ds <- make_test_matrix_dataset()
  good <- instantiate(tmpl, list(id = "ok", scans = ds$datamat, TR = 2,
                                 run_length = c(40L, 40L), events = ds$event_table))
  bad <- instantiate(tmpl, list(id = "bad", scans = ds$datamat, TR = -1,
                                run_length = c(40L, 40L), events = ds$event_table))
  rep <- preflight(list(good, bad), on_issue = "collect")
  expect_false(rep$ok)
  expect_equal(sort(unique(rep$issues$job_id)), "bad")
})
