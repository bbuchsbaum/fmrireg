make_inline_job <- function(id, template) {
  ds <- make_test_matrix_dataset()
  instantiate(template, list(id = id, scans = ds$datamat, TR = 2,
                             run_length = c(40L, 40L), events = ds$event_table))
}

test_that("run_job with no reducer returns a fitted fmri_lm", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
  fit <- run_job(make_inline_job("sub-01", tmpl))
  expect_s3_class(fit, "fmri_lm")
})

test_that("run_job applies the template reducer", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run, reducer = reduce_betas())
  out <- run_job(make_inline_job("sub-02", tmpl))
  expect_s3_class(out, "data.frame")
  expect_true(all(out$job_id == "sub-02"))
})

test_that("run_jobs runs a batch sequentially with per-job records", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run, reducer = reduce_betas())
  jobs <- lapply(c("sub-01", "sub-02", "sub-03"), make_inline_job, template = tmpl)
  res <- run_jobs(jobs)
  expect_s3_class(res, "fmri_batch_result")
  expect_true(all(res$ok))
  vals <- batch_values(res)
  expect_named(vals, c("sub-01", "sub-02", "sub-03"))
  expect_s3_class(vals[["sub-02"]], "data.frame")
})

test_that("run_jobs isolates a failing job", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run, reducer = reduce_betas())
  good <- make_inline_job("sub-ok", tmpl)
  # A deliberately broken job: run_length inconsistent with the data matrix.
  ds <- make_test_matrix_dataset()
  bad <- instantiate(tmpl, list(id = "sub-bad", scans = ds$datamat, TR = 2,
                                run_length = c(40L, 999L),
                                events = ds$event_table))
  res <- run_jobs(list(good, bad))
  expect_equal(sum(res$ok), 1L)
  expect_true("sub-bad" %in% names(batch_errors(res)))
  expect_named(batch_values(res), "sub-ok")

  # on_error = "stop" propagates
  expect_error(run_jobs(list(good, bad), on_error = "stop"))
})

test_that("parallel path via future is equivalent to sequential", {
  skip_if_not_installed("future.apply")
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run, reducer = reduce_betas())
  jobs <- lapply(c("sub-01", "sub-02"), make_inline_job, template = tmpl)

  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  future::plan(future::sequential)

  seq_res <- batch_values(run_jobs(jobs, parallel = FALSE))
  par_res <- batch_values(run_jobs(jobs, parallel = TRUE))
  expect_equal(par_res[["sub-01"]]$estimate, seq_res[["sub-01"]]$estimate)
  expect_equal(par_res[["sub-02"]]$estimate, seq_res[["sub-02"]]$estimate)
})

test_that("parallel multisession round-trips jobs across processes", {
  skip_on_cran()
  skip_if_not_installed("future")
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run, reducer = reduce_betas())
  jobs <- lapply(c("sub-01", "sub-02"), make_inline_job, template = tmpl)

  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  started <- tryCatch({ future::plan(future::multisession, workers = 2); TRUE },
                      error = function(e) FALSE)
  skip_if_not(started, "could not start a multisession plan")

  # Worker processes load the *installed* fmrireg, not the dev (load_all) code.
  # Skip when the installed build lacks this code (e.g. running under load_all);
  # the faithful cross-process check lives in the export_jobs subprocess test.
  worker_ready <- tryCatch(
    future::value(future::future(exists("run_job", where = asNamespace("fmrireg")))),
    error = function(e) FALSE
  )
  skip_if_not(isTRUE(worker_ready),
              "installed fmrireg on workers lacks run_job (load_all); see export_jobs test")

  res <- run_jobs(jobs, parallel = TRUE)
  expect_true(all(res$ok))
  expect_length(batch_values(res), 2)
})
