test_that("builtin backends are registered", {
  bks <- run_backends()
  expect_true(all(c("sequential", "future") %in% bks))
})

test_that("run_jobs dispatches to a named backend; unknown errors", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run, reducer = reduce_betas())
  ds <- make_test_matrix_dataset()
  jobs <- lapply(c("s1", "s2"), function(id)
    instantiate(tmpl, list(id = id, scans = ds$datamat, TR = 2,
                           run_length = c(40L, 40L), events = ds$event_table)))

  res_seq <- run_jobs(jobs, backend = "sequential")
  expect_true(all(res_seq$ok))
  expect_equal(res_seq$ids, c("s1", "s2"))

  expect_error(run_jobs(jobs, backend = "does-not-exist"), "unknown run backend")
})

test_that("custom backends can be registered and used", {
  register_run_backend("rev_order", function(jobs, run_one, ...) {
    rev(lapply(rev(jobs), run_one))
  }, overwrite = TRUE)
  expect_true("rev_order" %in% run_backends())

  tmpl <- fmri_template(onset ~ hrf(condition), ~ run, reducer = reduce_betas())
  ds <- make_test_matrix_dataset()
  jobs <- lapply(c("s1", "s2"), function(id)
    instantiate(tmpl, list(id = id, scans = ds$datamat, TR = 2,
                           run_length = c(40L, 40L), events = ds$event_table)))
  res <- run_jobs(jobs, backend = "rev_order")
  # order preserved by the backend wrapper
  expect_equal(res$ids, c("s1", "s2"))
  expect_true(all(res$ok))

  expect_error(register_run_backend("sequential", function(...) NULL),
               "already registered")
})
