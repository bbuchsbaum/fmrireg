# Review 2026-10, pass 3, stream R3-3: export_jobs() checked raw job ids for
# uniqueness, but run_one.R writes results to the *sanitised* id, so "sub 01"
# and "sub_01" (or "sub:01") wrote the same results file and one job's output
# silently overwrote the other's.

.r3_job <- function(id, template) {
  ds <- make_test_matrix_dataset()
  instantiate(template, list(id = id, scans = frame_data(ds), TR = 2,
                             run_length = c(40L, 40L), events = frame_events(ds)))
}

test_that("R3-3: export_jobs rejects ids that collide after sanitising", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run, reducer = reduce_betas())
  jobs <- lapply(c("sub 01", "sub_02", "sub_01", "sub:02", "sub-03"), .r3_job,
                 template = tmpl)
  dir <- withr::local_tempdir()
  err <- tryCatch(export_jobs(jobs, dir), error = conditionMessage)
  expect_true(is.character(err))
  expect_match(err, "sub 01", fixed = TRUE)
  expect_match(err, "sub_01", fixed = TRUE)
  # the second colliding group is reported too; the innocent id is not
  expect_match(err, "'sub_02', 'sub:02'", fixed = TRUE)
  expect_no_match(err, "sub-03", fixed = TRUE)
  # nothing is written when the export is refused
  expect_false(file.exists(file.path(dir, "manifest.rds")))
})

test_that("R3-3: the runner and the uniqueness check use the same file id", {
  ids <- c("sub 01", "a/b", "ok-1.x", "sé")
  expect_identical(fmrireg:::.job_file_id(ids),
                   gsub("[^A-Za-z0-9._-]", "_", ids))
  script <- fmrireg:::.run_one_script()
  expect_true(any(grepl(fmrireg:::.job_file_id_pattern, script, fixed = TRUE)))
})

test_that("R3-3: distinct sanitised ids still export", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run, reducer = reduce_betas())
  jobs <- lapply(c("sub 01", "sub-01"), .r3_job, template = tmpl)
  info <- export_jobs(jobs, withr::local_tempdir())
  expect_identical(info$ids, c("sub 01", "sub-01"))
})
