find_pkg_root <- function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "DESCRIPTION")) && dirname(d) != d) d <- dirname(d)
  d
}

make_inline_job2 <- function(id, template) {
  ds <- make_test_matrix_dataset()
  instantiate(template, list(id = id, scans = ds$datamat, TR = 2,
                             run_length = c(40L, 40L), events = ds$event_table))
}

test_that("export_jobs writes a manifest, ids, and a runner; read_jobs round-trips", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run, reducer = reduce_betas())
  jobs <- lapply(c("sub-01", "sub-02"), make_inline_job2, template = tmpl)

  dir <- withr::local_tempdir()
  info <- export_jobs(jobs, dir)
  expect_true(file.exists(info$manifest))
  expect_true(file.exists(info$runner))
  expect_true(file.exists(file.path(dir, "job_ids.txt")))
  expect_equal(info$ids, c("sub-01", "sub-02"))

  back <- read_jobs(dir)
  expect_length(back, 2)
  expect_equal(back[[2]]$id, "sub-02")
  expect_s3_class(back[[1]], "fmri_job")

  # No clobber without overwrite
  expect_error(export_jobs(jobs, dir), "already exists")
  expect_silent(export_jobs(jobs, dir, overwrite = TRUE))

  expect_error(export_jobs(jobs[c(1, 1)], withr::local_tempdir()), "unique")
})

test_that("run_one.R reconstructs and runs a job in a FRESH process (array-node sim)", {
  skip_on_cran()
  rscript <- file.path(R.home("bin"), "Rscript")
  skip_if_not(file.exists(rscript), "Rscript not found")
  skip_if_not_installed("devtools")

  pkg_root <- find_pkg_root()
  skip_if_not(file.exists(file.path(pkg_root, "DESCRIPTION")), "package root not found")

  tmpl <- fmri_template(onset ~ hrf(condition), ~ run, reducer = reduce_betas())
  jobs <- lapply(c("sub-01", "sub-02"), make_inline_job2, template = tmpl)

  dir <- withr::local_tempdir()
  # The child process loads the dev package (load_all) so it runs THIS code,
  # faithfully simulating a worker node that has fmrireg available.
  setup <- sprintf('suppressMessages(devtools::load_all(%s, quiet = TRUE))',
                   deparse(pkg_root))
  export_jobs(jobs, dir, setup = setup)

  runner <- file.path(dir, "run_one.R")
  status <- system2(rscript, c(shQuote(runner), "2", shQuote("results")),
                    stdout = TRUE, stderr = TRUE)
  out_rds <- file.path(dir, "results", "sub-02.rds")
  expect_true(file.exists(out_rds),
              info = paste(status, collapse = "\n"))

  value <- readRDS(out_rds)
  expect_s3_class(value, "data.frame")
  expect_true(all(value$job_id == "sub-02"))
})
