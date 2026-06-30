test_that("instantiate binds a single inline binding into a realizable job", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run,
                        baseline = baseline_spec(degree = 3))
  ds <- make_test_matrix_dataset()
  b <- list(id = "sub-01", scans = ds$datamat, TR = 2,
            run_length = c(40L, 40L), events = ds$event_table)

  job <- instantiate(tmpl, b)
  expect_s3_class(job, "fmri_job")
  expect_equal(job$id, "sub-01")
  expect_equal(job$dataset_spec$constructor, "matrix_dataset")
  expect_equal(job$dataset_spec$source, "inline")

  # Realize the dataset and assemble the model.
  rds <- realize_dataset(job)
  expect_s3_class(rds, "fmri_dataset")
  model <- build_model(job, rds)
  expect_s3_class(model, "fmri_model")
})

test_that("instantiate over a list of bindings returns one job each", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
  mk <- function(id) {
    ds <- make_test_matrix_dataset()
    list(id = id, scans = ds$datamat, TR = 2, run_length = c(40L, 40L),
         events = ds$event_table)
  }
  jobs <- instantiate(tmpl, list(mk("sub-01"), mk("sub-02"), mk("sub-03")))
  expect_length(jobs, 3)
  expect_equal(vapply(jobs, `[[`, character(1), "id"),
               c("sub-01", "sub-02", "sub-03"))
  expect_true(all(vapply(jobs, inherits, logical(1), "fmri_job")))
})

test_that("file-backed binding yields a tiny, path-only serializable job", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
  b <- list(id = "sub-09",
            scans = sprintf("sub-09/run-%d_bold.nii.gz", 1:2),
            TR = 2, run_length = c(200L, 200L),
            base_path = "/study/derivatives",
            meta = list(subject = "09", task = "stroop", space = "MNI152"))
  job <- instantiate(tmpl, b)
  expect_equal(job$dataset_spec$constructor, "fmri_dataset")
  expect_equal(job$dataset_spec$source, "file")
  expect_equal(job$meta$task, "stroop")

  f <- withr::local_tempfile(fileext = ".rds")
  saveRDS(job, f)
  expect_lt(file.info(f)$size, 50 * 1024)
  expect_equal(readRDS(f)$dataset_spec$args$scans, b$scans)
})

test_that("confounds are resolved into per-run nuisance regressors in the model", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run,
                        baseline = baseline_spec(degree = 3))
  ds <- make_test_matrix_dataset()
  conf <- matrix(rnorm(80 * 2), 80, 2,
                 dimnames = list(NULL, c("trans_x", "rot_z")))
  b <- list(id = "sub-05", scans = ds$datamat, TR = 2,
            run_length = c(40L, 40L), events = ds$event_table, confounds = conf)
  job <- instantiate(tmpl, b)

  nl <- fmrireg:::.resolve_nuisance(job$nuisance, c(40L, 40L))
  expect_length(nl, 2)
  expect_equal(nrow(nl[[1]]), 40)

  model <- build_model(job)
  bd <- fmridesign::design_matrix(model$baseline_model)
  expect_true(any(grepl("nuis", colnames(bd))))
})

test_that("instantiate over a manifest data.frame works", {
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
  ds <- make_test_matrix_dataset()
  mani <- data.frame(id = c("sub-01", "sub-02"), TR = c(2, 2),
                     stringsAsFactors = FALSE)
  mani$scans <- list(ds$datamat, ds$datamat)
  mani$run_length <- list(c(40L, 40L), c(40L, 40L))
  mani$events <- list(ds$event_table, ds$event_table)

  jobs <- instantiate(tmpl, mani)
  expect_length(jobs, 2)
  expect_equal(jobs[[2]]$id, "sub-02")
  expect_s3_class(build_model(jobs[[1]]), "fmri_model")
})

test_that(".as_binding_list keeps single-column data.frame list cells intact", {
  df <- data.frame(id = "s1", stringsAsFactors = FALSE)
  df$ev <- list(data.frame(onset = c(1, 2, 3)))  # single-column df in a list-column
  b <- fmrireg:::.as_binding_list(df)[[1]]
  expect_s3_class(b$ev, "data.frame")
  expect_equal(names(b$ev), "onset")
})
