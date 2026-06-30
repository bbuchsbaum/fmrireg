skip_if_not_installed("bidser")
skip_if_not_installed("RNifti")
skip_if_not_installed("jsonlite")

make_proj <- function() {
  root <- withr::local_tempdir(.local_envir = parent.frame())
  make_bids_fixture(root)
  bidser::bids_project(root, fmriprep = TRUE)
}

test_that("from_bids builds a per-subject manifest with scans/run_length/TR/events", {
  proj <- make_proj()
  mani <- from_bids(proj, task = "stroop", space = "MNI152NLin2009cAsym",
                    confounds = c("trans_x", "trans_y", "rot_z"))
  expect_s3_class(mani, "fmri_manifest")
  expect_length(mani, 2)

  b <- mani[[1]]
  expect_match(b$id, "^sub-0")
  expect_length(b$scans, 2)                      # two runs
  expect_true(all(file.exists(b$scans)))
  expect_equal(b$run_length, c(40L, 40L))        # nvols read from header
  expect_equal(b$TR, 2)
  expect_s3_class(b$events, "data.frame")
  expect_true(all(c("onset", "trial_type", "run") %in% names(b$events)))
  expect_setequal(unique(b$events$run), c(1L, 2L))

  # confounds: per-run list of matrices, restricted to the selected columns
  expect_length(b$confounds, 2)
  expect_equal(colnames(b$confounds[[1]]), c("trans_x", "trans_y", "rot_z"))
  expect_equal(nrow(b$confounds[[1]]), 40)

  expect_equal(b$meta$task, "stroop")
  expect_equal(b$meta$space, "MNI152NLin2009cAsym")
})

test_that("from_bids manifest instantiates, preflights, and runs end-to-end", {
  proj <- make_proj()
  mask <- file.path(proj$path, "brain_mask.nii.gz")
  mani <- from_bids(proj, task = "stroop", space = "MNI152NLin2009cAsym",
                    confounds = c("trans_x", "trans_y", "rot_z"), mask = mask)

  tmpl <- fmri_template(onset ~ hrf(trial_type), ~ run,
                        baseline = baseline_spec(degree = 3),
                        reducer = reduce_betas())
  jobs <- instantiate(tmpl, mani)
  expect_length(jobs, 2)
  expect_equal(jobs[[1]]$dataset_spec$constructor, "fmri_dataset")
  expect_equal(jobs[[1]]$dataset_spec$source, "file")

  expect_true(preflight(jobs, on_issue = "collect")$ok)

  res <- run_jobs(jobs)
  expect_true(all(res$ok), info = paste(batch_errors(res), collapse = "; "))
  vals <- batch_values(res)
  expect_named(vals, c("sub-01", "sub-02"))
  expect_s3_class(vals[["sub-01"]], "data.frame")
  expect_true(any(grepl("trial_type", vals[["sub-01"]]$term)))
})

test_that("space filter that matches nothing yields no usable subjects", {
  proj <- make_proj()
  expect_error(suppressWarnings(from_bids(proj, task = "stroop", space = "T1w")),
               "no usable subjects")
})

test_that("as_manifest coerces a plain binding list", {
  ds <- make_test_matrix_dataset()
  m <- as_manifest(list(
    list(id = "sub-01", scans = ds$datamat, TR = 2, run_length = c(40L, 40L),
         events = ds$event_table)))
  expect_s3_class(m, "fmri_manifest")
  expect_length(m, 1)
  expect_error(as_manifest(list(list(scans = 1))), "no 'id'")
})
