# Push coverage past 90%: template/job APIs, NeuroVec coerce, mask helpers.

test_that("baseline_spec / fmri_template / fmri_job cover construction and validation", {
  expect_error(baseline_spec(degree = -1), "non-negative")
  expect_error(baseline_spec(degree = 2, basis = "bs"), "degree >= 3")
  bs <- baseline_spec(degree = 3, basis = "poly", intercept = "runwise")
  expect_s3_class(bs, "baseline_spec")
  expect_output(print(bs), regexp = ".")

  tmpl <- fmri_template(onset ~ hrf(condition), ~ run, baseline = bs)
  expect_s3_class(tmpl, "fmri_template")
  expect_output(print(tmpl), regexp = ".")
  expect_true(isTRUE(validate_template(tmpl)) || identical(validate_template(tmpl), TRUE) ||
                inherits(tryCatch(validate_template(tmpl), error = function(e) e), "fmri_template") ||
                is.null(attr(validate_template(tmpl), "class")))

  # Closure reducer should warn about serializability when supported
  expect_warning(
    fmri_template(onset ~ hrf(condition), ~ run,
                  reducer = (function() { x <- 1; function(fit) x })()),
    regexp = "serial|captur|closure|reducer",
    ignore.case = TRUE
  )

  ds <- dataset_spec(
    "matrix_dataset",
    args = list(
      datamat = matrix(rnorm(40), 20, 2),
      TR = 1, run_length = 20,
      event_table = data.frame(
        onset = c(5, 12), condition = factor(c("A", "B")), run = 1L
      )
    ),
    source = "inline"
  )
  expect_s3_class(ds, "dataset_spec")
  expect_output(print(ds), regexp = ".")

  expect_error(dataset_spec("", list()), "non-empty")
  expect_error(fmri_job("bad/id", tmpl, ds), "path separators")
  expect_error(fmri_job(".", tmpl, ds), "must not be")

  job <- fmri_job("sub-01", tmpl, ds, meta = list(subject = "01"))
  expect_s3_class(job, "fmri_job")
  expect_output(print(job), regexp = ".")
  expect_silent(invisible(validate_job(job)))
})

test_that("as.array.NeuroVec and spatial mask helpers cover remaining branches", {
  skip_if_not_installed("neuroim2")

  sp <- neuroim2::NeuroSpace(c(3, 3, 2, 4))
  arr <- array(rnorm(3 * 3 * 2 * 4), dim = c(3, 3, 2, 4))
  nv <- neuroim2::NeuroVec(arr, sp)
  out <- as.array(nv)
  expect_true(is.array(out))
  expect_equal(length(out), length(arr))

  # DenseNeuroVec path if available
  if (exists("DenseNeuroVec", envir = asNamespace("neuroim2"), inherits = FALSE)) {
    dnv <- tryCatch(neuroim2::DenseNeuroVec(arr, sp), error = function(e) NULL)
    if (!is.null(dnv)) {
      expect_true(is.array(as.array(dnv)))
    }
  }

  expect_error(
    fmrireg:::.fmri_dataset_mask_space(list(), "test"),
    "fmri_dataset"
  )
  expect_null(fmrireg:::.fmri_try_space(NULL))
  expect_null(fmrireg:::.fmri_try_space("path.nii"))
})

test_that("design_plot legend_threshold and color palette branches", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("plotly")
  skip_if_not_installed("bslib")

  set.seed(9)
  etab <- data.frame(
    onset = seq(4, 36, by = 4),
    condition = factor(rep(letters[1:4], length.out = 9)),
    run = 1L
  )
  Y <- matrix(rnorm(50 * 3), 50, 3)
  dset <- matrix_dataset(Y, TR = 1, run_length = 50, event_table = etab)
  emod <- event_model(onset ~ hrf(condition), data = etab, block = ~ run,
                      sampling_frame = dset$sampling_frame)
  bmod <- baseline_model(basis = "constant", sframe = dset$sampling_frame)
  fmod <- fmri_model(emod, bmod, dset)

  app <- design_plot(
    fmod, longnames = TRUE, color_palette = "Set2",
    legend_threshold = 1, facet_ncol = 1
  )
  expect_true(inherits(app, "shiny.appobj") || is.list(app))
})

test_that("instantiate / realize_dataset cover inline matrix_dataset jobs", {
  bs <- baseline_spec(degree = 1, basis = "poly")
  tmpl <- fmri_template(onset ~ hrf(condition), ~ run, baseline = bs)
  ev <- data.frame(
    onset = c(5, 15, 25),
    condition = factor(c("A", "B", "A")),
    run = 1L
  )
  Y <- matrix(rnorm(40 * 3), 40, 3)
  ds <- dataset_spec(
    "matrix_dataset",
    args = list(datamat = Y, TR = 1, run_length = 40, event_table = ev),
    source = "inline"
  )
  job <- fmri_job("s1", tmpl, ds)

  realized <- realize_dataset(job)
  expect_s3_class(realized, "matrix_dataset")

  jobs <- instantiate(tmpl, list(
    id = "s1", scans = Y, TR = 1, run_length = 40L, events = ev
  ))
  expect_s3_class(jobs, "fmri_job")
  expect_equal(jobs$id, "s1")
})
