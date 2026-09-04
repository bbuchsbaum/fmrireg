# CLI benchmark metadata/help leftovers + prepare_fmri_lm_contrasts with real contrasts.

capture_cli <- function(args) {
  stdout <- character()
  stderr <- character()
  stdout_con <- textConnection("stdout", "w", local = TRUE)
  stderr_con <- textConnection("stderr", "w", local = TRUE)
  sink(stdout_con)
  sink(stderr_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(stdout_con)
    close(stderr_con)
  }, add = TRUE)
  status <- fmrireg_cli(args)
  list(status = status, stdout = paste(stdout, collapse = "\n"),
       stderr = paste(stderr, collapse = "\n"))
}

test_that("benchmark metadata/help and unknown subcommand CLI paths", {
  out_help <- capture_cli(c("benchmark", "--help"))
  expect_equal(out_help$status, 0L)
  expect_true(grepl("benchmark|list|summary|metadata", out_help$stdout, ignore.case = TRUE) ||
                grepl("benchmark|list|summary|metadata", out_help$stderr, ignore.case = TRUE))

  out_meta <- capture_cli(c("benchmark", "metadata"))
  expect_equal(out_meta$status, 0L)

  out_meta_json <- capture_cli(c("benchmark", "metadata", "--json"))
  expect_equal(out_meta_json$status, 0L)

  out_list_help <- capture_cli(c("benchmark", "list", "--help"))
  expect_equal(out_list_help$status, 0L)

  out_sum_help <- capture_cli(c("benchmark", "summary", "--help"))
  expect_equal(out_sum_help$status, 0L)

  out_meta_help <- capture_cli(c("benchmark", "metadata", "--help"))
  expect_equal(out_meta_help$status, 0L)

  bad <- capture_cli(c("benchmark", "not-a-subcommand"))
  expect_equal(bad$status, 2L)
})

test_that("prepare_fmri_lm_contrasts resolves model contrasts with col_indices", {
  set.seed(71)
  etab <- data.frame(
    onset = c(5, 15, 25, 35, 45, 55),
    condition = factor(rep(c("A", "B"), 3)),
    run = 1L
  )
  Y <- matrix(rnorm(70 * 3), 70, 3)
  dset <- matrix_dataset(Y, TR = 1, run_length = 70, event_table = etab)
  con <- contrast_set(pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B"))
  emod <- event_model(
    onset ~ hrf(condition, contrasts = con),
    data = etab, block = ~ run, sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = dset$sampling_frame)
  fmod <- fmri_model(emod, bmod, dset)

  prepared <- fmrireg:::prepare_fmri_lm_contrasts(fmod)
  expect_true(is.list(prepared))
  expect_true(length(prepared) >= 1L)

  # Strip col_indices -> error when contrasts exist
  dm <- fmod$event_model$design_matrix
  attr(dm, "col_indices") <- NULL
  fmod2 <- fmod
  fmod2$event_model$design_matrix <- dm
  cw <- tryCatch(contrast_weights(fmod2$event_model), error = function(e) list())
  if (length(cw) > 0L) {
    expect_error(fmrireg:::prepare_fmri_lm_contrasts(fmod2), "column-index metadata")
  }
})

test_that("fmri_lm_fit strategy dispatch and print with contrasts", {
  set.seed(72)
  etab <- data.frame(
    onset = c(5, 15, 25, 35),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  Y <- matrix(rnorm(60 * 2), 60, 2)
  dset <- matrix_dataset(Y, TR = 1, run_length = 60, event_table = etab)
  con <- contrast_set(pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B"))
  fit <- fmri_lm(
    onset ~ hrf(condition, contrasts = con),
    block = ~ run, dataset = dset,
    strategy = "chunkwise", nchunks = 2L
  )
  expect_s3_class(fit, "fmri_lm")
  out <- capture.output(print(fit))
  expect_true(any(grepl("Contrast|A_vs_B|Design|formula", out, ignore.case = TRUE)))

  # coef_names with contrasts present
  cn <- coef_names(fit)
  expect_true(is.character(cn) && length(cn) >= 1L)
})
