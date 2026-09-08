# Coverage-oriented tests for CLI metadata/help branches and group_data accessors.

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
  list(
    status = status,
    stdout = paste(stdout, collapse = "\n"),
    stderr = paste(stderr, collapse = "\n")
  )
}

test_that("CLI covers empty args, help topics, and metadata subcommand", {
  empty <- capture_cli(character())
  expect_equal(empty$status, 0L)
  expect_match(empty$stdout, "Usage:")

  help_bench <- capture_cli(c("help", "benchmark"))
  expect_equal(help_bench$status, 0L)
  expect_match(help_bench$stdout, "benchmark")

  bad_help <- capture_cli(c("help", "not-a-topic"))
  expect_equal(bad_help$status, 2L)

  unknown_cmd <- capture_cli("totally-unknown")
  expect_equal(unknown_cmd$status, 2L)
  expect_match(unknown_cmd$stderr, "Unknown command")

  bench_help <- capture_cli(c("benchmark", "--help"))
  expect_equal(bench_help$status, 0L)

  list_help <- capture_cli(c("benchmark", "list", "--help"))
  expect_equal(list_help$status, 0L)

  summary_help <- capture_cli(c("benchmark", "summary", "--help"))
  expect_equal(summary_help$status, 0L)

  meta_help <- capture_cli(c("benchmark", "metadata", "--help"))
  expect_equal(meta_help$status, 0L)

  meta_text <- capture_cli(c("benchmark", "metadata"))
  expect_equal(meta_text$status, 0L)
  expect_true(nzchar(meta_text$stdout))

  meta_json <- capture_cli(c("benchmark", "metadata", "--json"))
  expect_equal(meta_json$status, 0L)
  parsed <- jsonlite::fromJSON(meta_json$stdout, simplifyVector = FALSE)
  expect_true(is.list(parsed))

  unknown_sub <- capture_cli(c("benchmark", "explode"))
  expect_equal(unknown_sub$status, 2L)
  expect_match(unknown_sub$stderr, "Unknown benchmark subcommand")
})

test_that("install_cli rejects unknown commands and PATH messaging", {
  expect_error(install_cli(tempfile(), commands = "not-real"), "Unknown command")

  bin_dir <- tempfile("fmrireg-cli-path-")
  # Directory not on PATH -> message about PATH
  expect_message(
    install_cli(bin_dir, overwrite = TRUE),
    "not currently on PATH"
  )
})

test_that("detect_group_data_format and legacy group_data accessors", {
  expect_equal(fmrireg:::detect_group_data_format("a.h5"), "h5")
  expect_equal(fmrireg:::detect_group_data_format("a.hdf5"), "h5")
  expect_equal(fmrireg:::detect_group_data_format("a.nii.gz"), "nifti")
  expect_equal(fmrireg:::detect_group_data_format("a.csv"), "csv")
  expect_equal(
    fmrireg:::detect_group_data_format(list(beta = "b.nii", se = "s.nii")),
    "nifti"
  )
  expect_equal(
    fmrireg:::detect_group_data_format(list(beta = "b.nii", var = "v.nii")),
    "nifti"
  )
  expect_equal(fmrireg:::detect_group_data_format(data.frame(x = 1)), "csv")
  # bare numeric vector is ambiguous -> error
  expect_error(
    fmrireg:::detect_group_data_format(1:3),
    "Could not automatically detect"
  )

  # fmri_lm list detection
  fx_etab <- data.frame(onset = c(5, 20), condition = factor(c("A", "B")), run = 1L)
  Y <- matrix(rnorm(40 * 3), 40, 3)
  dset <- matrix_dataset(Y, TR = 1, run_length = 40, event_table = fx_etab)
  fit <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset)
  expect_equal(fmrireg:::detect_group_data_format(list(fit)), "fmrilm")

  # Synthetic legacy group_data object for print/summary/accessors
  gd <- structure(
    list(
      format = "csv",
      subjects = c("s1", "s2", "s1"),
      covariates = data.frame(
        age = c(20, 30, 25),
        group = factor(c("A", "B", "A")),
        stringsAsFactors = FALSE
      ),
      mask = "/tmp/mask.nii.gz",
      dim = c(10L, 10L, 10L)
    ),
    class = c("group_data_csv", "group_data")
  )

  expect_equal(n_subjects(gd), 2L)
  expect_equal(get_subjects(gd), c("s1", "s2"))
  expect_equal(names(get_covariates(gd)), c("age", "group"))

  out <- capture.output(print(gd))
  expect_true(any(grepl("Group Data Object", out)))
  expect_true(any(grepl("csv", out)))

  sout <- capture.output(summary(gd))
  expect_true(any(grepl("Group Data Summary", sout)))
  expect_true(any(grepl("age", sout)))

  # Validation helpers
  expect_true(fmrireg:::validate_group_data(structure(list(), class = "group_data_gds")))
  expect_error(fmrireg:::validate_group_data(list()), "group_data")
  expect_error(
    fmrireg:::validate_group_data(structure(list(format = "csv"), class = "group_data")),
    "Missing required fields"
  )
  expect_error(
    fmrireg:::validate_group_data(
      structure(list(format = "csv", subjects = 1:2), class = "group_data")
    ),
    "subjects"
  )
})

test_that("t_to_beta_se converts t statistics for meta inputs", {
  out <- t_to_beta_se(t = c(2, -1.5), df = 18)
  expect_equal(names(out), c("beta", "se"))
  expect_equal(length(out$beta), 2)
  expect_true(all(out$se > 0))

  out2 <- t_to_beta_se(matrix(c(1, 2, 3, 4), 2, 2), df = 30, n = 32)
  expect_equal(dim(out2$beta), c(2, 2))
})

test_that("soft subspace options, apply, and extract_nuisance_timeseries", {
  set.seed(3)
  N <- matrix(rnorm(60 * 8), 60, 8)
  Y <- matrix(rnorm(60 * 10), 60, 10)
  X <- cbind(1, rnorm(60))

  opts <- soft_subspace_options(enabled = FALSE)
  expect_s3_class(opts, "soft_subspace_options")
  expect_false(opts$enabled)

  expect_error(
    soft_subspace_options(enabled = TRUE),
    "nuisance_mask or nuisance_matrix"
  )

  expect_warning(
    soft_subspace_options(
      enabled = TRUE,
      nuisance_mask = rep(TRUE, 5),
      nuisance_matrix = N,
      lambda = 0.2
    ),
    "using nuisance_matrix"
  )

  proj <- soft_projection(N, lambda = 0.5)
  cleaned <- apply_soft_projection(proj, Y, X)
  expect_equal(dim(cleaned$Y), dim(Y))
  expect_equal(dim(cleaned$X), dim(X))
  expect_error(apply_soft_projection(list(), Y, X), "soft_projection")

  dset <- matrix_dataset(Y, TR = 1, run_length = c(30, 30))
  mask <- c(rep(TRUE, 3), rep(FALSE, 7))
  nuis <- extract_nuisance_timeseries(dset, mask)
  expect_equal(dim(nuis), c(60L, 3L))

  nuis_run <- extract_nuisance_timeseries(dset, mask, run = 2)
  expect_equal(nrow(nuis_run), 30L)

  expect_error(
    extract_nuisance_timeseries(dset, rep(TRUE, 3)),
    "Mask length"
  )

  # array mask path
  arr_mask <- array(FALSE, dim = c(2, 5, 1))
  arr_mask[1, 1:2, 1] <- TRUE
  # length after vectorize = 10, matches ncol(Y)
  nuis2 <- extract_nuisance_timeseries(dset, arr_mask)
  expect_equal(ncol(nuis2), 2L)

  # Redundant nuisance warning when baseline already has nuisance_term
  soft_opts <- soft_subspace_options(
    enabled = TRUE,
    nuisance_matrix = N,
    warn_redundant = TRUE
  )
  bmod_nuis <- list(nuisance_term = list(dummy = TRUE))
  expect_warning(
    fmrireg:::.check_redundant_nuisance(bmod_nuis, soft_opts),
    "Soft subspace projection is enabled"
  )
  expect_null(
    fmrireg:::.check_redundant_nuisance(
      bmod_nuis,
      soft_subspace_options(enabled = FALSE)
    )
  )

  # Pipeline path with precomputed nuisance matrix
  pipe <- fmrireg:::.apply_soft_projection_pipeline(
    Y, X, dset,
    soft_subspace_options(enabled = TRUE, nuisance_matrix = N[, 1:5], lambda = 0.3)
  )
  expect_equal(dim(pipe$Y), dim(Y))
  expect_true(inherits(pipe$projection, "soft_projection"))

  pipe_off <- fmrireg:::.apply_soft_projection_pipeline(
    Y, X, dset, soft_subspace_options(enabled = FALSE)
  )
  expect_identical(pipe_off$Y, Y)
})
