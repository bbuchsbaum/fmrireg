# CLI option parser branches + benchmark summary writer.

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

test_that(".parse_cli_options covers -- separator, --no-, flags, and errors", {
  spec <- list(
    help = list(type = "flag", default = FALSE),
    json = list(type = "flag", default = FALSE),
    dataset = list(type = "value", required = TRUE),
    limit = list(type = "value", default = "10")
  )

  # -- separator consumes remaining as positionals
  parsed <- fmrireg:::.parse_cli_options(
    c("--dataset", "BM_X", "--", "extra"),
    spec, min_positionals = 0L, max_positionals = 1L
  )
  expect_equal(parsed$options$dataset, "BM_X")
  expect_equal(parsed$positionals, "extra")

  # flag with --no-
  parsed2 <- fmrireg:::.parse_cli_options(
    c("--dataset=BM_X", "--no-json"),
    spec
  )
  expect_false(parsed2$options$json)

  # key=value form
  parsed3 <- fmrireg:::.parse_cli_options(c("--dataset=BM_Y", "--limit=5"), spec)
  expect_equal(parsed3$options$limit, "5")

  expect_error(
    fmrireg:::.parse_cli_options(c("--dataset", "BM_X", "--unknown"), spec),
    class = "fmrireg_cli_usage_error"
  )
  expect_error(
    fmrireg:::.parse_cli_options(c("--json=true", "--dataset", "BM_X"), spec),
    class = "fmrireg_cli_usage_error"
  )
  expect_error(
    fmrireg:::.parse_cli_options(c("--no-dataset"), spec),
    class = "fmrireg_cli_usage_error"
  )
  expect_error(
    fmrireg:::.parse_cli_options(c("--dataset"), spec),
    class = "fmrireg_cli_usage_error"
  )
  expect_error(
    fmrireg:::.parse_cli_options(c("--dataset", "BM_X", "a", "b"), spec,
                                 max_positionals = 0L),
    class = "fmrireg_cli_usage_error"
  )
  expect_error(
    fmrireg:::.parse_cli_options(character(), spec),
    class = "fmrireg_cli_usage_error"
  )
})

test_that("benchmark summary and list CLI paths write expected fields", {
  out <- capture_cli(c("benchmark", "list"))
  expect_equal(out$status, 0L)
  expect_true(nzchar(out$stdout))

  # Pick a real dataset name from list output when possible
  datasets <- tryCatch(list_benchmark_datasets(), error = function(e) NULL)
  skip_if(is.null(datasets) || nrow(as.data.frame(datasets)) == 0L, "no benchmark datasets")
  ds_df <- as.data.frame(datasets)
  ds <- if ("Dataset" %in% names(ds_df)) {
    as.character(ds_df$Dataset[[1]])
  } else {
    rownames(ds_df)[[1]]
  }
  skip_if(!nzchar(ds), "no benchmark dataset name")

  sum_txt <- capture_cli(c("benchmark", "summary", "--dataset", ds))
  expect_equal(sum_txt$status, 0L)
  expect_match(sum_txt$stdout, "Dataset:")
  expect_match(sum_txt$stdout, "Timepoints:")

  sum_json <- capture_cli(c("benchmark", "summary", "--dataset", ds, "--json"))
  expect_equal(sum_json$status, 0L)
  parsed <- jsonlite::fromJSON(sum_json$stdout, simplifyVector = FALSE)
  expect_true(is.list(parsed))

  # Direct writer
  info <- get_benchmark_summary(ds)
  expect_output(
    fmrireg:::.write_benchmark_summary(info, ds),
    "Condition labels"
  )
})

test_that("fmrireg_cli top-level -h/--help and domain error path", {
  h1 <- capture_cli("--help")
  expect_equal(h1$status, 0L)
  expect_match(h1$stdout, "Usage:")

  h2 <- capture_cli("-h")
  expect_equal(h2$status, 0L)

  # Force generic error handler via unknown nested path that throws non-usage
  # (unknown command already covered); exercise install overwrite refusal
  bin_dir <- tempfile("fmrireg-cli-ow-")
  dir.create(bin_dir)
  invisible(install_cli(bin_dir, overwrite = TRUE))
  expect_error(install_cli(bin_dir, overwrite = FALSE), "Refusing to overwrite")
})
