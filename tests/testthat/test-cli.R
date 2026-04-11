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

test_that("top-level CLI help succeeds", {
  result <- capture_cli("--help")

  expect_equal(result$status, 0L)
  expect_match(result$stdout, "Usage:")
  expect_match(result$stdout, "benchmark")
})

test_that("benchmark list supports text and json output", {
  text_result <- capture_cli(c("benchmark", "list"))
  json_result <- capture_cli(c("benchmark", "list", "--json"))

  expect_equal(text_result$status, 0L)
  expect_match(text_result$stdout, "Dataset")
  expect_match(text_result$stdout, "BM_Canonical_HighSNR")

  expect_equal(json_result$status, 0L)
  parsed <- jsonlite::fromJSON(json_result$stdout)
  expect_true(is.data.frame(parsed))
  expect_true("Dataset" %in% names(parsed))
})

test_that("benchmark summary enforces required arguments", {
  result <- capture_cli(c("benchmark", "summary"))

  expect_equal(result$status, 2L)
  expect_match(result$stderr, "Missing required option")
})

test_that("benchmark summary returns domain errors for invalid datasets", {
  result <- capture_cli(c("benchmark", "summary", "--dataset", "does-not-exist"))

  expect_equal(result$status, 1L)
  expect_match(result$stderr, "not found")
})

test_that("benchmark summary emits parseable json", {
  result <- capture_cli(c(
    "benchmark", "summary",
    "--dataset=BM_Canonical_HighSNR",
    "--json"
  ))

  expect_equal(result$status, 0L)
  parsed <- jsonlite::fromJSON(result$stdout, simplifyVector = FALSE)
  expect_equal(parsed$dimensions$n_voxels, 100)
  expect_equal(parsed$dimensions$n_conditions, 3)
})

test_that("unknown options and missing option values are usage errors", {
  unknown <- capture_cli(c("benchmark", "list", "--bogus"))
  missing <- capture_cli(c("benchmark", "summary", "--dataset"))

  expect_equal(unknown$status, 2L)
  expect_match(unknown$stderr, "Unknown option")
  expect_equal(missing$status, 2L)
  expect_match(missing$stderr, "Missing value for option '--dataset'")
})

test_that("install_cli copies wrapper and respects overwrite", {
  bin_dir <- tempfile("fmrireg-cli-")
  dir.create(bin_dir, recursive = TRUE)

  installed <- install_cli(bin_dir)
  target <- file.path(bin_dir, "fmrireg")

  expect_equal(length(installed), 1L)
  expect_true(file.exists(target))

  if (.Platform$OS.type != "windows") {
    expect_true(file.access(target, mode = 1) == 0)
  }

  expect_error(
    install_cli(bin_dir, overwrite = FALSE),
    "Refusing to overwrite"
  )

  expect_equal(length(install_cli(bin_dir, overwrite = TRUE)), 1L)
})

test_that("wrapper smoke test works through Rscript", {
  script <- file.path(getNamespaceInfo(asNamespace("fmrireg"), "path"), "exec", "fmrireg")
  skip_if_not(file.exists(script), "exec wrapper not available in source checkout")

  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(script, "benchmark", "list", "--json"),
    stdout = TRUE,
    stderr = TRUE
  )

  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }

  expect_equal(as.integer(status), 0L)
  start <- grep("^[[:space:]]*\\[", output)[1]
  end <- tail(grep("^[[:space:]]*\\]", output), 1)
  json_text <- paste(output[seq.int(start, end)], collapse = "\n")
  parsed <- jsonlite::fromJSON(json_text)
  expect_true(is.data.frame(parsed))
  expect_true(nrow(parsed) >= 1L)
})
