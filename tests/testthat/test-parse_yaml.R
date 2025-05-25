library(testthat)

context("parse_yaml_to_list")

test_that("valid YAML is parsed correctly", {
  tf <- tempfile(fileext = ".yml")
  writeLines(c("a: 1", "b:", "  - 2", "  - 3"), tf)
  on.exit(unlink(tf))

  res <- parse_yaml_to_list(tf)
  expect_equal(res$a, 1)
  expect_equal(res$b, c(2,3))
})

test_that("missing file throws error", {
  expect_error(parse_yaml_to_list("nonexistent.yml"), "not found")
})

test_that("malformed YAML throws error", {
  tf <- tempfile(fileext = ".yml")
  writeLines("a: [1,", tf)
  on.exit(unlink(tf))
  expect_error(parse_yaml_to_list(tf), "Failed to parse YAML")
})
