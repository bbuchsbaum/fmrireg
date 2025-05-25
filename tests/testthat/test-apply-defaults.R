library(testthat)

context("apply_defaults")

AD <- fmrireg:::apply_defaults

# basic merging

test_that("defaults are applied recursively", {
  data <- list(a = 1, b = list(x = 1))
  defaults <- list(a = 0, b = list(x = 0, y = 2), c = 3)
  res <- AD(data, defaults)
  expect_equal(res$a, 1)      # user value preserved
  expect_equal(res$b$x, 1)    # nested user value preserved
  expect_equal(res$b$y, 2)    # default added
  expect_equal(res$c, 3)      # default added
})

# null elements are left as-is

test_that("NULL entries are not overwritten", {
  data <- list(d = NULL)
  defaults <- list(d = list(e = 1))
  res <- AD(data, defaults)
  expect_null(res$d)
})
