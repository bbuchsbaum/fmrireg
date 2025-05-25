library(testthat)

context("schema helper functions")

make_errs <- function() ValidationErrors$new()

# check_required -----------------------------------------------------------

test_that("check_required detects missing fields", {
  errs <- make_errs()
  data <- list(a = 1)
  expect_true(fmrireg:::check_required(data, "a", "root", errs))
  expect_true(errs$is_valid())
  expect_false(fmrireg:::check_required(data, "b", "root", errs))
  expect_false(errs$is_valid())
  expect_true(any(vapply(errs$errors, function(e) grepl("root\\$b", e$path), logical(1))))
})

# check_type ---------------------------------------------------------------

test_that("check_type validates basic and array types", {
  errs <- make_errs()
  data <- list(
    s = "hi",
    i = 2L,
    n = 1.5,
    b = TRUE,
    arr = c("a", "b")
  )
  expect_true(fmrireg:::check_type(data, "s", "string", "p", errs))
  expect_true(fmrireg:::check_type(data, "i", "integer", "p", errs))
  expect_true(fmrireg:::check_type(data, "n", "number", "p", errs))
  expect_true(fmrireg:::check_type(data, "b", "boolean", "p", errs))
  expect_true(fmrireg:::check_type(data, "arr", "array[string]", "p", errs))

  errs2 <- make_errs()
  data2 <- list(x = "bad")
  expect_false(fmrireg:::check_type(data2, "x", "number", "p", errs2))
  expect_false(errs2$is_valid())
})

# check_enum ---------------------------------------------------------------

test_that("check_enum validates membership", {
  errs <- make_errs()
  data <- list(x = c("A", "B"))
  expect_true(fmrireg:::check_enum(data, "x", c("A", "B", "C"), "p", errs))

  errs2 <- make_errs()
  data2 <- list(x = "D")
  expect_false(fmrireg:::check_enum(data2, "x", c("A", "B"), "p", errs2))
  expect_false(errs2$is_valid())
})

# check_pattern ------------------------------------------------------------

test_that("check_pattern matches regex", {
  errs <- make_errs()
  data <- list(id = c("sub-01", "sub-02"))
  expect_true(fmrireg:::check_pattern(data, "id", "^sub-[0-9]+$", "p", errs))

  errs2 <- make_errs()
  data2 <- list(id = c("sub-1", "foo"))
  expect_false(fmrireg:::check_pattern(data2, "id", "^sub-[0-9]+$", "p", errs2))
  expect_false(errs2$is_valid())
})

