library(testthat)

context("ValidationErrors class")

test_that("validation error handling works", {
  errs <- ValidationErrors$new()
  expect_true(errs$is_valid())

  errs$add_error("a$b", "oops")
  expect_false(errs$is_valid())
  expect_true(grepl("oops", errs$format_errors()))
  expect_error(errs$stop_if_invalid("bad"), "bad")
})
