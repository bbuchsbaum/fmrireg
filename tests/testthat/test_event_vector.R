context("Event Vector Tests")

# Helper function to capture print output
capture_print <- function(x) {
  out <- capture.output(print(x))
  paste(out, collapse = "\n")
}

test_that("event_factor creation and validation works", {
  # Test basic creation
  fac <- factor(c("A", "B", "A", "C", "B"))
  onsets <- c(0, 10, 20, 30, 40)
  ef <- event_factor(fac, "condition", onsets, rep(1, 5))
  
  expect_s3_class(ef, "event_factor")
  expect_equal(ef$varname, "condition")
  expect_equal(length(ef$value), 5)
  expect_equal(levels(ef$value), c("A", "B", "C"))
  
  # Test error conditions
  expect_error(event_factor(fac, "condition", onsets[1:3], rep(1, 5)),
               "length.*not equal")
  expect_error(event_factor(fac, "condition", onsets, c(1, 2)),
               "length.*not equal")
  
  # Test print method
  printed <- capture_print(ef)
  expect_match(printed, "Event Factor", fixed = TRUE)
  expect_match(printed, "condition", fixed = TRUE)
  expect_match(printed, paste(levels(fac), collapse = ", "), fixed = TRUE)
})

test_that("event_variable creation and validation works", {
  # Test basic creation
  vals <- rnorm(5)
  onsets <- seq(0, 40, length.out = 5)
  ev <- event_variable(vals, "continuous", onsets, rep(1, 5))
  
  expect_s3_class(ev, "event_variable")
  expect_equal(ev$varname, "continuous")
  expect_equal(length(ev$value), 5)
  expect_equal(ev$value, vals)
  
  # Test error conditions
  expect_error(event_variable(vals, "continuous", onsets[1:3], rep(1, 5)),
               "length.*not equal")
  expect_error(event_variable(vals, "continuous", onsets, c(1, 2)),
               "length.*not equal")
  
  # Test print method
  printed <- capture_print(ev)
  expect_match(printed, "Event Variable", fixed = TRUE)
  expect_match(printed, "continuous", fixed = TRUE)
  expect_match(printed, sprintf("%.2f - %.2f", min(vals), max(vals)))
})

test_that("event_matrix creation and validation works", {
  # Test basic creation
  mat <- matrix(rnorm(15), ncol = 3)
  colnames(mat) <- c("V1", "V2", "V3")
  onsets <- seq(0, 40, length.out = 5)
  em <- event_matrix(mat, "matrix_var", onsets, rep(1, 5))
  
  expect_s3_class(em, "event_matrix")
  expect_equal(em$varname, "matrix_var")
  expect_equal(dim(em$value), c(5, 3))
  expect_equal(em$value, mat)
  
  # Test error conditions
  expect_error(event_matrix(mat, "matrix_var", onsets[1:3], rep(1, 5)),
               "length.*not equal")
  expect_error(event_matrix(mat, "matrix_var", onsets, c(1, 2)),
               "length.*not equal")
  
  # Test print method
  printed <- capture_print(em)
  expect_match(printed, "Event Matrix", fixed = TRUE)
  expect_match(printed, "matrix_var", fixed = TRUE)
  expect_match(printed, paste(dim(mat), collapse = " × "), fixed = TRUE)
})

test_that("event_term creation and validation works", {
  # Test basic creation
  x1 <- factor(rep(letters[1:3], 2))
  x2 <- rnorm(6)
  onsets <- seq(0, 50, length.out = 6)
  eterm <- event_term(list(x1 = x1, x2 = x2), onsets, rep(1, 6))
  
  expect_s3_class(eterm, "event_term")
  expect_equal(length(eterm$events), 2)
  expect_s3_class(eterm$events$x1, "event_factor")
  expect_s3_class(eterm$events$x2, "event_variable")
  
  # Test subset functionality
  subset_term <- event_term(list(x1 = x1, x2 = x2), onsets, rep(1, 6), 
                           subset = x2 > 0)
  expect_true(all(subset_term$events$x2$value > 0))
  
  # Test error conditions
  expect_error(event_term(list(x1 = x1[1:3], x2 = x2), onsets, rep(1, 6)),
               "all event variables must have the same length")
  
  # Test print method
  printed <- capture_print(eterm)
  expect_match(printed, "Event Term", fixed = TRUE)
  expect_match(printed, "x1", fixed = TRUE)
  expect_match(printed, "x2", fixed = TRUE)
})

test_that("event sequence validation works", {
  # Test increasing blockids validation
  fac <- factor(c("A", "B", "A"))
  onsets <- c(0, 10, 20)
  expect_error(event_factor(fac, "condition", onsets, c(2, 1, 3)))
  
  # Test onset ordering within blocks
  expect_error(event_factor(fac, "condition", c(10, 0, 20), c(1, 1, 2)))
  
  # Test NA handling
  expect_error(event_factor(fac, "condition", c(0, NA, 20), c(1, 1, 2)))
})

test_that("EV factory function works correctly", {
  # Test factor input
  fac_input <- factor(c("A", "B", "A"))
  fac_ev <- EV(fac_input, "fac", c(0, 10, 20), rep(1, 3))
  expect_s3_class(fac_ev, "event_factor")
  
  # Test numeric input
  num_input <- c(1, 2, 3)
  num_ev <- EV(num_input, "num", c(0, 10, 20), rep(1, 3))
  expect_s3_class(num_ev, "event_variable")
  
  # Test matrix input
  mat_input <- matrix(1:6, ncol = 2)
  mat_ev <- EV(mat_input, "mat", c(0, 10, 20), rep(1, 3))
  expect_s3_class(mat_ev, "event_matrix")
  
  # Test subset functionality
  subset_ev <- EV(num_input, "num", c(0, 10, 20), rep(1, 3), 
                  subset = c(TRUE, FALSE, TRUE))
  expect_equal(length(subset_ev$value), 2)
})