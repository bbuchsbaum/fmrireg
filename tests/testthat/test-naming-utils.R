context("Naming Utilities")

library(fmrireg)
source("helper-naming.R") # Load is_valid_heading

test_that("zeropad works correctly with min width 2", {
  expect_equal(fmrireg:::zeropad(3, 9), "03")  # Now pads to 2 digits even for n_total < 10
  expect_equal(fmrireg:::zeropad(9, 9), "09") 
  expect_equal(fmrireg:::zeropad(3, 10), "03") # Remains 2 digits for n_total = 10
  expect_equal(fmrireg:::zeropad(12, 99), "12")
  expect_equal(fmrireg:::zeropad(12, 100), "012") # Becomes 3 digits when n_total >= 100
  expect_equal(fmrireg:::zeropad(99, 100), "099")
  expect_equal(fmrireg:::zeropad(1:3, 10), c("01", "02", "03")) # Check vector input
  expect_equal(fmrireg:::zeropad(5, 0), "5") # handle n_total=0, pads to width 1 (not 2)
})

test_that("sanitize works correctly", {
  expect_equal(sanitize("a.b c"), "a.b.c")
  expect_equal(sanitize("a b.c", allow_dot = FALSE), "a_b_c")
  expect_equal(sanitize("1var"), "X1var")
  expect_equal(sanitize(c("a", "a")), c("a", "a")) # unique = FALSE
  expect_equal(sanitize("_start"), "X_start")
})

test_that("basis_suffix works correctly", {
  expect_equal(basis_suffix(1, 1), "_b1")  # Always adds suffix, even for single basis
  expect_equal(basis_suffix(1:3, 3), c("_b01", "_b02", "_b03")) # Now expects 2 digits
  expect_equal(basis_suffix(1:3, 15), c("_b01", "_b02", "_b03")) # Already 2 digits
  expect_equal(basis_suffix(1, 5), "_b01") # Now expects 2 digits
  expect_equal(basis_suffix(integer(0), 5), "_b") # Empty input returns "_b" due to paste0 behavior
})

test_that("make_unique_tags works correctly", {
  expect_equal(fmrireg:::make_unique_tags(c("a", "b")), c("a", "b"))
  expect_equal(fmrireg:::make_unique_tags(c("a", "a")), c("a", "a#1"))
  expect_equal(fmrireg:::make_unique_tags(c("a", "a", "a#1")), c("a", "a#2", "a#1"))
})

# Basic test structure for other helpers - expand later

test_that("make_term_tag generates tags", {
  # Requires hrfspec structure - placeholder
  # spec1 <- list(id = "my_term", vars = list())
  # spec2 <- list(vars = list(quote(RT), quote(cond)))
  # expect_equal(make_term_tag(spec1), "my_term")
  # expect_equal(make_term_tag(spec2), "RT_cond")
  # expect_equal(make_term_tag(spec2, existing_tags = "RT_cond"), "RT_cond#1")
  # expect_equal(make_term_tag(list(vars=list(quote(X.y))), existing_tags="X_y"), "X_y#1")
  expect_true(TRUE) # Placeholder
})

test_that("level_token creates Var.Level format", {
  expect_equal(fmrireg:::level_token("cond", "A"), "cond.A")
  expect_equal(fmrireg:::level_token("cond name", "Level 1"), "cond.name.Level.1")
  expect_equal(fmrireg:::level_token("Input", "20"), "Input.20")
  expect_equal(fmrireg:::level_token("Input", c("1", "2")), c("Input.1", "Input.2"))

})

test_that("continuous_token sanitizes", {
  expect_equal(fmrireg:::continuous_token("poly_RT_01"), "poly_RT_01")
  expect_equal(fmrireg:::continuous_token("z_RT by cond"), "z_RT.by.cond")
})

test_that("make_cond_tag combines with underscore", {
  expect_equal(fmrireg:::make_cond_tag(c("cond.A", "task.go")), "cond.A_task.go")
  expect_equal(fmrireg:::make_cond_tag("cond.A"), "cond.A")
})

test_that("add_basis expands tags correctly", {
  expect_equal(fmrireg:::add_basis("cond.A", 1), "cond.A")
  expect_equal(fmrireg:::add_basis("cond.A", 3), c("cond.A_b01", "cond.A_b02", "cond.A_b03"))
  expect_equal(fmrireg:::add_basis(c("t1", "t2"), 2), c("t1_b01", "t2_b01", "t1_b02", "t2_b02"))
})

test_that("make_column_names composes final names", {
  expect_equal(fmrireg:::make_column_names("term1", "cond.A", 1), "term1_cond.A")
  expect_equal(fmrireg:::make_column_names("term1", "cond.A", 3), c("term1_cond.A_b01", "term1_cond.A_b02", "term1_cond.A_b03"))
  expect_equal(fmrireg:::make_column_names("term1", c("c1", "c2"), 2), c("term1_c1_b01", "term1_c2_b01", "term1_c1_b02", "term1_c2_b02"))
  expect_error(fmrireg:::make_column_names("term__bad", "cond.A", 1)) # Double underscore guard
})

test_that("is_valid_heading works", {
    expect_true(is_valid_heading("term_cond.A"))
    expect_true(is_valid_heading("term_cond.A_b01"))
    expect_true(is_valid_heading("term#1_cond.A_b01"))
    expect_true(is_valid_heading(".internal"))
    expect_false(is_valid_heading("_term"))
    expect_true(is_valid_heading("term_"))
    expect_false(is_valid_heading("term name"))
    expect_true(is_valid_heading("term_cond.A_b##"))
}) 