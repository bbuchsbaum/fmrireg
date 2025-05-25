library(testthat)

# Ensure we use testthat edition 3
local_edition(3)

# Test that parse_event_formula returns hrfspec list intact

test_that("parse_event_formula returns hrfspec objects for multiple hrf calls", {
  events <- data.frame(
    onset = 1:4,
    condition = factor(c("a", "b", "a", "b")),
    modulator = rnorm(4)
  )

  form <- onset ~ hrf(condition) + hrf(modulator)

  parsed <- fmrireg:::parse_event_formula(form, events)

  expect_equal(nrow(parsed$spec_tbl), 2)
  expect_true(all(vapply(parsed$spec_tbl$spec, inherits, logical(1), "hrfspec")))
})
