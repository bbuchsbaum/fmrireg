test_that("mixed_betas rejects malformed dimensions before fitting", {
  X <- cbind(intercept = 1, random = seq_len(5))

  expect_error(
    fmrireg:::mixed_betas(X, seq_len(4), ran_ind = 2L, fixed_ind = 1L),
    "same number of rows"
  )
})

test_that("mixed_betas makes solver failures explicit and non-plausible", {
  X <- cbind(intercept = 1, random = seq_len(5))
  Y <- 2 + 0.5 * seq_len(5)
  failing_solver <- function(Y, Z, X) stop("forced solver failure")

  expect_warning(
    out <- fmrireg:::mixed_betas(
      X,
      Y,
      ran_ind = 2L,
      fixed_ind = 1L,
      solver = failing_solver
    ),
    "forced solver failure"
  )

  expect_length(out, 2L)
  expect_true(all(is.na(out)))
  expect_false(any(out == 0, na.rm = TRUE))
})
