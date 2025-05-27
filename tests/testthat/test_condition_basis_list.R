library(testthat)
library(fmrireg)

testthat::local_edition(3)

create_cbl_data <- function() {
  event_data <- data.frame(
    onset = c(0, 30, 60),
    cond  = factor(c("A", "B", "A")),
    run   = 1
  )
  sframe <- sampling_frame(blocklens = 90, TR = 1)
  emod <- event_model(onset ~ hrf(cond), data = event_data,
                      block = ~run, sampling_frame = sframe)
  list(term = emod$terms[[1]], sframe = sframe)
}

test_that("condition_basis_list splits matrix by condition", {
  dat <- create_cbl_data()
  term <- dat$term
  sf   <- dat$sframe

  res <- condition_basis_list(term, HRF_SPMG2, sf)

  expect_type(res, "list")
  expect_equal(length(res), 2)
  expect_named(res, c("cond.A", "cond.B"))

  nb <- nbasis(HRF_SPMG2)
  expect_true(all(vapply(res, ncol, integer(1)) == nb))
  expect_true(all(vapply(res, nrow, integer(1)) == sum(blocklens(sf))))
})

test_that("condition_basis_list matrix output matches convolve", {
  dat <- create_cbl_data()
  term <- dat$term
  sf   <- dat$sframe

  mat1 <- condition_basis_list(term, HRF_SPMG2, sf, output = "matrix")
  mat2 <- fmrireg::convolve.event_term(term, hrf = HRF_SPMG2,
                                       sampling_frame = sf)

  expect_equal(mat1, mat2)
})
