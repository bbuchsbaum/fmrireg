context("cfals design helpers")

test_that("reconstruction_matrix works", {
  sf <- sampling_frame(10, TR = 1)
  phi <- reconstruction_matrix(HRF_SPMG1, sf)
  expect_equal(ncol(phi), nbasis(HRF_SPMG1))
  expect_gt(nrow(phi), 1)
})

test_that("penalty_matrix defaults to identity", {
  Rm <- penalty_matrix(HRF_SPMG1)
  expect_equal(Rm, diag(nbasis(HRF_SPMG1)))
})

test_that("project_confounds projects via QR", {
  X <- matrix(rnorm(20), 5, 4)
  Y <- matrix(rnorm(10), 5, 2)
  Z <- matrix(seq_len(5), ncol = 1)
  res <- project_confounds(list(X), Y, Z)
  expect_equal(dim(res$X_list[[1]]), dim(X))
  expect_equal(dim(res$Y), dim(Y))
})
