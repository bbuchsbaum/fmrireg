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
  res <- project_confounds(Y, list(X), Z)
  expect_equal(dim(res$X_list[[1]]), dim(X))
  expect_equal(dim(res$Y), dim(Y))
})

test_that("create_fmri_design returns expected structure", {
  sf <- sampling_frame(20, TR = 1)
  events <- data.frame(onset = c(2, 6, 12),
                       condition = factor(c("A", "B", "A")),
                       block = 1)
  emod <- event_model(onset ~ hrf(condition), data = events,
                      block = ~ block, sampling_frame = sf)
  des <- create_fmri_design(emod, HRF_SPMG1)
  expect_type(des, "list")
  expect_equal(length(des$X_list), 2)
  expect_equal(des$d, nbasis(HRF_SPMG1))
  expect_equal(des$k, length(des$X_list))
  expect_true(is.matrix(des$Phi))
  expect_true(is.numeric(des$h_ref_shape_norm))
})
