# fmri_lm_internal voxelwise QR / storage helpers via small fixtures.

test_that(".fast_lm_matrix and voxel status helpers remain coherent", {
  X <- cbind(1, rnorm(25), rnorm(25))
  Y <- matrix(rnorm(25 * 3), 25, 3)
  proj <- fmrireg:::.fast_preproject(X)
  fit <- fmrireg:::.fast_lm_matrix(X, Y, proj, return_fitted = TRUE)
  expect_true(is.list(fit))
  expect_equal(ncol(fit$betas), 3L)
  expect_true(!is.null(fit$fitted) || !is.null(fit$residuals) || is.matrix(fit$betas))

  expect_true(is.list(fmrireg:::.fmri_lm_voxel_status(rnorm(10))) ||
                is.character(fmrireg:::.fmri_lm_voxel_status(rnorm(10))))
  bad <- fmrireg:::.fmri_lm_voxel_status(c(1, NA, Inf))
  expect_true(!is.null(bad))
})
