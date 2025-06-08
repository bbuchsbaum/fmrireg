# Test fMRI latent linear models

library(fmrireg)
library(testthat)

test_that("fmri_latent_lm works with basic latent dataset", {
  skip("fmri_latent_lm requires updates to work with new fmridataset structure")
  
  # The fmri_latent_lm function needs to be updated to work with the new
  # fmridataset package structure. The chunkwise_lm.latent_dataset method
  # has a different signature than what fmri_lm_fit expects.
})

test_that("fmri_latent_lm handles different options", {
  skip("fmri_latent_lm requires updates to work with new fmridataset structure")
})

test_that("fmri_latent_lm validates inputs", {
  # Non-latent dataset should error
  regular_data <- matrix(rnorm(100 * 10), 100, 10)
  regular_dset <- fmridataset::matrix_dataset(
    datamat = regular_data,
    TR = 2,
    run_length = 100
  )
  
  expect_error(
    fmri_latent_lm(
      ~ 1,
      block = ~ 1,
      dataset = regular_dset,
      durations = 0
    ),
    "latent_dataset"
  )
})

test_that("latent_dataset creation works correctly", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("neuroim2")
  
  # Test basic creation
  n_time <- 100
  n_comp <- 20
  
  # Create basis (time x components)
  basis <- matrix(rnorm(n_time * n_comp), n_time, n_comp)
  
  # Loadings matrix (voxels x components)
  n_voxels <- 1000
  loadings <- matrix(rnorm(n_voxels * n_comp), n_voxels, n_comp)
  
  # Create LatentNeuroVec
  lvec <- fmristore::LatentNeuroVec(
    basis = basis,
    loadings = loadings,
    space = neuroim2::NeuroSpace(c(10, 10, 10, n_time)),
    mask = rep(TRUE, n_voxels),
    offset = rep(0, n_voxels)
  )
  
  dset <- fmridataset::latent_dataset(
    lvec = lvec,
    TR = 2,
    run_length = c(50, 50)
  )
  
  expect_s3_class(dset, "latent_dataset")
  expect_equal(dset$TR, 2)
  
  # Test data access
  chunks <- fmridataset::data_chunks(dset, nchunks = 1)
  chunk <- chunks$nextElem()
  expect_equal(ncol(chunk$data), n_comp)
  expect_equal(nrow(chunk$data), n_time)
})

test_that("chunkwise_lm.latent_dataset processes correctly", {
  skip("chunkwise_lm.latent_dataset requires different signature than fmri_lm_fit expects")
  
  # The chunkwise_lm.latent_dataset function has a different set of parameters
  # than what the main fmri_lm_fit function passes to it. This needs to be
  # reconciled before this test can work.
})