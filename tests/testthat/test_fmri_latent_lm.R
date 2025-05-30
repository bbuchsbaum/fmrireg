# Test fMRI latent linear models

library(fmrireg)
library(testthat)

test_that("fmri_latent_lm works with basic latent dataset", {
  # Create test event data
  event_data <- data.frame(
    onsets = c(1, 15, 30, 45),
    condition = factor(c("A", "B", "A", "B")),
    run = c(1, 1, 1, 1)
  )
  
  # Create latent dataset
  # Simulating reduced data (e.g., from PCA)
  n_timepoints <- 50
  n_components <- 10
  latent_data <- matrix(rnorm(n_timepoints * n_components), n_timepoints, n_components)
  
  # Create mask and dataset
  mask <- array(1, dim = c(n_components, 1, 1))
  dset <- latent_dataset(
    latent_data,
    mask = mask,
    TR = 2,
    run_length = 50,
    event_table = event_data
  )
  
  # Fit latent model
  result <- fmri_latent_lm(
    onsets ~ hrf(condition),
    block = ~ run,
    dataset = dset,
    durations = 0,
    drop_empty = TRUE,
    robust = FALSE
  )
  
  expect_s3_class(result, "fmri_latent_lm")
  expect_s3_class(result, "fmri_lm")
  expect_equal(result$dataset, dset)
  
  # Check that results have expected structure
  expect_true(!is.null(result$result))
  expect_true(!is.null(result$result$betas))
})

test_that("fmri_latent_lm handles different options", {
  # Create simple test data
  event_data <- data.frame(
    onsets = seq(5, 45, by = 10),
    x = rnorm(5),
    run = rep(1, 5)
  )
  
  n_timepoints <- 50
  n_components <- 5
  latent_data <- matrix(rnorm(n_timepoints * n_components), n_timepoints, n_components)
  mask <- array(1, dim = c(n_components, 1, 1))
  
  dset <- latent_dataset(
    latent_data,
    mask = mask,
    TR = 2,
    run_length = 50,
    event_table = event_data
  )
  
  # Test with robust fitting
  result_robust <- fmri_latent_lm(
    onsets ~ hrf(x),
    block = ~ run,
    dataset = dset,
    durations = 0,
    robust = TRUE
  )
  
  expect_s3_class(result_robust, "fmri_latent_lm")
  
  # Test with baseline model
  baseline <- baseline_model(
    basis = "polynomial",
    degree = 2,
    sframe = sampling_frame(blocklens = 50, TR = 2)
  )
  
  result_baseline <- fmri_latent_lm(
    onsets ~ hrf(x),
    block = ~ run,
    baseline_model = baseline,
    dataset = dset,
    durations = 0
  )
  
  expect_s3_class(result_baseline, "fmri_latent_lm")
})

test_that("fmri_latent_lm validates inputs", {
  # Non-latent dataset should error
  regular_data <- matrix(rnorm(100 * 10), 100, 10)
  regular_dset <- matrix_dataset(
    regular_data,
    mask = array(1, dim = c(10, 1, 1)),
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
  # Test basic creation
  n_time <- 100
  n_comp <- 20
  latent_mat <- matrix(rnorm(n_time * n_comp), n_time, n_comp)
  
  # Loadings matrix (voxels x components)
  n_voxels <- 1000
  loadings <- matrix(rnorm(n_voxels * n_comp), n_voxels, n_comp)
  
  # Create 3D mask
  mask_3d <- array(0, dim = c(10, 10, 10))
  mask_3d[,,] <- 1
  
  dset <- latent_dataset(
    latent_data = latent_mat,
    loadings = loadings,
    mask = mask_3d,
    TR = 2,
    run_length = c(50, 50),
    preload = TRUE
  )
  
  expect_s3_class(dset, "latent_dataset")
  expect_equal(dim(dset), c(10, 10, 10, 100))
  expect_equal(dset@TR, 2)
  
  # Test data access
  chunk <- data_chunks(dset, nchunks = 1)[[1]]
  expect_equal(ncol(chunk$data), n_comp)
  expect_equal(nrow(chunk$data), n_time)
})

test_that("chunkwise_lm.latent_dataset processes correctly", {
  # Create test components
  n_time <- 60
  n_comp <- 8
  latent_data <- matrix(rnorm(n_time * n_comp), n_time, n_comp)
  
  event_data <- data.frame(
    onsets = c(5, 15, 25, 35, 45, 55),
    condition = factor(rep(c("A", "B"), 3)),
    run = rep(1, 6)
  )
  
  dset <- latent_dataset(
    latent_data,
    mask = array(1, dim = c(n_comp, 1, 1)),
    TR = 1,
    run_length = 60,
    event_table = event_data
  )
  
  # Create model
  sframe <- sampling_frame(blocklens = 60, TR = 1)
  baseline <- baseline_model(basis = "polynomial", degree = 1, sframe = sframe)
  
  model <- create_fmri_model(
    onsets ~ hrf(condition),
    block = ~ run,
    baseline_model = baseline,
    dataset = dset,
    drop_empty = TRUE
  )
  
  # Test chunking
  conlist <- list(
    AvsB = pair_contrast(~ condition == "A", ~ condition == "B", name = "AvsB")
  )
  
  result <- chunkwise_lm.latent_dataset(
    dset = dset,
    model = model,
    conlist = conlist,
    nchunks = 1,
    robust = FALSE,
    verbose = FALSE,
    autocor = "none",
    bootstrap = FALSE
  )
  
  expect_true(!is.null(result$bstats))
  expect_true(!is.null(result$contrasts))
  expect_equal(nrow(result$contrasts), 1)  # One contrast
})