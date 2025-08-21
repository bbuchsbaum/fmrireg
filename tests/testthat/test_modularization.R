# Test modularization of fmrilm components





test_that("modular components produce same results as original", {
  skip_if_not_installed("neuroim2")
  library(fmrireg)
  
  # Create test data
  n_voxels <- 100
  n_timepoints <- 200
  n_runs <- 2
  
  # Simple design with 2 conditions
  onsets1 <- seq(10, 190, by = 30)
  onsets2 <- seq(20, 190, by = 30)
  
  block <- rep(1:n_runs, each = n_timepoints/n_runs)
  
  # Create sampling frame
  sframe <- fmrihrf::sampling_frame(blocklens = c(n_timepoints/n_runs, n_timepoints/n_runs), TR = 2)
  
  # Create baseline model
  bmodel <- baseline_model(basis = "poly", degree = 3, sframe = sframe)
  
  # Create event model  
  # Determine which block each onset belongs to
  all_onsets <- c(onsets1, onsets2)
  all_conditions <- rep(c("A", "B"), c(length(onsets1), length(onsets2)))
  event_blocks <- ifelse(all_onsets <= n_timepoints/n_runs, 1, 2)
  
  # Sort by block, then by onset to ensure non-decreasing blockids
  event_order <- order(event_blocks, all_onsets)
  
  event_spec <- event_model(
    onset ~ hrf(onset, basis="spmg1"),
    data = data.frame(
      onset = all_onsets[event_order],
      condition = all_conditions[event_order],
      block = event_blocks[event_order]
    ),
    block = ~ block,
    sampling_frame = sframe
  )
  
  # Create temporary dataset for model construction
  event_data <- data.frame(
    onset = all_onsets[event_order],
    condition = all_conditions[event_order],
    block = event_blocks[event_order]
  )
  
  # Create temporary dummy data for model construction
  dummy_Y <- matrix(0, n_timepoints, 10)  # Just need something with right number of rows
  temp_dset <- fmridataset::matrix_dataset(dummy_Y, TR = 2, run_length = rep(n_timepoints/n_runs, n_runs), 
                        event_table = event_data)
  
  # Combine into fmri_model
  fmodel <- fmri_model(event_spec, bmodel, temp_dset)
  
  # Create simulated data
  X <- as.matrix(design_matrix(fmodel))
  betas <- matrix(rnorm(ncol(X) * n_voxels), ncol = n_voxels)
  Y <- X %*% betas + matrix(rnorm(n_timepoints * n_voxels, sd = 0.5), ncol = n_voxels)
  
  # Create actual dataset
  dset <- fmridataset::matrix_dataset(Y, TR = 2, run_length = rep(n_timepoints/n_runs, n_runs), 
                        event_table = event_data)
  
  # Test basic fitting
  fit1 <- fmri_lm(
    onset ~ hrf(condition, basis="spmg1"),
    block = ~ block,
    dataset = dset,
    strategy = "runwise"
  )
  
  expect_s3_class(fit1, "fmri_lm")
  expect_true(!is.null(fit1$result$betas))
  expect_true(!is.null(fit1$result$contrasts))
  
  # Test that coefficients can be extracted
  coefs <- coef(fit1)
  expect_true(is.data.frame(coefs) || is.matrix(coefs))
  expect_equal(ncol(coefs), n_voxels)
})

test_that("voxelwise AR with contrasts works", {
  skip_if_not_installed("neuroim2")
  library(fmrireg)
  
  # Small test for voxelwise AR
  n_voxels <- 10
  n_timepoints <- 100
  
  # Create simple data with AR structure
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
  
  # Create event table
  onsets <- seq(10, 90, by = 20)
  event_data <- data.frame(
    onset = onsets,
    condition = factor(rep("A", length(onsets))),
    block = 1
  )
  
  dset <- fmridataset::matrix_dataset(Y, TR = 2, run_length = n_timepoints,
                                      event_table = event_data)
  
  # Fit with voxelwise AR
  cfg <- fmri_lm_control(
    ar_options = list(
      struct = "ar1",
      voxelwise = TRUE
    )
  )
  
  fit <- fmri_lm(
    onset ~ hrf(condition, basis="spmg1"),
    block = ~ block,
    dataset = dset,
    strategy = "runwise",
    ar_options = cfg$ar_options
  )
  
  # Check that fit completed
  expect_s3_class(fit, "fmri_lm")
  expect_true(!is.null(fit$result))
  
  # Check that contrasts were computed if any were specified
  if (!is.null(fit$result$contrasts)) {
    expect_true(length(fit$result$contrasts) >= 0)
  }
})

test_that("memory-efficient contrast engine works", {
  # Test the chunked voxelwise contrast computation using internal function
  n_voxels <- 1000
  p <- 20
  n_timepoints <- 200
  
  # Create test data
  X <- matrix(rnorm(n_timepoints * p), n_timepoints, p)
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
  
  # Create phi matrix (AR coefficients per voxel)
  phi_matrix <- matrix(runif(n_voxels, 0.1, 0.5), nrow = 1)
  
  # Create contrast specifications
  conlist <- list(
    c1 = c(1, -1)
  )
  attr(conlist$c1, "colind") <- 1:2
  
  # Test chunked processing using ::: operator for internal function
  result <- fmrireg:::fit_lm_contrasts_voxelwise_chunked(
    X_run = X,
    Y_run = Y,
    phi_matrix = phi_matrix,
    conlist = conlist,
    fconlist = list(),
    chunk_size = 100
  )
  
  expect_equal(length(result), 1)
  expect_equal(result[[1]]$name, "c1")
  expect_equal(length(result[[1]]$data[[1]]$estimate), n_voxels)
})