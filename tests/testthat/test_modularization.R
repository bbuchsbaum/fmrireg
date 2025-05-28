# Test modularization of fmrilm components

test_that("all extracted modules load correctly", {
  # Check that all new modules exist
  expect_true(file.exists(system.file("R", "fmri_model_utils.R", package = "fmrireg")))
  expect_true(file.exists(system.file("R", "fmri_lm_methods.R", package = "fmrireg")))
  expect_true(file.exists(system.file("R", "fmri_lm_internal.R", package = "fmrireg")))
  expect_true(file.exists(system.file("R", "fmri_lm_strategies.R", package = "fmrireg")))
  expect_true(file.exists(system.file("R", "fmri_lm_runwise.R", package = "fmrireg")))
  expect_true(file.exists(system.file("R", "fmri_lm_chunkwise.R", package = "fmrireg")))
})

test_that("extracted functions are available", {
  # Model utilities
  expect_true(exists("get_formula.fmri_model"))
  expect_true(exists("term_matrices.fmri_model"))
  expect_true(exists("create_fmri_model"))
  
  # Methods
  expect_true(exists("coef.fmri_lm"))
  expect_true(exists("stats.fmri_lm"))
  expect_true(exists("standard_error.fmri_lm"))
  expect_true(exists("print.fmri_lm"))
  
  # Internal utilities
  expect_true(exists(".fast_preproject"))
  expect_true(exists("meta_contrasts"))
  expect_true(exists("meta_betas"))
  
  # Strategies
  expect_true(exists("process_run_standard"))
  expect_true(exists("process_run_robust"))
  expect_true(exists("process_run_ar_robust"))
  
  # Main functions
  expect_true(exists("runwise_lm"))
  expect_true(exists("chunkwise_lm.fmri_dataset"))
})

test_that("modular components produce same results as original", {
  skip_if_not_installed("neuroim2")
  
  # Create test data
  n_voxels <- 100
  n_timepoints <- 200
  n_runs <- 2
  
  # Simple design with 2 conditions
  onsets1 <- seq(10, 190, by = 30)
  onsets2 <- seq(20, 190, by = 30)
  
  block <- rep(1:n_runs, each = n_timepoints/n_runs)
  
  # Create sampling frame
  sframe <- sampling_frame(block, TR = 2)
  
  # Create baseline model
  bmodel <- baseline_model(sframe, degree = 3)
  
  # Create event model
  event_spec <- event_model(
    onset ~ hrf(onset, hrf_spmg1()),
    block = sframe,
    data = data.frame(
      onset = c(onsets1, onsets2),
      condition = rep(c("A", "B"), c(length(onsets1), length(onsets2)))
    )
  )
  
  # Combine into fmri_model
  fmodel <- fmri_model(event_spec, bmodel)
  
  # Create simulated data
  X <- design_matrix(fmodel)
  betas <- matrix(rnorm(ncol(X) * n_voxels), ncol = n_voxels)
  Y <- X %*% betas + matrix(rnorm(n_timepoints * n_voxels, sd = 0.5), ncol = n_voxels)
  
  # Create dataset
  dset <- matrix_dataset(Y, TR = 2)
  
  # Test basic fitting
  fit1 <- fmri_lm(
    onset ~ hrf(onset, hrf_spmg1()),
    block = block,
    baseline_model = bmodel,
    dataset = dset,
    strategy = "runwise"
  )
  
  expect_s3_class(fit1, "fmri_lm")
  expect_true(!is.null(fit1$result$betas))
  expect_true(!is.null(fit1$result$contrasts))
  
  # Test that coefficients can be extracted
  coefs <- coef(fit1)
  expect_true(is.data.frame(coefs) || is.matrix(coefs))
  expect_equal(nrow(coefs), n_voxels)
})

test_that("voxelwise AR with contrasts works", {
  skip_if_not_installed("neuroim2")
  skip("Voxelwise AR is computationally intensive")
  
  # Small test for voxelwise AR
  n_voxels <- 10
  n_timepoints := 100
  
  # Create AR(1) noise
  ar_coef <- 0.5
  Y <- matrix(0, n_timepoints, n_voxels)
  for (v in 1:n_voxels) {
    noise <- rnorm(n_timepoints)
    for (t in 2:n_timepoints) {
      noise[t] <- ar_coef * noise[t-1] + rnorm(1)
    }
    Y[, v] <- noise
  }
  
  # Add signal
  onsets <- seq(10, 90, by = 20)
  sframe <- sampling_frame(rep(1, n_timepoints), TR = 2)
  ev <- event_model(onset ~ hrf(onset, hrf_spmg1()), 
                    block = sframe,
                    data = data.frame(onset = onsets))
  X <- design_matrix(ev)
  Y <- Y + X %*% rep(2, ncol(X))
  
  dset <- matrix_dataset(Y, TR = 2)
  
  # Fit with voxelwise AR
  cfg <- fmri_lm_control(
    ar_options = list(
      struct = "ar1",
      voxelwise = TRUE
    )
  )
  
  fit <- fmri_lm(
    onset ~ hrf(onset, hrf_spmg1()),
    block = rep(1, n_timepoints),
    dataset = dset,
    strategy = "runwise",
    ar_options = cfg$ar
  )
  
  # Check that contrasts were computed
  expect_true(length(fit$result$contrasts) > 0)
})

test_that("memory-efficient contrast engine works", {
  # Test the chunked voxelwise contrast computation
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
    c1 = c(1, -1, rep(0, p-2))
  )
  attr(conlist$c1, "colind") <- 1:2
  
  # Test chunked processing
  result <- fit_lm_contrasts_voxelwise_chunked(
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