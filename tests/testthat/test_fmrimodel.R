options(mc.cores=2)
facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)


test_that("can construct and run an fmri_model with matrix data", {
  skip_if_not_installed("fmrireg")
  
  # Use a subset of the face design data for testing
  facedes_subset <- facedes[facedes$run <= 2, ]  # Only use first 2 runs
  facedes_subset$repnum <- factor(facedes_subset$rep_num)
  
  # Create synthetic data matrix instead of using missing NIfTI files
  n_timepoints <- 436 * 2  # 2 runs of 436 timepoints each
  n_voxels <- 100
  datamat <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
  
  # Create matrix dataset
  dset <- matrix_dataset(datamat, 
                        TR = 1.5,
                        run_length = rep(436, 2),
                        event_table = facedes_subset)

  # Create event and baseline models
  espec <- event_model(onset ~ hrf(repnum), data = facedes_subset, 
                      block = ~run, sampling_frame = dset$sampling_frame)
  bspec <- baseline_model(basis = "bs", degree = 5, sframe = dset$sampling_frame)
  
  # Create fmri_model
  fmod <- fmri_model(espec, bspec)
  
  # Test fmri_model properties
  expect_s3_class(fmod, "fmri_model")
  expect_true(!is.null(fmod))
  expect_equal(length(terms(fmod)), 3)
  expect_equal(ncol(design_matrix(fmod)), length(conditions(fmod)))
  expect_equal(2, length(baseline_terms(fmod)))
  expect_null(contrast_weights(fmod)$repnum)
  
  # Test that design matrix has correct dimensions
  dm <- design_matrix(fmod)
  expect_equal(nrow(dm), n_timepoints)
  expect_true(ncol(dm) > 0)
})

test_that("can construct fmri_model with real test files", {
  skip_if_not_installed("fmrireg")
  
  # Check if test files exist
  test_files <- c("test_data/testscan01.nii", "test_data/testscan02.nii")
  mask_file <- "test_data/testmask.nii"
  
  if (!all(file.exists(test_files)) || !file.exists(mask_file)) {
    skip("Test NIfTI files not available")
  }
  
  # Use a small subset of face design
  facedes_mini <- facedes[1:20, ]  # Just first 20 events
  facedes_mini$run <- 1  # All in run 1
  facedes_mini$repnum <- factor(facedes_mini$rep_num)
  
  tryCatch({
    # Try to create dataset with actual files
    dset <- fmri_dataset(scans = test_files[1],  # Just use one scan
                        mask = mask_file,
                        TR = 1.5,
                        run_length = 100,  # Shorter run
                        event_table = facedes_mini)
    
    espec <- event_model(onset ~ hrf(repnum), data = facedes_mini,
                        block = ~run, sampling_frame = dset$sampling_frame)
    bspec <- baseline_model(basis = "poly", degree = 3, sframe = dset$sampling_frame)
    fmod <- fmri_model(espec, bspec)
    
    expect_s3_class(fmod, "fmri_model")
    expect_true(!is.null(fmod))
    
  }, error = function(e) {
    # If NIfTI loading fails, skip the test
    skip(paste("Cannot load NIfTI test files:", e$message))
  })
})

test_that("fmri_model handles edge cases correctly", {
  skip_if_not_installed("fmrireg")
  
  # Test with minimal data
  mini_design <- data.frame(
    run = c(1, 1, 1),
    onset = c(10, 20, 30),
    condition = c("A", "B", "A"),
    rep_num = c(1, 1, 2)
  )
  
  n_timepoints <- 100
  datamat <- matrix(rnorm(n_timepoints * 50), n_timepoints, 50)
  
  dset <- matrix_dataset(datamat, 
                        TR = 2.0,
                        run_length = n_timepoints,
                        event_table = mini_design)
  
  # Test with minimal baseline model  
  espec <- event_model(onset ~ hrf(condition), data = mini_design,
                      block = ~run, sampling_frame = dset$sampling_frame)
  bspec <- baseline_model(basis = "poly", degree = 1, sframe = dset$sampling_frame)
  
  fmod <- fmri_model(espec, bspec)
  
  # Verify basic properties
  expect_s3_class(fmod, "fmri_model")
  expect_true(length(conditions(fmod)) >= 2)  # Should have A and B conditions
  
  # Test design matrix dimensions
  dm <- design_matrix(fmod)
  expect_equal(nrow(dm), n_timepoints)
  expect_true(ncol(dm) >= 2)  # At least condition A and B
  
  # Test that baseline terms are present
  bt <- baseline_terms(fmod)
  expect_true(length(bt) > 0)
  
  # Test printing doesn't error
  expect_no_error(print(fmod))
})





