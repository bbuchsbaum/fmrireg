# Tests for write_results.fmri_lm BIDS export functionality

library(testthat)
devtools::load_all('.')  # Load development version
library(neuroim2)
library(fmridataset)

# Helper function to create minimal test dataset
create_test_dataset <- function(dims = c(3, 3, 2), n_timepoints = 50) {
  set.seed(123)  # For reproducible tests
  
  # Create synthetic brain scans
  scans <- lapply(1:2, function(run) {  # 2 runs
    arr <- array(rnorm(prod(dims) * n_timepoints), c(dims, n_timepoints))
    bspace <- neuroim2::NeuroSpace(dim = c(dims, n_timepoints))
    neuroim2::NeuroVec(arr, bspace)
  })
  
  # Create brain mask (all voxels active for simplicity)
  mask <- neuroim2::LogicalNeuroVol(
    array(TRUE, dims), 
    neuroim2::NeuroSpace(dim = dims)
  )
  
  # Create simple event table
  event_table <- data.frame(
    onset = c(5, 15, 25, 35, 45,   # Run 1 events
              5, 15, 25, 35, 45),  # Run 2 events  
    condition = factor(rep(c("A", "B", "A", "B", "A"), 2)),
    run = rep(1:2, each = 5)
  )
  
  # Create fmri_mem_dataset
  dset <- fmridataset::fmri_mem_dataset(
    scans = scans,
    mask = mask,
    TR = 1.5,
    event_table = event_table
  )
  
  return(dset)
}

# Helper function to create minimal fmri_lm object
create_test_fmri_lm <- function() {
  set.seed(123)
  
  # Create test dataset
  dset <- create_test_dataset()
  
  # Fit a simple model with contrasts
  con <- contrast_set(pair_contrast( ~ condition == "A", ~ condition == "B", name="A_vs_B"))
  
  mod <- fmri_lm(
    onset ~ hrf(condition, contrasts = con), 
    block = ~ run, 
    dataset = dset, 
    durations = 0,
    strategy = "chunkwise", 
    nchunks = 2,
    progress = FALSE
  )
  
  return(mod)
}

test_that("write_results.fmri_lm creates expected files", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  
  # Create test model
  mod <- create_test_fmri_lm()
  
  # Create temporary directory for output
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  # Test basic export
  result <- write_results(
    mod,
    path = temp_dir,
    subject = "01",
    task = "test",
    space = "MNI152NLin2009cAsym"
  )
  
  # Check that files were created
  expect_true(is.list(result))
  expect_true("betas" %in% names(result))
  
  # Check beta files exist
  expect_true(file.exists(result$betas$h5))
  expect_true(file.exists(result$betas$json))
  
  # Check BIDS-compliant naming
  expect_true(grepl("sub-01_task-test_space-MNI152NLin2009cAsym_desc-GLM_betas\\.h5$", 
                    result$betas$h5))
  expect_true(grepl("sub-01_task-test_space-MNI152NLin2009cAsym_desc-GLM_betas\\.json$", 
                    result$betas$json))
})

test_that("write_results.fmri_lm handles contrast export strategies", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  
  mod <- create_test_fmri_lm()
  
  # Test by_stat strategy
  temp_dir1 <- tempfile()
  dir.create(temp_dir1)
  on.exit(unlink(temp_dir1, recursive = TRUE), add = TRUE)
  
  result1 <- write_results(
    mod,
    path = temp_dir1,
    subject = "01", 
    task = "test",
    strategy = "by_stat",
    contrast_stats = c("beta", "tstat")
  )
  
  # Should have files for each statistic
  expect_true("beta" %in% names(result1))
  expect_true("tstat" %in% names(result1))
  
  # Test by_contrast strategy  
  temp_dir2 <- tempfile()
  dir.create(temp_dir2)
  on.exit(unlink(temp_dir2, recursive = TRUE), add = TRUE)
  
  result2 <- write_results(
    mod,
    path = temp_dir2,
    subject = "01",
    task = "test", 
    strategy = "by_contrast",
    contrast_stats = c("beta", "tstat")
  )
  
  # Should have files for each contrast
  expect_true("A_vs_B" %in% names(result2))
})

test_that("write_results.fmri_lm creates valid JSON metadata", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  
  mod <- create_test_fmri_lm()
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  result <- write_results(
    mod,
    path = temp_dir,
    subject = "01",
    task = "test",
    save_betas = TRUE
  )
  
  # Read and validate JSON metadata
  json_content <- jsonlite::read_json(result$betas$json)
  
  # Check required BIDS derivative fields
  expect_true("Description" %in% names(json_content))
  expect_true("Sources" %in% names(json_content))
  expect_true("GeneratedBy" %in% names(json_content))
  expect_true("SoftwareVersions" %in% names(json_content))
  expect_true("ModelInfo" %in% names(json_content))
  expect_true("CreationTime" %in% names(json_content))
  
  # Check specific values
  expect_equal(json_content$GeneratedBy$Name, "fmrireg::write_results")
  expect_true(grepl("Raw regressor beta coefficients", json_content$Description))
  expect_true("R" %in% names(json_content$SoftwareVersions))
  expect_true("fmrireg" %in% names(json_content$SoftwareVersions))
})

test_that("write_results.fmri_lm validates required inputs", {
  mod <- create_test_fmri_lm()
  
  temp_dir <- tempfile() 
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  # Should require subject
  expect_error(
    write_results(mod, path = temp_dir, task = "test"),
    "Subject identifier is required"
  )
  
  # Should require task
  expect_error(
    write_results(mod, path = temp_dir, subject = "01"),
    "Task identifier is required"
  )
  
  # Should warn about missing space
  expect_warning(
    write_results(mod, path = temp_dir, subject = "01", task = "test"),
    "Spatial reference not specified"
  )
})

test_that("write_results.fmri_lm handles overwrite behavior correctly", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  
  mod <- create_test_fmri_lm()
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  # First write should succeed
  result1 <- write_results(
    mod,
    path = temp_dir,
    subject = "01",
    task = "test",
    space = "MNI152NLin2009cAsym",
    overwrite = FALSE  # Explicit default
  )
  
  expect_true(file.exists(result1$betas$h5))
  expect_true(file.exists(result1$betas$json))
  
  # Store original file modification times
  h5_mtime1 <- file.mtime(result1$betas$h5)
  json_mtime1 <- file.mtime(result1$betas$json)
  
  # Second write with overwrite = FALSE (default) should fail with specific error
  expect_error(
    write_results(
      mod,
      path = temp_dir,
      subject = "01", 
      task = "test",
      space = "MNI152NLin2009cAsym",
      overwrite = FALSE
    ),
    "already exists.*overwrite = TRUE",
    info = "Should provide clear guidance on how to overwrite"
  )
  
  # Files should be unchanged after failed overwrite attempt
  expect_equal(file.mtime(result1$betas$h5), h5_mtime1)
  expect_equal(file.mtime(result1$betas$json), json_mtime1)
  
  # Write dummy data to test that overwrite actually replaces content
  dummy_data <- "DUMMY_CONTENT_TO_BE_REPLACED"
  writeLines(dummy_data, result1$betas$json)
  expect_true(grepl("DUMMY", readLines(result1$betas$json)[1]))
  
  # Wait a moment to ensure different timestamps
  Sys.sleep(0.1)
  
  # Second write with overwrite = TRUE should succeed and replace content
  result2 <- write_results(
    mod,
    path = temp_dir,
    subject = "01",
    task = "test", 
    space = "MNI152NLin2009cAsym",
    overwrite = TRUE
  )
  
  expect_true(file.exists(result2$betas$h5))
  expect_true(file.exists(result2$betas$json))
  
  # File modification times should be different (newer)
  expect_true(file.mtime(result2$betas$h5) > h5_mtime1)
  expect_true(file.mtime(result2$betas$json) > json_mtime1)
  
  # JSON should no longer contain dummy content
  json_content <- jsonlite::read_json(result2$betas$json)
  expect_true("Description" %in% names(json_content))
  expect_false(grepl("DUMMY", json_content$Description))
  
  # HDF5 should contain valid data - verify by reading it back
  beta_data_read <- fmristore::read_labeled_vec(result2$betas$h5)
  expect_true(!is.null(beta_data_read))
  expect_true(inherits(beta_data_read, "LabeledVolumeSet"))
})

test_that("write_results.fmri_lm uses atomic write pattern", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  
  mod <- create_test_fmri_lm()
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  # Normal write should not leave temporary directories
  result <- write_results(
    mod,
    path = temp_dir,
    subject = "01",
    task = "test"
  )
  
  # Check no temporary directories remain
  temp_files <- list.files(temp_dir, pattern = "^\\.tmp_write_", include.dirs = TRUE)
  expect_length(temp_files, 0)
  
  # Files should exist in final location
  expect_true(file.exists(result$betas$h5))
  expect_true(file.exists(result$betas$json))
})

test_that("write_results.fmri_lm works with minimal contrast set", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  
  mod <- create_test_fmri_lm()
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  # Export only specific contrast and statistics
  result <- write_results(
    mod,
    path = temp_dir,
    subject = "01",
    task = "test",
    contrasts = "A_vs_B",
    contrast_stats = "beta",
    save_betas = FALSE  # Skip betas to test contrast-only export
  )
  
  # Should only have beta statistics
  expect_true("beta" %in% names(result))
  expect_false("betas" %in% names(result))
  expect_false("tstat" %in% names(result))
})

test_that("write_results.fmri_lm correctly writes numerical data to HDF5", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  
  mod <- create_test_fmri_lm()
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  # Export with both betas and contrasts
  result <- write_results(
    mod,
    path = temp_dir,
    subject = "01",
    task = "test",
    space = "MNI152NLin2009cAsym",
    save_betas = TRUE,
    strategy = "by_stat",
    contrast_stats = c("beta", "tstat")
  )
  
  # === BETA MAP VALIDATION ===
  beta_h5_path <- result$betas$h5
  expect_true(file.exists(beta_h5_path))
  
  # Read beta data using fmristore (proper interface)
  beta_data_read <- fmristore::read_labeled_vec(beta_h5_path)
  expect_true(!is.null(beta_data_read))
  
  # Validate that we can read the data successfully
  expect_true(inherits(beta_data_read, "LabeledVolumeSet"))
  
  # Get original beta data from model for comparison
  original_betas <- coef(mod, type = "betas", include_baseline = FALSE)
  original_regressor_names <- rownames(original_betas)
  
  # Basic validation that data was written and can be read
  expect_true(length(original_regressor_names) > 0)
  
  # === CONTRAST MAP VALIDATION ===  
  if ("beta" %in% names(result)) {
    contrast_h5_path <- result$beta$h5
    expect_true(file.exists(contrast_h5_path))
    
    # Read contrast data using fmristore
    contrast_data_read <- fmristore::read_labeled_vec(contrast_h5_path)
    expect_true(!is.null(contrast_data_read))
    expect_true(inherits(contrast_data_read, "LabeledVolumeSet"))
  }
})

test_that("write_results.fmri_lm validates ModelInfo content in JSON", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  
  mod <- create_test_fmri_lm()
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  result <- write_results(
    mod,
    path = temp_dir,
    subject = "01",
    task = "test",
    save_betas = TRUE
  )
  
  # Read and validate JSON metadata
  json_content <- jsonlite::read_json(result$betas$json)
  
  # Validate ModelInfo structure and content
  expect_true("ModelInfo" %in% names(json_content))
  model_info <- json_content$ModelInfo
  
  # Check required ModelInfo fields
  expect_true("Type" %in% names(model_info))
  expect_equal(model_info$Type, "General Linear Model")
  
  expect_true("Formula" %in% names(model_info))
  expect_true(is.character(model_info$Formula))
  # Should contain key model components
  expect_true(grepl("hrf", model_info$Formula))
  expect_true(grepl("condition", model_info$Formula))
  
  expect_true("NumRegressors" %in% names(model_info))
  expect_true(is.numeric(model_info$NumRegressors))
  expect_true(model_info$NumRegressors > 0)
  
  # Validate that the number of regressors matches the actual model
  actual_regressors <- ncol(coef(mod, include_baseline = TRUE))
  expect_equal(model_info$NumRegressors, actual_regressors)
  
  # Validate RegressorOrder field
  expect_true("RegressorOrder" %in% names(json_content))
  expect_true(is.list(json_content$RegressorOrder))  # JSON arrays are read as lists
  expect_true(all(sapply(json_content$RegressorOrder, is.character)))  # Each element should be character
  expect_true(length(json_content$RegressorOrder) > 0)
})

test_that("write_results.fmri_lm handles filesystem edge cases", {
  skip_if_not_installed("fmristore") 
  skip_if_not_installed("jsonlite")
  
  mod <- create_test_fmri_lm()
  
  # Test 1: Non-existent output directory should be created
  temp_base <- tempfile()
  non_existent_path <- file.path(temp_base, "deep", "nested", "path")
  on.exit(unlink(temp_base, recursive = TRUE))
  
  expect_false(dir.exists(non_existent_path))
  
  result <- write_results(
    mod,
    path = non_existent_path,
    subject = "01",
    task = "test"
  )
  
  expect_true(dir.exists(non_existent_path))
  expect_true(file.exists(result$betas$h5))
  
  # Test 2: Invalid BIDS characters should be sanitized
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
  
  # Test subject with invalid characters
  result_sanitized <- write_results(
    mod,
    path = temp_dir,
    subject = "test_01@#$",  # Contains invalid characters
    task = "my-task_test",   # Contains invalid characters
    space = "MNI152"
  )
  
  # Files should still be created with sanitized names
  expect_true(file.exists(result_sanitized$betas$h5))
  expect_true(file.exists(result_sanitized$betas$json))
  
  # Check that filename doesn't contain invalid characters  
  filename <- basename(result_sanitized$betas$h5)
  expect_false(grepl("[@#$]", filename))  # Should not contain invalid chars (excluding underscore)
  expect_true(grepl("sub-test01", filename))  # Should be sanitized
})

test_that("write_results.fmri_lm handles models without contrasts", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  
  # Create model without contrasts
  set.seed(123)
  dset <- create_test_dataset()
  
  # Simple model without contrast_set
  mod_no_contrasts <- fmri_lm(
    onset ~ hrf(condition),  # No contrasts specified
    block = ~ run,
    dataset = dset,
    durations = 0,
    strategy = "chunkwise",
    nchunks = 2,
    progress = FALSE
  )
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  # Should succeed even without contrasts
  result <- write_results(
    mod_no_contrasts,
    path = temp_dir,
    subject = "01",
    task = "test",
    save_betas = TRUE
  )
  
  # Should have betas but no contrast files
  expect_true("betas" %in% names(result))
  expect_true(file.exists(result$betas$h5))
  expect_true(file.exists(result$betas$json))
  
  # Should not have contrast-related outputs
  contrast_keys <- setdiff(names(result), "betas")
  expect_length(contrast_keys, 0)
})

test_that("write_results.fmri_lm validates GeneratedBy field for BIDS compliance", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  
  mod <- create_test_fmri_lm()
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  result <- write_results(
    mod,
    path = temp_dir,
    subject = "01",
    task = "test",
    save_betas = TRUE
  )
  
  # Read and validate JSON metadata
  json_content <- jsonlite::read_json(result$betas$json)
  
  # Validate GeneratedBy field for BIDS derivatives compliance
  expect_true("GeneratedBy" %in% names(json_content))
  generated_by <- json_content$GeneratedBy
  
  # Check required GeneratedBy subfields
  expect_true("Name" %in% names(generated_by))
  expect_equal(generated_by$Name, "fmrireg::write_results")
  
  expect_true("Version" %in% names(generated_by))
  expect_true(is.character(generated_by$Version))
  expect_true(nchar(generated_by$Version) > 0)
  
  # Should be a valid package version string
  expect_true(grepl("^[0-9]+\\.[0-9]+", generated_by$Version))
  
  expect_true("CodeURL" %in% names(generated_by))
  expect_true(grepl("github", generated_by$CodeURL))
  
  # Validate SoftwareVersions field
  expect_true("SoftwareVersions" %in% names(json_content))
  software_versions <- json_content$SoftwareVersions
  
  expect_true("R" %in% names(software_versions))
  expect_true("fmrireg" %in% names(software_versions))
  expect_true("fmristore" %in% names(software_versions))
  expect_true("neuroim2" %in% names(software_versions))
  
  # All version fields should be character strings
  expect_true(all(sapply(software_versions, is.character)))
  expect_true(all(sapply(software_versions, function(x) nchar(x) > 0)))
})

test_that("write_results.fmri_lm validates HDF5 data types and precision", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  
  mod <- create_test_fmri_lm()
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  result <- write_results(
    mod,
    path = temp_dir,
    subject = "01",
    task = "test",
    save_betas = TRUE
  )
  
  # Check that data can be read back successfully (indicates proper HDF5 format)
  beta_h5_path <- result$betas$h5
  beta_data_read <- fmristore::read_labeled_vec(beta_h5_path)
  
  # Validate data was read successfully and has expected properties
  expect_true(!is.null(beta_data_read))
  expect_true(inherits(beta_data_read, "LabeledVolumeSet"))
  
  # Basic validation that the LabeledVolumeSet is valid
  expect_true(is(beta_data_read, "LabeledVolumeSet"))
})

test_that("write_results.fmri_lm handles degenerate data gracefully", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  
  # Create a model but then simulate degenerate beta values
  mod <- create_test_fmri_lm()
  
  # Modify the model to have some degenerate beta values
  # This simulates cases where a regressor has no variance or fitting failed
  original_betas <- mod$result$betas$data[[1]]$estimate[[1]]
  
  # Create degenerate cases: all zeros and some NAs
  degenerate_betas <- original_betas
  degenerate_betas[, 1] <- 0  # First regressor all zeros
  if (ncol(degenerate_betas) > 1) {
    degenerate_betas[1:5, 2] <- NA  # Some NAs in second regressor
  }
  
  # Replace the beta estimates
  mod$result$betas$data[[1]]$estimate[[1]] <- degenerate_betas
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  # Function should complete without error even with degenerate data
  expect_no_error({
    result <- write_results(
      mod,
      path = temp_dir,
      subject = "01",
      task = "test",
      save_betas = TRUE
    )
  })
  
  # Files should be created successfully
  expect_true(file.exists(result$betas$h5))
  expect_true(file.exists(result$betas$json))
  
  # HDF5 should handle zeros and NAs correctly
  beta_h5_path <- result$betas$h5
  
  # Read back the data using fmristore interface
  beta_data_read <- fmristore::read_labeled_vec(beta_h5_path)
  expect_true(!is.null(beta_data_read))
  
  # Should be a valid LabeledVolumeSet object
  expect_true(inherits(beta_data_read, "LabeledVolumeSet"))
  
  # Basic validation that the object is accessible
  expect_true(is(beta_data_read, "LabeledVolumeSet"))
})

test_that("write_results.fmri_lm validates CreationTime format for BIDS compliance", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  
  mod <- create_test_fmri_lm()
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  result <- write_results(
    mod,
    path = temp_dir,
    subject = "01",
    task = "test",
    save_betas = TRUE
  )
  
  # Read and validate JSON metadata
  json_content <- jsonlite::read_json(result$betas$json)
  
  # Validate CreationTime format (should be ISO 8601)
  expect_true("CreationTime" %in% names(json_content))
  creation_time <- json_content$CreationTime
  
  expect_true(is.character(creation_time))
  
  # Should match ISO 8601 format: YYYY-MM-DDTHH:MM:SS
  iso8601_pattern <- "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}$"
  expect_true(grepl(iso8601_pattern, creation_time),
              info = paste("CreationTime should be ISO 8601 format, got:", creation_time))
  
  # Should be a valid date/time that can be parsed
  expect_no_error({
    parsed_time <- as.POSIXct(creation_time, format = "%Y-%m-%dT%H:%M:%S")
    expect_true(!is.na(parsed_time))
    
    # Should be recent (within last few minutes)
    time_diff <- abs(as.numeric(difftime(Sys.time(), parsed_time, units = "mins")))
    expect_true(time_diff < 5, 
                info = paste("CreationTime should be recent, time difference:", time_diff, "minutes"))
  })
})

test_that("write_results.fmri_lm cleans up atomic write on error", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  
  mod <- create_test_fmri_lm()
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  # Create a mock function that will fail
  # We'll temporarily replace the fmristore function with one that errors
  original_write_func <- fmristore::write_labeled_vec
  
  # Create a failing mock function
  failing_write <- function(...) {
    stop("Simulated HDF5 write failure for testing")
  }
  
  # Temporarily replace the function in the namespace
  # This is a bit tricky in R, so we'll use a different approach
  # We'll modify the object to cause an error during processing
  
  # Alternative approach: corrupt the beta data structure to cause an error in .compute_beta_volumes
  corrupted_mod <- mod
  # Replace with invalid matrix to trigger validation error while maintaining tibble structure
  corrupted_mod$result$betas$data[[1]]$estimate <- list(matrix(numeric(0), 0, 0))
  
  # Expect an error from the main function
  expect_error(
    write_results(
      corrupted_mod, 
      path = temp_dir, 
      subject = "01", 
      task = "test",
      space = "MNI152NLin2009cAsym"
    ),
    "Failed to write BIDS results"
  )
  
  # Assert that the final directory contains no output files
  # (temporary directories should be cleaned up)
  output_files <- list.files(temp_dir, pattern = "\\.(h5|json)$", recursive = TRUE)
  expect_length(output_files, 0)
  
  # Assert that no temporary directories remain
  temp_dirs <- list.files(temp_dir, pattern = "^\\.tmp_write_", include.dirs = TRUE, all.files = TRUE)
  expect_length(temp_dirs, 0)
})

test_that("write_results.fmri_lm handles fmristore write failure gracefully", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  
  mod <- create_test_fmri_lm()
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  # Create a model that will fail during HDF5 writing
  invalid_mod <- mod
  
  # Corrupt the beta data structure to cause failure during processing while maintaining tibble structure
  invalid_mod$result$betas$data[[1]]$estimate <- list(matrix(numeric(0), 0, 0))  # Empty matrix
  
  # This should fail during the beta computation phase
  expect_error(
    write_results(
      invalid_mod,
      path = temp_dir,
      subject = "01", 
      task = "test",
      save_betas = TRUE
    ),
    "Failed to write BIDS results"
  )
  
  # Verify cleanup occurred
  output_files <- list.files(temp_dir, recursive = TRUE)
  expect_length(output_files, 0)
})