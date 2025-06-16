# Simple test to debug write_results functionality

library(testthat)
library(fmrireg)

test_that("write_results basic debug test", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  
  # Create temporary directory for output
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  # Try to create the simplest possible real model by copying from the working test
  # but troubleshooting the data_chunks issue
  
  # Let's just test if we can load and run the basic infrastructure first
  cat("Testing basic fmri dataset creation...\n")
  
  # Use the exact same minimal structure from the working test
  set.seed(123)
  dims <- c(3, 3, 2)  # Very small
  n_timepoints <- 50
  
  scans <- lapply(1:2, function(run) {
    arr <- array(rnorm(prod(dims) * n_timepoints), c(dims, n_timepoints))
    bspace <- neuroim2::NeuroSpace(dim = c(dims, n_timepoints))
    neuroim2::NeuroVec(arr, bspace)
  })
  
  mask <- neuroim2::LogicalNeuroVol(
    array(TRUE, dims), 
    neuroim2::NeuroSpace(dim = dims)
  )
  
  event_table <- data.frame(
    onset = c(5, 15, 25, 35, 45,   # Run 1 events
              5, 15, 25, 35, 45),  # Run 2 events  
    condition = factor(rep(c("A", "B", "A", "B", "A"), 2)),
    run = rep(1:2, each = 5)
  )
  
  dset <- fmridataset::fmri_mem_dataset(
    scans = scans,
    mask = mask,
    TR = 1.5,
    event_table = event_table
  )
  
  cat("Dataset created successfully\n")
  
  # Now try to fit the model with different strategy
  cat("Trying to fit fmri_lm model...\n")
  
  con <- contrast_set(pair_contrast( ~ condition == "A", ~ condition == "B", name="A_vs_B"))
  
     # Try with simpler parameters
   mod <- fmri_lm(
     onset ~ hrf(condition, contrasts = con), 
     block = ~ run, 
     dataset = dset, 
     durations = 0,
     strategy = "runwise",  # Try runwise instead of chunkwise to avoid chunk issues
     progress = FALSE
   )
  
  cat("Model fitted successfully\n")
  cat("Model structure:\n")
  str(mod, max.level = 2)
  
  # Test basic export
  cat("Testing write_results...\n")
  
  # Debug the dimensions mismatch
  cat("=== DEBUGGING BETA/DESIGN MISMATCH ===\n")
  
  # Check beta data structure
  beta_data_raw <- mod$result$betas$data[[1]]$estimate[[1]]
  cat("Beta data dimensions:", dim(beta_data_raw), "\n")
  cat("Beta data colnames:", colnames(beta_data_raw), "\n")
  
  # Check design matrix
  design_mat <- design_matrix(mod$model)
  cat("Design matrix dimensions:", dim(design_mat), "\n")
  cat("Design matrix colnames:", colnames(design_mat), "\n")
  
  # Check event/baseline indices
  cat("Event indices:", mod$result$event_indices, "\n")
  cat("Baseline indices:", mod$result$baseline_indices, "\n")
  
  cat("=========================================\n")
  
  result <- write_results(
    mod,
    path = temp_dir,
    subject = "01",
    task = "test",
    space = "MNI152NLin2009cAsym",
    save_betas = TRUE
  )
  
  cat("write_results completed successfully!\n")
  cat("Result structure:\n")
  str(result)
  
  cat("Files in temp_dir:\n")
  print(list.files(temp_dir, recursive = TRUE))
  
  # Basic assertions
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
  
  cat("All tests passed!\n")
})