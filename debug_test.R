# Debug test in development environment
devtools::load_all('.')

# Check if write_results exists
cat("write_results function exists:", exists("write_results"), "\n")

# Try to call it
tryCatch({
  write_results
  cat("write_results function found\n")
}, error = function(e) {
  cat("Error finding write_results:", e$message, "\n")
})

# Test simple functionality
library(testthat)

test_that("write_results basic debug test", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  
  # Create temporary directory for output
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  # Create a very simple mock fmri_lm object
  mask_array <- array(TRUE, c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_array, neuroim2::NeuroSpace(c(2, 2, 2)))
  
  # Mock minimal dataset
  mock_dataset <- list(mask = mask)
  class(mock_dataset) <- c("fmri_mem_dataset", "fmri_dataset", "list")
  
  # Mock minimal fmri_lm structure
  mock_betas_data <- matrix(rnorm(8 * 3), nrow = 8, ncol = 3)
  colnames(mock_betas_data) <- c("Intercept", "CondA", "CondB")
  
  mock_result <- list(
    betas = list(
      data = list(
        list(
          estimate = list(mock_betas_data)
        )
      ),
      df.residual = c(50)
    ),
    contrasts = data.frame(),
    event_indices = c(2, 3),
    baseline_indices = c(1)
  )
  
  # Mock design matrix
  mock_design <- matrix(rnorm(50 * 3), nrow = 50, ncol = 3)
  colnames(mock_design) <- c("Intercept", "CondA", "CondB")
  
  mock_model <- list(
    event_model = list(),
    baseline_model = list()
  )
  attr(mock_model, "design_matrix") <- mock_design
  class(mock_model) <- "fmri_model"
  
  mock_fmri_lm <- list(
    result = mock_result,
    model = mock_model,
    dataset = mock_dataset
  )
  class(mock_fmri_lm) <- "fmri_lm"
  
  # Test basic export
  cat("Testing write_results with mock object...\n")
  
  result <- write_results(
    mock_fmri_lm,
    path = temp_dir,
    subject = "01",
    task = "test",
    space = "MNI152NLin2009cAsym",
    save_betas = TRUE
  )
  
  cat("Result structure:\n")
  str(result)
  
  # Basic assertions
  expect_true(is.list(result))
  if ("betas" %in% names(result)) {
    expect_true(file.exists(result$betas$h5))
    expect_true(file.exists(result$betas$json))
  }
})