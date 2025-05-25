test_that("benchmark dataset loading works", {
  skip_if_not_installed("fmrireg")
  
  # Test that we can list datasets
  datasets_info <- list_benchmark_datasets()
  expect_true(is.data.frame(datasets_info))
  expect_true(nrow(datasets_info) >= 5)
  expect_true(all(c("Dataset", "Description") %in% colnames(datasets_info)))
  
  # Test loading a specific dataset
  data <- load_benchmark_dataset("BM_Canonical_HighSNR")
  expect_true(is.list(data))
  expect_true("Y_noisy" %in% names(data))
  expect_true("event_onsets" %in% names(data))
  expect_true("condition_labels" %in% names(data))
  expect_true("true_betas_condition" %in% names(data))
  
  # Check dimensions are reasonable
  expect_true(is.matrix(data$Y_noisy))
  expect_true(ncol(data$Y_noisy) == 100)  # 100 voxels
  expect_true(nrow(data$Y_noisy) > 100)   # Should have many time points
  
  # Test metadata loading
  metadata <- load_benchmark_dataset("metadata")
  expect_true(is.list(metadata))
  expect_true("creation_date" %in% names(metadata))
  expect_true("description" %in% names(metadata))
})

test_that("benchmark dataset summary works", {
  skip_if_not_installed("fmrireg")
  
  summary_info <- get_benchmark_summary("BM_Canonical_HighSNR")
  expect_true(is.list(summary_info))
  expect_true("dimensions" %in% names(summary_info))
  expect_true("experimental_design" %in% names(summary_info))
  expect_true("hrf_information" %in% names(summary_info))
  
  # Check dimensions
  dims <- summary_info$dimensions
  expect_equal(dims$n_voxels, 100)
  expect_equal(dims$n_conditions, 3)
  expect_true(dims$n_events > 0)
})

test_that("design matrix creation works", {
  skip_if_not_installed("fmrireg")
  
  # Test design matrix creation
  X <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", HRF_SPMG1)
  expect_true(is.matrix(X))
  expect_true(ncol(X) == 4)  # 3 conditions + intercept
  expect_true(nrow(X) > 100)  # Should have many time points
  expect_true("Intercept" %in% colnames(X))
  
  # Test without intercept
  X_no_int <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", HRF_SPMG1, 
                                                  include_intercept = FALSE)
  expect_true(ncol(X_no_int) == 3)  # 3 conditions only
  expect_false("Intercept" %in% colnames(X_no_int))
})

test_that("performance evaluation works", {
  skip_if_not_installed("fmrireg")
  
  # Load dataset and create simple test
  data <- load_benchmark_dataset("BM_Canonical_HighSNR")
  true_betas <- data$true_betas_condition
  
  # Create some fake estimated betas (just add noise to true betas)
  set.seed(123)
  estimated_betas <- true_betas + matrix(rnorm(prod(dim(true_betas)), 0, 0.1), 
                                        nrow = nrow(true_betas))
  
  # Evaluate performance
  performance <- evaluate_method_performance("BM_Canonical_HighSNR", 
                                            estimated_betas, 
                                            "Test_Method")
  
  expect_true(is.list(performance))
  expect_true("overall_metrics" %in% names(performance))
  expect_true("condition_metrics" %in% names(performance))
  expect_true("voxel_metrics" %in% names(performance))
  
  # Check that correlation is high (since we added only small noise)
  expect_true(performance$overall_metrics$correlation > 0.9)
  expect_true(performance$overall_metrics$mse < 1.0)
})

test_that("all benchmark datasets can be loaded", {
  skip_if_not_installed("fmrireg")
  
  # Get list of all datasets
  datasets_info <- list_benchmark_datasets()
  dataset_names <- datasets_info$Dataset
  
  # Try to load each dataset
  for (dataset_name in dataset_names) {
    data <- load_benchmark_dataset(dataset_name)
    expect_true(is.list(data))
    expect_true("Y_noisy" %in% names(data))
    expect_true("event_onsets" %in% names(data))
    expect_true("condition_labels" %in% names(data))
    
    # Check that dimensions are consistent
    expect_equal(nrow(data$Y_noisy), length(seq(0, data$total_time, by = data$TR)))
    expect_equal(length(data$event_onsets), length(data$condition_labels))
  }
})

test_that("error handling works correctly", {
  skip_if_not_installed("fmrireg")
  
  # Test invalid dataset name
  expect_error(load_benchmark_dataset("NonExistent_Dataset"))
  expect_error(get_benchmark_summary("NonExistent_Dataset"))
  
  # Test performance evaluation with wrong dimensions
  data <- load_benchmark_dataset("BM_Canonical_HighSNR")
  wrong_betas <- matrix(1, nrow = 2, ncol = 50)  # Wrong dimensions
  expect_error(evaluate_method_performance("BM_Canonical_HighSNR", wrong_betas, "Test"))
}) 