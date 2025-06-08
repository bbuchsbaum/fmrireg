library(testthat)
library(fmrireg)

# Helper function to create test data
create_test_data <- function(seed = 123) {
  set.seed(seed)
  Y <- matrix(rnorm(1000), 100, 10)  # 100 timepoints, 10 voxels
  
  event_data <- data.frame(
    onset = c(10, 30, 50, 70),
    condition = factor(c('A', 'B', 'A', 'B')),
    run = rep(1, 4)
  )
  
  dset <- matrix_dataset(Y, TR = 2, run_length = 100, event_table = event_data)
  
  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
  model_obj <- event_model(onset ~ hrf(condition), 
                          data = event_data, 
                          block = ~ run, 
                          sampling_frame = sframe)
  
  list(dataset = dset, model = model_obj, event_data = event_data)
}

# Tests for glm_ols ----

test_that("glm_ols works with basic inputs", {
  test_data <- create_test_data()
  
  result <- glm_ols(test_data$dataset, test_data$model, fmrihrf::HRF_SPMG1)
  
  expect_s3_class(result, "fmri_betas")
  expect_true(!is.null(result$betas_ran))
  expect_true(!is.null(result$design_ran))
  
  # Should have 2 conditions (A, B) x 10 voxels
  expect_equal(dim(result$betas_ran), c(2, 10))
  expect_equal(dim(result$design_ran), c(100, 2))
})

test_that("glm_ols works with different HRF bases", {
  test_data <- create_test_data()
  
  # Test with SPMG1 (canonical)
  result_spmg1 <- glm_ols(test_data$dataset, test_data$model, fmrihrf::HRF_SPMG1)
  expect_equal(dim(result_spmg1$betas_ran), c(2, 10))
  
  # Test with SPMG2 (canonical + temporal derivative)
  result_spmg2 <- glm_ols(test_data$dataset, test_data$model, fmrihrf::HRF_SPMG2)
  expect_equal(dim(result_spmg2$betas_ran), c(4, 10))  # 2 conditions x 2 basis functions
  
  # Test with SPMG3 (canonical + both derivatives)
  result_spmg3 <- glm_ols(test_data$dataset, test_data$model, fmrihrf::HRF_SPMG3)
  expect_equal(dim(result_spmg3$betas_ran), c(6, 10))  # 2 conditions x 3 basis functions
})

test_that("glm_ols validates inputs correctly", {
  test_data <- create_test_data()
  
  # Test with invalid dataset
  expect_error(glm_ols("not_a_dataset", test_data$model, fmrihrf::HRF_SPMG1),
               "dataset must be a matrix_dataset object")
  
  # Test with invalid model
  expect_error(glm_ols(test_data$dataset, "not_a_model", fmrihrf::HRF_SPMG1),
               "model_obj must be an event_model object")
})

test_that("glm_ols matches traditional estimate_betas approach", {
  test_data <- create_test_data()
  
  # Convenience function approach
  result_convenience <- glm_ols(test_data$dataset, test_data$model, fmrihrf::HRF_SPMG1)
  
  # Traditional approach
  result_traditional <- estimate_betas(test_data$dataset, 
                                     fixed = NULL,
                                     ran = onset ~ hrf(condition, basis = fmrihrf::HRF_SPMG1), 
                                     block = ~ run,
                                     method = "ols")
  
  # Results should be identical
  expect_equal(result_convenience$betas_ran, result_traditional$betas_ran, tolerance = 1e-10)
  expect_equal(dim(result_convenience$design_ran), dim(result_traditional$design_ran))
})

# Tests for glm_lss ----

test_that("glm_lss works with basic inputs", {
  test_data <- create_test_data()
  
  result <- glm_lss(test_data$dataset, test_data$model, fmrihrf::HRF_SPMG1)
  
  expect_s3_class(result, "fmri_betas")
  expect_true(!is.null(result$betas_ran))
  expect_true(!is.null(result$design_ran))
  
  # Should have 2 trials (single trial estimation) x 10 voxels
  # Note: LSS estimates separate betas for each trial
  expect_equal(dim(result$betas_ran), c(2, 10))
  expect_equal(dim(result$design_ran), c(100, 2))
})

test_that("glm_lss works with different HRF bases", {
  test_data <- create_test_data()
  
  # Test with SPMG1 (canonical)
  result_spmg1 <- glm_lss(test_data$dataset, test_data$model, fmrihrf::HRF_SPMG1)
  expect_equal(dim(result_spmg1$betas_ran), c(2, 10))
  
  # Test with SPMG2 (canonical + temporal derivative)
  result_spmg2 <- glm_lss(test_data$dataset, test_data$model, fmrihrf::HRF_SPMG2)
  expect_equal(dim(result_spmg2$betas_ran), c(4, 10))  # 2 conditions x 2 basis functions
})

test_that("glm_lss validates inputs correctly", {
  test_data <- create_test_data()
  
  # Test with invalid dataset
  expect_error(glm_lss("not_a_dataset", test_data$model, fmrihrf::HRF_SPMG1),
               "dataset must be a matrix_dataset object")
  
  # Test with invalid model
  expect_error(glm_lss(test_data$dataset, "not_a_model", fmrihrf::HRF_SPMG1),
               "model_obj must be an event_model object")
})

test_that("glm_lss works with both C++ and R implementations", {
  test_data <- create_test_data()
  
  # Test with C++ implementation
  result_cpp <- glm_lss(test_data$dataset, test_data$model, fmrihrf::HRF_SPMG1, use_cpp = TRUE)
  expect_s3_class(result_cpp, "fmri_betas")
  
  # Test with R implementation
  result_r <- glm_lss(test_data$dataset, test_data$model, fmrihrf::HRF_SPMG1, use_cpp = FALSE)
  expect_s3_class(result_r, "fmri_betas")
  
  # Both should have same dimensions
  expect_equal(dim(result_cpp$betas_ran), dim(result_r$betas_ran))
})

test_that("glm_lss matches traditional estimate_betas LSS approach", {
  test_data <- create_test_data()
  
  # Convenience function approach
  result_convenience <- glm_lss(test_data$dataset, test_data$model, fmrihrf::HRF_SPMG1)
  
  # Traditional approach
  result_traditional <- estimate_betas(test_data$dataset, 
                                     fixed = NULL,
                                     ran = onset ~ hrf(condition, basis = fmrihrf::HRF_SPMG1), 
                                     block = ~ run,
                                     method = "lss")
  
  # Results should be identical
  expect_equal(result_convenience$betas_ran, result_traditional$betas_ran, tolerance = 1e-10)
  expect_equal(dim(result_convenience$design_ran), dim(result_traditional$design_ran))
})

# Comparison tests ----

test_that("glm_ols and glm_lss produce different results as expected", {
  test_data <- create_test_data()
  
  result_ols <- glm_ols(test_data$dataset, test_data$model, fmrihrf::HRF_SPMG1)
  result_lss <- glm_lss(test_data$dataset, test_data$model, fmrihrf::HRF_SPMG1)
  
  # Both should have same dimensions for this simple case
  expect_equal(dim(result_ols$betas_ran), dim(result_lss$betas_ran))
  
  # But the actual beta values should be different (OLS averages, LSS doesn't)
  expect_false(identical(result_ols$betas_ran, result_lss$betas_ran))
})

test_that("convenience functions work with baseline models", {
  test_data <- create_test_data()
  
  # Create a baseline model
  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
  basemod <- baseline_model("constant", sframe = sframe)
  
  # Test both functions with baseline model
  result_ols <- glm_ols(test_data$dataset, test_data$model, fmrihrf::HRF_SPMG1, basemod = basemod)
  result_lss <- glm_lss(test_data$dataset, test_data$model, fmrihrf::HRF_SPMG1, basemod = basemod)
  
  expect_s3_class(result_ols, "fmri_betas")
  expect_s3_class(result_lss, "fmri_betas")
  expect_true(!is.null(result_ols$design_base))
  expect_true(!is.null(result_lss$design_base))
})

test_that("convenience functions work with different block specifications", {
  test_data <- create_test_data()
  
  # Test with different block specification
  result_ols <- glm_ols(test_data$dataset, test_data$model, fmrihrf::HRF_SPMG1, block = ~ run)
  result_lss <- glm_lss(test_data$dataset, test_data$model, fmrihrf::HRF_SPMG1, block = ~ run)
  
  expect_s3_class(result_ols, "fmri_betas")
  expect_s3_class(result_lss, "fmri_betas")
}) 
