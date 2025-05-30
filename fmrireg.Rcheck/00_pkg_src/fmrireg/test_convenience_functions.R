# Test script for glm_ols and glm_lss convenience functions

library(fmrireg)

# Create test data
set.seed(123)
Y <- matrix(rnorm(1000), 100, 10)  # 100 timepoints, 10 voxels

# Create event model
event_data <- data.frame(
  onset = c(10, 30, 50, 70),
  condition = factor(c('A', 'B', 'A', 'B')),
  run = rep(1, 4)
)

# Create dataset with event table
dset <- matrix_dataset(Y, TR = 2, run_length = 100, event_table = event_data)

sframe <- sampling_frame(blocklens = 100, TR = 2)
model_obj <- event_model(onset ~ hrf(condition), 
                        data = event_data, 
                        block = ~ run, 
                        sampling_frame = sframe)

cat('Created test data and model successfully\n')

# Test glm_ols
cat('Testing glm_ols...\n')
tryCatch({
  fit_ols <- glm_ols(dset, model_obj, HRF_SPMG1)
  cat('glm_ols completed successfully\n')
  cat('OLS result class:', class(fit_ols), '\n')
  cat('OLS betas_ran dimensions:', dim(fit_ols$betas_ran), '\n')
  cat('OLS design matrix dimensions:', dim(fit_ols$design_ran), '\n')
}, error = function(e) {
  cat('glm_ols error:', e$message, '\n')
})

# Test glm_lss
cat('Testing glm_lss...\n')
tryCatch({
  fit_lss <- glm_lss(dset, model_obj, HRF_SPMG1)
  cat('glm_lss completed successfully\n')
  cat('LSS result class:', class(fit_lss), '\n')
  cat('LSS betas_ran dimensions:', dim(fit_lss$betas_ran), '\n')
  cat('LSS design matrix dimensions:', dim(fit_lss$design_ran), '\n')
}, error = function(e) {
  cat('glm_lss error:', e$message, '\n')
})

# Test with different HRF basis
cat('Testing with HRF_SPMG2...\n')
tryCatch({
  fit_ols_spmg2 <- glm_ols(dset, model_obj, HRF_SPMG2)
  cat('glm_ols with SPMG2 completed successfully\n')
  cat('SPMG2 betas_ran dimensions:', dim(fit_ols_spmg2$betas_ran), '\n')
}, error = function(e) {
  cat('glm_ols with SPMG2 error:', e$message, '\n')
})

# Compare with traditional approach
cat('Comparing with traditional estimate_betas approach...\n')
tryCatch({
  # Traditional approach
  fit_traditional <- estimate_betas(dset, 
                                   fixed = NULL,
                                   ran = onset ~ hrf(condition, basis = HRF_SPMG1), 
                                   block = ~ run,
                                   method = "ols")
  
  cat('Traditional approach completed successfully\n')
  cat('Traditional betas_ran dimensions:', dim(fit_traditional$betas_ran), '\n')
  
  # Check if results are similar (they should be identical)
  if (all.equal(fit_ols$betas_ran, fit_traditional$betas_ran, tolerance = 1e-10)) {
    cat('✓ Results match between convenience function and traditional approach\n')
  } else {
    cat('✗ Results differ between convenience function and traditional approach\n')
  }
}, error = function(e) {
  cat('Traditional approach error:', e$message, '\n')
})

cat('Test completed!\n') 