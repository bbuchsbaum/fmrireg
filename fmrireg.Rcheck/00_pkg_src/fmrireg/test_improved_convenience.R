# Test script for improved convenience functions
library(fmrireg)

cat("Testing convenience functions with proper matrix_dataset creation...\n")

# Create test data
Y_noisy <- matrix(rnorm(147 * 10), 147, 10)
TR <- 2

# Create sampling frame and event model
sf <- sampling_frame(blocklens = 147, TR = TR)
events_df <- data.frame(
  onset = c(10, 30, 50), 
  condition = factor(paste0('cond', 1:3)), 
  block = 1
)

model_obj <- event_model(onset ~ hrf(condition), 
                        data = events_df, 
                        block = ~ block, 
                        sampling_frame = sf)

cat("Model created successfully\n")

# Create matrix_dataset with event table (this is the correct approach)
# FIXED: Use event_table parameter instead of data
dset <- matrix_dataset(Y_noisy, TR = TR, run_length = 147, event_table = events_df)
cat("Matrix dataset created successfully\n")

# Test glm_ols with matrix_dataset
cat("Testing glm_ols with matrix_dataset...\n")
tryCatch({
  result1 <- fmrireg::glm_ols(dset, model_obj, HRF_SPMG1)
  cat("glm_ols SUCCESS! Result class:", class(result1), "\n")
  cat("Result dimensions:", dim(result1$betas_ran), "\n")
}, error = function(e) {
  cat("glm_ols ERROR:", e$message, "\n")
})

# Test glm_lss with matrix_dataset
cat("Testing glm_lss with matrix_dataset...\n")
tryCatch({
  result2 <- fmrireg::glm_lss(dset, model_obj, HRF_SPMG1)
  cat("glm_lss SUCCESS! Result class:", class(result2), "\n")
  cat("Result dimensions:", dim(result2$betas_ran), "\n")
}, error = function(e) {
  cat("glm_lss ERROR:", e$message, "\n")
})

cat("Test completed!\n") 