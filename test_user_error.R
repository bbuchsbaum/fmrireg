# Test script to reproduce the user's error and verify the fix
library(fmrireg)

set.seed(123)
Y_noisy <- matrix(rnorm(1000), 100, 10)
TR <- 2

event_data <- data.frame(
  onset = c(10, 30, 50, 70), 
  condition = factor(c('A', 'B', 'A', 'B')), 
  run = rep(1, 4)
)

sframe <- sampling_frame(blocklens = 100, TR = 2)
model_obj <- event_model(onset ~ hrf(condition), 
                        data = event_data, 
                        block = ~ run, 
                        sampling_frame = sframe)

basis_cfals <- HRF_SPMG1

# This was the line causing the error
result <- fmrireg::glm_ols(matrix_dataset(Y_noisy, TR = TR, run_length = nrow(Y_noisy), event_table = event_data), 
                          model_obj, basis_cfals)

cat('Success! Result class:', class(result), '\n')
cat('Result dimensions:', dim(result$betas_ran), '\n') 