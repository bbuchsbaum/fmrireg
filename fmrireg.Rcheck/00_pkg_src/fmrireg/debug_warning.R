# Debug script to understand the duplicate column names warning
library(fmrireg)

# Create simple test data
set.seed(123)
Y <- matrix(rnorm(1000), 100, 10)
TR <- 2

event_data <- data.frame(
  onset = c(10, 30, 50, 70), 
  condition = factor(c('A', 'B', 'A', 'B')), 
  run = rep(1, 4)
)

# Create the event model step by step to see where the warning comes from
cat("Creating sampling frame...\n")
sframe <- sampling_frame(blocklens = 100, TR = 2)

cat("Creating event model...\n")
# This is where the warning is likely coming from
model_obj <- event_model(onset ~ hrf(condition), 
                        data = event_data, 
                        block = ~ run, 
                        sampling_frame = sframe)

cat("Event model created successfully\n")
cat("Design matrix column names:", colnames(design_matrix(model_obj)), "\n")

# Check if there are any empty column names
dm <- design_matrix(model_obj)
empty_names <- which(colnames(dm) == "")
if (length(empty_names) > 0) {
  cat("Found empty column names at positions:", empty_names, "\n")
} else {
  cat("No empty column names found\n")
}

# Check the terms
cat("Model terms:\n")
print(names(model_obj$terms)) 