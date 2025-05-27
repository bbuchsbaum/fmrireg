# Debug script to trace empty column names
library(fmrireg)

# Temporarily modify the warning code to add debugging
# Let's trace the issue by adding some debug output

# Create simple test data
set.seed(123)
event_data <- data.frame(
  onset = c(10, 30, 50, 70), 
  condition = factor(c('A', 'B', 'A', 'B')), 
  run = rep(1, 4)
)

sframe <- sampling_frame(blocklens = 100, TR = 2)

# Let's check what make.names does with empty strings
cat("Testing make.names with empty strings:\n")
test_names <- c("", "", "condition_A", "condition_B")
cat("Original:", paste0("'", test_names, "'", collapse = ", "), "\n")
unique_names <- make.names(test_names, unique = TRUE)
cat("After make.names:", paste0("'", unique_names, "'", collapse = ", "), "\n")
cat("Are they identical?", identical(test_names, unique_names), "\n")

# Check which indices changed
changed_indices <- which(test_names != unique_names)
cat("Changed indices:", changed_indices, "\n")
if (length(changed_indices) > 0) {
  for (i in changed_indices) {
    cat(sprintf("  Index %d: '%s' -> '%s'\n", i, test_names[i], unique_names[i]))
  }
} 