# Test script to verify that real duplicate warnings still work
library(fmrireg)

# Create test data with variables that will cause real duplicates
set.seed(123)
event_data <- data.frame(
  onset = c(10, 30, 50, 70, 90, 110), 
  RT1 = rnorm(6),
  RT2 = rnorm(6),
  Age = rnorm(6),
  block = rep(1, 6)
)

sframe <- sampling_frame(blocklens = 120, TR = 2)

cat("Testing case that should produce a warning (real duplicates)...\n")

# This should still produce a warning because RT1 appears in both terms
tryCatch({
  model_with_duplicates <- event_model(onset ~ hrf(Ident(RT1, Age)) + hrf(Ident(RT1, RT2)), 
                                       data = event_data, 
                                       block = ~ block, 
                                       sampling_frame = sframe)
  cat("Model created successfully\n")
  cat("Column names:", colnames(design_matrix(model_with_duplicates)), "\n")
}, warning = function(w) {
  cat("Warning (as expected):", w$message, "\n")
}, error = function(e) {
  cat("Error:", e$message, "\n")
}) 