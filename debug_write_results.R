# Quick debug script
library(fmrireg)

# Check if write_results exists
cat("write_results function exists:", exists("write_results"), "\n")

# Try to call it
tryCatch({
  write_results
  cat("write_results function found\n")
}, error = function(e) {
  cat("Error finding write_results:", e$message, "\n")
})

# Check namespace
cat("write_results in namespace:", "write_results" %in% ls("package:fmrireg"), "\n")

# Check S3 methods
cat("S3 methods for write_results:", methods("write_results"), "\n")