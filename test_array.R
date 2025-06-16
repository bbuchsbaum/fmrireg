# Test as.array.NeuroVec fix
devtools::load_all('.')

# Test the as.array method
library(neuroim2)

# Create a simple NeuroVec (needs to be 4D)
arr <- array(1:24, c(2, 3, 2, 2))
space <- NeuroSpace(c(2, 3, 2, 2))
nvec <- NeuroVec(arr, space)

cat("Testing as.array on NeuroVec:\n")
cat("Original array class:", class(arr), "\n")
cat("NeuroVec class:", class(nvec), "\n")

# Test the conversion
tryCatch({
  result <- as.array(nvec)
  cat("as.array() succeeded!\n")
  cat("Result class:", class(result), "\n")
  cat("Result dimensions:", dim(result), "\n")
  cat("Values match original:", identical(as.vector(result), as.vector(arr)), "\n")
}, error = function(e) {
  cat("as.array() failed:", e$message, "\n")
})