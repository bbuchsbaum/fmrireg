# Additional S3 methods for neuroim2 objects

#' Coerce NeuroVec to base array
#' 
#' This method allows NeuroVec objects from neuroim2 to be converted to base R arrays
#' using the standard as.array() function. This is particularly useful in testing
#' and data manipulation contexts.
#'
#' @param x A NeuroVec object from neuroim2
#' @param ... Additional arguments (currently unused)
#' @return A base R array containing the data from the NeuroVec object
#' @export
#' @method as.array NeuroVec
as.array.NeuroVec <- function(x, ...) {
  # Extract the underlying array data from the NeuroVec object
  # Try different approaches based on the object structure
  if (inherits(x, "DenseNeuroVec")) {
    # For DenseNeuroVec, access the .Data slot directly
    return(x@.Data)
  } else {
    # Fallback: try to access as an array directly
    return(as.array(as(x, "array")))
  }
}