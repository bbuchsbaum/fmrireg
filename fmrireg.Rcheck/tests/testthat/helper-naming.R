# Helper functions for naming tests

#' Check if a String is a Valid fmrireg Column Heading
#'
#' Based on the agreed grammar, allowing letters, numbers, dot, underscore, and hash.
#' Starts with a letter or dot.
#'
#' @param x Character vector of potential headings.
#' @return Logical vector.
is_valid_heading <- function(x) {
  # Allows letters, numbers, ., _, and # (for unique tags)
  # Must start with a letter or dot.
  grepl("^[A-Za-z\\.][A-Za-z0-9\\._#]*$", x)
} 