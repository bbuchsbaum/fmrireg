# Script to fix S3 method documentation for CRAN compliance

# List of files to check and fix
files_to_check <- list.files("R/", pattern = "\\.R$", full.names = TRUE)

# Common S3 generics that have main documentation in all_generic.R
generics_in_all_generic <- c(
  "design_matrix", "term_matrices", "conditions", "elements", 
  "construct", "convolve", "levels", "onsets", "durations",
  "blockids", "samples", "global_onsets", "cells", "split_onsets",
  "event_terms", "baseline_terms", "columns", "nbasis", "term_names",
  "shortnames", "longnames", "parent_terms", "shift", "p_values",
  "standard_error", "stats", "formula", "sub_basis", "is_continuous",
  "is_categorical", "terms", "coef", "get_formula", "get_mask",
  "get_data", "coefficients"
)

# Read each file and check for S3 methods
for (file in files_to_check) {
  if (file == "R/all_generic.R") next  # Skip the generics file
  
  lines <- readLines(file)
  in_doc_block <- FALSE
  doc_start <- NULL
  method_line <- NULL
  
  for (i in seq_along(lines)) {
    line <- lines[i]
    
    # Check if we're starting a documentation block
    if (grepl("^#'", line) && !in_doc_block) {
      in_doc_block <- TRUE
      doc_start <- i
    }
    
    # Check if we're ending a documentation block and hitting a function
    if (in_doc_block && !grepl("^#'", line)) {
      # Check if this is an S3 method definition
      if (grepl("^[a-zA-Z_]+\\.[a-zA-Z_].*<-.*function", line)) {
        method_line <- i
        
        # Extract method name
        method_name <- sub("^([a-zA-Z_]+\\.[a-zA-Z_]+).*", "\\1", line)
        generic_name <- sub("^([a-zA-Z_]+)\\..*", "\\1", method_name)
        
        # Check if this generic is documented in all_generic.R
        if (generic_name %in% generics_in_all_generic) {
          # Check if documentation has more than just @rdname and @export
          doc_lines <- lines[doc_start:(i-1)]
          
          has_redundant_docs <- any(grepl("^#' @(param|return|description|details|examples|seealso)", doc_lines))
          has_rdname <- any(grepl("^#' @rdname", doc_lines))
          has_export <- any(grepl("^#' @export", doc_lines))
          
          if (has_redundant_docs) {
            cat("\n=== Found redundant docs in", file, "===\n")
            cat("Method:", method_name, "\n")
            cat("Generic:", generic_name, "\n")
            cat("Line:", method_line, "\n")
            
            # Show current documentation
            cat("\nCurrent documentation:\n")
            cat(doc_lines, sep = "\n")
            cat("\n")
          }
        }
      }
      in_doc_block <- FALSE
    }
  }
}

# Also check for \dontrun{} usage
cat("\n\n=== Checking for \\dontrun{} usage ===\n")
for (file in files_to_check) {
  lines <- readLines(file)
  dontrun_lines <- grep("\\\\dontrun\\{", lines)
  if (length(dontrun_lines) > 0) {
    cat("\nFound \\dontrun{} in", file, "at lines:", dontrun_lines, "\n")
    for (line_num in dontrun_lines) {
      # Show context
      start <- max(1, line_num - 2)
      end <- min(length(lines), line_num + 2)
      cat("Context:\n")
      cat(paste(start:end, lines[start:end]), sep = "\n")
    }
  }
}