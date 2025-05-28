# Script to automatically fix S3 method documentation for CRAN compliance

library(stringr)

# List of generics documented in all_generic.R
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

# Function to fix a single file
fix_file <- function(file_path) {
  if (basename(file_path) == "all_generic.R") return(NULL)
  
  lines <- readLines(file_path)
  changes_made <- FALSE
  output_lines <- character()
  i <- 1
  
  while (i <= length(lines)) {
    line <- lines[i]
    
    # Check if this is the start of a documentation block
    if (grepl("^#'", line)) {
      doc_start <- i
      doc_lines <- character()
      
      # Collect all documentation lines
      while (i <= length(lines) && grepl("^#'", lines[i])) {
        doc_lines <- c(doc_lines, lines[i])
        i <- i + 1
      }
      
      # Check if the next line is an S3 method definition
      if (i <= length(lines) && grepl("^[a-zA-Z_]+\\.[a-zA-Z_].*<-.*function", lines[i])) {
        method_line <- lines[i]
        method_name <- sub("^([a-zA-Z_]+\\.[a-zA-Z_]+).*", "\\1", method_line)
        generic_name <- sub("^([a-zA-Z_]+)\\..*", "\\1", method_name)
        
        # Check if this generic is documented in all_generic.R
        if (generic_name %in% generics_in_all_generic) {
          # Check if documentation has redundant tags
          has_redundant <- any(grepl("^#' @(param|return|description|details|examples|seealso)", doc_lines))
          has_rdname <- any(grepl("^#' @rdname", doc_lines))
          has_export <- any(grepl("^#' @export", doc_lines))
          has_method <- any(grepl("^#' @method", doc_lines))
          
          if (has_redundant) {
            cat("Fixing", method_name, "in", file_path, "\n")
            changes_made <- TRUE
            
            # Keep only essential lines
            new_doc <- character()
            
            # Add a simple title if there isn't one
            if (!has_rdname) {
              new_doc <- c(new_doc, paste0("#' @rdname ", generic_name))
            }
            
            # Keep @method if present
            method_lines <- doc_lines[grepl("^#' @method", doc_lines)]
            if (length(method_lines) > 0) {
              new_doc <- c(new_doc, method_lines)
            } else if (!has_method) {
              # Add @method if missing
              class_name <- sub("^[a-zA-Z_]+\\.([a-zA-Z_]+).*", "\\1", method_name)
              new_doc <- c(new_doc, paste0("#' @method ", generic_name, " ", class_name))
            }
            
            # Keep @rdname
            rdname_lines <- doc_lines[grepl("^#' @rdname", doc_lines)]
            if (length(rdname_lines) > 0) {
              new_doc <- c(new_doc, rdname_lines)
            } else {
              new_doc <- c(new_doc, paste0("#' @rdname ", generic_name))
            }
            
            # Keep @export
            if (has_export) {
              new_doc <- c(new_doc, "#' @export")
            }
            
            # Keep any @importFrom lines
            import_lines <- doc_lines[grepl("^#' @importFrom", doc_lines)]
            if (length(import_lines) > 0) {
              new_doc <- c(new_doc, import_lines)
            }
            
            # Keep @keywords internal if present
            internal_lines <- doc_lines[grepl("^#' @keywords internal", doc_lines)]
            if (length(internal_lines) > 0) {
              new_doc <- c(new_doc, internal_lines)
            }
            
            # Keep @noRd if present
            nord_lines <- doc_lines[grepl("^#' @noRd", doc_lines)]
            if (length(nord_lines) > 0) {
              new_doc <- c(new_doc, nord_lines)
            }
            
            output_lines <- c(output_lines, new_doc)
          } else {
            # Keep original documentation
            output_lines <- c(output_lines, doc_lines)
          }
        } else {
          # Not a generic in all_generic.R, keep original
          output_lines <- c(output_lines, doc_lines)
        }
        
        # Add the method definition line
        output_lines <- c(output_lines, method_line)
        i <- i + 1
      } else {
        # Not followed by a method definition, keep original
        output_lines <- c(output_lines, doc_lines)
        if (i <= length(lines)) {
          output_lines <- c(output_lines, lines[i])
          i <- i + 1
        }
      }
    } else {
      # Not a documentation line
      output_lines <- c(output_lines, line)
      i <- i + 1
    }
  }
  
  if (changes_made) {
    writeLines(output_lines, file_path)
    return(file_path)
  } else {
    return(NULL)
  }
}

# Process all R files
files_to_check <- list.files("R/", pattern = "\\.R$", full.names = TRUE)
fixed_files <- character()

for (file in files_to_check) {
  result <- fix_file(file)
  if (!is.null(result)) {
    fixed_files <- c(fixed_files, result)
  }
}

if (length(fixed_files) > 0) {
  cat("\n\nFixed", length(fixed_files), "files:\n")
  cat(paste("-", fixed_files), sep = "\n")
} else {
  cat("\n\nNo files needed fixing.\n")
}