# Fix duplicate @rdname tags and malformed @importFrom

files_to_fix <- list.files("R/", pattern = "\\.R$", full.names = TRUE)

for (file in files_to_fix) {
  lines <- readLines(file)
  fixed_lines <- character()
  i <- 1
  
  while (i <= length(lines)) {
    line <- lines[i]
    
    # Check for duplicate @rdname on consecutive lines
    if (i < length(lines) && 
        grepl("^#' @rdname", line) && 
        grepl("^#' @rdname", lines[i+1])) {
      # Skip the duplicate
      fixed_lines <- c(fixed_lines, line)
      i <- i + 2  # Skip the duplicate line
    } else {
      fixed_lines <- c(fixed_lines, line)
      i <- i + 1
    }
  }
  
  # Write back if changes were made
  if (length(fixed_lines) != length(lines)) {
    cat("Fixed duplicate @rdname in", file, "\n")
    writeLines(fixed_lines, file)
  }
}

# Also fix the malformed @importFrom directives
fix_importfrom <- function(file) {
  content <- readLines(file)
  
  # Fix malformed @importFrom in file headers
  for (i in seq_along(content)) {
    if (grepl("^#' @importFrom .+ [A-Z]", content[i])) {
      # This looks like a description accidentally in @importFrom
      content[i] <- "#' @description"
    }
  }
  
  writeLines(content, file)
}

# Fix specific files with malformed @importFrom
problem_files <- c(
  "R/fmri_lm_chunkwise.R",
  "R/fmri_lm_internal.R", 
  "R/fmri_lm_methods.R",
  "R/fmri_lm_runwise.R",
  "R/fmri_lm_strategies.R",
  "R/fmrilm_new.R"
)

for (file in problem_files) {
  if (file.exists(file)) {
    cat("Checking", file, "for malformed @importFrom\n")
    fix_importfrom(file)
  }
}