#!/usr/bin/env Rscript

# Automated migration script for fmrireg to fmrihrf refactoring
# Run this script from the fmrireg package root directory

library(stringr)
library(purrr)

# Configuration
DRY_RUN <- TRUE  # Set to FALSE to actually make changes

# Files to remove (already in fmrihrf)
files_to_remove <- c(
  "R/hrf.R",
  "R/hrf-functions.R",
  "R/hrf-afni.R",
  "R/hrf_decorators.R",
  "R/hrf_from_coefficients.R",
  "R/reg-constructor.R",
  "R/reg-methods.R",
  "R/regressor.R",
  "R/sampling_frame.R",
  "R/penalty_matrix.R",
  "R/reconstruction_matrix.R",
  "R/evaluate-helpers.R"
)

# Files to keep but update
files_to_update <- c(
  # Critical files
  "R/hrf-formula.R",
  "R/event_model.R",
  "R/event_model_helpers.R",
  "R/all_generic.R",
  
  # Core functionality
  "R/afni.R",
  "R/fmri_model.R",
  "R/baseline_model.R",
  "R/fmri_betas.R",
  "R/fmri_lm_methods.R",
  "R/event_vector.R",
  
  # Supporting functions
  "R/simulate.R",
  "R/benchmark_datasets.R",
  "R/data_fmri_benchmark_datasets.R",
  "R/design_plot.R",
  "R/con_stats.R",
  "R/lss.R",
  "R/basis.R",
  "R/effective_df.R",
  "R/fmri_lm_config.R"
)

# Function to update imports in a file
update_file_imports <- function(file_path, dry_run = TRUE) {
  if (!file.exists(file_path)) {
    cat("File not found:", file_path, "\n")
    return(invisible(NULL))
  }
  
  cat("Processing:", file_path, "\n")
  content <- readLines(file_path, warn = FALSE)
  original_content <- content
  
  # Define replacement patterns
  # Use word boundaries to avoid partial matches
  replacements <- list(
    # Functions that need fmrihrf:: prefix
    "\\bas_hrf\\(" = "fmrihrf::as_hrf(",
    "\\bgen_hrf\\(" = "fmrihrf::gen_hrf(",
    "\\bbind_basis\\(" = "fmrihrf::bind_basis(",
    "\\bgetHRF\\(" = "fmrihrf::getHRF(",
    
    # HRF constructors
    "\\bhrf_spmg1\\(" = "fmrihrf::hrf_spmg1(",
    "\\bhrf_spmg2\\(" = "fmrihrf::hrf_spmg2(",
    "\\bhrf_spmg3\\(" = "fmrihrf::hrf_spmg3(",
    "\\bhrf_gamma\\(" = "fmrihrf::hrf_gamma(",
    "\\bhrf_gaussian\\(" = "fmrihrf::hrf_gaussian(",
    "\\bhrf_bspline\\(" = "fmrihrf::hrf_bspline(",
    "\\bhrf_mexhat\\(" = "fmrihrf::hrf_mexhat(",
    "\\bhrf_sine\\(" = "fmrihrf::hrf_sine(",
    "\\bhrf_inv_logit\\(" = "fmrihrf::hrf_inv_logit(",
    "\\bhrf_half_cosine\\(" = "fmrihrf::hrf_half_cosine(",
    
    # Decorators
    "\\blag_hrf\\(" = "fmrihrf::lag_hrf(",
    "\\bblock_hrf\\(" = "fmrihrf::block_hrf(",
    "\\bnormalise_hrf\\(" = "fmrihrf::normalise_hrf(",
    
    # Regressor functions
    "\\bregressor\\(" = "fmrihrf::regressor(",
    "\\bsampling_frame\\(" = "fmrihrf::sampling_frame(",
    
    # Constants (only in non-quoted contexts)
    "([^\"'])\\bHRF_SPMG1\\b" = "\\1fmrihrf::HRF_SPMG1",
    "([^\"'])\\bHRF_SPMG2\\b" = "\\1fmrihrf::HRF_SPMG2",
    "([^\"'])\\bHRF_SPMG3\\b" = "\\1fmrihrf::HRF_SPMG3",
    "([^\"'])\\bHRF_GAMMA\\b" = "\\1fmrihrf::HRF_GAMMA",
    "([^\"'])\\bHRF_GAUSSIAN\\b" = "\\1fmrihrf::HRF_GAUSSIAN"
  )
  
  # Apply replacements
  changes_made <- FALSE
  for (pattern in names(replacements)) {
    new_content <- gsub(pattern, replacements[[pattern]], content, perl = TRUE)
    if (!identical(new_content, content)) {
      changes_made <- TRUE
      content <- new_content
    }
  }
  
  # Special handling for all_generic.R
  if (basename(file_path) == "all_generic.R") {
    # Skip replacements for this file - just check if we need to remove generic definitions
    content <- original_content
    
    # Look for HRF-specific generic definitions that should be removed
    hrf_generics <- c("evaluate", "nbasis", "penalty_matrix", "reconstruction_matrix")
    generic_pattern <- paste0("^(", paste(hrf_generics, collapse="|"), ") <- function")
    
    lines_to_remove <- grep(generic_pattern, content)
    if (length(lines_to_remove) > 0) {
      content <- content[-lines_to_remove]
      changes_made <- TRUE
      cat("  Will remove", length(lines_to_remove), "HRF generic definitions\n")
    }
    
    # Skip the normal replacement process for this file
    return(invisible(changes_made))
  }
  
  # Write changes
  if (changes_made) {
    if (dry_run) {
      cat("  Would update", sum(content != original_content), "lines\n")
      # Show a few example changes
      diffs <- which(content != original_content)
      if (length(diffs) > 0) {
        for (i in head(diffs, 3)) {
          cat("  Line", i, ":\n")
          cat("    Old:", original_content[i], "\n")
          cat("    New:", content[i], "\n")
        }
        if (length(diffs) > 3) {
          cat("  ... and", length(diffs) - 3, "more changes\n")
        }
      }
    } else {
      writeLines(content, file_path)
      cat("  Updated", sum(content != original_content), "lines\n")
    }
  } else {
    cat("  No changes needed\n")
  }
  
  invisible(changes_made)
}

# Function to update DESCRIPTION file
update_description <- function(dry_run = TRUE) {
  cat("\nUpdating DESCRIPTION file...\n")
  
  desc_lines <- readLines("DESCRIPTION", warn = FALSE)
  
  # Add fmrihrf to Imports (after the Imports: line)
  imports_line <- which(grepl("^Imports:", desc_lines))
  if (length(imports_line) > 0) {
    # Check if fmrihrf already added
    if (!any(grepl("fmrihrf", desc_lines))) {
      new_line <- "    fmrihrf (>= 0.1.0),"
      desc_lines <- append(desc_lines, new_line, after = imports_line)
      
      if (dry_run) {
        cat("  Would add fmrihrf to Imports\n")
      } else {
        writeLines(desc_lines, "DESCRIPTION")
        cat("  Added fmrihrf to Imports\n")
      }
    }
  }
  
  # Add to Remotes
  remotes_line <- which(grepl("^Remotes:", desc_lines))
  if (length(remotes_line) > 0) {
    if (!any(grepl("bbuchsbaum/fmrihrf", desc_lines))) {
      # Find the last Remote entry
      remote_lines <- which(grepl("^    \\S+", desc_lines) & seq_along(desc_lines) > remotes_line)
      last_remote <- if (length(remote_lines) > 0) max(remote_lines) else remotes_line
      
      new_line <- "    bbuchsbaum/fmrihrf,"
      desc_lines <- append(desc_lines, new_line, after = last_remote)
      
      if (dry_run) {
        cat("  Would add fmrihrf to Remotes\n")
      } else {
        writeLines(desc_lines, "DESCRIPTION")
        cat("  Added fmrihrf to Remotes\n")
      }
    }
  }
  
  # Update Collate to remove deleted files
  # This is complex, so just remind user
  cat("  Remember to manually update Collate section to remove deleted files\n")
}

# Function to create imports file
create_imports_file <- function(dry_run = TRUE) {
  imports_content <- '# Auto-generated file for fmrihrf imports

#\' @import fmrihrf
#\' @importFrom fmrihrf HRF as_hrf gen_hrf bind_basis getHRF
#\' @importFrom fmrihrf hrf_spmg1 hrf_spmg2 hrf_spmg3 hrf_gamma hrf_gaussian
#\' @importFrom fmrihrf hrf_bspline hrf_mexhat hrf_sine hrf_inv_logit hrf_half_cosine
#\' @importFrom fmrihrf lag_hrf block_hrf normalise_hrf
#\' @importFrom fmrihrf regressor sampling_frame
#\' @importFrom fmrihrf evaluate nbasis penalty_matrix reconstruction_matrix
#\' @importFrom fmrihrf HRF_SPMG1 HRF_SPMG2 HRF_SPMG3 HRF_GAMMA HRF_GAUSSIAN
NULL

# Re-export commonly used functions for backward compatibility
#\' @export
fmrihrf::HRF

#\' @export
fmrihrf::as_hrf

#\' @export
fmrihrf::gen_hrf

#\' @export
fmrihrf::sampling_frame

#\' @export
fmrihrf::hrf_spmg1

#\' @export
fmrihrf::regressor
'
  
  if (dry_run) {
    cat("\nWould create R/fmrihrf-imports.R with import statements\n")
  } else {
    writeLines(imports_content, "R/fmrihrf-imports.R")
    cat("\nCreated R/fmrihrf-imports.R\n")
  }
}

# Main migration function
migrate_to_fmrihrf <- function(dry_run = TRUE) {
  cat("Starting fmrireg to fmrihrf migration...\n")
  cat("DRY RUN:", dry_run, "\n\n")
  
  # Step 1: Update DESCRIPTION
  update_description(dry_run)
  
  # Step 2: Create imports file
  create_imports_file(dry_run)
  
  # Step 3: Update existing files
  cat("\nUpdating imports in R files...\n")
  walk(files_to_update, ~ update_file_imports(.x, dry_run))
  
  # Step 4: Remove redundant files
  cat("\nRemoving redundant files...\n")
  for (file in files_to_remove) {
    if (file.exists(file)) {
      if (dry_run) {
        cat("  Would remove:", file, "\n")
      } else {
        unlink(file)
        cat("  Removed:", file, "\n")
      }
    }
  }
  
  # Step 5: Update test files
  cat("\nUpdating test files...\n")
  test_files <- list.files("tests/testthat", pattern = "test.*\\.R$", full.names = TRUE)
  walk(test_files, ~ {
    if (any(grepl("hrf|regressor|sampling_frame", readLines(.x, warn = FALSE), ignore.case = TRUE))) {
      update_file_imports(.x, dry_run)
    }
  })
  
  cat("\nMigration", ifelse(dry_run, "simulation", ""), "complete!\n")
  
  if (dry_run) {
    cat("\nTo actually perform the migration, run:\n")
    cat("  source('migrate_to_fmrihrf.R')\n")
    cat("  migrate_to_fmrihrf(dry_run = FALSE)\n")
  } else {
    cat("\nNext steps:\n")
    cat("1. Run devtools::document() to update NAMESPACE\n")
    cat("2. Run devtools::test() to check for issues\n")
    cat("3. Update Collate in DESCRIPTION to remove deleted files\n")
    cat("4. Run devtools::check() for full validation\n")
  }
}

# Run the migration
if (interactive()) {
  cat("This script will migrate fmrireg to use fmrihrf.\n")
  cat("Currently in DRY RUN mode - no changes will be made.\n\n")
  migrate_to_fmrihrf(dry_run = DRY_RUN)
} else {
  # If sourced, just define the function
  cat("Migration function loaded. Run migrate_to_fmrihrf() to start.\n")
}