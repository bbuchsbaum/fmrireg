# Group Data Structures for Meta-Analysis
# Generic constructor and format dispatch for group-level fMRI data

#' Create Group Dataset for Meta-Analysis
#'
#' Generic constructor that creates a group dataset from various input formats
#' for use in group-level meta-analysis with \code{\link{fmri_meta}}.
#'
#' @param data Input data. Format depends on the \code{format} argument:
#'   \itemize{
#'     \item For "h5": Character vector of HDF5 file paths
#'     \item For "nifti": List or data frame with beta/SE/variance paths
#'     \item For "csv": Path to CSV file or data frame
#'     \item For "fmrilm": List of fmri_lm objects
#'   }
#' @param format Character string specifying the input format. One of "auto" (default),
#'   "h5", "nifti", "csv", or "fmrilm". If "auto", attempts to detect format.
#' @param ... Additional arguments passed to format-specific constructors
#'
#' @return A group_data object suitable for meta-analysis
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # From HDF5 files created by write_results.fmri_lm
#' gd <- group_data(
#'   h5_paths,
#'   format = "h5",
#'   subjects = subject_ids,
#'   covariates = covariates_df,
#'   contrast = "FaceVsPlace"
#' )
#'
#' # From NIfTI files
#' gd <- group_data(
#'   list(beta = beta_paths, se = se_paths),
#'   format = "nifti",
#'   subjects = subject_ids,
#'   mask = "group_mask.nii.gz"
#' )
#' }
group_data <- function(data, format = c("auto", "h5", "nifti", "csv", "fmrilm"), ...) {
  format <- match.arg(format)
  
  if (format == "auto") {
    format <- detect_group_data_format(data)
  }
  
  switch(format,
    h5 = group_data_from_h5(data, ...),
    nifti = group_data_from_nifti(data, ...),
    csv = group_data_from_csv(data, ...),
    fmrilm = group_data_from_fmrilm(data, ...),
    fmrireg_abort(
      c(
        "Unsupported data format",
        "x" = "Format {.val {format}} is not supported",
        "i" = "Supported formats: {.val h5}, {.val nifti}, {.val csv}, {.val fmrilm}"
      )
    )
  )
}

#' Detect Input Format for Group Data
#'
#' @param data Input data to detect format from
#' @return Character string indicating detected format
#' @keywords internal
detect_group_data_format <- function(data) {
  if (is.character(data)) {
    # Check file extensions
    if (length(data) > 0) {
      first_file <- data[1]
      if (grepl("\\.h5$|\\.hdf5$", first_file, ignore.case = TRUE)) {
        return("h5")
      } else if (grepl("\\.nii$|\\.nii\\.gz$", first_file, ignore.case = TRUE)) {
        return("nifti")
      } else if (grepl("\\.csv$", first_file, ignore.case = TRUE) && length(data) == 1) {
        return("csv")
      }
    }
  } else if (is.list(data)) {
    # Check for list of fmri_lm objects
    if (length(data) > 0 && inherits(data[[1]], "fmri_lm")) {
      return("fmrilm")
    }
    # Check for named list with beta/se paths
    if (all(c("beta", "se") %in% names(data)) || all(c("beta", "var") %in% names(data))) {
      return("nifti")
    }
  } else if (is.data.frame(data)) {
    return("csv")
  }
  
  fmrireg_abort(
    c(
      "Could not automatically detect data format",
      "i" = "Unable to infer format from input data",
      ">" = "Specify format explicitly: group_data(data, format = 'h5')"
    )
  )
}

#' Validate Group Data Object
#'
#' @param x A group_data object
#' @return Invisible TRUE if valid, error otherwise
#' @keywords internal
#' @noRd
validate_group_data <- function(x) {
  if (!inherits(x, "group_data")) {
    check_inherits(x, "group_data", arg = "x")
  }

  # Check required fields
  required_fields <- c("format", "subjects")
  missing_fields <- setdiff(required_fields, names(x))
  if (length(missing_fields) > 0) {
    fmrireg_abort(
      c(
        "Invalid group_data object",
        "x" = "Missing required field{?s}: {.field {missing_fields}}"
      )
    )
  }

  # Validate subjects
  if (!is.character(x$subjects) && !is.factor(x$subjects)) {
    fmrireg_abort_input("subjects", "a character vector or factor", class(x$subjects)[1])
  }
  
  # Format-specific validation
  if (inherits(x, "group_data_h5")) {
    validate_group_data_h5(x)
  } else if (inherits(x, "group_data_nifti")) {
    validate_group_data_nifti(x)
  } else if (inherits(x, "group_data_csv")) {
    validate_group_data_csv(x)
  }
  
  invisible(TRUE)
}

#' Print Group Data Object
#'
#' @param x A group_data object
#' @param ... Additional print arguments
#' @return Invisibly returns the input object x
#' @export
print.group_data <- function(x, ...) {
  cat("Group Data Object\n")
  cat("Format:", x$format, "\n")
  cat("Subjects:", length(unique(x$subjects)), "\n")
  
  if (!is.null(x$covariates)) {
    cat("Covariates:", paste(names(x$covariates), collapse = ", "), "\n")
  }
  
  if (!is.null(x$mask)) {
    cat("Mask:", x$mask, "\n")
  }
  
  # Format-specific info
  if (inherits(x, "group_data_h5")) {
    if (!is.null(x$contrast)) {
      cat("Contrast:", x$contrast, "\n")
    }
    if (!is.null(x$stat)) {
      cat("Statistics:", paste(x$stat, collapse = ", "), "\n")
    }
  } else if (inherits(x, "group_data_nifti")) {
    available_data <- c()
    if (!is.null(x$beta_paths)) available_data <- c(available_data, "beta")
    if (!is.null(x$se_paths)) available_data <- c(available_data, "SE")
    if (!is.null(x$var_paths)) available_data <- c(available_data, "variance")
    if (!is.null(x$t_paths)) available_data <- c(available_data, "t-stat")
    cat("Available data:", paste(available_data, collapse = ", "), "\n")
  }
  
  invisible(x)
}

#' Summary of Group Data Object
#'
#' @param object A group_data object
#' @param ... Additional summary arguments
#' @return Invisibly returns the input object
#' @export
summary.group_data <- function(object, ...) {
  cat("Group Data Summary\n")
  cat("==================\n\n")
  
  # Basic info
  cat("Format:", object$format, "\n")
  cat("Number of subjects:", length(unique(object$subjects)), "\n")
  
  # Subject breakdown
  if (!is.null(object$covariates)) {
    cat("\nCovariates:\n")
    for (cov_name in names(object$covariates)) {
      cov_data <- object$covariates[[cov_name]]
      if (is.factor(cov_data) || is.character(cov_data)) {
        cat("  ", cov_name, ": ", paste(unique(cov_data), collapse = ", "), "\n", sep = "")
      } else if (is.numeric(cov_data)) {
        cat("  ", cov_name, ": range [", min(cov_data, na.rm = TRUE), 
            ", ", max(cov_data, na.rm = TRUE), "]\n", sep = "")
      }
    }
  }
  
  # Data dimensions if available
  if (!is.null(object$dim)) {
    cat("\nData dimensions:", paste(object$dim, collapse = " x "), "\n")
  }
  
  # Missing data check
  if (inherits(object, "group_data_h5") || inherits(object, "group_data_nifti")) {
    paths <- if (inherits(object, "group_data_h5")) {
      object$paths
    } else {
      object$beta_paths
    }
    
    missing_files <- paths[!file.exists(paths)]
    if (length(missing_files) > 0) {
      cat("\nWarning: ", length(missing_files), " file(s) not found\n", sep = "")
    }
  }
  
  invisible(object)
}

#' Extract Number of Subjects
#'
#' @param x A group_data object
#' @return Integer number of unique subjects
#' @examples
#' gd <- fmrireg:::.demo_group_data_csv()
#' n_subjects(gd)
#' @export
n_subjects <- function(x) {
  UseMethod("n_subjects")
}

#' @export
n_subjects.group_data <- function(x) {
  length(unique(x$subjects))
}

#' Get Subject IDs
#'
#' @param x A group_data object
#' @return Character vector of unique subject IDs
#' @examples
#' gd <- fmrireg:::.demo_group_data_csv()
#' get_subjects(gd)
#' @export
get_subjects <- function(x) {
  UseMethod("get_subjects")
}

#' @export
get_subjects.group_data <- function(x) {
  unique(as.character(x$subjects))
}

#' Get Covariates
#'
#' @param x A group_data object
#' @return Data frame of covariates or NULL
#' @examples
#' gd <- fmrireg:::.demo_group_data_csv()
#' get_covariates(gd)
#' @export
get_covariates <- function(x) {
  UseMethod("get_covariates")
}

#' @export
get_covariates.group_data <- function(x) {
  x$covariates
}
