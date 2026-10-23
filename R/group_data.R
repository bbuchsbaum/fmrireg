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
#'     \item For "fmrilm": not supported; see below
#'   }
#' @param format Character string specifying the input format. One of "auto" (default),
#'   "h5", "nifti", "csv", or "fmrilm". If "auto", attempts to detect format.
#'   \code{format = "fmrilm"} (and a list of \code{fmri_lm} fits) is not
#'   supported and fails immediately: write each subject's fit with
#'   \code{\link{write_results}} and load the resulting beta/SE maps with
#'   \code{format = "nifti"} (or HDF5 files with \code{format = "h5"}).
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
  .group_data_gds(data, format = match.arg(format), ...)
}

# LEGACY-GROUP-DELETE (next minor): reads legacy group_data_* fields ($format, $subjects, $covariates, $paths).
#' Validate Group Data Object
#'
#' @param x A group_data object
#' @return Invisible TRUE if valid, error otherwise
#' @keywords internal
#' @noRd
validate_group_data <- function(x) {
  # Skip validation for gds-backed objects (new format)
  if (inherits(x, "group_data_gds") || inherits(x, "gds_plan") || inherits(x, "gds")) {
    return(invisible(TRUE))
  }

  if (!inherits(x, "group_data")) {
    stop("Object must be of class 'group_data'", call. = FALSE)
  }

  # Check required fields (only for legacy formats)
  required_fields <- c("format", "subjects")
  missing_fields <- setdiff(required_fields, names(x))
  if (length(missing_fields) > 0) {
    stop("Missing required fields: ", paste(missing_fields, collapse = ", "),
         call. = FALSE)
  }
  
  # Validate subjects
  if (!is.character(x$subjects) && !is.factor(x$subjects)) {
    stop("'subjects' must be a character vector or factor", call. = FALSE)
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

# LEGACY-GROUP-DELETE (next minor): reads legacy group_data_* fields ($format, $subjects, $covariates, $paths).
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

#' @export
print.group_data_gds <- function(x, ...) {
  .print_group_data_gds(x, ...)
}

# LEGACY-GROUP-DELETE (next minor): reads legacy group_data_* fields ($format, $subjects, $covariates, $paths). NB: gds objects have no summary method of their own and dispatch here today; add summary.group_data_gds or accept fmrigds's.
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
#' gd <- fmrireg:::.demo_group_data()
#' n_subjects(gd)
#' @export
n_subjects <- function(x) {
  UseMethod("n_subjects")
}

# LEGACY-GROUP-DELETE (next minor): reads legacy group_data_* fields ($format, $subjects, $covariates, $paths).
#' @export
n_subjects.group_data <- function(x) {
  length(unique(x$subjects))
}

#' @export
n_subjects.group_data_gds <- function(x) {
  .n_subjects_group_data_gds(x)
}

#' Get Subject IDs
#'
#' @param x A group_data object
#' @return Character vector of unique subject IDs
#' @examples
#' gd <- fmrireg:::.demo_group_data()
#' get_subjects(gd)
#' @export
get_subjects <- function(x) {
  UseMethod("get_subjects")
}

# LEGACY-GROUP-DELETE (next minor): reads legacy group_data_* fields ($format, $subjects, $covariates, $paths).
#' @export
get_subjects.group_data <- function(x) {
  unique(as.character(x$subjects))
}

#' @export
get_subjects.group_data_gds <- function(x) {
  .get_subjects_group_data_gds(x)
}

#' Get Covariates
#'
#' @param x A group_data object
#' @return Data frame of covariates or NULL
#' @examples
#' gd <- fmrireg:::.demo_group_data()
#' get_covariates(gd)
#' @export
get_covariates <- function(x) {
  UseMethod("get_covariates")
}

# LEGACY-GROUP-DELETE (next minor): reads legacy group_data_* fields ($format, $subjects, $covariates, $paths).
#' @export
get_covariates.group_data <- function(x) {
  x$covariates
}

#' @export
get_covariates.group_data_gds <- function(x) {
  .get_covariates_group_data_gds(x)
}

# LEGACY-GROUP-DELETE (next minor): deprecation plumbing for the legacy group-data stack.
# .legacy_group_deprecated() warns once per call of a public legacy entry point
# (silenced by options(fmrireg.suppress_deprecation = TRUE), like the
# deprecated constructors); .legacy_group_quiet() evaluates an internal call to
# a deprecated accessor without a second warning.
.legacy_group_deprecated <- function(what, instead) {
  if (!isTRUE(getOption("fmrireg.suppress_deprecation", FALSE))) {
    .Deprecated(msg = paste0(what, " is deprecated and will be removed in the next minor release. ",
                             instead))
  }
}

.legacy_group_quiet <- function(expr) {
  op <- options(fmrireg.suppress_deprecation = TRUE)
  on.exit(options(op), add = TRUE)
  expr
}
