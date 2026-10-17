# Group Data Support for HDF5/LabeledVolumeSet Files
# Reads fmristore HDF5 files created by write_results.fmri_lm

#' Create Group Dataset from HDF5 Files
#'
#' Creates a group dataset from HDF5 files produced by \code{\link{write_results.fmri_lm}}.
#' These files use the fmristore LabeledVolumeSet format for efficient storage of
#' multiple statistical maps.
#'
#' @param paths Character vector of HDF5 file paths, one per subject
#' @param subjects Character vector of subject identifiers. If NULL, extracted from file paths.
#' @param covariates Data frame of subject-level covariates (optional)
#' @param mask Path to mask file or mask object (optional)
#' @param contrast Character string specifying which contrast to extract (for multi-contrast files)
#' @param stat Character vector of statistics to extract (e.g., c("beta", "se", "tstat"))
#' @param validate Logical. Validate that all files exist and contain expected data (default: TRUE)
#'
#' @return A group_data_h5 object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Read HDF5 files from write_results.fmri_lm
#' subjects <- data.frame(
#'   subject = sprintf("sub-%02d", 1:20),
#'   group = rep(c("young", "old"), each = 10),
#'   age = c(rnorm(10, 25, 3), rnorm(10, 70, 5))
#' )
#'
#' h5_paths <- sprintf("derivatives/sub-%02d_task-nback_desc-GLMstatmap_bold.h5", 1:20)
#'
#' gd <- group_data_from_h5(
#'   h5_paths,
#'   subjects = subjects$subject,
#'   covariates = subjects[c("group", "age")],
#'   contrast = "FaceVsPlace",
#'   stat = c("beta", "se")
#' )
#' }
group_data_from_h5 <- function(paths,
                               subjects = NULL,
                               covariates = NULL,
                               mask = NULL,
                               contrast = NULL,
                               stat = c("beta", "se"),
                               validate = TRUE) {
  if (!isTRUE(getOption("fmrireg.suppress_deprecation", FALSE))) {
    .Deprecated("group_data", msg = "group_data_from_h5() is deprecated. Use group_data(paths, format='h5', ...) instead.")
  }

  # Validate paths
  if (!is.character(paths)) {
    stop("'paths' must be a character vector of file paths", call. = FALSE)
  }
  
  if (validate) {
    missing_files <- paths[!file.exists(paths)]
    if (length(missing_files) > 0) {
      stop("The following files do not exist:\n  ", 
           paste(head(missing_files, 5), collapse = "\n  "),
           if (length(missing_files) > 5) paste0("\n  ... and ", length(missing_files) - 5, " more"),
           call. = FALSE)
    }
  }
  
  # Extract or validate subjects
  if (is.null(subjects)) {
    subjects <- extract_subject_ids_from_paths(paths)
  } else {
    if (length(subjects) != length(paths)) {
      stop("Length of 'subjects' must match length of 'paths'", call. = FALSE)
    }
  }
  
  # Process covariates
  if (!is.null(covariates)) {
    if (!is.data.frame(covariates)) {
      stop("'covariates' must be a data frame", call. = FALSE)
    }
    if (nrow(covariates) != length(paths)) {
      stop("Number of rows in 'covariates' must match length of 'paths'", call. = FALSE)
    }
  }
  
  # Read metadata from first file to get dimensions and available labels
  meta <- read_h5_metadata(paths[1])
  
  # Validate requested statistics are available
  if (!is.null(contrast)) {
    if (!contrast %in% meta$labels) {
      stop("Contrast '", contrast, "' not found in HDF5 file. Available: ",
           paste(meta$labels, collapse = ", "), call. = FALSE)
    }
  }
  
  # Check if this is a multi-statistic file (e.g., from save_contrasts_by_stat)
  # or a multi-contrast file (e.g., from save_contrasts_by_contrast)
  file_type <- detect_h5_file_type(meta)
  
  # Validate requested statistics based on file type
  if (file_type == "by_stat") {
    available_stats <- meta$labels
    missing_stats <- setdiff(stat, available_stats)
    if (length(missing_stats) > 0) {
      stop("The following statistics are not available: ", 
           paste(missing_stats, collapse = ", "),
           "\nAvailable: ", paste(available_stats, collapse = ", "),
           call. = FALSE)
    }
  } else if (file_type == "by_contrast") {
    # For by_contrast files, stats are the labels within each contrast
    if (!is.null(contrast)) {
      # Will validate when reading specific contrast
    }
  }
  
  # Create group_data_h5 object with lazy loading
  structure(
    list(
      paths = paths,
      subjects = subjects,
      covariates = covariates,
      mask = mask,
      contrast = contrast,
      stat = stat,
      dim = meta$dim,
      labels = meta$labels,
      file_type = file_type,
      format = "h5",
      n_voxels = prod(meta$dim[1:3]),  # Spatial dimensions
      n_volumes = meta$dim[4]           # Number of labels/statistics
    ),
    class = c("group_data_h5", "group_data")
  )
}

#' Read Metadata from HDF5 File
#'
#' @param path Path to HDF5 file
#' @return List with metadata including dimensions and labels
#' @keywords internal
read_h5_metadata <- function(path) {
  # Use fmristore to read metadata without fully loading data
  tryCatch({
    h5_handle <- fmristore::read_labeled_vec(path)

    # Best-effort close on exit (support both list-based and S4 handles)
    on.exit({
      try(h5_handle$close_all(), silent = TRUE)
      try(h5_handle$close(), silent = TRUE)
      try(close(h5_handle), silent = TRUE)
      try(close(h5_handle@obj), silent = TRUE)
      try(h5_handle@obj$close_all(), silent = TRUE)
    }, add = TRUE)

    # Legacy list-like handle with $ fields
    if (is.list(h5_handle) && all(c("dim","labels") %in% names(h5_handle))) {
      return(list(
        dim = h5_handle$dim,
        labels = h5_handle$labels,
        mask_dim = h5_handle$mask_dim %||% NULL,
        has_mask = !is.null(h5_handle$mask)
      ))
    }

    # S4 handle (current fmristore LabeledVolumeSet)
    if (isS4(h5_handle)) {
      sn <- try(methods::slotNames(h5_handle), silent = TRUE)
      labels_val <- tryCatch(as.character(methods::slot(h5_handle, "labels")), error = function(e) NULL)
      mask_val <- tryCatch(methods::slot(h5_handle, "mask"), error = function(e) NULL)
      mask_dim <- if (!is.null(mask_val)) tryCatch(dim(mask_val), error = function(e) NULL) else NULL
      dim_val <- NULL

      if (!is.null(mask_dim) && length(mask_dim) >= 3 && !is.null(labels_val)) {
        dim_val <- c(as.integer(mask_dim[1:3]), length(labels_val))
      } else if (!inherits(sn, "try-error") && "dim" %in% sn) {
        dim_val <- tryCatch(methods::slot(h5_handle, "dim"), error = function(e) NULL)
      }

      if (!is.null(dim_val) && !is.null(labels_val)) {
        return(list(
          dim = as.integer(dim_val),
          labels = labels_val,
          mask_dim = if (!is.null(mask_dim)) as.integer(mask_dim) else NULL,
          has_mask = !is.null(mask_val)
        ))
      }

      stop("Failed to introspect fmristore S4 HDF5 handle")
    }

    # Unknown handle type
    stop("Unsupported fmristore handle type")
  }, error = function(e) {
    stop("Failed to read HDF5 metadata from ", path, ": ", e$message, call. = FALSE)
  })
}

#' Detect HDF5 File Type
#'
#' Determines if HDF5 file is organized by statistic or by contrast
#'
#' @param meta Metadata list from read_h5_metadata
#' @return Character string: "by_stat", "by_contrast", or "betas"
#' @keywords internal
detect_h5_file_type <- function(meta) {
  labels <- meta$labels
  
  # Check if labels match known statistic names
  stat_names <- c("beta", "se", "tstat", "pval", "fstat", "sigma", "estimate", "stat", "prob")
  if (all(labels %in% stat_names)) {
    return("by_stat")
  }
  
  # Check if this looks like raw betas (numeric or regressor names)
  if (all(grepl("^beta_|^reg_|Intercept", labels))) {
    return("betas")
  }
  
  # Otherwise assume it's contrasts
  return("by_contrast")
}

#' Extract Subject IDs from File Paths
#'
#' @param paths Character vector of file paths
#' @return Character vector of subject IDs
#' @keywords internal
extract_subject_ids_from_paths <- function(paths) {
  # Try to extract BIDS-style subject IDs
  subjects <- gsub(".*sub-([^_/]+).*", "sub-\\1", paths)
  
  # If that didn't work, just use file basenames
  if (all(subjects == paths)) {
    subjects <- tools::file_path_sans_ext(basename(paths))
    # Remove double extension for .nii.gz
    subjects <- sub("\\.nii$", "", subjects)
  }
  
  return(subjects)
}

#' Read Data Chunk from HDF5 Files
#'
#' Reads a subset of voxels from all subjects' HDF5 files
#'
#' @param gd A group_data_h5 object
#' @param voxel_indices Integer vector of voxel indices to read
#' @param stat Character vector of statistics to extract
#' @return List with one element per subject, each containing extracted data
#' @keywords internal
read_h5_chunk <- function(gd, voxel_indices, stat = NULL) {
  if (is.null(stat)) {
    stat <- gd$stat
  }
  if (gd$file_type != "by_stat") {
    stop("Only by_stat HDF5 layout is supported by read_h5_chunk", call. = FALSE)
  }

  n_subjects <- length(gd$paths)
  n_voxels <- length(voxel_indices)

  want_beta <- "beta" %in% stat
  want_se <- any(c("se", "var") %in% stat)
  want_t <- any(c("t", "tstat") %in% stat)
  beta_data <- if (want_beta) matrix(NA_real_, nrow = n_subjects, ncol = n_voxels) else NULL
  se_data <- if (want_se) matrix(NA_real_, nrow = n_subjects, ncol = n_voxels) else NULL
  t_data <- if (want_t) matrix(NA_real_, nrow = n_subjects, ncol = n_voxels) else NULL

  for (i in seq_along(gd$paths)) {
    vals <- .read_h5_subject(gd$paths[i], voxel_indices, want_beta, want_se, want_t)
    if (want_beta) beta_data[i, ] <- vals$beta
    if (want_se) se_data[i, ] <- vals$se
    if (want_t) t_data[i, ] <- vals$t
  }

  result <- list()
  if (!is.null(beta_data)) result$beta <- beta_data
  if (!is.null(se_data)) result$se <- se_data
  if (!is.null(t_data)) result$t <- t_data
  if (!is.null(gd$df)) result$df <- gd$df

  result
}

# Internal: close an fmristore handle (list-based or S4) exactly once.
.close_h5_handle <- function(h5_handle) {
  if (is.list(h5_handle)) {
    if (is.function(h5_handle$close_all)) return(invisible(try(h5_handle$close_all(), silent = TRUE)))
    if (is.function(h5_handle$close)) return(invisible(try(h5_handle$close(), silent = TRUE)))
    return(invisible(NULL))
  }
  if (isS4(h5_handle)) {
    obj <- tryCatch(methods::slot(h5_handle, "obj"), error = function(e) NULL)
    if (!is.null(obj) && is.function(obj$close_all)) {
      return(invisible(try(obj$close_all(), silent = TRUE)))
    }
  }
  invisible(try(close(h5_handle), silent = TRUE))
}

# Internal: read the requested statistics for one subject's HDF5 file. The
# handle is closed when this function exits (also on error), and a file that
# cannot be read, or lacks a requested statistic, is an error naming the file
# rather than a silent row of NA.
.read_h5_subject <- function(path, voxel_indices, want_beta, want_se, want_t) {
  fail <- function(msg) {
    stop("Failed to read HDF5 data from ", path, ": ", msg, call. = FALSE)
  }
  h5_handle <- tryCatch(fmristore::read_labeled_vec(path),
                        error = function(e) fail(conditionMessage(e)))
  on.exit(.close_h5_handle(h5_handle), add = TRUE)

  if (is.list(h5_handle) && !is.null(h5_handle$read_voxels)) {
    data_list <- tryCatch(
      h5_handle$read_voxels(voxel_indices, labels = c("beta", "se", "var", "t", "tstat")),
      error = function(e) fail(conditionMessage(e))
    )
    read_label <- function(lbl) {
      v <- data_list[[lbl]]
      if (is.null(v)) NULL else as.numeric(v)
    }
  } else if (isS4(h5_handle)) {
    labels_available <- tryCatch(as.character(methods::slot(h5_handle, "labels")),
                                 error = function(e) character(0))
    read_label <- function(lbl) {
      if (!(lbl %in% labels_available)) return(NULL)
      tryCatch(as.numeric(h5_handle@obj[["data"]][[lbl]][voxel_indices]),
               error = function(e) fail(conditionMessage(e)))
    }
  } else {
    fail("unsupported fmristore handle type")
  }

  out <- list()
  if (want_beta) {
    out$beta <- read_label("beta")
    if (is.null(out$beta)) fail("statistic 'beta' not found")
  }
  if (want_se) {
    se <- read_label("se")
    if (is.null(se)) {
      v <- read_label("var")
      se <- if (is.null(v)) NULL else sqrt(v)
    }
    if (is.null(se)) fail("neither 'se' nor 'var' found")
    out$se <- se
  }
  if (want_t) {
    tv <- read_label("tstat")
    if (is.null(tv)) tv <- read_label("t")
    if (is.null(tv)) fail("neither 't' nor 'tstat' found")
    out$t <- tv
  }
  for (nm in names(out)) {
    if (length(out[[nm]]) != length(voxel_indices)) {
      fail(sprintf("statistic '%s' returned %d values for %d voxels",
                   nm, length(out[[nm]]), length(voxel_indices)))
    }
  }
  out
}

#' Read All Data from HDF5 Files
#'
#' Reads complete data from all subjects' HDF5 files.
#' Warning: This can use a lot of memory for whole-brain data.
#'
#' @param gd A group_data_h5 object
#' @param stat Character vector of statistics to extract
#' @return Array with dimensions (voxels, subjects, statistics)
#' @export
read_h5_full <- function(gd, stat = NULL) {
  if (!inherits(gd, "group_data_h5")) {
    stop("Input must be a group_data_h5 object", call. = FALSE)
  }
  
  if (is.null(stat)) {
    stat <- gd$stat
  }
  
  n_subjects <- length(gd$paths)
  n_stats <- length(stat)
  n_voxels <- gd$n_voxels
  data_array <- array(NA_real_, dim = c(n_voxels, n_subjects, n_stats))

  chunk <- read_h5_chunk(gd, voxel_indices = seq_len(n_voxels), stat = stat)
  for (j in seq_along(stat)) {
    sname <- stat[j]
    if (sname == "var") {
      if (!is.null(chunk$se)) {
        data_array[, , j] <- t(chunk$se^2)
      }
    } else if (sname %in% names(chunk)) {
      data_array[, , j] <- t(chunk[[sname]])
    } else if (sname == "tstat" && !is.null(chunk$t)) {
      data_array[, , j] <- t(chunk$t)
    } else if (sname == "t" && !is.null(chunk$t)) {
      data_array[, , j] <- t(chunk$t)
    }
  }
  
  # Add dimension names
  dimnames(data_array) <- list(
    voxel = NULL,
    subject = gd$subjects,
    statistic = stat
  )
  
  return(data_array)
}

#' Validate group_data_h5 Object
#'
#' @param x A group_data_h5 object
#' @keywords internal
#' @noRd
validate_group_data_h5 <- function(x) {
  # Check for required fields
  required_fields <- c("paths", "subjects", "format", "dim", "labels")
  missing_fields <- setdiff(required_fields, names(x))
  if (length(missing_fields) > 0) {
    stop("Missing required fields for group_data_h5: ", 
         paste(missing_fields, collapse = ", "), call. = FALSE)
  }
  
  # Validate paths and subjects have same length
  if (length(x$paths) != length(x$subjects)) {
    stop("Length of paths and subjects must match", call. = FALSE)
  }
  
  invisible(TRUE)
}

#' Create Group Dataset from fmri_lm Objects
#'
#' Creates a group dataset directly from a list of fitted fmri_lm objects
#'
#' @param lm_list List of fmri_lm objects
#' @param contrast Character string specifying which contrast to extract
#' @param stat Character vector of statistics to extract
#' @param subjects Character vector of subject identifiers
#' @param covariates Data frame of subject-level covariates
#'
#' @return A group_data_h5 object (in-memory variant)
#' @export
group_data_from_fmrilm <- function(lm_list,
                                   contrast = NULL,
                                   stat = c("beta", "se"),
                                   subjects = NULL,
                                   covariates = NULL) {
  if (!isTRUE(getOption("fmrireg.suppress_deprecation", FALSE))) {
    .Deprecated("group_data", msg = paste("group_data_from_fmrilm() is deprecated.",
                                          "Save each fit with write_results() and use",
                                          "group_data(..., format = 'nifti') or format = 'h5' instead."))
  }

  if (!is.list(lm_list)) {
    stop("'lm_list' must be a list", call. = FALSE)
  }
  
  if (!all(sapply(lm_list, inherits, "fmri_lm"))) {
    stop("All elements of 'lm_list' must be fmri_lm objects", call. = FALSE)
  }
  
  # Extract subjects if not provided
  if (is.null(subjects)) {
    subjects <- paste0("subject_", seq_along(lm_list))
  }
  
  # Get dimensions from first object
  first_lm <- lm_list[[1]]
  spatial <- .fmri_dataset_mask_space(first_lm$dataset, "group data reconstruction")
  mask <- spatial$mask_array
  space <- spatial$space
  dim <- c(dim(space), length(stat))
  
  # Create in-memory group_data object
  structure(
    list(
      lm_list = lm_list,
      subjects = subjects,
      covariates = covariates,
      contrast = contrast,
      stat = stat,
      dim = dim,
      format = "fmrilm",
      n_voxels = length(mask),
      mask = mask,
      space = space
    ),
    class = c("group_data_fmrilm", "group_data_h5", "group_data")
  )
}
