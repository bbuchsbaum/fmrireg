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

    # S4 or other handle type: attempt to retrieve via slots or accessors
    if (methods::is(h5_handle, "S4")) {
      sn <- try(methods::slotNames(h5_handle), silent = TRUE)
      get_slot <- function(nm) try(methods::slot(h5_handle, nm), silent = TRUE)
      dim_val <- NULL; labels_val <- NULL; mask_val <- NULL; mask_dim <- NULL
      if (!inherits(sn, "try-error")) {
        if ("dim" %in% sn) dim_val <- get_slot("dim")
        if ("labels" %in% sn) labels_val <- get_slot("labels")
        if ("mask" %in% sn) mask_val <- get_slot("mask")
        if ("mask_dim" %in% sn) mask_dim <- get_slot("mask_dim")
      }
      # Fallback to potential accessor generics if available
      if (is.null(dim_val)) dim_val <- try(dim(h5_handle), silent = TRUE)
      if (is.null(labels_val)) labels_val <- try(labels(h5_handle), silent = TRUE)
      if (is.null(mask_val)) mask_val <- tryCatch(get_mask(h5_handle), error = function(e) NULL)

      if (!inherits(dim_val, "try-error") && !inherits(labels_val, "try-error") && !is.null(dim_val) && !is.null(labels_val)) {
        return(list(
          dim = dim_val,
          labels = labels_val,
          mask_dim = if (!inherits(mask_dim, "try-error")) mask_dim else NULL,
          has_mask = isTRUE(!inherits(mask_val, "try-error") && !is.null(mask_val))
        ))
      }

      # Could not safely introspect S4 handle—informative error for tests to skip
      stop("Legacy HDF5 reader cannot introspect fmristore S4 handle; use group_data(..., format='h5') via fmrigds instead.")
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
  
  # Read from each subject's file
  data_list <- lapply(seq_along(gd$paths), function(i) {
    path <- gd$paths[i]
    
    tryCatch({
      h5_handle <- fmristore::read_labeled_vec(path)
      on.exit(h5_handle$close_all())
      
      # Read specific voxels and statistics
      if (gd$file_type == "by_stat") {
        # File has statistics as labels
        data <- h5_handle$read_voxels(voxel_indices, labels = stat)
      } else if (gd$file_type == "by_contrast") {
        # File has contrasts as labels, need to extract specific contrast
        if (!is.null(gd$contrast)) {
          data <- h5_handle$read_voxels(voxel_indices, labels = gd$contrast)
          # Now extract the statistics from within the contrast
          # This depends on how the data is structured
        } else {
          stop("Contrast must be specified for by_contrast files", call. = FALSE)
        }
      } else {
        # Raw betas file
        data <- h5_handle$read_voxels(voxel_indices)
      }
      
      return(data)
    }, error = function(e) {
      warning("Failed to read from ", path, ": ", e$message)
      return(NULL)
    })
  })
  
  names(data_list) <- gd$subjects
  return(data_list)
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
  
  # Pre-allocate array
  data_array <- array(NA_real_, dim = c(n_voxels, n_subjects, n_stats))
  
  # Read each subject's data
  for (i in seq_along(gd$paths)) {
    path <- gd$paths[i]
    
    tryCatch({
      h5_handle <- fmristore::read_labeled_vec(path)
      on.exit(h5_handle$close_all())
      
      # Read all voxels for specified statistics
      if (gd$file_type == "by_stat") {
        for (j in seq_along(stat)) {
          data_array[, i, j] <- h5_handle$read_data(labels = stat[j])
        }
      } else {
        # Handle other file types
        # Implementation depends on exact structure
      }
    }, error = function(e) {
      warning("Failed to read from ", path, ": ", e$message)
    })
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
    .Deprecated("group_data", msg = "group_data_from_fmrilm() is deprecated. Use group_data(lm_list, format='fmrilm', ...) instead.")
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
  mask <- fmridataset::get_mask(first_lm$dataset)
  space <- neuroim2::space(mask)
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
