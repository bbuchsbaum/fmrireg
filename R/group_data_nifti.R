# Group Data Support for NIfTI Files
# Handles standard NIfTI neuroimaging files for group analysis

#' Create Group Dataset from NIfTI Files
#'
#' Creates a group dataset from NIfTI files containing effect sizes and their
#' standard errors or variances. Supports various input configurations including
#' beta/SE pairs, beta/variance pairs, or t-statistics with degrees of freedom.
#'
#' @param beta_paths Character vector of paths to beta/effect size NIfTI files
#' @param se_paths Character vector of paths to standard error NIfTI files
#' @param var_paths Character vector of paths to variance NIfTI files (alternative to se_paths)
#' @param t_paths Character vector of paths to t-statistic NIfTI files
#' @param df Degrees of freedom (scalar or vector). Required if using t_paths.
#' @param subjects Character vector of subject identifiers. If NULL, extracted from file paths.
#' @param covariates Data frame of subject-level covariates (optional)
#' @param mask Path to mask NIfTI file or mask object (optional but recommended)
#' @param target_space Path to template NIfTI for spatial alignment checking
#' @param validate Logical. Validate that all files exist and have matching dimensions (default: TRUE)
#'
#' @return A group_data_nifti object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # From FSL FEAT output (COPE and VARCOPE files)
#' gd <- group_data_from_nifti(
#'   beta_paths = Sys.glob("feat_output/sub-*/cope1.nii.gz"),
#'   var_paths = Sys.glob("feat_output/sub-*/varcope1.nii.gz"),
#'   subjects = sprintf("sub-%02d", 1:20),
#'   mask = "group_mask.nii.gz"
#' )
#'
#' # From SPM contrast images
#' gd <- group_data_from_nifti(
#'   beta_paths = Sys.glob("SPM/sub*/con_0001.nii"),
#'   se_paths = Sys.glob("SPM/sub*/se_0001.nii"),
#'   mask = "SPM/mask.nii"
#' )
#'
#' # From t-statistics only (for Stouffer's Z combination)
#' gd <- group_data_from_nifti(
#'   t_paths = Sys.glob("stats/sub-*/tstat1.nii.gz"),
#'   df = 100,  # Or vector of per-subject df
#'   mask = "group_mask.nii.gz"
#' )
#' }
group_data_from_nifti <- function(beta_paths = NULL,
                                  se_paths = NULL,
                                  var_paths = NULL,
                                  t_paths = NULL,
                                  df = NULL,
                                  subjects = NULL,
                                  covariates = NULL,
                                  mask = NULL,
                                  target_space = NULL,
                                  validate = TRUE) {
  
  # Validate input combinations
  has_beta <- !is.null(beta_paths)
  has_se <- !is.null(se_paths)
  has_var <- !is.null(var_paths)
  has_t <- !is.null(t_paths)
  
  if (!has_beta && !has_t) {
    stop("Must provide either 'beta_paths' or 't_paths'", call. = FALSE)
  }
  
  if (has_beta && !has_se && !has_var) {
    stop("When providing 'beta_paths', must also provide either 'se_paths' or 'var_paths'", 
         call. = FALSE)
  }
  
  if (has_se && has_var) {
    stop("Provide either 'se_paths' or 'var_paths', not both", call. = FALSE)
  }
  
  if (has_t && is.null(df)) {
    stop("When providing 't_paths', must also provide 'df' (degrees of freedom)", 
         call. = FALSE)
  }
  
  # Determine primary paths for validation
  primary_paths <- if (has_beta) beta_paths else t_paths
  n_subjects <- length(primary_paths)
  
  # Validate path lengths match
  if (has_se && length(se_paths) != n_subjects) {
    stop("Length of 'se_paths' must match 'beta_paths'", call. = FALSE)
  }
  if (has_var && length(var_paths) != n_subjects) {
    stop("Length of 'var_paths' must match 'beta_paths'", call. = FALSE)
  }
  if (has_t && length(t_paths) != n_subjects) {
    stop("All path vectors must have the same length", call. = FALSE)
  }
  
  # Validate files exist
  if (validate) {
    all_paths <- c(beta_paths, se_paths, var_paths, t_paths)
    all_paths <- all_paths[!is.null(all_paths)]
    missing_files <- all_paths[!file.exists(all_paths)]
    
    if (length(missing_files) > 0) {
      stop("The following files do not exist:\n  ",
           paste(head(missing_files, 5), collapse = "\n  "),
           if (length(missing_files) > 5) paste0("\n  ... and ", length(missing_files) - 5, " more"),
           call. = FALSE)
    }
  }
  
  # Extract or validate subjects
  if (is.null(subjects)) {
    subjects <- extract_subject_ids_from_paths(primary_paths)
  } else {
    if (length(subjects) != n_subjects) {
      stop("Length of 'subjects' must match number of input files", call. = FALSE)
    }
  }
  
  # Process covariates
  if (!is.null(covariates)) {
    if (!is.data.frame(covariates)) {
      stop("'covariates' must be a data frame", call. = FALSE)
    }
    if (nrow(covariates) != n_subjects) {
      stop("Number of rows in 'covariates' must match number of subjects", call. = FALSE)
    }
  }
  
  # Process degrees of freedom
  if (!is.null(df)) {
    if (length(df) == 1) {
      df <- rep(df, n_subjects)
    } else if (length(df) != n_subjects) {
      stop("'df' must be a scalar or have length equal to number of subjects", call. = FALSE)
    }
  }
  
  # Get dimensions from first file
  first_file <- primary_paths[1]
  img_info <- read_nifti_header(first_file)
  
  # Validate all files have matching dimensions
  if (validate) {
    validate_nifti_dimensions(c(beta_paths, se_paths, var_paths, t_paths), img_info$dim)
  }
  
  # Load or validate mask
  mask_data <- NULL
  if (!is.null(mask)) {
    if (is.character(mask)) {
      mask_data <- read_nifti_data(mask)
      # Validate mask dimensions match data
      if (!all(dim(mask_data) == img_info$dim[1:3])) {
        stop("Mask dimensions do not match data dimensions", call. = FALSE)
      }
    } else if (inherits(mask, "NeuroVol")) {
      mask_data <- mask
    }
  }
  
  # Create group_data_nifti object
  structure(
    list(
      beta_paths = beta_paths,
      se_paths = se_paths,
      var_paths = var_paths,
      t_paths = t_paths,
      df = df,
      subjects = subjects,
      covariates = covariates,
      mask = mask,
      mask_data = mask_data,
      dim = img_info$dim,
      voxel_size = img_info$voxel_size,
      format = "nifti",
      n_subjects = n_subjects,
      n_voxels = prod(img_info$dim[1:3])
    ),
    class = c("group_data_nifti", "group_data")
  )
}

#' Read NIfTI Header Information
#'
#' @param path Path to NIfTI file
#' @return List with header information
#' @keywords internal
read_nifti_header <- function(path) {
  tryCatch({
    # Try RNifti first if available
    if (requireNamespace("RNifti", quietly = TRUE)) {
      hdr <- RNifti::niftiHeader(path)
      return(list(
        dim = hdr$dim[2:4],  # Spatial dimensions
        voxel_size = hdr$pixdim[2:4],
        datatype = hdr$datatype
      ))
    }
    
    # Fall back to neuroim2
    img <- neuroim2::read_vol(path)
    return(list(
      dim = dim(img),
      voxel_size = neuroim2::spacing(img),
      datatype = NULL
    ))
  }, error = function(e) {
    stop("Failed to read NIfTI header from ", path, ": ", e$message, call. = FALSE)
  })
}

#' Read NIfTI Data
#'
#' @param path Path to NIfTI file
#' @param mask Optional mask to apply
#' @return Array or vector of data
#' @keywords internal
read_nifti_data <- function(path, mask = NULL) {
  # Robust NIfTI reader: try neuroim2, then RNifti
  # Returns a full array or masked vector (if mask provided)
  #
  # path: character path to NIfTI
  # mask: optional array/NeuroVol; if provided, returns vector of masked values
  
  # Helper to apply mask or return full array
  apply_mask <- function(arr, mask) {
    if (!is.null(mask)) {
      return(as.vector(arr[mask > 0]))
    }
    arr
  }
  
  # Attempt read with neuroim2
  res <- try(suppressWarnings(neuroim2::read_vol(path)), silent = TRUE)
  if (!inherits(res, "try-error")) {
    arr <- try(as.array(res), silent = TRUE)
    if (!inherits(arr, "try-error") && is.atomic(arr)) {
      return(apply_mask(arr, mask))
    }
  }
  
  # Fallback to RNifti if available
  if (requireNamespace("RNifti", quietly = TRUE)) {
    arr <- try(RNifti::readNifti(path, internal = FALSE), silent = TRUE)
    if (!inherits(arr, "try-error")) {
      # Ensure base type and drop singleton fourth dim if present
      arr <- as.array(arr)
      if (length(dim(arr)) == 4 && dim(arr)[4] == 1) {
        dim(arr) <- dim(arr)[1:3]
      }
      return(apply_mask(arr, mask))
    }
  }
  
  stop("Failed to read NIfTI data from ", path, ": could not parse with neuroim2 or RNifti", call. = FALSE)
}

#' Validate NIfTI Dimensions
#'
#' @param paths Character vector of NIfTI file paths
#' @param expected_dim Expected dimensions
#' @keywords internal
#' @noRd
validate_nifti_dimensions <- function(paths, expected_dim) {
  paths <- paths[!is.null(paths)]
  
  for (path in paths) {
    if (file.exists(path)) {
      hdr <- read_nifti_header(path)
      if (!all(hdr$dim == expected_dim)) {
        stop("Dimension mismatch in file ", path, 
             "\n  Expected: ", paste(expected_dim, collapse = " x "),
             "\n  Found: ", paste(hdr$dim, collapse = " x "),
             call. = FALSE)
      }
    }
  }
  
  invisible(TRUE)
}

#' Read Data Chunk from NIfTI Files
#'
#' Reads a subset of voxels from all subjects' NIfTI files
#'
#' @param gd A group_data_nifti object
#' @param voxel_indices Integer vector of voxel indices to read
#' @return Matrix with dimensions (subjects, voxels)
#' @keywords internal
read_nifti_chunk <- function(gd, voxel_indices) {
  n_subjects <- gd$n_subjects
  n_voxels <- length(voxel_indices)
  
  # Initialize output matrices
  beta_data <- NULL
  se_data <- NULL
  var_data <- NULL
  t_data <- NULL
  
  # Read beta data
  if (!is.null(gd$beta_paths)) {
    beta_data <- matrix(NA_real_, nrow = n_subjects, ncol = n_voxels)
    for (i in seq_len(n_subjects)) {
      arr <- read_nifti_data(gd$beta_paths[i])
      beta_data[i, ] <- as.vector(arr)[voxel_indices]
    }
  }
  
  # Read SE data
  if (!is.null(gd$se_paths)) {
    se_data <- matrix(NA_real_, nrow = n_subjects, ncol = n_voxels)
    for (i in seq_len(n_subjects)) {
      arr <- read_nifti_data(gd$se_paths[i])
      se_data[i, ] <- as.vector(arr)[voxel_indices]
    }
  }
  
  # Read variance data and convert to SE
  if (!is.null(gd$var_paths)) {
    se_data <- matrix(NA_real_, nrow = n_subjects, ncol = n_voxels)
    for (i in seq_len(n_subjects)) {
      arr <- read_nifti_data(gd$var_paths[i])
      var_vals <- as.vector(arr)[voxel_indices]
      se_data[i, ] <- sqrt(var_vals)
    }
  }
  
  # Read t-statistic data
  if (!is.null(gd$t_paths)) {
    t_data <- matrix(NA_real_, nrow = n_subjects, ncol = n_voxels)
    for (i in seq_len(n_subjects)) {
      arr <- read_nifti_data(gd$t_paths[i])
      t_data[i, ] <- as.vector(arr)[voxel_indices]
    }
  }
  
  # Return list of matrices
  result <- list()
  if (!is.null(beta_data)) result$beta <- beta_data
  if (!is.null(se_data)) result$se <- se_data
  if (!is.null(t_data)) result$t <- t_data
  if (!is.null(gd$df)) result$df <- gd$df
  
  return(result)
}

#' Read All Data from NIfTI Files
#'
#' Reads complete data from all subjects' NIfTI files using memory mapping when possible.
#'
#' @param gd A group_data_nifti object
#' @param use_mask Logical. Apply mask to data (default: TRUE if mask exists)
#' @return List with data matrices
#' @export
read_nifti_full <- function(gd, use_mask = NULL) {
  if (!inherits(gd, "group_data_nifti")) {
    stop("Input must be a group_data_nifti object", call. = FALSE)
  }
  
  # Determine whether to use mask
  if (is.null(use_mask)) {
    use_mask <- !is.null(gd$mask_data)
  }
  
  # Get voxel indices to read
  if (use_mask && !is.null(gd$mask_data)) {
    voxel_indices <- which(as.vector(gd$mask_data) > 0)
  } else {
    voxel_indices <- seq_len(gd$n_voxels)
  }
  
  # Read all data
  return(read_nifti_chunk(gd, voxel_indices))
}

#' Validate group_data_nifti Object
#'
#' @param x A group_data_nifti object
#' @keywords internal
#' @noRd
validate_group_data_nifti <- function(x) {
  # Check for at least one data source
  if (is.null(x$beta_paths) && is.null(x$t_paths)) {
    stop("group_data_nifti must have either beta_paths or t_paths", call. = FALSE)
  }
  
  # If beta paths exist, must have SE or variance
  if (!is.null(x$beta_paths)) {
    if (is.null(x$se_paths) && is.null(x$var_paths)) {
      stop("When beta_paths are provided, must also have se_paths or var_paths", 
           call. = FALSE)
    }
  }
  
  # If t paths exist, must have df
  if (!is.null(x$t_paths) && is.null(x$df)) {
    stop("When t_paths are provided, must also have df", call. = FALSE)
  }
  
  invisible(TRUE)
}

#' Convert t-statistics to Effect Sizes
#'
#' Helper function to convert t-statistics and df to beta and SE estimates
#'
#' @param t T-statistic values
#' @param df Degrees of freedom
#' @param n Sample size (optional, improves SE estimation)
#' @return List with beta and se estimates
#' @export
t_to_beta_se <- function(t, df, n = NULL) {
  # Standard error of the t-statistic
  se_t <- 1
  
  # If sample size is known, can get better SE estimate
  if (!is.null(n)) {
    se_t <- sqrt(1/n)
  }
  
  # Convert t to effect size (Cohen's d approximation)
  # beta = t * SE
  beta <- t * se_t
  
  # Standard error from t-distribution
  se <- se_t
  
  return(list(beta = beta, se = se))
}
