# BIDS Export Functions using LabeledVolumeSet
# Functions for saving fMRI analysis results in BIDS-compliant format
#' @importFrom jsonlite write_json
NULL

#' Write Results from fMRI Linear Model
#'
#' Exports statistical maps from an fmri_lm object to HDF5 files with BIDS-compliant
#' naming and JSON metadata sidecars using the fmristore LabeledVolumeSet infrastructure.
#'
#' @param x An fmri_lm object containing fitted model results
#' @param path Output directory path. If NULL, uses current working directory
#' @param subject Subject identifier (e.g., "01", "1001"). Required.
#' @param task Task identifier (e.g., "nback", "rest"). Required for BIDS compliance.
#' @param space Spatial reference (e.g., "MNI152NLin2009cAsym"). Optional but recommended.
#' @param desc Description of the analysis (default: "GLM")
#' @param format Character vector specifying output formats. Supported values
#'   are `"h5"` (default BIDS/HDF5 export) and `"gds"` (fmrigds-compatible
#'   assays plus an `.rds` plan describing the dataset). You may request both.
#' @param strategy Storage strategy: "by_stat" (group contrasts by statistic) or "by_contrast" (separate files)
#' @param save_betas Logical. Save raw regressor betas (default: TRUE)
#' @param contrasts Character vector of contrast names to save. NULL saves all contrasts
#' @param contrast_stats Character vector of contrast statistics to save (default: c("beta", "tstat", "pval", "se"))
#' @param overwrite Logical. Overwrite existing files (default: FALSE)
#' @param validate_inputs Logical. Validate fmrilm object structure (default: TRUE)
#' @param ... Additional arguments passed to internal functions
#'
#' @return Invisible list of file paths created
#' @export
#'
#' @examples
#' \dontrun{
#' # Save all results using default settings
#' write_results(fitted_model, subject = "01", task = "nback")
#' 
#' # Save only specific contrasts and statistics  
#' write_results(fitted_model, 
#'               subject = "01", task = "nback", space = "MNI152NLin2009cAsym",
#'               contrasts = c("FacesVsPlaces", "GoVsNoGo"),
#'               contrast_stats = c("beta", "tstat"))
#' }
#' @method write_results fmri_lm
write_results.fmri_lm <- function(x,
                                  path = NULL,
                                  subject = NULL,
                                  task = NULL, 
                                  space = NULL,
                                  desc = "GLM",
                                  format = c("h5"),
                                  strategy = c("by_stat", "by_contrast"),
                                  save_betas = TRUE,
                                  contrasts = NULL,
                                  contrast_stats = c("beta", "tstat", "pval", "se"),
                                  overwrite = FALSE,
                                  validate_inputs = TRUE,
                                  ...) {
  
  strategy <- match.arg(strategy)
  format <- match.arg(format, c("h5", "gds"), several.ok = TRUE)
  if (length(format) == 0) format <- "h5"
  produce_h5 <- "h5" %in% format
  produce_gds <- "gds" %in% format
  
  # Input validation
  if (validate_inputs) {
    .validate_fmrilm_object(x)
  }
  
  # Set default path
  if (is.null(path)) {
    path <- getwd()
  }
  
  # Create BIDS entities list
  entities <- .create_bids_entities(subject, task, space)
  
  # Validate required entities
  .validate_required_entities(entities)
  
  predicted_files <- .predict_output_files(
    path = path,
    entities = entities,
    desc = desc,
    format = format,
    strategy = strategy,
    save_betas = save_betas,
    contrasts = contrasts,
    contrast_stats = contrast_stats,
    fmrilm_obj = x
  )
  
  if (!overwrite) {
    existing_files <- predicted_files[file.exists(predicted_files)]
    if (length(existing_files) > 0) {
      stop(
        "Output files already exist: ",
        paste(existing_files, collapse = ", "),
        ". Set overwrite = TRUE to replace them.",
        call. = FALSE
      )
    }
  }
  
  # Create output directory if needed
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  # Use atomic write pattern with temporary directory
  temp_dir <- .create_temp_write_dir(path)
  finalized <- NULL
  
  tryCatch({
    created_files <- list()
    
    # Save raw regressor betas
    if (produce_h5 && save_betas) {
      beta_files <- .save_regressor_betas(x, temp_dir, entities, desc, overwrite)
      created_files <- c(created_files, list(betas = beta_files))
    }
    
    # Save contrast results
    if (produce_h5) {
      if (strategy == "by_stat") {
        contrast_files <- .save_contrasts_by_stat(x, temp_dir, entities, desc, 
                                                  contrasts, contrast_stats, overwrite)
      } else {
        contrast_files <- .save_contrasts_by_contrast(x, temp_dir, entities, desc,
                                                      contrasts, contrast_stats, overwrite)
      }
      created_files <- c(created_files, contrast_files)
    }

    if (produce_gds) {
      if (!requireNamespace("fmrigds", quietly = TRUE)) {
        stop("fmrigds package is required for GDS export. Please install fmrigds.", call. = FALSE)
      }
      gds_files <- .save_gds_outputs(x, temp_dir, entities, desc)
      created_files <- c(created_files, list(gds = gds_files))
    }
    
    # Atomic move: if all writes succeeded, move temp_dir contents to final location
    finalized <- .finalize_atomic_write(temp_dir, path, created_files, overwrite = overwrite)
    created_files <- finalized$files
    
    if (produce_gds && "gds" %in% names(created_files)) {
      .write_gds_plan_file(
        gds_h5_final = created_files$gds$h5,
        gds_plan_final = created_files$gds$plan,
        overwrite = overwrite,
        backup_root = file.path(temp_dir, ".overwrite_backup")
      )
    }
    
    .cleanup_temp_write_dir(temp_dir)
    
    invisible(created_files)
    
  }, error = function(e) {
    if (!is.null(finalized)) {
      .rollback_finalized_write(finalized$transaction)
    }
    .cleanup_temp_write_dir(temp_dir)
    stop("Failed to write BIDS results: ", e$message, call. = FALSE)
  })
}

#' Predict Output Files for Overwrite Check
#' @keywords internal
#' @noRd
.predict_output_files <- function(path, entities, desc, format, strategy, save_betas, contrasts, contrast_stats, fmrilm_obj) {
  predicted_files <- c()
  produce_h5 <- "h5" %in% format
  produce_gds <- "gds" %in% format
  
  # Predict beta files
  if (produce_h5 && save_betas) {
    beta_h5 <- .generate_bids_filename(entities, desc = desc, suffix = "betas", extension = "h5")
    beta_json <- .generate_bids_filename(entities, desc = desc, suffix = "betas", extension = "json")
    predicted_files <- c(predicted_files, 
                        file.path(path, beta_h5), 
                        file.path(path, beta_json))
  }
  
  # Predict contrast files
  if (produce_h5 && !is.null(fmrilm_obj$result$contrasts) && nrow(fmrilm_obj$result$contrasts) > 0) {
    available_contrasts <- fmrilm_obj$result$contrasts
    if (!is.null(contrasts)) {
      available_contrasts <- available_contrasts[available_contrasts$name %in% contrasts, ]
    }
    
    if (nrow(available_contrasts) > 0) {
      if (strategy == "by_stat") {
        # Files grouped by statistic
        for (stat in contrast_stats) {
          has_stat <- any(vapply(available_contrasts$data, function(entry) {
            .map_stat_name(stat) %in% names(entry)
          }, logical(1)))
          
          if (has_stat) {
            stat_desc <- paste0(desc, "statmap")
            filename_h5 <- .generate_bids_filename(entities, desc = stat_desc, stat = stat, suffix = "bold", extension = "h5")
            filename_json <- .generate_bids_filename(entities, desc = stat_desc, stat = stat, suffix = "bold", extension = "json")
            predicted_files <- c(predicted_files, 
                                file.path(path, filename_h5),
                                file.path(path, filename_json))
          }
        }
      } else {
        # Files grouped by contrast
        for (i in seq_len(nrow(available_contrasts))) {
          available_stats <- vapply(contrast_stats, function(stat) {
            .map_stat_name(stat) %in% names(available_contrasts$data[[i]])
          }, logical(1))
          
          if (any(available_stats)) {
            contrast_name <- available_contrasts$name[i]
            stat_desc <- paste0(desc, "statmap")
            filename_h5 <- .generate_bids_filename(entities, desc = stat_desc, contrast = contrast_name, suffix = "bold", extension = "h5")
            filename_json <- .generate_bids_filename(entities, desc = stat_desc, contrast = contrast_name, suffix = "bold", extension = "json")
            predicted_files <- c(predicted_files, 
                                file.path(path, filename_h5),
                                file.path(path, filename_json))
          }
        }
      }
    }
  }
  
  if (produce_gds) {
    gds_h5 <- .generate_bids_filename(entities, desc = desc, suffix = "gds", extension = "h5")
    gds_plan <- .generate_bids_filename(entities, desc = desc, suffix = "gds", extension = "rds")
    predicted_files <- c(predicted_files, file.path(path, gds_h5), file.path(path, gds_plan))
  }
  
  unique(predicted_files)
}

#' Validate fmri_lm Object Structure
#' @keywords internal
#' @noRd
.validate_fmrilm_object <- function(x) {
  if (!inherits(x, "fmri_lm")) {
    stop("Input must be an 'fmri_lm' object", call. = FALSE)
  }
  
  if (is.null(x$result)) {
    stop("fmri_lm object has not been fitted (missing $result component)", call. = FALSE)
  }
  
  if (is.null(x$result$betas)) {
    stop("fmri_lm object missing beta estimates", call. = FALSE)
  }
  
  if (is.null(x$model)) {
    stop("fmri_lm object missing model specification", call. = FALSE)
  }
  
  if (is.null(x$dataset)) {
    stop("fmri_lm object missing dataset information", call. = FALSE)
  }
}

#' Create BIDS Entities List
#' @keywords internal
#' @noRd
.create_bids_entities <- function(subject, task, space) {
  
  entities <- list()
  
  # Subject: sanitize and store
  if (!is.null(subject)) {
    entities$subject <- .sanitize_label(subject)
  }
  
  # Task: sanitize and store
  if (!is.null(task)) {
    entities$task <- .sanitize_label(task)
  }
  
  # Space: sanitize and store
  if (!is.null(space)) {
    entities$space <- .sanitize_label(space)
  }
  
  return(entities)
}

#' Sanitize Labels for BIDS Compliance
#' @keywords internal
#' @noRd
.sanitize_label <- function(label) {
  if (is.null(label) || is.na(label)) return(NULL)
  
  # Convert to character and remove invalid characters
  label <- as.character(label)
  # BIDS labels must be alphanumeric only (no hyphens allowed in labels)
  label <- gsub("[^a-zA-Z0-9]", "", label)
  
  if (nchar(label) == 0) {
    warning("Label sanitization resulted in empty string")
    return(NULL)  # Return NULL instead of empty string to properly trigger validation
  }
  
  return(label)
}

#' Validate Required BIDS Entities
#' @keywords internal
#' @noRd
.validate_required_entities <- function(entities) {
  if (is.null(entities$subject)) {
    stop("Subject identifier is required. Please provide 'subject' parameter.", 
         call. = FALSE)
  }
  
  if (is.null(entities$task)) {
    stop("Task identifier is required for BIDS compliance. Please provide 'task' parameter.",
         call. = FALSE)
  }
  
  # Space is optional but recommended
  if (is.null(entities$space)) {
    warning("Spatial reference not specified. Consider providing 'space' parameter for better BIDS compliance.")
  }
}

#' Generate BIDS-Compliant Filename
#' @keywords internal
#' @noRd
.generate_bids_filename <- function(entities, desc = NULL, contrast = NULL, 
                                   stat = NULL, suffix = NULL, extension = "h5") {
  
  parts <- c()
  
  # Required/standard entities in BIDS order
  if (!is.null(entities$subject)) {
    parts <- c(parts, paste0("sub-", entities$subject))
  }
  
  if (!is.null(entities$task)) {
    parts <- c(parts, paste0("task-", entities$task))
  }
  
  if (!is.null(entities$space)) {
    parts <- c(parts, paste0("space-", entities$space))
  }
  
  if (!is.null(contrast)) {
    parts <- c(parts, paste0("contrast-", .sanitize_label(contrast)))
  }
  
  if (!is.null(stat)) {
    parts <- c(parts, paste0("stat-", stat))
  }
  
  if (!is.null(desc)) {
    parts <- c(parts, paste0("desc-", desc))
  }
  
  # Add suffix
  if (!is.null(suffix)) {
    parts <- c(parts, suffix)
  }
  
  # Combine parts
  filename <- paste(parts, collapse = "_")
  
  # Add extension
  if (!is.null(extension)) {
    filename <- paste0(filename, ".", extension)
  }
  
  return(filename)
}

#' Save Raw Regressor Betas using LabeledVolumeSet
#' @keywords internal
#' @noRd
.save_regressor_betas <- function(fmrilm_obj, path, entities, desc, overwrite) {
  
  # Compute beta volumes
  beta_data <- .compute_beta_volumes(fmrilm_obj)
  
  # Get regressor names from the actual beta data to ensure they match
  # First get the raw beta matrix from coef() to access column names
  beta_matrix <- coef(fmrilm_obj, type = "betas", include_baseline = TRUE)
  
  # Beta data should now always contain all regressors (event + baseline)
  n_beta_cols <- dim(beta_data)[4]  # 4th dimension is number of regressors
  
  # Use the column names from the beta matrix if available
  if (!is.null(colnames(beta_matrix))) {
    regressor_names <- colnames(beta_matrix)
  } else {
    # Fallback: generate generic names
    regressor_names <- paste0("beta_", seq_len(n_beta_cols))
  }
  
  # Ensure we have the right number of names
  if (length(regressor_names) != n_beta_cols) {
    warning("Regressor names don't match beta dimensions. Using generic names.")
    regressor_names <- paste0("beta_", seq_len(n_beta_cols))
  }
  
  spatial <- .fmri_dataset_mask_space(fmrilm_obj$dataset, "raw beta export")
  mask_array <- spatial$mask_array
  space <- spatial$space
  
  # Generate filename for raw betas  
  filename <- .generate_bids_filename(entities, desc = desc, suffix = "betas", extension = "h5")
  filepath <- file.path(path, filename)
  
  # Save using fmristore infrastructure
  tryCatch({
    # Create 3D mask volume for fmristore
    mask_vol <- neuroim2::LogicalNeuroVol(mask_array, space)
    h5_handle <- fmristore::write_labeled_vec(beta_data, mask_vol, regressor_names, file = filepath)
    h5_handle$close_all()
  }, error = function(e) {
    stop("Failed to write LabeledVolumeSet: ", e$message, call. = FALSE)
  })
  
  # Generate and save JSON metadata
  json_filepath <- .save_betas_json_metadata(fmrilm_obj, path, entities, desc, regressor_names, filepath)
  
  return(list(h5 = filepath, json = json_filepath))
}

#' Save Contrasts by Statistic Type using LabeledVolumeSet
#' @keywords internal
#' @noRd
.save_contrasts_by_stat <- function(fmrilm_obj, path, entities, desc, contrasts, contrast_stats, overwrite) {
  
  if (is.null(fmrilm_obj$result$contrasts) || nrow(fmrilm_obj$result$contrasts) == 0) {
    warning("No contrasts found in fmrilm object")
    return(list())
  }
  
  # Filter contrasts if specified
  available_contrasts <- fmrilm_obj$result$contrasts
  if (!is.null(contrasts)) {
    available_contrasts <- available_contrasts[available_contrasts$name %in% contrasts, ]
    if (nrow(available_contrasts) == 0) {
      warning("None of the specified contrasts found in model")
      return(list())
    }
  }
  
  # Get spatial information
  spatial <- .fmri_dataset_mask_space(fmrilm_obj$dataset, "contrast export")
  mask_array <- spatial$mask_array
  space <- spatial$space
  brain_dims <- dim(mask_array)
  
  # Validate mask dimensionality
  if (length(mask_array) != prod(brain_dims)) {
    stop("Mask length (", length(mask_array), 
         ") does not match product of space dims (", prod(brain_dims), ").", 
         call. = FALSE)
  }
  
  mask_vol <- neuroim2::LogicalNeuroVol(mask_array, space)
  
  created_files <- list()
  
  # Save each statistic type separately
  for (stat in contrast_stats) {
    stat_data <- .compute_statistical_volumes(available_contrasts, stat, brain_dims, mask_array, space)
    
    if (!is.null(stat_data)) {
      # Generate filename with desc-statmap suffix (BIDS compliant)
      stat_desc <- paste0(desc, "statmap")
      filename <- .generate_bids_filename(entities, desc = stat_desc, stat = stat, suffix = "bold")
      filepath <- file.path(path, filename)
      
      h5_handle <- fmristore::write_labeled_vec(stat_data$neurovec, mask_vol, stat_data$contrast_names, file = filepath)
      h5_handle$close_all()
      
      # Save JSON metadata
      json_filepath <- .save_contrasts_json_metadata(fmrilm_obj, path, entities, desc, 
                                                     stat, stat_data$contrast_names, available_contrasts, filepath)
      
      created_files[[stat]] <- list(h5 = filepath, json = json_filepath)
    }
  }
  
  return(created_files)
}

#' Save Contrasts by Individual Contrast using LabeledVolumeSet
#' @keywords internal  
#' @noRd
.save_contrasts_by_contrast <- function(fmrilm_obj, path, entities, desc, contrasts, contrast_stats, overwrite) {
  
  if (is.null(fmrilm_obj$result$contrasts) || nrow(fmrilm_obj$result$contrasts) == 0) {
    warning("No contrasts found in fmrilm object")
    return(list())
  }
  
  # Filter contrasts if specified
  available_contrasts <- fmrilm_obj$result$contrasts
  if (!is.null(contrasts)) {
    available_contrasts <- available_contrasts[available_contrasts$name %in% contrasts, ]
    if (nrow(available_contrasts) == 0) {
      warning("None of the specified contrasts found in model")
      return(list())
    }
  }
  
  # Get spatial information
  spatial <- .fmri_dataset_mask_space(fmrilm_obj$dataset, "contrast export")
  mask_array <- spatial$mask_array
  space <- spatial$space
  brain_dims <- dim(mask_array)
  
  # Validate mask dimensionality
  if (length(mask_array) != prod(brain_dims)) {
    stop("Mask length (", length(mask_array), 
         ") does not match product of space dims (", prod(brain_dims), ").", 
         call. = FALSE)
  }
  
  mask_vol <- neuroim2::LogicalNeuroVol(mask_array, space)
  
  created_files <- list()
  
  # Save each contrast separately with all its statistics
  for (i in seq_len(nrow(available_contrasts))) {
    contrast_name <- available_contrasts$name[i]
    contrast_data <- available_contrasts$data[[i]]
    
    # Map user statistic names to internal field names and collect available statistics
    available_stats <- c()
    mapped_stats <- c()
    for (stat in contrast_stats) {
      internal_stat_name <- .map_stat_name(stat)
      if (internal_stat_name %in% names(contrast_data)) {
        available_stats <- c(available_stats, stat)  # Keep user-friendly name for output
        mapped_stats <- c(mapped_stats, internal_stat_name)  # Store internal name for data access
      }
    }
    
    if (length(available_stats) > 0) {
      # Create 4D array for this contrast's statistics - optimized version
      stat_array <- array(0, dim = c(brain_dims, length(available_stats)))
      mask_indices <- which(mask_array, arr.ind = TRUE)
      n_mask_voxels <- sum(mask_array)
      
      # Vectorized assignment using mask indices
      for (j in seq_along(available_stats)) {
        internal_stat_name <- mapped_stats[j]
        stat_values <- contrast_data[[internal_stat_name]][[1]]
        if (length(stat_values) != n_mask_voxels) {
          stop(
            "Statistic '", available_stats[[j]], "' for contrast '", contrast_name,
            "' has ", length(stat_values),
            " values but mask contains ", n_mask_voxels, " in-mask voxels.",
            call. = FALSE
          )
        }
        stat_array[cbind(mask_indices, j)] <- stat_values
      }
      
      # Create NeuroVec using add_dim
      ospace_stat <- neuroim2::add_dim(space, length(available_stats))
      stat_vec <- neuroim2::NeuroVec(stat_array, ospace_stat)
      
      # Generate filename with desc-statmap suffix (BIDS compliant)
      stat_desc <- paste0(desc, "statmap")
      filename <- .generate_bids_filename(entities, desc = stat_desc, contrast = contrast_name, suffix = "bold")
      filepath <- file.path(path, filename)
      
      h5_handle <- fmristore::write_labeled_vec(stat_vec, mask_vol, available_stats, file = filepath)
      h5_handle$close_all()
      
      # Save JSON metadata
      json_filepath <- .save_single_contrast_json_metadata(fmrilm_obj, path, entities, desc,
                                                          contrast_name, available_stats, available_contrasts[i, ], filepath)
      
      created_files[[contrast_name]] <- list(h5 = filepath, json = json_filepath)
    }
  }
  
  return(created_files)
}


#' Save GDS representation of subject results
#' @keywords internal
#' @noRd
.save_gds_outputs <- function(fmrilm_obj, path, entities, desc) {
  # Collect beta-level statistics
  betas_tbl <- fmrilm_obj$result$betas$data[[1]]
  extract_component <- function(col) {
    if (is.null(col)) return(NULL)
    if (is.list(col)) {
      col <- col[[1]]
    }
    col
  }
  beta_mat <- extract_component(betas_tbl$estimate)
  se_mat   <- extract_component(betas_tbl$se)
  stat_mat <- extract_component(betas_tbl$stat)
  prob_mat <- extract_component(betas_tbl$prob)

  if (is.null(beta_mat) || is.null(se_mat)) {
    stop("Unable to extract beta estimates or standard errors for GDS export.", call. = FALSE)
  }

  # Determine dimensions and labels
  n_vox <- nrow(beta_mat)
  n_coef <- ncol(beta_mat)
  contrast_names <- colnames(beta_mat) %||% paste0("beta", seq_len(n_coef))
  subject_label <- entities$subject %||% "subject1"
  subject_dimnames <- subject_label

  sample_labels <- rownames(beta_mat)
  if (is.null(sample_labels)) {
    sample_labels <- sprintf("voxel_%0*d", ceiling(log10(n_vox + 1)), seq_len(n_vox))
  }

  # Assemble assays
  beta_array <- array(beta_mat, dim = c(n_vox, 1, n_coef),
                      dimnames = list(sample_labels, subject_dimnames, contrast_names))
  se_array   <- array(se_mat, dim = c(n_vox, 1, n_coef),
                      dimnames = list(sample_labels, subject_dimnames, contrast_names))
  var_array  <- array(se_mat^2, dim = c(n_vox, 1, n_coef),
                      dimnames = list(sample_labels, subject_dimnames, contrast_names))

  assays <- list(beta = beta_array, var = var_array, se = se_array)

  if (!is.null(stat_mat)) {
    t_array <- array(stat_mat, dim = c(n_vox, 1, n_coef),
                     dimnames = list(sample_labels, subject_dimnames, contrast_names))
    assays$t <- t_array
  }
  if (!is.null(prob_mat)) {
    p_array <- array(prob_mat, dim = c(n_vox, 1, n_coef),
                     dimnames = list(sample_labels, subject_dimnames, contrast_names))
    assays$p <- p_array
  }
  df_resid <- .extract_degrees_of_freedom(fmrilm_obj)

  # Build spatial description from mask
  spatial <- .fmri_dataset_mask_space(fmrilm_obj$dataset, "GDS export")
  mask_array <- spatial$mask_array
  neuro_space <- spatial$space
  affine <- neuroim2::trans(neuro_space)
  voxel_space <- fmrigds::space_voxel(dim = dim(mask_array), affine = affine, mask_bitmap = mask_array)

  # Assemble column metadata
  col_df <- data.frame(subject = subject_label, stringsAsFactors = FALSE)
  rownames(col_df) <- subject_label
  if (!is.null(entities$task)) col_df$task <- entities$task
  if (!is.null(entities$space)) col_df$space <- entities$space

  # Metadata about the model
  model_formula <- tryCatch({
    form <- get_formula(fmrilm_obj$model)
    if (is.null(form)) NULL else as.character(form)
  }, error = function(e) NULL)

  metadata <- list(
    schema_version = "0.1.0",
    model = list(
      formula = model_formula,
    regressors = contrast_names,
    df_residual = df_resid
  ),
    generated_by = list(
      package = "fmrireg",
      version = as.character(utils::packageVersion("fmrireg"))
    ),
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  )

  gds_obj <- fmrigds::as_gds(
    assays,
    space = voxel_space,
    subjects = subject_label,
    contrasts = contrast_names,
    col_data = col_df,
    metadata = metadata
  )

  gds_h5_filename <- .generate_bids_filename(entities, desc = desc, suffix = "gds", extension = "h5")
  gds_h5_path <- file.path(path, gds_h5_filename)

  write_plan <- fmrigds::write_out(gds_obj, path = gds_h5_path, format = "h5")
  fmrigds::compute(write_plan)

  plan_filename <- .generate_bids_filename(entities, desc = desc, suffix = "gds", extension = "rds")
  plan_path <- file.path(path, plan_filename)

  list(h5 = gds_h5_path, plan = plan_path)
}


#' Save JSON Metadata for Betas
#' @keywords internal
#' @noRd
.save_betas_json_metadata <- function(fmrilm_obj, path, entities, desc, regressor_names, h5_filepath) {
  
  # Generate JSON filename
  filename <- .generate_bids_filename(entities, desc = desc, suffix = "betas", extension = "json")
  filepath <- file.path(path, filename)
  
  # Extract source information
  sources <- .extract_source_files(fmrilm_obj)
  
  # Create metadata structure
  metadata <- list(
    Description = "Raw regressor beta coefficients from GLM analysis",
    Sources = sources,
    SoftwareVersions = list(
      R = paste(R.Version()$major, R.Version()$minor, sep = "."),
      fmrireg = as.character(utils::packageVersion("fmrireg")),
      fmristore = as.character(utils::packageVersion("fmristore")),
      neuroim2 = as.character(utils::packageVersion("neuroim2"))
    ),
    ModelInfo = list(
      Type = "General Linear Model",
      Formula = as.character(deparse(fmrilm_obj$model$event_model$model_spec$formula_or_list)),
      NumRegressors = as.integer(length(regressor_names))
    ),
    RegressorOrder = as.character(regressor_names),
    DataInfo = list(
      FileFormat = "HDF5",
      Storage = "LabeledVolumeSet",
      Units = "arbitrary (depends on input data scaling)"
    ),
    GeneratedBy = list(
      Name = "fmrireg::write_results",
      Version = as.character(utils::packageVersion("fmrireg")),
      CodeURL = "https://github.com/bbuchsbaum/fmrireg"
    ),
    CreationTime = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  )
  
  # Add spatial reference if available
  if (!is.null(entities$space)) {
    metadata$SpatialReference = entities$space
  }
  
  # Save JSON
  jsonlite::write_json(metadata, filepath, pretty = TRUE, auto_unbox = TRUE)
  
  return(filepath)
}

#' Save JSON Metadata for Contrasts (by stat)
#' @keywords internal
#' @noRd
.save_contrasts_json_metadata <- function(fmrilm_obj, path, entities, desc, stat, contrast_names, contrasts, h5_filepath) {
  
  # Generate JSON filename
  stat_desc <- paste0(desc, "statmap")
  filename <- .generate_bids_filename(entities, desc = stat_desc, stat = stat, suffix = "bold", extension = "json")
  filepath <- file.path(path, filename)
  
  # Extract contrast definitions
  contrast_definitions <- list()
  for (i in seq_len(nrow(contrasts))) {
    contrast_name <- contrasts$name[i]
    if (contrast_name %in% contrast_names) {
      contrast_type <- .extract_contrast_type(contrasts[i, ], fmrilm_obj)
      contrast_definitions[[contrast_name]] <- list(
        Name = contrast_name,
        Type = contrast_type
      )
    }
  }
  
  # Extract source information and degrees of freedom
  sources <- .extract_source_files(fmrilm_obj)
  df_resid <- .extract_degrees_of_freedom(fmrilm_obj)
  
  # Create metadata structure
  metadata <- list(
    Description = paste("Contrast", stat, "statistics from GLM analysis"),
    StatisticType = toupper(stat),
    DegreesOfFreedom = df_resid,
    Sources = sources,
    SoftwareVersions = list(
      R = paste(R.Version()$major, R.Version()$minor, sep = "."),
      fmrireg = as.character(utils::packageVersion("fmrireg")),
      fmristore = as.character(utils::packageVersion("fmristore")),
      neuroim2 = as.character(utils::packageVersion("neuroim2"))
    ),
    ModelInfo = list(
      Type = "General Linear Model", 
      Formula = as.character(deparse(get_formula(fmrilm_obj$model)))
    ),
    ContrastOrder = contrast_names,
    ContrastDefinitions = contrast_definitions,
    DataInfo = list(
      FileFormat = "HDF5",
      Storage = "LabeledVolumeSet",
      Units = .get_stat_units(stat)
    ),
    GeneratedBy = list(
      Name = "fmrireg::write_results",
      Version = as.character(utils::packageVersion("fmrireg")),
      CodeURL = "https://github.com/bbuchsbaum/fmrireg"
    ),
    CreationTime = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  )
  
  # Add spatial reference if available
  if (!is.null(entities$space)) {
    metadata$SpatialReference = entities$space
  }
  
  # Save JSON
  jsonlite::write_json(metadata, filepath, pretty = TRUE, auto_unbox = TRUE)
  
  return(filepath)
}

#' Save JSON Metadata for Single Contrast
#' @keywords internal
#' @noRd
.save_single_contrast_json_metadata <- function(fmrilm_obj, path, entities, desc, contrast_name, stat_names, contrast_row, h5_filepath) {
  
  # Generate JSON filename
  stat_desc <- paste0(desc, "statmap")
  filename <- .generate_bids_filename(entities, desc = stat_desc, contrast = contrast_name, suffix = "bold", extension = "json")
  filepath <- file.path(path, filename)
  
  # Extract source information and degrees of freedom
  sources <- .extract_source_files(fmrilm_obj)
  df_resid <- .extract_degrees_of_freedom(fmrilm_obj)
  
  # Create metadata structure
  metadata <- list(
    Description = paste("Statistical maps for contrast", contrast_name, "from GLM analysis"),
    ContrastName = contrast_name,
    DegreesOfFreedom = df_resid,
    Sources = sources,
    SoftwareVersions = list(
      R = paste(R.Version()$major, R.Version()$minor, sep = "."),
      fmrireg = as.character(utils::packageVersion("fmrireg")),
      fmristore = as.character(utils::packageVersion("fmristore")),
      neuroim2 = as.character(utils::packageVersion("neuroim2"))
    ),
    ModelInfo = list(
      Type = "General Linear Model",
      Formula = as.character(deparse(get_formula(fmrilm_obj$model)))
    ),
    ContrastInfo = list(
      Name = contrast_name,
      Type = .extract_contrast_type(contrast_row, fmrilm_obj)
    ),
    StatisticOrder = stat_names,
    DataInfo = list(
      FileFormat = "HDF5",
      Storage = "LabeledVolumeSet",
      Units = "various (see StatisticOrder for mapping)"
    ),
    GeneratedBy = list(
      Name = "fmrireg::write_results",
      Version = as.character(utils::packageVersion("fmrireg")),
      CodeURL = "https://github.com/bbuchsbaum/fmrireg"
    ),
    CreationTime = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  )
  
  # Add spatial reference if available
  if (!is.null(entities$space)) {
    metadata$SpatialReference = entities$space
  }
  
  # Save JSON
  jsonlite::write_json(metadata, filepath, pretty = TRUE, auto_unbox = TRUE)
  
  return(filepath)
}

#' Get Units for Statistic Type
#' @keywords internal
#' @noRd
.get_stat_units <- function(stat) {
  switch(stat,
    "beta" = "arbitrary (depends on input data scaling)",
    "tstat" = "t-statistic",
    "pval" = "probability",
    "se" = "standard error",
    "z" = "z-score",
    "f" = "F-statistic",
    "arbitrary"
  )
}

#' Extract Source Files from fmri_lm Object
#' @keywords internal
.extract_source_files <- function(fmrilm_obj) {
  # Try to extract source file information from dataset
  sources <- character(0)
  
  tryCatch({
    # Check if dataset has source information
    if (!is.null(fmrilm_obj$dataset) && !is.null(attr(fmrilm_obj$dataset, "source_files"))) {
      sources <- attr(fmrilm_obj$dataset, "source_files")
    } else if (!is.null(fmrilm_obj$dataset) && inherits(fmrilm_obj$dataset, "fmri_dataset")) {
      # Try to get file path if available
      dataset_path <- attr(fmrilm_obj$dataset, "file_path")
      if (!is.null(dataset_path)) {
        # Convert to relative BIDS URI if possible
        sources <- paste0("bids::", basename(dataset_path))
      }
    }
  }, error = function(e) {
    # If extraction fails, leave empty
    sources <<- character(0)
  })
  
  return(sources)
}

#' Extract Contrast Type Information
#' @keywords internal
#' @noRd
.extract_contrast_type <- function(contrast_row, fmrilm_obj) {
  # Try to determine contrast type from available information
  tryCatch({
    contrast_data <- contrast_row$data[[1]]
    
    # Check what statistics are available to infer type
    available_stats <- names(contrast_data)
    
    # If F-statistic is present, it's likely an F-contrast
    if ("fstat" %in% available_stats || "f" %in% available_stats) {
      return("F-contrast")
    }
    
    # If only t-statistic, beta, and p-values, likely a t-contrast
    if (any(c("tstat", "t") %in% available_stats)) {
      return("t-contrast")
    }
    
    # Check if contrast has matrix information (safer approach)
    # Use safe column access to avoid warnings
    has_contrast_matrix <- "contrast_matrix" %in% names(contrast_row)
    has_weights <- "weights" %in% names(contrast_row)
    
    if (has_contrast_matrix && !is.null(contrast_row[["contrast_matrix"]])) {
      # If we can access the contrast matrix, check dimensions
      contrast_matrix <- contrast_row[["contrast_matrix"]]
      if (is.matrix(contrast_matrix) && nrow(contrast_matrix) > 1) {
        return("F-contrast")
      } else {
        return("t-contrast")
      }
    } else if (has_weights && !is.null(contrast_row[["weights"]])) {
      # Check weights if available
      return("t-contrast")
    }
    
    # Default fallback based on naming patterns
    contrast_name <- contrast_row$name
    if (grepl("vs|VS|_vs_|-vs-", contrast_name)) {
      return("t-contrast")  # Pairwise comparisons are usually t-contrasts
    }
    
    # Default to generic contrast
    return("contrast")
    
  }, error = function(e) {
    # If extraction fails, return generic type
    return("contrast")
  })
}

#' Extract Degrees of Freedom from fmri_lm Object
#' @keywords internal
#' @noRd
.extract_degrees_of_freedom <- function(fmrilm_obj) {
  # Extract residual degrees of freedom
  df_resid <- NULL
  
  tryCatch({
    # Check if df.residual is directly in the result structure
    if (!is.null(fmrilm_obj$result$df.residual)) {
      df_resid <- as.integer(fmrilm_obj$result$df.residual)
    } else if (!is.null(fmrilm_obj$result$betas) && "df.residual" %in% names(fmrilm_obj$result$betas)) {
      # Check if it's in the betas data frame
      df_resid <- as.integer(fmrilm_obj$result$betas$df.residual[1])
    } else if (!is.null(fmrilm_obj$result$contrasts) && "df.residual" %in% names(fmrilm_obj$result$contrasts)) {
      df_resid <- as.integer(fmrilm_obj$result$contrasts$df.residual[1])
    }
  }, error = function(e) {
    # If extraction fails, leave NULL
    df_resid <<- NULL
  })
  
  return(df_resid)
}

#' Create Temporary Directory for Atomic Writes
#' @keywords internal
#' @noRd
.create_temp_write_dir <- function(base_path) {
  # Create unique temporary directory name
  temp_name <- paste0(".tmp_write_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", 
                      sample(10000:99999, 1))
  temp_dir <- file.path(base_path, temp_name)
  
  # Create temporary directory
  if (!dir.create(temp_dir, recursive = TRUE)) {
    stop("Failed to create temporary directory: ", temp_dir, call. = FALSE)
  }
  
  return(temp_dir)
}

#' Stage an Existing File for Rollback
#' @keywords internal
#' @noRd
.stash_existing_file <- function(file, backup_root, overwrite = FALSE, relative_name = basename(file)) {
  if (!file.exists(file)) {
    return(NULL)
  }
  if (!overwrite) {
    stop("Output file already exists: ", file, call. = FALSE)
  }
  
  backup_file <- file.path(backup_root, relative_name)
  if (!dir.exists(dirname(backup_file))) {
    dir.create(dirname(backup_file), recursive = TRUE)
  }
  
  if (!file.rename(file, backup_file)) {
    if (!file.copy(file, backup_file, overwrite = TRUE)) {
      stop("Failed to back up existing file: ", file, call. = FALSE)
    }
    if (!file.remove(file)) {
      stop("Failed to remove existing file after backup: ", file, call. = FALSE)
    }
  }
  
  backup_file
}

#' Finalize Atomic Write by Moving Files
#' @keywords internal
.finalize_atomic_write <- function(temp_dir, final_path, created_files, overwrite = FALSE) {
  # Move all files from temp_dir to final_path
  temp_files <- list.files(temp_dir, full.names = TRUE, recursive = TRUE)
  rel_paths <- sub(paste0("^", temp_dir, "/"), "", temp_files)
  final_files <- file.path(final_path, rel_paths)
  existing_targets <- final_files[file.exists(final_files)]
  
  if (!overwrite && length(existing_targets) > 0) {
    stop(
      "Output files already exist: ",
      paste(existing_targets, collapse = ", "),
      ". Set overwrite = TRUE to replace them.",
      call. = FALSE
    )
  }
  
  backup_root <- file.path(temp_dir, ".overwrite_backup")
  backup_map <- character(0)
  moved_files <- character(0)
  success <- FALSE
  
  on.exit({
    if (!success) {
      .rollback_finalized_write(list(
        moved_files = moved_files,
        backup_map = backup_map
      ))
    }
  }, add = TRUE)
  
  for (i in seq_along(temp_files)) {
    temp_file <- temp_files[[i]]
    # Get relative path within temp directory
    rel_path <- rel_paths[[i]]
    final_file <- final_files[[i]]
    
    # Ensure destination directory exists
    final_dir <- dirname(final_file)
    if (!dir.exists(final_dir)) {
      dir.create(final_dir, recursive = TRUE)
    }
    
    backup_file <- .stash_existing_file(
      final_file,
      backup_root = backup_root,
      overwrite = overwrite,
      relative_name = rel_path
    )
    if (!is.null(backup_file)) {
      backup_map[[final_file]] <- backup_file
    }
    
    # Move file
    if (!file.rename(temp_file, final_file)) {
      # If rename fails, try copy + delete
      if (!file.copy(temp_file, final_file, overwrite = FALSE)) {
        stop("Failed to move file from ", temp_file, " to ", final_file, call. = FALSE)
      }
      file.remove(temp_file)
    }
    moved_files <- c(moved_files, final_file)
  }
  
  # Update file paths in created_files to point to final locations
  updated_files <- .update_file_paths_in_results(created_files, temp_dir, final_path)
  success <- TRUE
  
  list(
    files = updated_files,
    transaction = list(
      moved_files = moved_files,
      backup_map = backup_map
    )
  )
}

#' Update File Paths in Results After Atomic Move
#' @keywords internal
#' @noRd
.update_file_paths_in_results <- function(created_files, temp_dir, final_path) {
  # Recursively update file paths in the results structure
  for (i in seq_along(created_files)) {
    if (is.list(created_files[[i]])) {
      created_files[[i]] <- .update_file_paths_in_results(created_files[[i]], temp_dir, final_path)
    } else if (is.character(created_files[[i]]) && length(created_files[[i]]) == 1) {
      # Replace temp_dir with final_path in file paths
      created_files[[i]] <- gsub(paste0("^", temp_dir), final_path, created_files[[i]])
    }
  }
  return(created_files)
}

#' Clean Up Temporary Write Directory
#' @keywords internal
#' @noRd
.cleanup_temp_write_dir <- function(temp_dir) {
  if (dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE)
  }
}

#' Roll Back Finalized Writes
#' @keywords internal
#' @noRd
.rollback_finalized_write <- function(transaction) {
  moved_files <- rev(transaction$moved_files %||% character(0))
  backup_map <- transaction$backup_map %||% character(0)
  
  for (file in moved_files) {
    if (file.exists(file)) {
      file.remove(file)
    }
  }
  
  if (length(backup_map) > 0) {
    for (final_file in names(backup_map)) {
      backup_file <- backup_map[[final_file]]
      if (file.exists(backup_file)) {
        if (!dir.exists(dirname(final_file))) {
          dir.create(dirname(final_file), recursive = TRUE)
        }
        if (!file.rename(backup_file, final_file)) {
          file.copy(backup_file, final_file, overwrite = TRUE)
          file.remove(backup_file)
        }
      }
    }
  }
}

#' Write GDS Plan File
#' @keywords internal
#' @noRd
.write_gds_plan_file <- function(gds_h5_final, gds_plan_final, overwrite = FALSE, backup_root = NULL) {
  backup_file <- NULL
  if (!is.null(backup_root)) {
    backup_file <- .stash_existing_file(
      gds_plan_final,
      backup_root = backup_root,
      overwrite = overwrite,
      relative_name = basename(gds_plan_final)
    )
  }
  
  success <- FALSE
  on.exit({
    if (!success && !is.null(backup_file) && file.exists(backup_file)) {
      if (file.exists(gds_plan_final)) {
        file.remove(gds_plan_final)
      }
      if (!file.rename(backup_file, gds_plan_final)) {
        file.copy(backup_file, gds_plan_final, overwrite = TRUE)
        file.remove(backup_file)
      }
    }
  }, add = TRUE)
  
  disk_plan <- fmrigds::gds(source = gds_h5_final, format = "h5")
  saveRDS(disk_plan, gds_plan_final)
  success <- TRUE
  
  invisible(gds_plan_final)
}

#' Compute Beta Volumes from fmri_lm Object
#' @keywords internal
#' @noRd
.compute_beta_volumes <- function(fmrilm_obj) {
  # Input validation
  if (is.null(fmrilm_obj$result$betas$data) || 
      length(fmrilm_obj$result$betas$data) == 0 || 
      is.null(fmrilm_obj$result$betas$data[[1]]$estimate) ||
      length(fmrilm_obj$result$betas$data[[1]]$estimate) == 0) {
    stop("Beta data structure is invalid or empty", call. = FALSE)
  }
  
  # Extract beta data - use coef() to ensure we get all betas including baseline
  beta_data <- coef(fmrilm_obj, type = "betas", include_baseline = TRUE)
  
  # Additional validation for beta data structure
  if (!is.matrix(beta_data) || nrow(beta_data) == 0 || ncol(beta_data) == 0) {
    stop("Beta data structure is invalid or empty", call. = FALSE)
  }
  
  # Get spatial information from dataset
  spatial <- .fmri_dataset_mask_space(fmrilm_obj$dataset, "beta volume reconstruction")
  mask_array <- spatial$mask_array
  space <- spatial$space
  
  # Reshape beta data to 4D array [X, Y, Z, #regressors]
  brain_dims <- dim(mask_array)
  n_regressors <- ncol(beta_data)
  n_mask_voxels <- sum(mask_array)
  
  # Validate mask dimensionality
  if (length(mask_array) != prod(brain_dims)) {
    stop("Mask length (", length(mask_array), 
         ") does not match product of space dims (", prod(brain_dims), ").", 
         call. = FALSE)
  }
  
  # Validate beta data dimensions match mask
  if (nrow(beta_data) != n_mask_voxels) {
    stop(
      "Beta data has ", nrow(beta_data),
      " rows but mask contains ", n_mask_voxels, " in-mask voxels.",
      call. = FALSE
    )
  }
  
  # Create 4D array - optimized vectorized version
  beta_array <- array(0, dim = c(brain_dims, n_regressors))
  
  # Vectorized assignment: create mask indices once, then assign all regressors
  mask_indices <- which(mask_array, arr.ind = TRUE)
  for (i in seq_len(n_regressors)) {
    beta_array[cbind(mask_indices, i)] <- beta_data[, i]
  }
  
  # Create NeuroVec for the betas
  # Use add_dim to properly extend 3D space to 4D
  ospace_beta <- neuroim2::add_dim(space, n_regressors)
  beta_vec <- neuroim2::NeuroVec(beta_array, ospace_beta)
  
  return(beta_vec)
}

#' Map User Statistic Names to Internal Field Names
#' @keywords internal
#' @noRd
.map_stat_name <- function(stat) {
  # Map user-friendly names to internal contrast data field names
  stat_mapping <- c(
    "beta" = "estimate",
    "tstat" = "stat", 
    "se" = "se",
    "pval" = "prob",
    "sigma" = "sigma"
  )
  
  if (stat %in% names(stat_mapping)) {
    return(stat_mapping[[stat]])
  } else {
    return(stat)  # Return as-is if no mapping found
  }
}

#' Compute Statistical Map Volumes from Contrast Data
#' @keywords internal
#' @noRd
.compute_statistical_volumes <- function(contrast_data, stat, brain_dims, mask, space) {
  stat_matrices <- list()
  contrast_names <- c()
  
  # Map user statistic name to internal field name
  internal_stat_name <- .map_stat_name(stat)
  
  for (i in seq_len(nrow(contrast_data))) {
    contrast_entry <- contrast_data$data[[i]]
    if (internal_stat_name %in% names(contrast_entry)) {
      # Handle both indexed and named access patterns
      stat_value <- if (is.list(contrast_entry[[internal_stat_name]]) && length(contrast_entry[[internal_stat_name]]) > 0) {
        contrast_entry[[internal_stat_name]][[1]]
      } else {
        contrast_entry[[internal_stat_name]]
      }
      stat_matrices[[i]] <- stat_value
      contrast_names[i] <- contrast_data$name[i]
    }
  }
  
  # Remove NULL entries
  valid_indices <- !sapply(stat_matrices, is.null)
  stat_matrices <- stat_matrices[valid_indices]
  contrast_names <- contrast_names[valid_indices]
  
  if (length(stat_matrices) == 0) {
    return(NULL)
  }
  
  # Create 4D array [X, Y, Z, #contrasts] - optimized version
  n_contrasts <- length(stat_matrices)
  stat_array <- array(0, dim = c(brain_dims, n_contrasts))
  mask_array <- as.logical(mask)
  mask_indices <- which(mask_array, arr.ind = TRUE)
  n_mask_voxels <- sum(mask_array)
  
  # Vectorized assignment using mask indices
  for (i in seq_len(n_contrasts)) {
    if (length(stat_matrices[[i]]) != n_mask_voxels) {
      stop(
        "Statistic '", stat, "' for contrast '", contrast_names[[i]],
        "' has ", length(stat_matrices[[i]]),
        " values but mask contains ", n_mask_voxels, " in-mask voxels.",
        call. = FALSE
      )
    }
    stat_array[cbind(mask_indices, i)] <- stat_matrices[[i]]
  }
  
  # Create NeuroVec using add_dim to properly extend 3D space to 4D
  ospace_contrast <- neuroim2::add_dim(space, n_contrasts)
  neurovec <- neuroim2::NeuroVec(stat_array, ospace_contrast)
  
  return(list(
    neurovec = neurovec,
    contrast_names = contrast_names
  ))
}
