# S3 Methods for fmri_meta Objects
# Methods for extracting and manipulating meta-analysis results

#' Extract Coefficients from Meta-Analysis
#'
#' @param object An fmri_meta object
#' @param ... Additional arguments
#' @return Matrix of coefficients
#' @export
coef.fmri_meta <- function(object, ...) {
  object$coefficients
}

#' Extract Standard Errors from Meta-Analysis
#'
#' @param object An fmri_meta object
#' @return Matrix of standard errors
#' @export
se <- function(object) {
  UseMethod("se")
}

#' @export
se.fmri_meta <- function(object) {
  object$se
}

#' Compute Z-scores from Meta-Analysis
#'
#' @param object An fmri_meta object
#' @return Matrix of z-scores
#' @export
zscores <- function(object) {
  UseMethod("zscores")
}

#' @export
zscores.fmri_meta <- function(object) {
  object$coefficients / object$se
}

#' Compute P-values from Meta-Analysis
#'
#' @param object An fmri_meta object
#' @param two_tailed Logical. Use two-tailed test (default: TRUE)
#' @return Matrix of p-values
#' @export
pvalues <- function(object, two_tailed = TRUE) {
  UseMethod("pvalues")
}

#' @export
pvalues.fmri_meta <- function(object, two_tailed = TRUE) {
  z <- zscores(object)
  if (two_tailed) {
    2 * pnorm(-abs(z))
  } else {
    pnorm(-z)
  }
}

#' Apply Contrast to Meta-Analysis Results
#'
#' @param x An fmri_meta object
#' @param contrast Contrast specification. Can be:
#'   \itemize{
#'     \item A numeric vector of contrast weights
#'     \item A formula (e.g., ~ groupold - groupyoung)
#'     \item A named vector (e.g., c("groupold" = 1, "groupyoung" = -1))
#'   }
#' @param ... Additional arguments
#' @return An fmri_meta_contrast object with contrast results
#' @export
contrast.fmri_meta <- function(x, contrast, ...) {
  # Parse contrast specification
  if (is.numeric(contrast) && is.null(names(contrast))) {
    # Numeric vector - must match number of coefficients
    if (length(contrast) != ncol(x$coefficients)) {
      stop("Contrast vector length must match number of coefficients", call. = FALSE)
    }
    contrast_weights <- contrast
    contrast_name <- "custom"
  } else if (inherits(contrast, "formula")) {
    # Formula specification
    contrast_weights <- parse_contrast_formula(contrast, x)
    contrast_name <- deparse(contrast)
  } else if (is.numeric(contrast) && !is.null(names(contrast))) {
    # Named numeric vector
    contrast_weights <- build_contrast_from_names(contrast, x)
    contrast_name <- paste(names(contrast)[contrast != 0], collapse = "_vs_")
  } else {
    stop("Invalid contrast specification", call. = FALSE)
  }
  
  # Apply contrast
  contrast_est <- x$coefficients %*% contrast_weights
  
  # Compute standard error of contrast
  # SE = sqrt(c' * V * c) where V is diagonal matrix of variances
  contrast_var <- rowSums((x$se^2) * matrix(contrast_weights^2, 
                                                  nrow = nrow(x$se), 
                                                  ncol = length(contrast_weights), 
                                                  byrow = TRUE))
  contrast_se <- sqrt(contrast_var)
  
  # Compute z-scores and p-values
  contrast_z <- contrast_est / contrast_se
  contrast_p <- 2 * pnorm(-abs(contrast_z))
  
  # Create contrast object
  structure(
    list(
      estimate = as.vector(contrast_est),
      se = contrast_se,
      z = as.vector(contrast_z),
      p = contrast_p,
      weights = contrast_weights,
      name = contrast_name,
      parent = x
    ),
    class = "fmri_meta_contrast"
  )
}

#' Parse Contrast Formula
#'
#' @param formula Contrast formula
#' @param object fmri_meta object
#' @return Numeric vector of contrast weights
#' @keywords internal
parse_contrast_formula <- function(formula, object) {
  # Get coefficient names
  coef_names <- colnames(object$coefficients)
  
  # Simple implementation - would expand for complex formulas
  formula_str <- deparse(formula)
  
  # Initialize weights
  weights <- rep(0, length(coef_names))
  names(weights) <- coef_names
  
  # Parse simple contrasts like ~ groupold - groupyoung
  # This is a placeholder - full implementation would use proper formula parsing
  if (grepl("-", formula_str)) {
    parts <- strsplit(formula_str, "-")[[1]]
    parts <- trimws(gsub("~", "", parts))
    
    if (length(parts) == 2) {
      if (parts[1] %in% coef_names) weights[parts[1]] <- 1
      if (parts[2] %in% coef_names) weights[parts[2]] <- -1
    }
  }
  
  return(weights)
}

#' Build Contrast from Named Vector
#'
#' @param contrast Named numeric vector
#' @param object fmri_meta object
#' @return Numeric vector of contrast weights
#' @keywords internal
build_contrast_from_names <- function(contrast, object) {
  coef_names <- colnames(object$coefficients)
  weights <- rep(0, length(coef_names))
  names(weights) <- coef_names
  
  for (name in names(contrast)) {
    if (name %in% coef_names) {
      weights[name] <- contrast[name]
    } else {
      warning("Coefficient '", name, "' not found in model", call. = FALSE)
    }
  }
  
  return(weights)
}

#' Print Meta-Analysis Results
#'
#' @param x An fmri_meta object
#' @param ... Additional print arguments
#' @return Invisibly returns the input object x
#' @export
print.fmri_meta <- function(x, ...) {
  cat("fMRI Meta-Analysis Results\n")
  cat("==========================\n\n")
  cat("Method:", x$method, "\n")
  cat("Robust:", x$robust, "\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Subjects:", x$n_subjects, "\n")
  
  if (!is.null(x$n_voxels)) {
    cat("Voxels analyzed:", x$n_voxels, "\n")
  } else if (!is.null(x$n_rois)) {
    cat("ROIs analyzed:", x$n_rois, "\n")
  }
  
  # Show mean heterogeneity
  if (!is.null(x$tau2)) {
    cat("\nHeterogeneity:\n")
    cat("  Mean tau²:", mean(x$tau2, na.rm = TRUE), "\n")
    cat("  Mean I²:", mean(x$I2, na.rm = TRUE), "%\n")
  }
  
  invisible(x)
}

#' Summary of Meta-Analysis Results
#'
#' @param object An fmri_meta object
#' @param threshold P-value threshold for significance (default: 0.05)
#' @param ... Additional summary arguments
#' @return A list containing summary statistics invisibly
#' @export
summary.fmri_meta <- function(object, threshold = 0.05, ...) {
  cat("fMRI Meta-Analysis Summary\n")
  cat("==========================\n\n")
  
  # Basic info
  print(object)
  
  # Coefficient summary
  cat("\nCoefficients:\n")
  coef_names <- colnames(object$coefficients)
  
  for (i in seq_along(coef_names)) {
    coef_vals <- object$coefficients[, i]
    se_vals <- object$se[, i]
    z_vals <- coef_vals / se_vals
    p_vals <- 2 * pnorm(-abs(z_vals))
    
    # Count significant voxels/ROIs
    n_sig <- sum(p_vals < threshold, na.rm = TRUE)
    n_total <- sum(!is.na(p_vals))
    
    cat("  ", coef_names[i], ":\n", sep = "")
    cat("    Mean effect:", mean(coef_vals, na.rm = TRUE), "\n")
    cat("    Mean SE:", mean(se_vals, na.rm = TRUE), "\n")
    cat("    Significant:", n_sig, "/", n_total, 
        " (", round(100 * n_sig / n_total, 1), "%)\n", sep = "")
  }
  
  invisible(object)
}

#' Extract Image/Volume for Coefficient
#'
#' @param object An fmri_meta object
#' @param coef Coefficient name or index
#' @param statistic Type of statistic to extract ("estimate", "se", "z", "p")
#' @return NeuroVol object or matrix
#' @export
coef_image <- function(object, coef = 1, statistic = c("estimate", "se", "z", "p")) {
  UseMethod("coef_image")
}

#' @export
coef_image.fmri_meta <- function(object, coef = 1, statistic = c("estimate", "se", "z", "p")) {
  statistic <- match.arg(statistic)
  
  # Get coefficient index
  if (is.character(coef)) {
    coef_idx <- which(colnames(object$coefficients) == coef)
    if (length(coef_idx) == 0) {
      stop("Coefficient '", coef, "' not found", call. = FALSE)
    }
  } else {
    coef_idx <- coef
  }
  
  # Extract values
  values <- switch(statistic,
    estimate = object$coefficients[, coef_idx],
    se = object$se[, coef_idx],
    z = object$coefficients[, coef_idx] / object$se[, coef_idx],
    p = 2 * pnorm(-abs(object$coefficients[, coef_idx] / object$se[, coef_idx]))
  )
  
  # Convert to image if spatial data
  if (inherits(object$data, "group_data_h5") || inherits(object$data, "group_data_nifti")) {
    values <- reconstruct_image(values, object)
  }
  
  return(values)
}

#' Reconstruct 3D Image from Vector
#'
#' @param values Numeric vector of values
#' @param object fmri_meta object with spatial information
#' @return NeuroVol object
#' @keywords internal
reconstruct_image <- function(values, object) {
  if (inherits(object$data, "group_data_h5")) {
    # Use space information from group_data_h5
    if (!is.null(object$data$space)) {
      space <- object$data$space
    } else {
      # Create space from dimensions
      space <- neuroim2::NeuroSpace(object$data$dim[1:3])
    }
    
    # Create full volume
    vol_data <- array(NA_real_, dim = object$data$dim[1:3])
    
    # Fill in values
    if (!is.null(object$data$mask)) {
      mask_indices <- which(object$data$mask > 0)
      vol_data[mask_indices] <- values
    } else {
      vol_data[] <- values
    }
    
    # Create NeuroVol
    return(neuroim2::NeuroVol(vol_data, space))
    
  } else if (inherits(object$data, "group_data_nifti")) {
    # Similar for NIfTI
    vol_data <- array(NA_real_, dim = object$data$dim[1:3])
    
    if (!is.null(object$voxel_indices)) {
      vol_data[object$voxel_indices] <- values
    } else {
      vol_data[] <- values
    }
    
    space <- neuroim2::NeuroSpace(object$data$dim[1:3], 
                                  spacing = object$data$voxel_size)
    
    return(neuroim2::NeuroVol(vol_data, space))
  }
  
  # For non-spatial data, just return values
  return(values)
}

#' Write Meta-Analysis Results
#'
#' @param x An fmri_meta object
#' @param path Output directory
#' @param prefix File name prefix
#' @param format Output format ("nifti" or "h5")
#' @param ... Additional arguments passed to write_results.fmri_lm
#' @return Invisible list of created files
#' @export
write_results.fmri_meta <- function(x, path = ".", prefix = "meta", 
                                    format = c("nifti", "h5"), ...) {
  format <- match.arg(format)
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  created_files <- list()
  
  # Write each coefficient
  coef_names <- colnames(x$coefficients)
  
  for (i in seq_along(coef_names)) {
    coef_name <- make.names(coef_names[i])
    
    # Extract images
    est_img <- coef_image(x, i, "estimate")
    se_img <- coef_image(x, i, "se")
    z_img <- coef_image(x, i, "z")
    p_img <- coef_image(x, i, "p")
    
    if (format == "nifti") {
      # Write NIfTI files
      est_file <- file.path(path, paste0(prefix, "_", coef_name, "_beta.nii.gz"))
      se_file <- file.path(path, paste0(prefix, "_", coef_name, "_se.nii.gz"))
      z_file <- file.path(path, paste0(prefix, "_", coef_name, "_z.nii.gz"))
      p_file <- file.path(path, paste0(prefix, "_", coef_name, "_p.nii.gz"))
      
      neuroim2::write_vol(est_img, est_file)
      neuroim2::write_vol(se_img, se_file)
      neuroim2::write_vol(z_img, z_file)
      neuroim2::write_vol(p_img, p_file)
      
      created_files[[coef_name]] <- list(
        beta = est_file,
        se = se_file,
        z = z_file,
        p = p_file
      )
    } else {
      # Write HDF5 using fmristore
      # This would use the LabeledVolumeSet format
      # Implementation would be similar to write_results.fmri_lm
    }
  }
  
  # Write heterogeneity maps
  if (!is.null(x$tau2)) {
    tau2_img <- reconstruct_image(x$tau2, x)
    I2_img <- reconstruct_image(x$I2, x)
    
    if (format == "nifti") {
      tau2_file <- file.path(path, paste0(prefix, "_tau2.nii.gz"))
      I2_file <- file.path(path, paste0(prefix, "_I2.nii.gz"))
      
      neuroim2::write_vol(tau2_img, tau2_file)
      neuroim2::write_vol(I2_img, I2_file)
      
      created_files$heterogeneity <- list(
        tau2 = tau2_file,
        I2 = I2_file
      )
    }
  }
  
  invisible(created_files)
}

#' Internal handler for spatial FDR on fmri_meta objects
#'
#' @param object An fmri_meta object
#' @param coef Character or integer; coefficient name or index to test
#' @param group Grouping variable for spatial FDR. Can be:
#'   \itemize{
#'     \item NULL: Auto-create blocks (voxelwise) or use ROIs (ROI-wise)
#'     \item Integer/factor vector: Custom grouping
#'     \item "blocks": Create 3D blocks (voxelwise only)
#'     \item "parcels": Use existing parcellation
#'   }
#' @param ... Additional arguments passed to \code{spatial_fdr}
#' @keywords internal
spatial_fdr_fmri_meta <- function(object, coef = 1, group = NULL,
                                 alpha = 0.05, tau = 0.5, lambda = 1.0,
                                 neighbors = NULL, min_pi0 = 0.05,
                                 empirical_null = TRUE, verbose = FALSE) {
  # Get coefficient index
  if (is.character(coef)) {
    coef_idx <- which(colnames(object$coefficients) == coef)
    if (length(coef_idx) == 0) {
      stop("Coefficient '", coef, "' not found", call. = FALSE)
    }
  } else {
    coef_idx <- coef
  }
  
  # Extract z-scores for this coefficient
  z <- object$coefficients[, coef_idx] / object$se[, coef_idx]
  
  # Determine grouping
  if (is.null(group)) {
    if (inherits(object, "fmri_meta_roi")) {
      # ROI analysis - each ROI is its own group
      group <- seq_len(nrow(object$coefficients))
    } else {
      # Voxelwise - create 3D blocks
      if (!is.null(object$data$mask)) {
        blocks <- create_3d_blocks(object$data$mask)
        group <- blocks$group_id
        # Store for potential neighbor use
        attr(group, "neighbors") <- blocks$neighbors
      } else {
        # No mask - treat each voxel as its own group
        group <- seq_len(length(z))
      }
    }
  } else if (is.character(group) && group == "blocks") {
    # Explicitly requested blocks
    if (!is.null(object$data$mask)) {
      blocks <- create_3d_blocks(object$data$mask)
      group <- blocks$group_id
      attr(group, "neighbors") <- blocks$neighbors
    } else {
      stop("Cannot create blocks without mask", call. = FALSE)
    }
  }
  
  # Get neighbors if attached to group
  if (is.null(neighbors)) {
    neighbors <- attr(group, "neighbors", exact = TRUE)
  }
  
  # Apply spatial FDR
  result <- spatial_fdr(
    z = z,
    group = group,
    alpha = alpha,
    tau = tau,
    lambda = lambda,
    neighbors = neighbors,
    min_pi0 = min_pi0,
    empirical_null = empirical_null,
    verbose = verbose
  )
  
  # Add coefficient info
  result$coef_name <- colnames(object$coefficients)[coef_idx]
  result$coef_idx <- coef_idx
  
  return(result)
}

#' Tidy Meta-Analysis Results
#'
#' @param x An fmri_meta object
#' @param conf.int Logical. Include confidence intervals (default: FALSE)
#' @param conf.level Confidence level (default: 0.95)
#' @param ... Additional arguments
#' @return A tibble with tidy results
#' @export
tidy.fmri_meta <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required for tidy() method", call. = FALSE)
  }
  
  # For ROI analysis, return tidy data frame
  if (inherits(x, "fmri_meta_roi")) {
    result <- tibble::tibble(
      roi = rep(x$roi_names, ncol(x$coefficients)),
      term = rep(colnames(x$coefficients), each = length(x$roi_names)),
      estimate = as.vector(x$coefficients),
      std.error = as.vector(x$se),
      statistic = as.vector(x$coefficients / x$se),
      p.value = as.vector(2 * pnorm(-abs(x$coefficients / x$se))),
      tau2 = rep(x$tau2, ncol(x$coefficients)),
      I2 = rep(x$I2, ncol(x$coefficients))
    )
    
    if (conf.int) {
      z_crit <- qnorm((1 + conf.level) / 2)
      result$conf.low <- result$estimate - z_crit * result$std.error
      result$conf.high <- result$estimate + z_crit * result$std.error
    }
    
    return(result)
  }
  
  # For voxelwise analysis, return summary statistics
  coef_names <- colnames(x$coefficients)
  
  result <- tibble::tibble(
    term = coef_names,
    mean_estimate = colMeans(x$coefficients, na.rm = TRUE),
    mean_std.error = colMeans(x$se, na.rm = TRUE),
    n_significant = sapply(1:ncol(x$coefficients), function(i) {
      p <- 2 * pnorm(-abs(x$coefficients[, i] / x$se[, i]))
      sum(p < 0.05, na.rm = TRUE)
    }),
    prop_significant = n_significant / x$n_voxels,
    mean_tau2 = mean(x$tau2, na.rm = TRUE),
    mean_I2 = mean(x$I2, na.rm = TRUE)
  )
  
  return(result)
}
