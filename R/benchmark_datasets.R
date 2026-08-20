#' Load fMRI Benchmark Datasets
#'
#' This function provides easy access to the benchmark datasets included with the fmrireg package.
#' These datasets are designed for testing HRF fitting, beta estimation, and other fMRI analysis methods.
#'
#' @param dataset_name Character string specifying which dataset to load. Options include:
#'   \itemize{
#'     \item \code{"BM_Canonical_HighSNR"}: Canonical HRF with high SNR (3 conditions)
#'     \item \code{"BM_Canonical_LowSNR"}: Canonical HRF with low SNR (3 conditions)
#'     \item \code{"BM_HRF_Variability_AcrossVoxels"}: HRF varies across voxel groups (2 conditions)
#'     \item \code{"BM_Trial_Amplitude_Variability"}: Trial-to-trial amplitude variability (1 condition)
#'     \item \code{"BM_Complex_Realistic"}: Complex scenario with multiple factors (3 conditions)
#'     \item \code{"all"}: Returns all datasets as a list
#'     \item \code{"metadata"}: Returns metadata about the datasets
#'   }
#'
#' @return A list containing the specified benchmark dataset(s) with the following components:
#'   \itemize{
#'     \item \code{description}: Text description of the dataset
#'     \item \code{Y_noisy}: Matrix of noisy BOLD time series (time x voxels)
#'     \item \code{Y_clean}: Matrix of clean BOLD time series (when available)
#'     \item \code{X_list_true_hrf}: List of design matrices convolved with true HRF
#'     \item \code{true_hrf_parameters}: Reconstructed true HRF object(s) plus
#'       complete, versioned generation recipes
#'     \item \code{oracle_contract}: Scope and tolerance of the valid stored
#'       numerical oracle, including limitations for partial-truth datasets
#'     \item \code{event_onsets}: Vector of event onset times
#'     \item \code{condition_labels}: Vector of condition labels for each event
#'     \item \code{true_betas_condition}: Matrix of true condition-level beta values
#'     \item \code{true_amplitudes_trial}: Matrix of true trial-level amplitudes
#'     \item \code{TR}: Repetition time
#'     \item \code{total_time}: Total scan duration
#'     \item \code{run_length}: Number of acquired time points
#'     \item \code{noise_parameters}: Information about noise generation
#'     \item \code{simulation_seed}: Random seed used for generation
#'     \item \code{target_snr}: Target signal-to-noise ratio
#'   }
#'
#' @examples
#' # Load a specific dataset
#' high_snr_data <- load_benchmark_dataset("BM_Canonical_HighSNR")
#' 
#' # Get information about all available datasets
#' metadata <- load_benchmark_dataset("metadata")
#' 
#' # Load all datasets
#' all_data <- load_benchmark_dataset("all")
#' 
#' # Access the BOLD data
#' Y <- high_snr_data$Y_noisy
#' 
#' # Get event information
#' onsets <- high_snr_data$event_onsets
#' conditions <- high_snr_data$condition_labels
#'
#' @export
load_benchmark_dataset <- function(dataset_name = "BM_Canonical_HighSNR") {
  all_processed_data <- .ensure_benchmark_data_loaded()
  
  if (dataset_name == "all") {
    # Return all datasets except metadata
    all_names <- names(all_processed_data)
    all_names <- all_names[all_names != "metadata"]
    return(all_processed_data[all_names])
  } else if (dataset_name == "metadata") {
    return(all_processed_data$metadata)
  } else if (dataset_name %in% names(all_processed_data)) {
    return(all_processed_data[[dataset_name]])
  } else {
    available_names <- names(all_processed_data)
    available_names <- available_names[available_names != "metadata"]
    stop("Dataset '", dataset_name, "' not found. Available datasets: ", 
         paste(available_names, collapse = ", "))
  }
}

#' List Available Benchmark Datasets
#'
#' Returns a summary of all available benchmark datasets with their descriptions.
#'
#' @return A data.frame with dataset names and descriptions
#'
#' @examples
#' # See what benchmark datasets are available
#' list_benchmark_datasets()
#'
#' @export
list_benchmark_datasets <- function() {
  all_processed_data <- .ensure_benchmark_data_loaded()
  
  # Get all dataset names except metadata
  all_names <- names(all_processed_data)
  dataset_names <- all_names[all_names != "metadata"]
  
  # Extract descriptions
  descriptions <- sapply(dataset_names, function(name) {
    all_processed_data[[name]]$description %||% "No description available"
  })
  
  data.frame(
    Dataset = dataset_names,
    Description = descriptions,
    stringsAsFactors = FALSE
  )
}

#' Get Benchmark Dataset Summary
#'
#' Provides a detailed summary of a specific benchmark dataset including
#' dimensions, experimental design, and ground truth information.
#'
#' @param dataset_name Character string specifying which dataset to summarize
#'
#' @return A list with summary information about the dataset
#'
#' @examples
#' # Get summary of a specific dataset
#' summary_info <- get_benchmark_summary("BM_Canonical_HighSNR")
#' print(summary_info)
#'
#' @export
get_benchmark_summary <- function(dataset_name) {
  all_processed_data <- .ensure_benchmark_data_loaded()
  
  if (!dataset_name %in% names(all_processed_data)) {
    available_names <- names(all_processed_data)
    available_names <- available_names[available_names != "metadata"]
    stop("Dataset '", dataset_name, "' not found. Available datasets: ", 
         paste(available_names, collapse = ", "))
  }
  
  dataset <- all_processed_data[[dataset_name]]
  
  # Basic dimensions
  n_timepoints <- nrow(dataset$Y_noisy)
  n_voxels <- ncol(dataset$Y_noisy)
  n_events <- length(dataset$event_onsets)
  unique_conditions <- unique(dataset$condition_labels)
  n_conditions <- length(unique_conditions)
  
  # Experimental design
  condition_counts <- table(dataset$condition_labels)
  
  # HRF information
  hrf_info <- dataset$true_hrf_parameters
  
  list(
    description = dataset$description,
    dimensions = list(
      n_timepoints = n_timepoints,
      n_voxels = n_voxels,
      n_events = n_events,
      n_conditions = n_conditions
    ),
    experimental_design = list(
      conditions = unique_conditions,
      events_per_condition = as.list(condition_counts),
      TR = dataset$TR,
      total_time = dataset$total_time,
      run_length = dataset$run_length,
      target_snr = dataset$target_snr %||% "Not specified"
    ),
    hrf_information = hrf_info,
    noise_information = dataset$noise_parameters %||% "Not available"
  )
}

#' Create Design Matrix from Benchmark Dataset
#'
#' Helper function to create a design matrix from a benchmark dataset using
#' a specified HRF. This is useful for testing different HRF assumptions
#' against the ground truth.
#'
#' @param dataset_name Character string specifying which dataset to use
#' @param hrf HRF object to use for convolution (e.g., fmrihrf::HRF_SPMG1)
#' @param include_intercept Logical, whether to include an intercept column
#'
#' @return A matrix with the design matrix (time x conditions)
#'
#' @examples
#' # Create design matrix using canonical HRF
#' X <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", fmrihrf::HRF_SPMG1)
#' 
#' # Test with a different HRF
#' X_wrong <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", fmrihrf::HRF_SPMG2)
#'
#' @export
create_design_matrix_from_benchmark <- function(dataset_name, hrf, include_intercept = TRUE) {
  dataset <- load_benchmark_dataset(dataset_name)
  
  # Create time grid matching the actual data dimensions
  n_timepoints <- nrow(dataset$Y_noisy)
  time_grid <- seq(0, by = dataset$TR, length.out = n_timepoints)
  
  # Get unique conditions
  unique_conditions <- unique(dataset$condition_labels)
  n_conditions <- length(unique_conditions)
  
  # Create design matrix
  X_list <- list()
  
  for (i in seq_along(unique_conditions)) {
    cond_name <- unique_conditions[i]
    cond_onsets <- dataset$event_onsets[dataset$condition_labels == cond_name]
    
    if (length(cond_onsets) > 0) {
      # Use first onset to get duration if available from core_data_args
      duration <- 0.5  # Default duration
      if (!is.null(dataset$core_data_args) && !is.null(dataset$core_data_args$event_table)) {
        first_event_idx <- which(dataset$condition_labels == cond_name)[1]
        if (first_event_idx <= nrow(dataset$core_data_args$event_table)) {
          duration <- dataset$core_data_args$event_table$duration[first_event_idx]
        }
      }
      
      reg <- fmrihrf::regressor(cond_onsets, hrf, duration = duration, amplitude = 1)
      X_list[[cond_name]] <- fmrihrf::evaluate(reg, time_grid)
    }
  }
  
  # Combine into matrix - handle both vectors and matrices
  X <- do.call(cbind, X_list)
  
  # Create appropriate column names
  col_names <- c()
  for (name in names(X_list)) {
    element <- X_list[[name]]
    if (is.matrix(element) && ncol(element) > 1) {
      # If HRF has multiple basis functions (e.g., with derivatives)
      col_names <- c(col_names, paste0(name, "_", 1:ncol(element)))
    } else {
      col_names <- c(col_names, name)
    }
  }
  colnames(X) <- col_names
  
  # Add intercept if requested
  if (include_intercept) {
    intercept <- rep(1, nrow(X))
    X <- cbind(Intercept = intercept, X)
  }
  
  return(X)
}

# Helper function for safe correlation calculation
.safe_cor <- function(x, y) {
  if (sd(x) == 0 || sd(y) == 0) {
    # Return 1 if both are constant and equal, 0 if different constants, NA otherwise
    if (sd(x) == 0 && sd(y) == 0) {
      return(ifelse(all(x[1] == y), 1, 0))
    } else {
      return(NA)
    }
  } else {
    return(cor(x, y))
  }
}

#' Evaluate Method Performance on Benchmark Dataset
#'
#' Helper function to evaluate the performance of beta estimation methods
#' on benchmark datasets by comparing estimated betas to ground truth.
#'
#' @param dataset_name Character string specifying which dataset to use
#' @param estimated_betas Matrix of estimated beta values (conditions x voxels)
#' @param method_name Character string describing the method used
#'
#' @return A list with performance metrics
#'
#' @examples
#' \dontrun{
#' # Load dataset and create design matrix
#' dataset <- load_benchmark_dataset("BM_Canonical_HighSNR")
#' X <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", fmrihrf::HRF_SPMG1)
#' 
#' # Fit simple linear model
#' betas <- solve(t(X) %*% X) %*% t(X) %*% dataset$Y_noisy
#' 
#' # Evaluate performance
#' performance <- evaluate_method_performance("BM_Canonical_HighSNR", 
#'                                           betas[-1, ], # Remove intercept
#'                                           "OLS")
#' }
#'
#' @export
evaluate_method_performance <- function(dataset_name, estimated_betas, method_name = "Unknown") {
  dataset <- load_benchmark_dataset(dataset_name)
  true_betas <- dataset$true_betas_condition
  
  # Check dimensions
  if (nrow(estimated_betas) != nrow(true_betas) || ncol(estimated_betas) != ncol(true_betas)) {
    stop("Estimated betas dimensions (", nrow(estimated_betas), " x ", ncol(estimated_betas), 
         ") don't match true betas dimensions (", nrow(true_betas), " x ", ncol(true_betas), ")")
  }
  
  # Overall metrics
  correlation <- .safe_cor(as.vector(true_betas), as.vector(estimated_betas))
  mse <- mean((true_betas - estimated_betas)^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(true_betas - estimated_betas))
  
  # Condition-specific metrics
  condition_metrics <- list()
  for (i in 1:nrow(true_betas)) {
    cond_corr <- .safe_cor(true_betas[i, ], estimated_betas[i, ])
    cond_mse <- mean((true_betas[i, ] - estimated_betas[i, ])^2)
    condition_metrics[[paste0("condition_", i)]] <- list(
      correlation = cond_corr,
      mse = cond_mse,
      rmse = sqrt(cond_mse)
    )
  }
  
  # Voxel-specific metrics (for each voxel across conditions)
  voxel_correlations <- numeric(ncol(true_betas))
  voxel_mse <- numeric(ncol(true_betas))
  for (j in 1:ncol(true_betas)) {
    voxel_correlations[j] <- .safe_cor(true_betas[, j], estimated_betas[, j])
    voxel_mse[j] <- mean((true_betas[, j] - estimated_betas[, j])^2)
  }
  
  list(
    method_name = method_name,
    dataset_name = dataset_name,
    overall_metrics = list(
      correlation = correlation,
      mse = mse,
      rmse = rmse,
      mae = mae
    ),
    condition_metrics = condition_metrics,
    voxel_metrics = list(
      correlations = voxel_correlations,
      mse_values = voxel_mse,
      mean_correlation = mean(voxel_correlations, na.rm = TRUE),
      mean_mse = mean(voxel_mse, na.rm = TRUE)
    )
  )
}

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# Environment to store the loaded and processed benchmark datasets
.benchmark_data_env <- new.env(parent = emptyenv())

# Complete recipes for every HRF used by the packaged benchmark artifact.
.benchmark_hrf_recipe_catalog <- function() {
  make_recipe <- function(recipe_name, base_hrf_name, lag = 0, width = 0,
                          normalize = FALSE) {
    list(
      schema_version = 1L,
      recipe_name = recipe_name,
      generator = "fmrihrf::gen_hrf",
      base_hrf_name = base_hrf_name,
      lag = as.numeric(lag),
      width = as.numeric(width),
      precision = 0.1,
      half_life = Inf,
      summate = TRUE,
      normalize = isTRUE(normalize),
      name_override = NULL,
      span_override = NULL,
      fmrihrf_version = as.character(utils::packageVersion("fmrihrf"))
    )
  }

  list(
    HRF_SPMG1 = make_recipe("HRF_SPMG1", "HRF_SPMG1"),
    HRF_SPMG2 = make_recipe("HRF_SPMG2", "HRF_SPMG2"),
    HRF_SPMG3 = make_recipe("HRF_SPMG3", "HRF_SPMG3"),
    variant1 = make_recipe(
      "variant1", "HRF_SPMG1", lag = 1, width = 1.2, normalize = TRUE
    ),
    variant2 = make_recipe(
      "variant2", "HRF_SPMG1", lag = -0.5, width = 0.8, normalize = TRUE
    )
  )
}

.benchmark_base_hrfs <- function() {
  list(
    HRF_SPMG1 = fmrihrf::HRF_SPMG1,
    HRF_SPMG2 = fmrihrf::HRF_SPMG2,
    HRF_SPMG3 = fmrihrf::HRF_SPMG3
  )
}

.same_recipe_value <- function(x, y) {
  isTRUE(all.equal(x, y, tolerance = 0, check.attributes = TRUE))
}

# Helper function to reconstruct HRF objects from a fail-closed recipe.
.reconstruct_hrf_object <- function(recipe) {
  if (!is.list(recipe) || is.null(names(recipe))) {
    stop("HRF recipe must be a named list.", call. = FALSE)
  }

  required <- c(
    "schema_version", "recipe_name", "generator", "base_hrf_name",
    "lag", "width", "precision", "half_life", "summate", "normalize",
    "name_override", "span_override", "fmrihrf_version"
  )
  missing_fields <- setdiff(required, names(recipe))
  if (length(missing_fields)) {
    stop(
      "HRF recipe is missing required fields: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }
  if (!identical(recipe$schema_version, 1L)) {
    stop("Unknown benchmark HRF recipe schema version.", call. = FALSE)
  }
  if (!identical(recipe$generator, "fmrihrf::gen_hrf")) {
    stop("Unknown benchmark HRF recipe generator.", call. = FALSE)
  }

  catalog <- .benchmark_hrf_recipe_catalog()
  if (!is.character(recipe$recipe_name) || length(recipe$recipe_name) != 1L ||
      !recipe$recipe_name %in% names(catalog)) {
    stop("Unknown benchmark HRF recipe name.", call. = FALSE)
  }

  base_hrfs <- .benchmark_base_hrfs()
  if (!is.character(recipe$base_hrf_name) ||
      length(recipe$base_hrf_name) != 1L ||
      !recipe$base_hrf_name %in% names(base_hrfs)) {
    stop("Unknown base HRF in benchmark recipe.", call. = FALSE)
  }

  expected <- catalog[[recipe$recipe_name]]
  semantic_fields <- setdiff(required, "fmrihrf_version")
  mismatched <- semantic_fields[!vapply(
    semantic_fields,
    function(field) .same_recipe_value(recipe[[field]], expected[[field]]),
    logical(1)
  )]
  if (length(mismatched)) {
    stop(
      "Benchmark HRF recipe '", recipe$recipe_name,
      "' has inconsistent fields: ", paste(mismatched, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.character(recipe$fmrihrf_version) ||
      length(recipe$fmrihrf_version) != 1L ||
      !nzchar(recipe$fmrihrf_version)) {
    stop("Benchmark HRF recipe has an invalid fmrihrf version.", call. = FALSE)
  }
  current_fmrihrf_version <- as.character(utils::packageVersion("fmrihrf"))
  if (!identical(recipe$fmrihrf_version, current_fmrihrf_version)) {
    stop(
      "Benchmark HRF recipe requires fmrihrf ", recipe$fmrihrf_version,
      "; installed version is ", current_fmrihrf_version, ".",
      call. = FALSE
    )
  }

  base_hrf <- base_hrfs[[recipe$base_hrf_name]]
  no_decoration <- recipe$lag == 0 && recipe$width == 0 &&
    !recipe$normalize && is.null(recipe$name_override) &&
    is.null(recipe$span_override)
  if (no_decoration) {
    return(base_hrf)
  }

  # fmrihrf 0.4.0's block decorator forwards the base HRF's P1/P2/A1
  # metadata to a wrapper that accepts only `t`, producing a misleading warning.
  # Rewrap the identical evaluator without those stale parameter attributes;
  # dense-grid tests verify that this warning-free route is numerically exact.
  clean_base <- fmrihrf::as_hrf(
    function(time) base_hrf(time),
    name = attr(base_hrf, "name"),
    nbasis = fmrihrf::nbasis(base_hrf),
    span = attr(base_hrf, "span")
  )

  withCallingHandlers(
    fmrihrf::gen_hrf(
      clean_base,
      lag = recipe$lag,
      width = recipe$width,
      precision = recipe$precision,
      half_life = recipe$half_life,
      summate = recipe$summate,
      normalize = recipe$normalize,
      name = recipe$name_override,
      span = recipe$span_override
    ),
    warning = function(condition) {
      stop(
        "Unexpected warning while reconstructing benchmark HRF recipe '",
        recipe$recipe_name, "': ", conditionMessage(condition),
        call. = FALSE
      )
    }
  )
}

# Recursive helper to process true_hrf_parameters
.process_hrf_params <- function(params_list) {
  if (is.list(params_list) && "hrf_object_name" %in% names(params_list)) {
    # Single HRF parameter set
    if (is.null(params_list$hrf_recipe)) {
      stop(
        "Missing reconstructable HRF recipe for '",
        params_list$hrf_object_name, "'.",
        call. = FALSE
      )
    }
    if (!identical(
      params_list$hrf_object_name,
      params_list$hrf_recipe$recipe_name
    )) {
      stop(
        "HRF object name does not match its stored recipe name.",
        call. = FALSE
      )
    }
    params_list$hrf_object <- .reconstruct_hrf_object(params_list$hrf_recipe)
    return(params_list)
  } else if (is.list(params_list)) {
    # Multiple HRF parameter sets (e.g., for different voxel groups)
    return(lapply(params_list, .process_hrf_params))
  } else {
    return(params_list)
  }
}

.load_and_reconstruct_raw_dataset <- function(ds_raw) {
  # Process the raw dataset to reconstruct HRF objects
  ds_processed <- ds_raw
  
  # Reconstruct HRF parameters if present
  if (!is.null(ds_raw$true_hrf_parameters)) {
    ds_processed$true_hrf_parameters <- .process_hrf_params(ds_raw$true_hrf_parameters)
  }
  
  return(ds_processed)
}

.ensure_benchmark_data_loaded <- function() {
  if (!exists("fmri_benchmark_datasets_processed", envir = .benchmark_data_env)) {
    # Load the raw data
    data("fmri_benchmark_datasets", package = "fmrireg", envir = environment())
    
    # Process each dataset
    processed_data <- list()
    for (name in names(fmri_benchmark_datasets)) {
      if (name == "metadata") {
        processed_data[[name]] <- fmri_benchmark_datasets[[name]]
      } else {
        processed_data[[name]] <- .load_and_reconstruct_raw_dataset(fmri_benchmark_datasets[[name]])
      }
    }
    
    # Store in environment
    .benchmark_data_env$fmri_benchmark_datasets_processed <- processed_data
  }
  
  return(.benchmark_data_env$fmri_benchmark_datasets_processed)
} 
