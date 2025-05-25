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
#'     \item \code{true_hrf_parameters}: Information about the true HRF(s) used
#'     \item \code{event_onsets}: Vector of event onset times
#'     \item \code{condition_labels}: Vector of condition labels for each event
#'     \item \code{true_betas_condition}: Matrix of true condition-level beta values
#'     \item \code{true_amplitudes_trial}: Matrix of true trial-level amplitudes
#'     \item \code{TR}: Repetition time
#'     \item \code{total_time}: Total scan duration
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
  
  # Load the benchmark datasets (this will be available after package installation)
  if (!exists("fmri_benchmark_datasets")) {
    data("fmri_benchmark_datasets", package = "fmrireg", envir = environment())
  }
  
  available_datasets <- c("BM_Canonical_HighSNR", "BM_Canonical_LowSNR", 
                         "BM_HRF_Variability_AcrossVoxels", "BM_Trial_Amplitude_Variability",
                         "BM_Complex_Realistic", "all", "metadata")
  
  if (!dataset_name %in% available_datasets) {
    stop("Dataset '", dataset_name, "' not found. Available datasets: ", 
         paste(available_datasets, collapse = ", "))
  }
  
  if (dataset_name == "all") {
    # Return all datasets except metadata
    datasets <- fmri_benchmark_datasets
    datasets$metadata <- NULL
    return(datasets)
  } else if (dataset_name == "metadata") {
    return(fmri_benchmark_datasets$metadata)
  } else {
    return(fmri_benchmark_datasets[[dataset_name]])
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
  
  # Load the benchmark datasets
  if (!exists("fmri_benchmark_datasets")) {
    data("fmri_benchmark_datasets", package = "fmrireg", envir = environment())
  }
  
  dataset_names <- names(fmri_benchmark_datasets)
  dataset_names <- dataset_names[dataset_names != "metadata"]
  
  descriptions <- sapply(dataset_names, function(name) {
    fmri_benchmark_datasets[[name]]$description
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
  
  dataset <- load_benchmark_dataset(dataset_name)
  
  if (is.null(dataset)) {
    stop("Dataset not found: ", dataset_name)
  }
  
  # Basic dimensions
  n_timepoints <- nrow(dataset$Y_noisy)
  n_voxels <- ncol(dataset$Y_noisy)
  n_events <- length(dataset$event_onsets)
  
  # Experimental design
  unique_conditions <- unique(dataset$condition_labels)
  n_conditions <- length(unique_conditions)
  events_per_condition <- table(dataset$condition_labels)
  
  # HRF information
  hrf_info <- dataset$true_hrf_parameters
  
  # Noise information
  noise_info <- dataset$noise_parameters
  
  summary_list <- list(
    dataset_name = dataset_name,
    description = dataset$description,
    dimensions = list(
      n_timepoints = n_timepoints,
      n_voxels = n_voxels,
      n_events = n_events,
      n_conditions = n_conditions
    ),
    experimental_design = list(
      conditions = unique_conditions,
      events_per_condition = events_per_condition,
      TR = dataset$TR,
      total_time = dataset$total_time,
      scan_duration_minutes = dataset$total_time / 60
    ),
    hrf_information = hrf_info,
    noise_information = noise_info,
    target_snr = dataset$target_snr %||% "Not specified",
    simulation_seed = dataset$simulation_seed
  )
  
  return(summary_list)
}

#' Create Design Matrix from Benchmark Dataset
#'
#' Helper function to create a design matrix from a benchmark dataset using
#' a specified HRF. This is useful for testing different HRF assumptions
#' against the ground truth.
#'
#' @param dataset_name Character string specifying which dataset to use
#' @param hrf HRF object to use for convolution (e.g., HRF_SPMG1)
#' @param include_intercept Logical, whether to include an intercept column
#'
#' @return A matrix with the design matrix (time x conditions)
#'
#' @examples
#' # Create design matrix using canonical HRF
#' X <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", HRF_SPMG1)
#' 
#' # Test with a different HRF
#' X_wrong <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", HRF_SPMG2)
#'
#' @export
create_design_matrix_from_benchmark <- function(dataset_name, hrf, include_intercept = TRUE) {
  
  dataset <- load_benchmark_dataset(dataset_name)
  
  # Create time grid
  time_grid <- seq(0, dataset$total_time, by = dataset$TR)
  
  # Get unique conditions
  unique_conditions <- unique(dataset$condition_labels)
  n_conditions <- length(unique_conditions)
  
  # Create design matrix
  X_list <- list()
  
  for (i in seq_along(unique_conditions)) {
    cond <- unique_conditions[i]
    cond_onsets <- dataset$event_onsets[dataset$condition_labels == cond]
    
    if (length(cond_onsets) > 0) {
      # Assume duration of 0.5 if not specified in dataset
      duration <- if ("true_durations_trial" %in% names(dataset)) {
        mean(dataset$true_durations_trial[dataset$condition_labels == cond, 1])
      } else {
        0.5
      }
      
      reg <- regressor(cond_onsets, hrf, duration = duration, amplitude = 1)
      X_list[[cond]] <- evaluate(reg, time_grid)
    }
  }
  
  # Combine into matrix
  X <- do.call(cbind, X_list)
  colnames(X) <- names(X_list)
  
  # Add intercept if requested
  if (include_intercept) {
    X <- cbind(Intercept = 1, X)
  }
  
  return(X)
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
#' X <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", HRF_SPMG1)
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
  
  # Ensure dimensions match
  if (!all(dim(estimated_betas) == dim(true_betas))) {
    stop("Dimensions of estimated_betas (", paste(dim(estimated_betas), collapse = "x"), 
         ") do not match true_betas (", paste(dim(true_betas), collapse = "x"), ")")
  }
  
  # Calculate performance metrics
  mse <- mean((estimated_betas - true_betas)^2)
  mae <- mean(abs(estimated_betas - true_betas))
  correlation <- cor(as.vector(estimated_betas), as.vector(true_betas))
  
  # Per-condition metrics
  condition_mse <- apply((estimated_betas - true_betas)^2, 1, mean)
  condition_correlation <- sapply(1:nrow(true_betas), function(i) {
    cor(estimated_betas[i, ], true_betas[i, ])
  })
  
  # Per-voxel metrics
  voxel_mse <- apply((estimated_betas - true_betas)^2, 2, mean)
  voxel_correlation <- sapply(1:ncol(true_betas), function(i) {
    cor(estimated_betas[, i], true_betas[, i])
  })
  
  performance <- list(
    dataset_name = dataset_name,
    method_name = method_name,
    overall_metrics = list(
      mse = mse,
      mae = mae,
      correlation = correlation,
      rmse = sqrt(mse)
    ),
    condition_metrics = list(
      mse = condition_mse,
      correlation = condition_correlation
    ),
    voxel_metrics = list(
      mse = voxel_mse,
      correlation = voxel_correlation
    ),
    summary_stats = list(
      mean_condition_correlation = mean(condition_correlation, na.rm = TRUE),
      mean_voxel_correlation = mean(voxel_correlation, na.rm = TRUE),
      min_correlation = min(c(condition_correlation, voxel_correlation), na.rm = TRUE),
      max_correlation = max(c(condition_correlation, voxel_correlation), na.rm = TRUE)
    )
  )
  
  return(performance)
}

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x 