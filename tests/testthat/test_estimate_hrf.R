# Helper function to generate event table
generate_event_table <- function(onsets, duration=0, level=1) {
  data.frame(
    onset = onsets,
    duration = rep(duration, length(onsets)),
    stem = rep(level, length(onsets)),
    run = rep(1, length(onsets))
  )
}

# Helper function to simulate fMRI data
simulate_fmri_data <- function(onsets, n_voxels, n_timepoints, hrf_means, noise_sd=0.1) {
  timepoints <- seq(0, n_timepoints - 1)
  
  # Generate HRFs for each onset with different means
  hrfs <- lapply(hrf_means, function(mean) {
    h <- gen_hrf(hrf_gaussian, mean=mean)
    reg <- regressor(onsets = onsets, hrf = h, duration = 0, amplitude = 1)
    evaluate(reg, timepoints)
  })
  
  # Generate matrix dataset
  matrix_dataset <- sapply(1:n_voxels, function(i) {
    signal <- hrfs[[i]]
    noise <- rnorm(n_timepoints, mean=0, sd=noise_sd)
    signal + noise
  })
  
  colnames(matrix_dataset) <- paste0("Voxel_", 1:n_voxels)
  rownames(matrix_dataset) <- paste0("Time_", 1:n_timepoints)
  
  matrix_dataset
}

# Helper function to create matrix dataset object
create_matrix_dataset <- function(datamat, TR=1, run_length, event_table, sampling_frame) {
  matrix_dataset(
    datamat = datamat,
    TR = TR,
    run_length = run_length,
    event_table = event_table
  )
}

# Test function for single HRF
# test_find_best_hrf <- function(n_voxels = 100, n_timepoints = 100, onsets = seq(10, 100, by=10), hrf_fun = hrf_gaussian, hrflib=NULL) {
#   # Generate event table
#   event_table <- generate_event_table(onsets)
#   
#   # Simulate fMRI data
#   hrf_means <- rep(5, n_voxels)  # All voxels have the same HRF mean for test 1
#   datamat <- simulate_fmri_data(onsets, n_voxels, n_timepoints, hrf_means)
#   
#   # Create sampling frame
#   sampling_frame <- sampling_frame(blocklens=c(n_timepoints), TR=1)
#   
#   # Combine into a dataset
#   dataset <- create_matrix_dataset(datamat, TR=1, run_length=n_timepoints, event_table=event_table, sampling_frame=sampling_frame)
#   
#   # Run the find_best_hrf function
#   result <- find_best_hrf.matrix_dataset(dataset, fac_var="stem", onset_var="onset", hrflib=hrflib, block=~run)
#   
#   expect_true(!is.null(result))
# }

# Test function for multiple HRFs
# test_find_best_hrf_multiple_means <- function(n_voxels = 100, n_timepoints = 100, onsets = seq(10, 100, by=10), 
#                                               hrf_means = c(5, 6, 7), hrflib=NULL) {
#   # Generate event table
#   event_table <- generate_event_table(onsets)
#   
#   # Adjust hrf_means to ensure 33, 33, 34 distribution
#   hrf_means <- rep(hrf_means, length.out=n_voxels)
#   
#   # Simulate fMRI data
#   datamat <- simulate_fmri_data(onsets, n_voxels, n_timepoints, hrf_means)
#   
#   # Create sampling frame
#   sampling_frame <- sampling_frame(blocklens=c(n_timepoints), TR=1)
#   
#   # Combine into a dataset
#   dataset <- create_matrix_dataset(datamat, TR=1, run_length=n_timepoints, event_table=event_table, sampling_frame=sampling_frame)
#   
#   # Run the find_best_hrf function
#   result <- find_best_hrf.matrix_dataset(dataset, fac_var="stem", onset_var="onset", hrflib=hrflib, block=~run)
#   
#   expect_true(!is.null(result))
# }

# Define tests
# test_that("find_best_hrf.matrix_dataset works with a single HRF", {
#   result <- test_find_best_hrf()
#   expect_true(!is.null(result))
# })

# test_that("find_best_hrf.matrix_dataset works with multiple HRFs", {
#   result <- test_find_best_hrf_multiple_means()
#   expect_true(!is.null(result))
# })