## Generate Benchmark Datasets for fmrireg Package
## This script creates a comprehensive set of benchmark datasets for testing
## HRF fitting, beta estimation, and other fMRI analysis methods.

library(fmrireg)

# Set master seed for reproducibility
set.seed(12345)

# Common parameters across all datasets
COMMON_PARAMS <- list(
  total_time = 300,  # 5 minutes
  TR = 2.0,
  n_voxels = 100,
  buffer = 16
)

# Helper function to calculate SNR-appropriate noise_sd
calculate_noise_sd <- function(signal_matrix, target_snr) {
  signal_sd <- sd(as.vector(signal_matrix))
  return(signal_sd / target_snr)
}

# Helper function to create HRF variants
create_hrf_variants <- function() {
  list(
    canonical = HRF_SPMG1,
    with_temporal_deriv = HRF_SPMG2,
    with_both_derivs = HRF_SPMG3,
    # Create a modified canonical HRF (shifted and scaled)
    variant1 = gen_hrf(HRF_SPMG1, lag = 1, width = 1.2, normalize = TRUE),
    variant2 = gen_hrf(HRF_SPMG1, lag = -0.5, width = 0.8, normalize = TRUE)
  )
}

# Get HRF variants
hrf_variants <- create_hrf_variants()

# Initialize benchmark datasets list
fmri_benchmark_datasets <- list()

#' Generate Benchmark Dataset 1: Canonical HRF, High SNR
#' Goal: Test basic beta estimation accuracy when HRF is known and SNR is high
generate_BM_Canonical_HighSNR <- function() {
  cat("Generating BM_Canonical_HighSNR...\n")
  
  # Generate clean signal first to calculate appropriate noise level
  clean_sim <- simulate_fmri_matrix(
    n = COMMON_PARAMS$n_voxels,
    total_time = COMMON_PARAMS$total_time,
    TR = COMMON_PARAMS$TR,
    hrf = hrf_variants$canonical,
    n_events = 45,  # 15 events per condition * 3 conditions
    isi_dist = "uniform",
    isi_min = 4,
    isi_max = 8,
    durations = 0.5,
    duration_sd = 0,
    amplitudes = rep(c(1.0, 1.5, 0.8), each = 15),  # 3 conditions
    amplitude_sd = 0,
    noise_type = "none",
    random_seed = 12345
  )
  
  # Calculate noise for high SNR (signal/noise ratio = 4)
  noise_sd <- calculate_noise_sd(clean_sim$time_series$datamat, target_snr = 4)
  
  # Generate final dataset with noise
  final_sim <- simulate_fmri_matrix(
    n = COMMON_PARAMS$n_voxels,
    total_time = COMMON_PARAMS$total_time,
    TR = COMMON_PARAMS$TR,
    hrf = hrf_variants$canonical,
    n_events = 45,
    isi_dist = "uniform",
    isi_min = 4,
    isi_max = 8,
    durations = 0.5,
    duration_sd = 0,
    amplitudes = rep(c(1.0, 1.5, 0.8), each = 15),
    amplitude_sd = 0,
    noise_type = "ar1",
    noise_ar = 0.4,
    noise_sd = noise_sd,
    random_seed = 12345
  )
  
  # Create condition labels
  condition_labels <- rep(c("Cond1", "Cond2", "Cond3"), each = 15)
  
  # Create design matrix with true HRF
  event_onsets <- final_sim$time_series$event_table$onset
  time_grid <- seq(0, COMMON_PARAMS$total_time, by = COMMON_PARAMS$TR)
  
  # Create condition-specific regressors
  X_list_true_hrf <- list()
  for (i in 1:3) {
    cond_onsets <- event_onsets[condition_labels == paste0("Cond", i)]
    if (length(cond_onsets) > 0) {
      reg <- regressor(cond_onsets, hrf_variants$canonical, duration = 0.5, amplitude = 1)
      X_list_true_hrf[[paste0("Cond", i)]] <- evaluate(reg, time_grid)
    }
  }
  
  # Create event_table for matrix_dataset
  benchmark_event_table <- data.frame(
    onset = event_onsets,
    condition = condition_labels,
    duration = 0.5 # Assuming fixed duration for this benchmark
  )
  
  # Create matrix_dataset
  core_data <- matrix_dataset(
    datamat = final_sim$time_series$datamat,
    TR = COMMON_PARAMS$TR,
    run_length = COMMON_PARAMS$total_time / COMMON_PARAMS$TR,
    event_table = benchmark_event_table
  )
  
  list(
    description = "Canonical HRF (SPMG1), high SNR, 3 conditions, fixed amplitudes per condition",
    core_data = core_data, # Encapsulated data
    Y_noisy = final_sim$time_series$datamat, # Still available for direct access
    Y_clean = clean_sim$time_series$datamat,
    X_list_true_hrf = X_list_true_hrf,
    true_hrf_parameters = list(type = "SPMG1", hrf_object = hrf_variants$canonical),
    event_onsets = event_onsets,
    condition_labels = condition_labels,
    true_betas_condition = matrix(rep(c(1.0, 1.5, 0.8), each = COMMON_PARAMS$n_voxels), 
                                 nrow = 3, byrow = TRUE),
    true_amplitudes_trial = final_sim$ampmat,
    TR = COMMON_PARAMS$TR,
    total_time = COMMON_PARAMS$total_time,
    noise_parameters = final_sim$noise_params,
    simulation_seed = 12345,
    target_snr = 4
  )
}

#' Generate Benchmark Dataset 2: Canonical HRF, Low SNR
#' Goal: Test beta estimation robustness with known HRF but low SNR
generate_BM_Canonical_LowSNR <- function() {
  cat("Generating BM_Canonical_LowSNR...\n")
  
  # Generate clean signal first
  clean_sim <- simulate_fmri_matrix(
    n = COMMON_PARAMS$n_voxels,
    total_time = COMMON_PARAMS$total_time,
    TR = COMMON_PARAMS$TR,
    hrf = hrf_variants$canonical,
    n_events = 45,
    isi_dist = "uniform",
    isi_min = 4,
    isi_max = 8,
    durations = 0.5,
    duration_sd = 0,
    amplitudes = rep(c(1.0, 1.5, 0.8), each = 15),
    amplitude_sd = 0,
    noise_type = "none",
    random_seed = 12346
  )
  
  # Calculate noise for low SNR (signal/noise ratio = 0.5)
  noise_sd <- calculate_noise_sd(clean_sim$time_series$datamat, target_snr = 0.5)
  
  # Generate final dataset with noise
  final_sim <- simulate_fmri_matrix(
    n = COMMON_PARAMS$n_voxels,
    total_time = COMMON_PARAMS$total_time,
    TR = COMMON_PARAMS$TR,
    hrf = hrf_variants$canonical,
    n_events = 45,
    isi_dist = "uniform",
    isi_min = 4,
    isi_max = 8,
    durations = 0.5,
    duration_sd = 0,
    amplitudes = rep(c(1.0, 1.5, 0.8), each = 15),
    amplitude_sd = 0,
    noise_type = "ar1",
    noise_ar = 0.4,
    noise_sd = noise_sd,
    random_seed = 12346
  )
  
  condition_labels <- rep(c("Cond1", "Cond2", "Cond3"), each = 15)
  event_onsets <- final_sim$time_series$event_table$onset
  time_grid <- seq(0, COMMON_PARAMS$total_time, by = COMMON_PARAMS$TR)
  
  X_list_true_hrf <- list()
  for (i in 1:3) {
    cond_onsets <- event_onsets[condition_labels == paste0("Cond", i)]
    if (length(cond_onsets) > 0) {
      reg <- regressor(cond_onsets, hrf_variants$canonical, duration = 0.5, amplitude = 1)
      X_list_true_hrf[[paste0("Cond", i)]] <- evaluate(reg, time_grid)
    }
  }
  
  # Create event_table for matrix_dataset
  benchmark_event_table <- data.frame(
    onset = event_onsets,
    condition = condition_labels,
    duration = 0.5 # Assuming fixed duration for this benchmark
  )
  
  # Create matrix_dataset
  core_data <- matrix_dataset(
    datamat = final_sim$time_series$datamat,
    TR = COMMON_PARAMS$TR,
    run_length = COMMON_PARAMS$total_time / COMMON_PARAMS$TR,
    event_table = benchmark_event_table
  )
  
  list(
    description = "Canonical HRF (SPMG1), low SNR, 3 conditions, fixed amplitudes per condition",
    core_data = core_data,
    Y_noisy = final_sim$time_series$datamat,
    Y_clean = clean_sim$time_series$datamat,
    X_list_true_hrf = X_list_true_hrf,
    true_hrf_parameters = list(type = "SPMG1", hrf_object = hrf_variants$canonical),
    event_onsets = event_onsets,
    condition_labels = condition_labels,
    true_betas_condition = matrix(rep(c(1.0, 1.5, 0.8), each = COMMON_PARAMS$n_voxels), 
                                 nrow = 3, byrow = TRUE),
    true_amplitudes_trial = final_sim$ampmat,
    TR = COMMON_PARAMS$TR,
    total_time = COMMON_PARAMS$total_time,
    noise_parameters = final_sim$noise_params,
    simulation_seed = 12346,
    target_snr = 0.5
  )
}

#' Generate Benchmark Dataset 3: HRF Variability Across Voxels
#' Goal: Test HRF fitting when true HRF varies across voxels
generate_BM_HRF_Variability_AcrossVoxels <- function() {
  cat("Generating BM_HRF_Variability_AcrossVoxels...\n")
  
  n_voxels_per_group <- COMMON_PARAMS$n_voxels / 2
  
  # Generate data for first group of voxels (canonical HRF)
  sim_group1 <- simulate_fmri_matrix(
    n = n_voxels_per_group,
    total_time = COMMON_PARAMS$total_time,
    TR = COMMON_PARAMS$TR,
    hrf = hrf_variants$canonical,
    n_events = 30,  # 15 events per condition * 2 conditions
    isi_dist = "uniform",
    isi_min = 4,
    isi_max = 8,
    durations = 0.5,
    duration_sd = 0,
    amplitudes = rep(c(1.2, 0.9), each = 15),
    amplitude_sd = 0,
    noise_type = "ar1",
    noise_ar = 0.4,
    noise_sd = 1.0,
    random_seed = 12347
  )
  
  # Generate data for second group of voxels (variant HRF)
  sim_group2 <- simulate_fmri_matrix(
    n = n_voxels_per_group,
    total_time = COMMON_PARAMS$total_time,
    TR = COMMON_PARAMS$TR,
    hrf = hrf_variants$variant1,
    n_events = 30,
    isi_dist = "uniform",
    isi_min = 4,
    isi_max = 8,
    durations = 0.5,
    duration_sd = 0,
    amplitudes = rep(c(1.2, 0.9), each = 15),
    amplitude_sd = 0,
    noise_type = "ar1",
    noise_ar = 0.4,
    noise_sd = 1.0,
    random_seed = 12347  # Same seed for same event timing
  )
  
  # Combine the data
  Y_combined <- cbind(sim_group1$time_series$datamat, sim_group2$time_series$datamat)
  ampmat_combined <- cbind(sim_group1$ampmat, sim_group2$ampmat)
  
  condition_labels <- rep(c("Cond1", "Cond2"), each = 15)
  event_onsets <- sim_group1$time_series$event_table$onset
  
  # Create HRF group assignment
  hrf_group_assignment <- c(rep("canonical", n_voxels_per_group), 
                           rep("variant1", n_voxels_per_group))
  
  # Create event_table for matrix_dataset
  benchmark_event_table <- data.frame(
    onset = event_onsets,
    condition = condition_labels,
    duration = 0.5 # Assuming fixed duration for this benchmark
  )
  
  # Create matrix_dataset
  core_data <- matrix_dataset(
    datamat = Y_combined,
    TR = COMMON_PARAMS$TR,
    run_length = COMMON_PARAMS$total_time / COMMON_PARAMS$TR,
    event_table = benchmark_event_table
  )
  
  list(
    description = "HRF varies across voxel groups, 2 conditions, moderate SNR",
    core_data = core_data,
    Y_noisy = Y_combined,
    true_hrf_parameters = list(
      canonical = list(type = "SPMG1", hrf_object = hrf_variants$canonical),
      variant1 = list(type = "SPMG1_modified", hrf_object = hrf_variants$variant1)
    ),
    event_onsets = event_onsets,
    condition_labels = condition_labels,
    true_betas_condition = matrix(rep(c(1.2, 0.9), each = COMMON_PARAMS$n_voxels), 
                                 nrow = 2, byrow = TRUE),
    true_amplitudes_trial = ampmat_combined,
    true_hrf_group_assignment = hrf_group_assignment,
    TR = COMMON_PARAMS$TR,
    total_time = COMMON_PARAMS$total_time,
    noise_parameters = sim_group1$noise_params,
    simulation_seed = 12347,
    target_snr = 1.0
  )
}

#' Generate Benchmark Dataset 4: Trial Amplitude Variability
#' Goal: Test LSS-style beta estimation with trial-to-trial amplitude changes
generate_BM_Trial_Amplitude_Variability <- function() {
  cat("Generating BM_Trial_Amplitude_Variability...\n")
  
  final_sim <- simulate_fmri_matrix(
    n = COMMON_PARAMS$n_voxels,
    total_time = COMMON_PARAMS$total_time,
    TR = COMMON_PARAMS$TR,
    hrf = hrf_variants$canonical,
    n_events = 20,
    isi_dist = "uniform",
    isi_min = 6,
    isi_max = 12,
    durations = 0.5,
    duration_sd = 0,
    amplitudes = 1.0,  # Base amplitude
    amplitude_sd = 0.3,  # Significant trial-to-trial variability
    amplitude_dist = "gaussian",
    noise_type = "ar1",
    noise_ar = 0.4,
    noise_sd = 1.0,
    random_seed = 12348
  )
  
  condition_labels <- rep("Cond1", 20)
  event_onsets <- final_sim$time_series$event_table$onset
  
  # Create event_table for matrix_dataset
  benchmark_event_table <- data.frame(
    onset = event_onsets,
    condition = condition_labels,
    duration = 0.5 # Assuming fixed duration for this benchmark
  )
  
  # Create matrix_dataset
  core_data <- matrix_dataset(
    datamat = final_sim$time_series$datamat,
    TR = COMMON_PARAMS$TR,
    run_length = COMMON_PARAMS$total_time / COMMON_PARAMS$TR,
    event_table = benchmark_event_table
  )
  
  list(
    description = "Single condition with significant trial-to-trial amplitude variability",
    core_data = core_data,
    Y_noisy = final_sim$time_series$datamat,
    true_hrf_parameters = list(type = "SPMG1", hrf_object = hrf_variants$canonical),
    event_onsets = event_onsets,
    condition_labels = condition_labels,
    true_betas_condition = matrix(1.0, nrow = 1, ncol = COMMON_PARAMS$n_voxels),
    true_amplitudes_trial = final_sim$ampmat,  # This is the key ground truth for LSS
    TR = COMMON_PARAMS$TR,
    total_time = COMMON_PARAMS$total_time,
    noise_parameters = final_sim$noise_params,
    simulation_seed = 12348,
    target_snr = 1.0
  )
}

#' Generate Benchmark Dataset 5: Complex Realistic Scenario
#' Goal: A challenging dataset combining multiple factors
generate_BM_Complex_Realistic <- function() {
  cat("Generating BM_Complex_Realistic...\n")
  
  # Divide voxels into 3 groups with different HRFs
  n_voxels_per_group <- COMMON_PARAMS$n_voxels / 3
  
  # Group 1: Canonical HRF
  sim_group1 <- simulate_fmri_matrix(
    n = n_voxels_per_group,
    total_time = COMMON_PARAMS$total_time,
    TR = COMMON_PARAMS$TR,
    hrf = hrf_variants$canonical,
    n_events = 36,  # 12 events per condition * 3 conditions
    isi_dist = "exponential",
    isi_min = 2,
    isi_rate = 0.2,
    durations = 0.2,
    duration_sd = 0.1,
    amplitudes = rep(c(1.0, 1.2, 0.9), each = 12),
    amplitude_sd = 0.15,
    amplitude_dist = "gaussian",
    noise_type = "ar2",
    noise_ar = c(0.3, 0.15),
    noise_sd = 1.2,
    random_seed = 12349
  )
  
  # Group 2: Variant1 HRF
  sim_group2 <- simulate_fmri_matrix(
    n = n_voxels_per_group,
    total_time = COMMON_PARAMS$total_time,
    TR = COMMON_PARAMS$TR,
    hrf = hrf_variants$variant1,
    n_events = 36,
    isi_dist = "exponential",
    isi_min = 2,
    isi_rate = 0.2,
    durations = 0.2,
    duration_sd = 0.1,
    amplitudes = rep(c(1.0, 1.2, 0.9), each = 12),
    amplitude_sd = 0.15,
    amplitude_dist = "gaussian",
    noise_type = "ar2",
    noise_ar = c(0.3, 0.15),
    noise_sd = 1.2,
    random_seed = 12349  # Same seed for consistent event timing
  )
  
  # Group 3: Variant2 HRF
  sim_group3 <- simulate_fmri_matrix(
    n = n_voxels_per_group,
    total_time = COMMON_PARAMS$total_time,
    TR = COMMON_PARAMS$TR,
    hrf = hrf_variants$variant2,
    n_events = 36,
    isi_dist = "exponential",
    isi_min = 2,
    isi_rate = 0.2,
    durations = 0.2,
    duration_sd = 0.1,
    amplitudes = rep(c(1.0, 1.2, 0.9), each = 12),
    amplitude_sd = 0.15,
    amplitude_dist = "gaussian",
    noise_type = "ar2",
    noise_ar = c(0.3, 0.15),
    noise_sd = 1.2,
    random_seed = 12349
  )
  
  # Combine all groups
  Y_combined <- cbind(sim_group1$time_series$datamat, 
                     sim_group2$time_series$datamat,
                     sim_group3$time_series$datamat)
  ampmat_combined <- cbind(sim_group1$ampmat, sim_group2$ampmat, sim_group3$ampmat)
  durmat_combined <- cbind(sim_group1$durmat, sim_group2$durmat, sim_group3$durmat)
  
  condition_labels <- rep(c("Cond1", "Cond2", "Cond3"), each = 12)
  event_onsets <- sim_group1$time_series$event_table$onset
  
  hrf_group_assignment <- c(rep("canonical", n_voxels_per_group),
                           rep("variant1", n_voxels_per_group),
                           rep("variant2", n_voxels_per_group))
  
  # Create event_table for matrix_dataset
  # Use actual event durations from sim_group1 (they are the same across groups for this sim)
  benchmark_event_table <- data.frame(
    onset = event_onsets,
    condition = condition_labels,
    duration = sim_group1$time_series$event_table$duration 
  )
  
  # Create matrix_dataset
  core_data <- matrix_dataset(
    datamat = Y_combined,
    TR = COMMON_PARAMS$TR,
    run_length = COMMON_PARAMS$total_time / COMMON_PARAMS$TR,
    event_table = benchmark_event_table
  )
  
  list(
    description = "Complex realistic scenario: 3 HRF groups, 3 conditions, variable durations/amplitudes, AR(2) noise",
    core_data = core_data,
    Y_noisy = Y_combined,
    true_hrf_parameters = list(
      canonical = list(type = "SPMG1", hrf_object = hrf_variants$canonical),
      variant1 = list(type = "SPMG1_modified1", hrf_object = hrf_variants$variant1),
      variant2 = list(type = "SPMG1_modified2", hrf_object = hrf_variants$variant2)
    ),
    event_onsets = event_onsets,
    condition_labels = condition_labels,
    true_betas_condition = matrix(rep(c(1.0, 1.2, 0.9), each = COMMON_PARAMS$n_voxels), 
                                 nrow = 3, byrow = TRUE),
    true_amplitudes_trial = ampmat_combined,
    true_durations_trial = durmat_combined,
    true_hrf_group_assignment = hrf_group_assignment,
    TR = COMMON_PARAMS$TR,
    total_time = COMMON_PARAMS$total_time,
    noise_parameters = sim_group1$noise_params,
    simulation_seed = 12349,
    target_snr = 0.8
  )
}

# Generate all benchmark datasets
cat("Starting benchmark dataset generation...\n")

fmri_benchmark_datasets$BM_Canonical_HighSNR <- generate_BM_Canonical_HighSNR()
fmri_benchmark_datasets$BM_Canonical_LowSNR <- generate_BM_Canonical_LowSNR()
fmri_benchmark_datasets$BM_HRF_Variability_AcrossVoxels <- generate_BM_HRF_Variability_AcrossVoxels()
fmri_benchmark_datasets$BM_Trial_Amplitude_Variability <- generate_BM_Trial_Amplitude_Variability()
fmri_benchmark_datasets$BM_Complex_Realistic <- generate_BM_Complex_Realistic()

# Add metadata
fmri_benchmark_datasets$metadata <- list(
  creation_date = Sys.Date(),
  fmrireg_version = packageVersion("fmrireg"),
  r_version = R.version.string,
  description = "Benchmark datasets for testing HRF fitting and beta estimation methods",
  common_parameters = COMMON_PARAMS,
  hrf_variants_used = names(hrf_variants)
)

cat("Benchmark dataset generation complete!\n")
cat("Generated", length(fmri_benchmark_datasets) - 1, "benchmark datasets\n")

# Save the datasets
usethis::use_data(fmri_benchmark_datasets, overwrite = TRUE)

cat("Datasets saved to data/fmri_benchmark_datasets.rda\n") 