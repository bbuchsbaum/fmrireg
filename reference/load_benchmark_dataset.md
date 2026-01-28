# Load fMRI Benchmark Datasets

This function provides easy access to the benchmark datasets included
with the fmrireg package. These datasets are designed for testing HRF
fitting, beta estimation, and other fMRI analysis methods.

## Usage

``` r
load_benchmark_dataset(dataset_name = "BM_Canonical_HighSNR")
```

## Arguments

- dataset_name:

  Character string specifying which dataset to load. Options include:

  - `"BM_Canonical_HighSNR"`: Canonical HRF with high SNR (3 conditions)

  - `"BM_Canonical_LowSNR"`: Canonical HRF with low SNR (3 conditions)

  - `"BM_HRF_Variability_AcrossVoxels"`: HRF varies across voxel groups
    (2 conditions)

  - `"BM_Trial_Amplitude_Variability"`: Trial-to-trial amplitude
    variability (1 condition)

  - `"BM_Complex_Realistic"`: Complex scenario with multiple factors (3
    conditions)

  - `"all"`: Returns all datasets as a list

  - `"metadata"`: Returns metadata about the datasets

## Value

A list containing the specified benchmark dataset(s) with the following
components:

- `description`: Text description of the dataset

- `Y_noisy`: Matrix of noisy BOLD time series (time x voxels)

- `Y_clean`: Matrix of clean BOLD time series (when available)

- `X_list_true_hrf`: List of design matrices convolved with true HRF

- `true_hrf_parameters`: Information about the true HRF(s) used

- `event_onsets`: Vector of event onset times

- `condition_labels`: Vector of condition labels for each event

- `true_betas_condition`: Matrix of true condition-level beta values

- `true_amplitudes_trial`: Matrix of true trial-level amplitudes

- `TR`: Repetition time

- `total_time`: Total scan duration

- `noise_parameters`: Information about noise generation

- `simulation_seed`: Random seed used for generation

- `target_snr`: Target signal-to-noise ratio

## Examples

``` r
# Load a specific dataset
high_snr_data <- load_benchmark_dataset("BM_Canonical_HighSNR")

# Get information about all available datasets
metadata <- load_benchmark_dataset("metadata")

# Load all datasets
all_data <- load_benchmark_dataset("all")

# Access the BOLD data
Y <- high_snr_data$Y_noisy

# Get event information
onsets <- high_snr_data$event_onsets
conditions <- high_snr_data$condition_labels
```
