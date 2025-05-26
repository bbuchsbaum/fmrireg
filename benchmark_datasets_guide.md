# fMRI Benchmark Datasets Guide

## Overview

The `fmrireg` package includes a comprehensive set of benchmark datasets designed for testing and evaluating HRF fitting, beta estimation, and other fMRI analysis methods. These datasets provide known ground truth for various challenging scenarios commonly encountered in fMRI analysis.

## Dataset Collection

The benchmark collection includes 5 carefully designed datasets, each targeting specific aspects of fMRI analysis:

### 1. BM_Canonical_HighSNR
**Purpose**: Test basic beta estimation accuracy when HRF is known and SNR is high

**Characteristics**:
- 3 experimental conditions with fixed amplitudes (1.0, 1.5, 0.8)
- 15 events per condition (45 total events)
- High signal-to-noise ratio (SNR = 4)
- Canonical SPM HRF (HRF_SPMG1)
- AR(1) noise with coefficient 0.4
- 5-minute scan duration, TR = 2s

**Use Cases**:
- Baseline performance testing
- Method validation under ideal conditions
- Comparing different beta estimation approaches

### 2. BM_Canonical_LowSNR
**Purpose**: Test beta estimation robustness with known HRF but low SNR

**Characteristics**:
- Same experimental design as BM_Canonical_HighSNR
- Low signal-to-noise ratio (SNR = 0.5)
- Tests method robustness to noise

**Use Cases**:
- Evaluating noise robustness
- Testing regularization methods
- Comparing denoising approaches

### 3. BM_HRF_Variability_AcrossVoxels
**Purpose**: Test HRF fitting when true HRF varies across voxels

**Characteristics**:
- 2 experimental conditions with amplitudes (1.2, 0.9)
- 15 events per condition (30 total events)
- 50% of voxels use canonical HRF, 50% use modified HRF
- Moderate SNR (SNR = 1.0)
- Tests spatial HRF heterogeneity

**Use Cases**:
- Testing HRF estimation methods
- Evaluating spatial modeling approaches
- Comparing single vs. multiple HRF models

### 4. BM_Trial_Amplitude_Variability
**Purpose**: Test LSS-style beta estimation with trial-to-trial amplitude changes

**Characteristics**:
- Single experimental condition
- 20 events with significant trial-to-trial amplitude variability (SD = 0.3)
- Base amplitude = 1.0 with Gaussian variability
- Canonical HRF
- Ideal for testing single-trial estimation methods

**Use Cases**:
- Testing Least Squares Separate (LSS) methods
- Evaluating trial-wise beta estimation
- Comparing single-trial vs. condition-level approaches

### 5. BM_Complex_Realistic
**Purpose**: A challenging dataset combining multiple factors

**Characteristics**:
- 3 experimental conditions with base amplitudes (1.0, 1.2, 0.9)
- 12 events per condition (36 total events)
- 3 different HRF variants across voxel groups
- Variable event durations (mean = 0.2s, SD = 0.1s)
- Trial-to-trial amplitude variability (SD = 0.15)
- AR(2) noise model
- Exponential inter-stimulus intervals

**Use Cases**:
- Comprehensive method evaluation
- Testing complex experimental designs
- Evaluating multi-factor models

## Data Structure

Each benchmark dataset contains the following components:

### Core Data Object (New)
- `core_data`: A `matrix_dataset` object (see `fmrireg::matrix_dataset`) encapsulating the primary BOLD data and experimental timing. It contains:
    - `datamat`: The noisy BOLD time series (time × voxels), identical to `Y_noisy`.
    - `TR`: Repetition time.
    - `run_length`: Length of the run.
    - `event_table`: A `data.frame` with `onset`, `condition`, and `duration` for events.
    - `sampling_frame`: An object detailing the time sampling for each run.

### Core Data (Direct Access - for convenience and backward compatibility)
- `Y_noisy`: Matrix of noisy BOLD time series (time × voxels)
- `Y_clean`: Matrix of clean BOLD time series (when available)
- `event_onsets`: Vector of event onset times (seconds)
- `condition_labels`: Vector of condition labels for each event

### Ground Truth Information
- `true_betas_condition`: Matrix of true condition-level beta values (conditions × voxels)
- `true_amplitudes_trial`: Matrix of true trial-level amplitudes (events × voxels)
- `true_hrf_parameters`: Information about the true HRF(s) used
- `X_list_true_hrf`: List of design matrices convolved with true HRF

### Experimental Parameters (Direct Access)
- `TR`: Repetition time (seconds) - also in `core_data`
- `total_time`: Total scan duration (seconds)
- `simulation_seed`: Random seed used for reproducibility
- `target_snr`: Target signal-to-noise ratio

### Noise Information
- `noise_parameters`: Details about noise generation (type, AR coefficients, SD)

## Usage Examples

### Basic Dataset Loading

```r
library(fmrireg)

# Load a specific dataset
data <- load_benchmark_dataset("BM_Canonical_HighSNR")

# Get BOLD time series
Y <- data$Y_noisy

# Get experimental design
onsets <- data$event_onsets
conditions <- data$condition_labels

# Get ground truth betas
true_betas <- data$true_betas_condition
```

### Exploring Available Datasets

```r
# List all available datasets
list_benchmark_datasets()

# Get detailed summary of a dataset
summary <- get_benchmark_summary("BM_Canonical_HighSNR")
print(summary)

# Load metadata
metadata <- load_benchmark_dataset("metadata")
```

### Creating Design Matrices

```r
# Create design matrix with canonical HRF
X_canonical <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", HRF_SPMG1)

# Test with different HRF
X_derivative <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", HRF_SPMG2)
```

### Method Evaluation

```r
# Load dataset
data <- load_benchmark_dataset("BM_Canonical_HighSNR")
X <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", HRF_SPMG1)

# Fit ordinary least squares
betas_ols <- solve(t(X) %*% X) %*% t(X) %*% data$Y_noisy

# Evaluate performance (remove intercept)
performance <- evaluate_method_performance("BM_Canonical_HighSNR", 
                                          betas_ols[-1, ], 
                                          "OLS")

# View results
print(performance$overall_metrics)
```

### HRF Estimation Testing

```r
# Load dataset with HRF variability
data <- load_benchmark_dataset("BM_HRF_Variability_AcrossVoxels")

# Get true HRF assignments
hrf_groups <- data$true_hrf_group_assignment

# Test HRF estimation method
# (Your HRF estimation code here)

# Compare estimated vs. true HRFs
true_canonical <- data$true_hrf_parameters$canonical$hrf_object
true_variant <- data$true_hrf_parameters$variant1$hrf_object
```

### Trial-wise Analysis

```r
# Load trial variability dataset
data <- load_benchmark_dataset("BM_Trial_Amplitude_Variability")

# Get true trial-wise amplitudes
true_trial_amps <- data$true_amplitudes_trial

# Implement LSS or other single-trial method
# (Your single-trial estimation code here)

# Compare estimated vs. true trial amplitudes
```

## Performance Metrics

The `evaluate_method_performance()` function provides comprehensive metrics:

### Overall Metrics
- **MSE**: Mean squared error
- **MAE**: Mean absolute error  
- **RMSE**: Root mean squared error
- **Correlation**: Overall correlation between estimated and true betas

### Condition-wise Metrics
- MSE and correlation for each experimental condition

### Voxel-wise Metrics
- MSE and correlation for each voxel

### Summary Statistics
- Mean correlations across conditions and voxels
- Min/max correlations

## Best Practices

### 1. Start Simple
Begin with `BM_Canonical_HighSNR` to establish baseline performance before testing more challenging scenarios.

### 2. Test Incrementally
Progress through datasets in order of complexity:
1. BM_Canonical_HighSNR (baseline)
2. BM_Canonical_LowSNR (noise robustness)
3. BM_HRF_Variability_AcrossVoxels (HRF estimation)
4. BM_Trial_Amplitude_Variability (single-trial methods)
5. BM_Complex_Realistic (comprehensive evaluation)

### 3. Use Multiple Metrics
Don't rely on a single performance metric. Consider MSE, correlation, and condition-specific performance.

### 4. Test HRF Assumptions
Use `create_design_matrix_from_benchmark()` with different HRFs to test sensitivity to HRF misspecification.

### 5. Validate Reproducibility
All datasets use fixed random seeds. Verify your results are reproducible by running analyses multiple times.

## Extending the Benchmark Suite

The benchmark generation framework is extensible. To add new datasets:

1. Modify `data-raw/generate_benchmark_datasets.R`
2. Add new generation functions following existing patterns
3. Update `load_benchmark_dataset()` and related functions
4. Regenerate the data with `source("data-raw/generate_benchmark_datasets.R")`

## Technical Details

### Simulation Parameters
- **Time Grid**: 0 to 300 seconds with TR = 2s (150 time points)
- **Voxels**: 100 simulated voxels per dataset
- **Buffer**: 16-second buffer to avoid edge effects
- **Noise Models**: AR(1) and AR(2) with realistic coefficients
- **HRF Variants**: Based on SPM canonical HRF with modifications

### Ground Truth Construction
- True betas represent the amplitude scaling factors for each condition
- Trial-wise amplitudes include both condition effects and trial-to-trial variability
- Design matrices are constructed using the exact HRF used for simulation
- Noise parameters are carefully calibrated to achieve target SNR levels

### Reproducibility
- All datasets use fixed random seeds
- Simulation parameters are stored with each dataset
- Package version and R version are recorded in metadata

## Citation

When using these benchmark datasets in publications, please cite the fmrireg package and mention the specific datasets used.

## Support

For questions about the benchmark datasets or suggestions for new scenarios, please open an issue on the fmrireg GitHub repository. 