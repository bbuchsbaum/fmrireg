# fMRI Benchmark Datasets Integration

## Overview

This proposal integrates a comprehensive benchmark dataset system into the `fmrireg` package for testing HRF fitting and beta estimation methods. The system leverages the existing excellent simulation functions (`simulate_fmri_matrix`, `simulate_simple_dataset`) to create standardized, reproducible benchmark scenarios.

## Key Features

### 🎯 **5 Carefully Designed Benchmark Scenarios**
1. **BM_Canonical_HighSNR**: Baseline testing with known HRF and high SNR
2. **BM_Canonical_LowSNR**: Noise robustness testing with low SNR
3. **BM_HRF_Variability_AcrossVoxels**: HRF estimation with spatial heterogeneity
4. **BM_Trial_Amplitude_Variability**: Single-trial methods with trial-to-trial variability
5. **BM_Complex_Realistic**: Comprehensive testing with multiple challenging factors

### 🔧 **Easy-to-Use API**
```r
# Load any benchmark dataset
data <- load_benchmark_dataset("BM_Canonical_HighSNR")

# List all available datasets
list_benchmark_datasets()

# Get detailed summary
get_benchmark_summary("BM_Canonical_HighSNR")

# Create design matrices with different HRF assumptions
X <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", HRF_SPMG1)

# Evaluate method performance
performance <- evaluate_method_performance("BM_Canonical_HighSNR", estimated_betas, "MyMethod")
```

### 📊 **Comprehensive Ground Truth**
Each dataset includes:
- **BOLD time series** (noisy and clean versions)
- **True HRF parameters** and objects
- **True beta values** (condition-level and trial-level)
- **Event timing** and condition labels
- **Noise parameters** and simulation settings
- **Design matrices** convolved with true HRFs

### 🔬 **Rigorous Evaluation Framework**
- Multiple performance metrics (MSE, MAE, correlation, RMSE)
- Condition-wise and voxel-wise analysis
- Built-in comparison tools
- Standardized evaluation procedures

## Implementation Structure

### Files Created/Modified

1. **`data-raw/generate_benchmark_datasets.R`**
   - Comprehensive dataset generation script
   - Uses existing `simulate_fmri_matrix` function
   - Creates 5 benchmark scenarios with different challenges
   - Includes proper SNR calibration and ground truth recording

2. **`R/benchmark_datasets.R`**
   - User-facing API functions
   - `load_benchmark_dataset()`, `list_benchmark_datasets()`
   - `get_benchmark_summary()`, `create_design_matrix_from_benchmark()`
   - `evaluate_method_performance()`

3. **`vignettes/benchmark_datasets.Rmd`**
   - Comprehensive tutorial and examples
   - Demonstrates all major use cases
   - Shows method comparison workflows

4. **`inst/doc/benchmark_datasets_guide.md`**
   - Detailed documentation
   - Technical specifications
   - Best practices and guidelines

5. **`tests/testthat/test_benchmark_datasets.R`**
   - Comprehensive test suite
   - Ensures all functions work correctly
   - Tests error handling and edge cases

6. **`data-raw/DATASET.R`** (updated)
   - Sources the generation script
   - Integrates with standard R package data workflow

## Technical Specifications

### Common Parameters
- **Duration**: 5 minutes (300 seconds)
- **TR**: 2.0 seconds (150 time points)
- **Voxels**: 100 simulated voxels
- **Noise Models**: AR(1) and AR(2) with realistic coefficients
- **HRF Variants**: Based on SPM canonical with modifications
- **Reproducibility**: Fixed random seeds for all datasets

### Dataset-Specific Details

| Dataset | Conditions | Events | SNR | HRF Variants | Special Features |
|---------|------------|--------|-----|--------------|------------------|
| BM_Canonical_HighSNR | 3 | 45 | 4.0 | 1 (canonical) | Baseline testing |
| BM_Canonical_LowSNR | 3 | 45 | 0.5 | 1 (canonical) | Noise robustness |
| BM_HRF_Variability | 2 | 30 | 1.0 | 2 (spatial) | HRF estimation |
| BM_Trial_Variability | 1 | 20 | 1.0 | 1 (canonical) | Trial-wise analysis |
| BM_Complex_Realistic | 3 | 36 | 0.8 | 3 (spatial) | Multi-factor challenge |

## Integration Benefits

### 1. **Leverages Existing Infrastructure**
- Uses proven `simulate_fmri_matrix` function
- Integrates seamlessly with existing HRF objects (`HRF_SPMG1`, `HRF_SPMG2`, etc.)
- Follows established package conventions

### 2. **Standardized Evaluation**
- Consistent benchmarking across methods
- Fair comparison framework
- Reproducible results

### 3. **Educational Value**
- Clear examples of proper fMRI simulation
- Demonstrates package capabilities
- Teaching tool for fMRI analysis concepts

### 4. **Research Enablement**
- Facilitates method development
- Enables rigorous validation
- Supports publication-quality evaluations

## Usage Examples

### Basic Method Testing
```r
# Load dataset
data <- load_benchmark_dataset("BM_Canonical_HighSNR")

# Create design matrix
X <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", HRF_SPMG1)

# Fit your method
betas <- your_method(data$Y_noisy, X)

# Evaluate performance
performance <- evaluate_method_performance("BM_Canonical_HighSNR", betas, "YourMethod")
print(performance$overall_metrics)
```

### HRF Sensitivity Testing
```r
# Test sensitivity to HRF misspecification
X_correct <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", HRF_SPMG1)
X_wrong <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", HRF_SPMG2)

# Compare performance
perf_correct <- evaluate_method_performance("BM_Canonical_HighSNR", 
                                           fit_model(X_correct), "Correct_HRF")
perf_wrong <- evaluate_method_performance("BM_Canonical_HighSNR", 
                                         fit_model(X_wrong), "Wrong_HRF")
```

### Comprehensive Method Evaluation
```r
# Test across all benchmark scenarios
datasets <- c("BM_Canonical_HighSNR", "BM_Canonical_LowSNR", 
              "BM_HRF_Variability_AcrossVoxels", "BM_Trial_Amplitude_Variability",
              "BM_Complex_Realistic")

results <- lapply(datasets, function(ds) {
  data <- load_benchmark_dataset(ds)
  X <- create_design_matrix_from_benchmark(ds, HRF_SPMG1)
  betas <- your_method(data$Y_noisy, X)
  evaluate_method_performance(ds, betas, "YourMethod")
})
```

## Future Extensions

The framework is designed to be easily extensible:

1. **Additional Scenarios**: New benchmark datasets can be added by extending the generation script
2. **Different HRF Models**: Support for other HRF families (Glover, AFNI, etc.)
3. **Multi-run Datasets**: Extension to multi-session scenarios
4. **Real Data Integration**: Hybrid simulated/real data benchmarks
5. **Specialized Challenges**: Task-specific or population-specific scenarios

## Conclusion

This benchmark dataset integration provides the `fmrireg` package with a powerful, standardized framework for method evaluation. It leverages existing simulation capabilities while adding comprehensive ground truth tracking and evaluation tools. The system enables rigorous, reproducible testing of HRF fitting and beta estimation methods across a range of challenging but realistic scenarios.

The implementation follows R package best practices, integrates seamlessly with existing code, and provides both novice-friendly and advanced-user interfaces. This positions `fmrireg` as not just a method implementation package, but as a comprehensive platform for fMRI analysis method development and validation. 