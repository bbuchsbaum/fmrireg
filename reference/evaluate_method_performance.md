# Evaluate Method Performance on Benchmark Dataset

Helper function to evaluate the performance of beta estimation methods
on benchmark datasets by comparing estimated betas to ground truth.

## Usage

``` r
evaluate_method_performance(
  dataset_name,
  estimated_betas,
  method_name = "Unknown"
)
```

## Arguments

- dataset_name:

  Character string specifying which dataset to use

- estimated_betas:

  Matrix of estimated beta values (conditions x voxels)

- method_name:

  Character string describing the method used

## Value

A list with performance metrics

## Examples

``` r
if (FALSE) { # \dontrun{
# Load dataset and create design matrix
dataset <- load_benchmark_dataset("BM_Canonical_HighSNR")
X <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", fmrihrf::HRF_SPMG1)

# Fit simple linear model
betas <- solve(t(X) %*% X) %*% t(X) %*% dataset$Y_noisy

# Evaluate performance
performance <- evaluate_method_performance("BM_Canonical_HighSNR", 
                                          betas[-1, ], # Remove intercept
                                          "OLS")
} # }
```
