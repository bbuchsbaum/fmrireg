# Create Design Matrix from Benchmark Dataset

Helper function to create a design matrix from a benchmark dataset using
a specified HRF. This is useful for testing different HRF assumptions
against the ground truth.

## Usage

``` r
create_design_matrix_from_benchmark(
  dataset_name,
  hrf,
  include_intercept = TRUE
)
```

## Arguments

- dataset_name:

  Character string specifying which dataset to use

- hrf:

  HRF object to use for convolution (e.g., fmrihrf::HRF_SPMG1)

- include_intercept:

  Logical, whether to include an intercept column

## Value

A matrix with the design matrix (time x conditions)

## Examples

``` r
# Create design matrix using canonical HRF
X <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", fmrihrf::HRF_SPMG1)

# Test with a different HRF
X_wrong <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", fmrihrf::HRF_SPMG2)
```
