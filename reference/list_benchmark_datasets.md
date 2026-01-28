# List Available Benchmark Datasets

Returns a summary of all available benchmark datasets with their
descriptions.

## Usage

``` r
list_benchmark_datasets()
```

## Value

A data.frame with dataset names and descriptions

## Examples

``` r
# See what benchmark datasets are available
list_benchmark_datasets()
#>                                                         Dataset
#> BM_Canonical_HighSNR                       BM_Canonical_HighSNR
#> BM_Canonical_LowSNR                         BM_Canonical_LowSNR
#> BM_HRF_Variability_AcrossVoxels BM_HRF_Variability_AcrossVoxels
#> BM_Trial_Amplitude_Variability   BM_Trial_Amplitude_Variability
#> BM_Complex_Realistic                       BM_Complex_Realistic
#>                                                                                                                        Description
#> BM_Canonical_HighSNR                                 Canonical HRF (SPMG1), high SNR, 3 conditions, fixed amplitudes per condition
#> BM_Canonical_LowSNR                                   Canonical HRF (SPMG1), low SNR, 3 conditions, fixed amplitudes per condition
#> BM_HRF_Variability_AcrossVoxels                                         HRF varies across voxel groups, 2 conditions, moderate SNR
#> BM_Trial_Amplitude_Variability                              Single condition with significant trial-to-trial amplitude variability
#> BM_Complex_Realistic            Complex realistic scenario: 3 HRF groups, 3 conditions, variable durations/amplitudes, AR(2) noise
```
