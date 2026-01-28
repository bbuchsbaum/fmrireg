# Get Benchmark Dataset Summary

Provides a detailed summary of a specific benchmark dataset including
dimensions, experimental design, and ground truth information.

## Usage

``` r
get_benchmark_summary(dataset_name)
```

## Arguments

- dataset_name:

  Character string specifying which dataset to summarize

## Value

A list with summary information about the dataset

## Examples

``` r
# Get summary of a specific dataset
summary_info <- get_benchmark_summary("BM_Canonical_HighSNR")
print(summary_info)
#> $description
#> [1] "Canonical HRF (SPMG1), high SNR, 3 conditions, fixed amplitudes per condition"
#> 
#> $dimensions
#> $dimensions$n_timepoints
#> [1] 150
#> 
#> $dimensions$n_voxels
#> [1] 100
#> 
#> $dimensions$n_events
#> [1] 45
#> 
#> $dimensions$n_conditions
#> [1] 3
#> 
#> 
#> $experimental_design
#> $experimental_design$conditions
#> [1] "Cond1" "Cond2" "Cond3"
#> 
#> $experimental_design$events_per_condition
#> $experimental_design$events_per_condition$Cond1
#> [1] 15
#> 
#> $experimental_design$events_per_condition$Cond2
#> [1] 15
#> 
#> $experimental_design$events_per_condition$Cond3
#> [1] 15
#> 
#> 
#> $experimental_design$TR
#> [1] 2
#> 
#> $experimental_design$total_time
#> [1] 300
#> 
#> $experimental_design$target_snr
#> [1] 4
#> 
#> 
#> $hrf_information
#> $hrf_information$type
#> [1] "SPMG1"
#> 
#> $hrf_information$hrf_object_name
#> [1] "HRF_SPMG1"
#> 
#> $hrf_information$hrf_object
#> function (t) 
#> {
#>     do.call(orig_f, c(list(t = t), callable_params_list))
#> }
#> <bytecode: 0x555a52022688>
#> <environment: 0x555a5202ad18>
#> attr(,"class")
#> [1] "HRF"      "function"
#> attr(,"name")
#> [1] "SPMG1"
#> attr(,"nbasis")
#> [1] 1
#> attr(,"span")
#> [1] 24
#> attr(,"param_names")
#> [1] "P1" "P2" "A1"
#> attr(,"params")
#> attr(,"params")$P1
#> [1] 5
#> 
#> attr(,"params")$P2
#> [1] 15
#> 
#> attr(,"params")$A1
#> [1] 0.0833
#> 
#> 
#> 
#> $noise_information
#> $noise_information$noise_type
#> [1] "ar1"
#> 
#> $noise_information$noise_ar
#> [1] 0.4
#> 
#> $noise_information$noise_sd
#> [1] 0.5217219
#> 
#> 
```
