# Run the fmrireg command line interface

Runs the package-owned command line interface for `fmrireg`. The current
public command surface focuses on bundled benchmark dataset inspection
through the `benchmark` subcommands.

## Usage

``` r
fmrireg_cli(args = commandArgs(trailingOnly = TRUE))
```

## Arguments

- args:

  Character vector of command line arguments.

## Value

Integer exit status. Returns `0` on success, `1` for command-level
validation failures, and `2` for usage or runtime errors.

## Examples

``` r
fmrireg_cli(c("benchmark", "list"))
#> Dataset  Description
#> BM_Canonical_HighSNR Canonical HRF (SPMG1), high SNR, 3 conditions, fixed amplitudes per condition
#> BM_Canonical_LowSNR  Canonical HRF (SPMG1), low SNR, 3 conditions, fixed amplitudes per condition
#> BM_HRF_Variability_AcrossVoxels  HRF varies across voxel groups, 2 conditions, moderate SNR
#> BM_Trial_Amplitude_Variability   Single condition with significant trial-to-trial amplitude variability
#> BM_Complex_Realistic Complex realistic scenario: 3 HRF groups, 3 conditions, variable durations/amplitudes, AR(2) noise
#> [1] 0
fmrireg_cli(c(
  "benchmark", "summary",
  "--dataset", "BM_Canonical_HighSNR",
  "--json"
))
#> {
#>   "description": "Canonical HRF (SPMG1), high SNR, 3 conditions, fixed amplitudes per condition",
#>   "dimensions": {
#>     "n_timepoints": 150,
#>     "n_voxels": 100,
#>     "n_events": 45,
#>     "n_conditions": 3
#>   },
#>   "experimental_design": {
#>     "conditions": ["Cond1", "Cond2", "Cond3"],
#>     "events_per_condition": {
#>       "Cond1": 15,
#>       "Cond2": 15,
#>       "Cond3": 15
#>     },
#>     "TR": 2,
#>     "total_time": 300,
#>     "target_snr": 4
#>   },
#>   "hrf_information": {
#>     "type": "SPMG1",
#>     "hrf_object_name": "HRF_SPMG1",
#>     "hrf_object": ["structure(function (t) ", "{", "    do.call(orig_f, c(list(t = t), callable_params_list))", "}, class = \"function\", name = \"SPMG1\", nbasis = 1L, span = 24, param_names = c(\"P1\", ", "\"P2\", \"A1\"), params = list(P1 = 5, P2 = 15, A1 = 0.0833))"]
#>   },
#>   "noise_information": {
#>     "noise_type": "ar1",
#>     "noise_ar": 0.4,
#>     "noise_sd": 0.5217
#>   }
#> }
#> [1] 0
```
