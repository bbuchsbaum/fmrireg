# Extract per-job errors from a batch result

Extract per-job errors from a batch result

## Usage

``` r
batch_errors(x)
```

## Arguments

- x:

  An `fmri_batch_result`.

## Value

A named character vector (by job id) of error messages for jobs that
failed; empty if all succeeded.

## See also

[`run_jobs()`](https://bbuchsbaum.github.io/fmrireg/reference/run_jobs.md),
[`batch_values()`](https://bbuchsbaum.github.io/fmrireg/reference/batch_values.md)
