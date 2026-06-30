# Reducer: tidy beta-estimate table

Extracts the fitted condition estimates as a long data frame (one row
per term x voxel), with standard error and a Wald statistic when
available. Built on [`coef()`](https://rdrr.io/r/stats/coef.html) /
[`standard_error()`](https://bbuchsbaum.github.io/fmrireg/reference/standard_error.md)
for robustness.

## Usage

``` r
reduce_betas(add_id = TRUE, include_baseline = FALSE)
```

## Arguments

- add_id:

  Logical; prepend a `job_id` column. Default `TRUE`.

- include_baseline:

  Logical; include baseline/nuisance terms. Default `FALSE`.

## Value

A reducer function suitable for
[`fmri_template()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_template.md).

## See also

Other reducers:
[`reduce_contrasts()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_contrasts.md),
[`reduce_identity()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_identity.md),
[`reduce_write_results()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_write_results.md)
