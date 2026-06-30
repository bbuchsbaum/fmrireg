# Reducer: tidy contrast table

Extracts the fitted contrasts as a long data frame (one row per contrast
x voxel), optionally prefixed with the job id so results stack cleanly
across subjects. Surfaces template-level contrasts (those passed to
[`fmri_template()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_template.md),
applied by
[`run_job()`](https://bbuchsbaum.github.io/fmrireg/reference/run_job.md))
when present, otherwise the model's formula-embedded contrasts via
[`coef()`](https://rdrr.io/r/stats/coef.html). Returns a zero-row frame
when no contrasts are defined.

## Usage

``` r
reduce_contrasts(add_id = TRUE)
```

## Arguments

- add_id:

  Logical; prepend a `job_id` column. Default `TRUE`.

## Value

A reducer function suitable for
[`fmri_template()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_template.md).

## See also

Other reducers:
[`reduce_betas()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_betas.md),
[`reduce_identity()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_identity.md),
[`reduce_write_results()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_write_results.md)
