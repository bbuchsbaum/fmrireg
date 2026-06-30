# Reducer: write per-subject results to disk (BIDS-keyed)

Calls
[`write_results()`](https://bbuchsbaum.github.io/fmrireg/reference/write_results.md)
on the fitted model, keying the output filenames by the job's `meta`
(`subject`, `task`, `space`). Returns the vector of created file paths.
This is the default reducer for large studies: workers write compact
statistical maps that
[`group_data()`](https://bbuchsbaum.github.io/fmrireg/reference/group_data.md)
/
[`fmri_meta()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_meta.md)
can pick up for the group level.

## Usage

``` r
reduce_write_results(
  format = c("h5", "nifti", "gds"),
  stats = c("beta", "tstat"),
  contrasts = NULL,
  path = NULL,
  desc = "GLM",
  overwrite = FALSE
)
```

## Arguments

- format:

  Output format passed to
  [`write_results()`](https://bbuchsbaum.github.io/fmrireg/reference/write_results.md).

- stats:

  Which contrast statistics to write (mapped to `contrast_stats`).

- contrasts:

  Optional restriction to specific contrasts.

- path:

  Output directory (passed to
  [`write_results()`](https://bbuchsbaum.github.io/fmrireg/reference/write_results.md)).

- desc:

  BIDS `desc-` label.

- overwrite:

  Overwrite existing files.

## Value

A reducer function suitable for
[`fmri_template()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_template.md).

## Note

[`write_results()`](https://bbuchsbaum.github.io/fmrireg/reference/write_results.md)
requires BIDS entities: each job's `meta` should carry `subject`
(defaults to the job id), `task`, and `space`. A missing `task`/`space`
surfaces as an error when the reducer runs.

## See also

Other reducers:
[`reduce_betas()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_betas.md),
[`reduce_contrasts()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_contrasts.md),
[`reduce_identity()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_identity.md)
