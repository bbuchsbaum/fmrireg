# Instantiate a template into per-subject jobs

Binds an
[fmri_template](https://bbuchsbaum.github.io/fmrireg/reference/fmri_template.md)
to one or more per-subject data bindings, producing serializable
[fmri_job](https://bbuchsbaum.github.io/fmrireg/reference/fmri_job.md)
recipes. No data is loaded: file-backed scans are kept as paths and
realized lazily by `run()`.

## Usage

``` r
instantiate(template, x, ...)
```

## Arguments

- template:

  An
  [fmri_template](https://bbuchsbaum.github.io/fmrireg/reference/fmri_template.md).

- x:

  A single binding (a named list with at least `id`), a list of
  bindings, or a manifest `data.frame` (one row per analysis unit, with
  list-columns for `scans`/`events`/`run_length`/etc.).

- ...:

  Unused.

## Value

A single
[fmri_job](https://bbuchsbaum.github.io/fmrireg/reference/fmri_job.md)
if `x` is one binding, otherwise a list of
[fmri_job](https://bbuchsbaum.github.io/fmrireg/reference/fmri_job.md)s.

## See also

[`fmri_job()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_job.md),
`run()`,
[`from_bids()`](https://bbuchsbaum.github.io/fmrireg/reference/from_bids.md)

## Examples

``` r
tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
b <- list(id = "sub-01", scans = matrix(rnorm(80 * 3), 80, 3),
          TR = 2, run_length = c(40, 40),
          events = data.frame(onset = c(5, 25, 45, 65),
                              condition = factor(c("A", "B", "A", "B")),
                              run = c(1, 1, 2, 2)))
job <- instantiate(tmpl, b)
```
