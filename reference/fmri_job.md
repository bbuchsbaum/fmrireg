# Construct a per-subject job recipe

Binds an
[fmri_template](https://bbuchsbaum.github.io/fmrireg/reference/fmri_template.md)
to one subject's
[dataset_spec](https://bbuchsbaum.github.io/fmrireg/reference/dataset_spec.md).
The result is a fully serializable recipe (no voxel data, no captured
environment) that `run()` turns into a fitted model and reduced output.
Build jobs with
[`instantiate()`](https://bbuchsbaum.github.io/fmrireg/reference/instantiate.md)
rather than by hand in the common case.

## Usage

``` r
fmri_job(id, template, dataset_spec, meta = list(), nuisance = NULL)
```

## Arguments

- id:

  A unique job identifier (e.g. a subject label).

- template:

  An
  [fmri_template](https://bbuchsbaum.github.io/fmrireg/reference/fmri_template.md).

- dataset_spec:

  A
  [dataset_spec](https://bbuchsbaum.github.io/fmrireg/reference/dataset_spec.md)
  describing this subject's data.

- meta:

  Optional named list of metadata used for output keying (e.g.
  `list(subject = "01", task = "stroop", space = "MNI152...")`).

- nuisance:

  Optional per-subject nuisance / confound regressors fed to the
  baseline model at run time: `NULL`, a numeric matrix with one row per
  scan (split across runs), or a list of per-run matrices. (A deferred
  file-read spec is also accepted by `run()`.)

## Value

An object of class `fmri_job`.

## See also

[`instantiate()`](https://bbuchsbaum.github.io/fmrireg/reference/instantiate.md),
`run()`,
[`export_jobs()`](https://bbuchsbaum.github.io/fmrireg/reference/export_jobs.md)

## Examples

``` r
tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
ds <- dataset_spec("fmri_dataset",
                   args = list(scans = "run-1_bold.nii.gz", TR = 2,
                               run_length = 200), source = "file")
fmri_job("sub-01", tmpl, ds, meta = list(subject = "01"))
#> <fmri_job> sub-01
#>   dataset: fmri_dataset() [file]
#>   meta:    subject=01 
#>   formula: onset ~ hrf(condition) 
```
