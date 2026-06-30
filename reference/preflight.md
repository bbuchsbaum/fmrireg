# Preflight-check jobs before fan-out

Validates one or more
[fmri_job](https://bbuchsbaum.github.io/fmrireg/reference/fmri_job.md)s
against the structural contract their template implies, collecting
issues per job. Intended to run on the driver before
[`run_jobs()`](https://bbuchsbaum.github.io/fmrireg/reference/run_jobs.md)
/
[`export_jobs()`](https://bbuchsbaum.github.io/fmrireg/reference/export_jobs.md).

## Usage

``` r
preflight(x, check_files = FALSE, on_issue = c("warn", "error", "collect"))
```

## Arguments

- x:

  A single
  [fmri_job](https://bbuchsbaum.github.io/fmrireg/reference/fmri_job.md)
  or a list of them.

- check_files:

  Logical; for file-backed jobs, verify scan paths exist.

- on_issue:

  One of `"warn"` (default), `"error"` (stop if any issue), or
  `"collect"` (return silently).

## Value

Invisibly, an object of class `fmri_preflight` with an `$issues` data
frame (`job_id`, `message`) and `$ok`.

## Details

Checks performed per job: template validity; every variable referenced
by the design `formula` and `block` is a column of that job's event
table; `TR` is positive; run lengths are consistent with the data
(`matrix_dataset`: rows match `sum(run_length)`; `fmri_dataset`: one
`run_length` per scan file); nuisance regressor rows match the total
number of scans; and, when `check_files = TRUE`, that file-backed scans
exist on disk.

The design-column check uses
[`all.vars()`](https://rdrr.io/r/base/allnames.html) on the formula, so
it is deliberately conservative: a formula with a variable-valued HRF
argument (e.g. `hrf(x, basis = my_basis)`) may flag `my_basis` as a
missing column. Event tables supplied as file paths are not yet
validated here.

## See also

[`validate_template()`](https://bbuchsbaum.github.io/fmrireg/reference/validate_template.md),
[`run_jobs()`](https://bbuchsbaum.github.io/fmrireg/reference/run_jobs.md),
[`export_jobs()`](https://bbuchsbaum.github.io/fmrireg/reference/export_jobs.md)

## Examples

``` r
tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
job <- instantiate(tmpl, list(id = "sub-01",
                              scans = matrix(rnorm(80 * 2), 80, 2), TR = 2,
                              run_length = c(40, 40),
                              events = data.frame(onset = c(5, 45),
                                                  condition = factor(c("A", "B")),
                                                  run = c(1, 2))))
preflight(job)
```
