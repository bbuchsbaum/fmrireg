# Export jobs for external / array execution

Writes a serialized job manifest plus a minimal, backend-agnostic
`run_one.R` runner to `dir`. Deliberately emits no scheduler code: drive
it with whatever array system you have, passing the 1-based job index as
the first argument, e.g. `Rscript run_one.R $SLURM_ARRAY_TASK_ID`.

## Usage

``` r
export_jobs(
  jobs,
  dir,
  overwrite = FALSE,
  setup = "library(fmrireg)",
  results_dir = "results"
)
```

## Arguments

- jobs:

  A single
  [fmri_job](https://bbuchsbaum.github.io/fmrireg/reference/fmri_job.md)
  or a list of them (e.g. from
  [`instantiate()`](https://bbuchsbaum.github.io/fmrireg/reference/instantiate.md)).

- dir:

  Output directory (created if needed).

- overwrite:

  Overwrite an existing manifest. Default `FALSE`.

- setup:

  Character vector of R lines run at the top of `run_one.R` before jobs
  are executed (load packages, set threads, etc.). Default loads
  fmrireg.

- results_dir:

  Default output subdirectory written by `run_one.R`.

## Value

Invisibly, a list describing what was written (`dir`, `manifest`,
`runner`, `n`, `ids`).

## See also

[`read_jobs()`](https://bbuchsbaum.github.io/fmrireg/reference/read_jobs.md),
[`run_jobs()`](https://bbuchsbaum.github.io/fmrireg/reference/run_jobs.md),
[`instantiate()`](https://bbuchsbaum.github.io/fmrireg/reference/instantiate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
jobs <- instantiate(tmpl, manifest)
export_jobs(jobs, "study/jobs")
# then, per array task:  Rscript study/jobs/run_one.R $SLURM_ARRAY_TASK_ID
} # }
```
