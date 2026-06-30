# Run a batch of jobs

Executes a list of
[fmri_job](https://bbuchsbaum.github.io/fmrireg/reference/fmri_job.md)s,
returning an `fmri_batch_result`. By default each job is isolated: a
failure is captured rather than aborting the batch. Set
`parallel = TRUE` to dispatch through the active
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
(e.g. `multisession`, `cluster`, or a `future.batchtools` cluster
backend).

## Usage

``` r
run_jobs(
  jobs,
  parallel = FALSE,
  progress = FALSE,
  on_error = c("isolate", "stop"),
  backend = NULL
)
```

## Arguments

- jobs:

  A single
  [fmri_job](https://bbuchsbaum.github.io/fmrireg/reference/fmri_job.md)
  or a list of them.

- parallel:

  Logical shorthand; `TRUE` selects the `"future"` backend, `FALSE` the
  `"sequential"` one. Ignored if `backend` is given.

- progress:

  Logical; per-job fitting progress.

- on_error:

  Either `"isolate"` (default; capture per-job errors) or `"stop"` (fail
  the batch on the first error).

- backend:

  Optional execution backend name (see
  [`run_backends()`](https://bbuchsbaum.github.io/fmrireg/reference/run_backends.md)).
  Overrides `parallel`. The `"future"` backend honours the active
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  (including `future.batchtools` cluster plans).

## Value

An object of class `fmri_batch_result`.

## See also

[`run_job()`](https://bbuchsbaum.github.io/fmrireg/reference/run_job.md),
[`batch_values()`](https://bbuchsbaum.github.io/fmrireg/reference/batch_values.md),
[`batch_errors()`](https://bbuchsbaum.github.io/fmrireg/reference/batch_errors.md),
[`register_run_backend()`](https://bbuchsbaum.github.io/fmrireg/reference/register_run_backend.md)

## Examples

``` r
if (FALSE) { # \dontrun{
jobs <- instantiate(tmpl, manifest)
future::plan(future::multisession, workers = 4)
res <- run_jobs(jobs, parallel = TRUE)
values <- batch_values(res)
} # }
```
