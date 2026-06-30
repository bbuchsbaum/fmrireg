# Register an execution backend for run_jobs()

Adds a named backend that
[`run_jobs()`](https://bbuchsbaum.github.io/fmrireg/reference/run_jobs.md)
can dispatch to. This is how to plug in a custom scheduler driver
without modifying fmrireg. The builtins `"sequential"` and `"future"`
are always available.

## Usage

``` r
register_run_backend(name, fn, overwrite = FALSE)
```

## Arguments

- name:

  Backend name (a string).

- fn:

  A function `function(jobs, run_one, ...)` returning an
  order-preserving list of per-job results.

- overwrite:

  Allow replacing an existing backend of the same name.

## Value

The backend name, invisibly.

## See also

[`run_jobs()`](https://bbuchsbaum.github.io/fmrireg/reference/run_jobs.md),
[`run_backends()`](https://bbuchsbaum.github.io/fmrireg/reference/run_backends.md)

## Examples

``` r
register_run_backend("first_only", function(jobs, run_one, ...) {
  list(run_one(jobs[[1]]))
})
```
