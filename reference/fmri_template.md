# Define a subject-invariant model template

Captures a complete model specification that is constant across
subjects. Per-subject data (BOLD paths, event table, run lengths,
confounds) is bound later with
[`instantiate()`](https://bbuchsbaum.github.io/fmrireg/reference/instantiate.md)
to produce serializable
[fmri_job](https://bbuchsbaum.github.io/fmrireg/reference/fmri_job.md)
recipes that can be executed locally, in parallel via future, or on a
cluster.

## Usage

``` r
fmri_template(
  formula,
  block,
  baseline = baseline_spec(),
  durations = 0,
  contrasts = NULL,
  control = fmri_lm_control(),
  strategy = c("runwise", "chunkwise"),
  engine = NULL,
  engine_args = list(),
  reducer = NULL
)
```

## Arguments

- formula:

  Event model formula, e.g. `onset ~ hrf(condition)`.

- block:

  Block / run-structure formula, e.g. `~ run`.

- baseline:

  A
  [`baseline_spec()`](https://bbuchsbaum.github.io/fmrireg/reference/baseline_spec.md)
  describing the nuisance model.

- durations:

  Event durations passed through to the event model.

- contrasts:

  Optional contrast specification (e.g. from
  [`contrast_set()`](https://bbuchsbaum.github.io/fmrireg/reference/contrast_set.md)).

- control:

  An `fmri_lm_config` from
  [`fmri_lm_control()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm_control.md)
  holding the robust / AR / preprocessing options for fitting.

- strategy:

  Fitting strategy: `"runwise"` or `"chunkwise"`.

- engine:

  Optional fitting engine name (see
  [`register_engine()`](https://bbuchsbaum.github.io/fmrireg/reference/register_engine.md)).

- engine_args:

  Optional list of engine arguments.

- reducer:

  Optional function `function(fit, job)` that turns a fitted `fmri_lm`
  into the per-subject output (written to disk and/or returned). `NULL`
  (the default) returns the fitted object unchanged. See builtins such
  as
  [`reduce_write_results()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_write_results.md).

## Value

An object of class `fmri_template`.

## Details

The template is pure, serializable data: it holds no voxel data and (by
design) no captured execution environment. The `reducer` is the one
field that can accidentally capture state; see **Reducer
serializability** below.

## Reducer serializability

A `reducer` should be a top-level / package function (or a closure that
captures nothing), so it survives serialization to a worker node. A
closure that captures local variables will drag those bindings along
when the job is written with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html); `fmri_template()`
warns in that case.

## See also

[`baseline_spec()`](https://bbuchsbaum.github.io/fmrireg/reference/baseline_spec.md),
[`instantiate()`](https://bbuchsbaum.github.io/fmrireg/reference/instantiate.md),
[`fmri_lm_control()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm_control.md)

## Examples

``` r
tmpl <- fmri_template(onset ~ hrf(condition), ~ run,
                      baseline = baseline_spec(degree = 3))
tmpl
#> <fmri_template>
#>   formula:   onset ~ hrf(condition) 
#>   block:     ~run 
#>   strategy: runwise
#>   baseline:  bs(degree=3) 
#>   contrasts: none 
#>   reducer:   none (returns fitted object) 
```
