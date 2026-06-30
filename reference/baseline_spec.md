# Baseline specification (subject-invariant)

Describes how the nuisance / baseline model should be built for every
subject, without binding any particular subject's confound values. The
concrete `baseline_model` is assembled per subject at instantiation time
from this spec, the subject's `sampling_frame`, and the subject's
confound matrix.

## Usage

``` r
baseline_spec(
  degree = 3,
  basis = c("bs", "poly", "ns", "constant"),
  confounds = NULL,
  intercept = c("runwise", "global", "none"),
  nuisance_check = c("warn", "error", "drop", "none")
)
```

## Arguments

- degree:

  Integer drift degree passed to
  [`baseline_model()`](https://bbuchsbaum.github.io/fmridesign/reference/baseline_model.html).

- basis:

  Drift basis: one of `"bs"`, `"poly"`, `"ns"`, `"constant"`.

- confounds:

  Optional confound selection used to populate the per-subject nuisance
  regressors. Either `NULL` (no confounds), a character vector of
  confound column names / patterns, or a `bidser` confound-set object.
  Resolved to actual values by
  [`from_bids()`](https://bbuchsbaum.github.io/fmrireg/reference/from_bids.md)
  /
  [`instantiate()`](https://bbuchsbaum.github.io/fmrireg/reference/instantiate.md).

- intercept:

  Intercept handling passed to
  [`baseline_model()`](https://bbuchsbaum.github.io/fmridesign/reference/baseline_model.html).

- nuisance_check:

  How to handle problematic nuisance regressors, passed to
  [`baseline_model()`](https://bbuchsbaum.github.io/fmridesign/reference/baseline_model.html).

## Value

An object of class `baseline_spec`.

## See also

[`fmri_template()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_template.md),
[`baseline_model()`](https://bbuchsbaum.github.io/fmridesign/reference/baseline_model.html)

## Examples

``` r
baseline_spec(degree = 3, basis = "bs")
#> <baseline_spec>
#>   basis: bs  degree: 3  intercept: runwise
#>   confounds: none
```
