# Configure an fMRI linear model

`fmri_lm_control()` is the single statistical configuration boundary for
[`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md).
Each section is a validated typed specification; mechanical choices such
as chunking and parallelism belong in
[`compute_spec()`](https://bbuchsbaum.github.io/fmrireg/reference/compute_spec.md).

## Usage

``` r
fmri_lm_control(
  estimation = NULL,
  noise = NULL,
  robust = NULL,
  variance = NULL,
  weights = NULL,
  projection = NULL,
  na_action = c("error", "propagate"),
  ...
)
```

## Arguments

- estimation:

  An
  [`estimation_spec()`](https://bbuchsbaum.github.io/fmrireg/reference/estimation_spec.md).

- noise:

  A
  [`noise_spec()`](https://bbuchsbaum.github.io/fmrireg/reference/noise_spec.md).

- robust:

  A
  [`robust_spec()`](https://bbuchsbaum.github.io/fmrireg/reference/robust_spec.md).

- variance:

  A
  [`variance_spec()`](https://bbuchsbaum.github.io/fmrireg/reference/variance_spec.md).

- weights:

  A
  [`weights_spec()`](https://bbuchsbaum.github.io/fmrireg/reference/weights_spec.md).

- projection:

  A
  [`projection_spec()`](https://bbuchsbaum.github.io/fmrireg/reference/projection_spec.md).

- na_action:

  Either `"error"` or `"propagate"`.

- ...:

  Transitional legacy option lists retained for one compatibility
  window. New code should use the typed sections above.

## Value

An object of class `fmri_lm_control`.
