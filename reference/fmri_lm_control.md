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
  When omitted (`NULL`), the control object uses
  `estimation_spec("runwise_meta")`, the same default as
  [`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md),
  so supplying a control object only to change another section (for
  example `robust`) does not change the estimation scope. Note that this
  differs from calling
  [`estimation_spec()`](https://bbuchsbaum.github.io/fmrireg/reference/estimation_spec.md)
  with no arguments, which returns the joint scope. There are two
  exceptions: HAC, sandwich, or Satterthwaite inference requested
  through `variance` is implemented only for a joint fit, so an omitted
  `estimation` then resolves to `estimation_spec("joint")` with a
  message saying so; and with a registered `engine`, an omitted
  `estimation` resolves to `estimation_spec("joint")` (see
  [`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md)).

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
