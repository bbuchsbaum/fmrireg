# Register a plugin engine for `fmri_lm`

Register a plugin engine for `fmri_lm`

## Usage

``` r
register_engine(name, fit, preflight = NULL, capabilities = list())
```

## Arguments

- name:

  Character scalar identifier advertised to users (e.g. "friman").

- fit:

  Function invoked as `fit(model, dataset, args, cfg)` and expected to
  return an `fmri_lm` object.

- preflight:

  Optional function invoked before fitting; receives the same arguments
  as `fit` and can signal errors early.

- capabilities:

  Optional named list describing engine support for global
  [`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md)
  options. Recognized fields currently include `robust`,
  `preprocessing`, `ar_voxelwise`, `ar_by_cluster`, plus contextual
  rules such as `requires_event_regressors`,
  `requires_parcels_for_by_cluster`, and
  `forbid_by_cluster_dataset_classes`.

## Value

Invisibly, `TRUE`.
