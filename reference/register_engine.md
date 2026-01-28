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

  Optional named list describing the engine (for future use).

## Value

Invisibly, `TRUE`.
