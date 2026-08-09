# Nuisance-projection specification

Nuisance-projection specification

## Usage

``` r
projection_spec(
  method = c("none", "soft_subspace"),
  nuisance_mask = NULL,
  nuisance_matrix = NULL,
  lambda = "auto",
  warn_redundant = TRUE
)
```

## Arguments

- method:

  Projection rule.

- nuisance_mask:

  Optional logical or numeric selector for nuisance features.

- nuisance_matrix:

  Optional explicit nuisance-feature matrix.

- lambda:

  Nonnegative regularization strength, `"auto"`, or `"gcv"`.

- warn_redundant:

  Warn when projection overlaps modeled nuisance terms.
