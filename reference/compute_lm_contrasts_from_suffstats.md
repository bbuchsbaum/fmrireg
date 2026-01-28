# Compute contrast statistics from sufficient statistics (public API)

Convenience wrapper that accepts design/data sufficient statistics,
computes betas and residual variance, and delegates to
[`compute_lm_contrasts()`](https://bbuchsbaum.github.io/fmrireg/reference/compute_lm_contrasts.md).

## Usage

``` r
compute_lm_contrasts_from_suffstats(
  XtX,
  XtS,
  StS,
  df,
  sigma = NULL,
  sigma2 = NULL,
  contrasts = NULL,
  t_contrasts = NULL,
  f_contrasts = NULL,
  columns = NULL,
  output = c("stacked", "list"),
  robust_weights = NULL,
  ar_order = 0,
  drop_failed = TRUE
)
```

## Arguments

- XtX:

  Numeric (p x p) cross-product of the design.

- XtS:

  Numeric (p x V) cross-product of design with data.

- StS:

  Numeric length-V vector of sum of squares per voxel.

- df:

  Residual degrees of freedom.

- sigma:

  Optional numeric vector (length V) of residual std. dev.; ignored if
  `sigma2` provided.

- sigma2:

  Optional numeric vector (length V) of residual variances.

- contrasts:

  Optional named list mixing t- and F-contrasts; vectors are t, matrices
  are F.

- t_contrasts:

  Optional named list of numeric vectors (t-contrasts).

- f_contrasts:

  Optional named list of numeric matrices (F-contrasts).

- columns:

  Optional character vector (length p) naming coefficients; used for
  name matching.

- output:

  Either "stacked" (default; tibble) or "list" (raw list of tibbles).

- robust_weights:

  Optional numeric vector of robust weights or NULL. Used only for df
  adjustments.

- ar_order:

  Integer AR order; used only for effective df adjustments.

- drop_failed:

  Logical; drop contrasts that fail validation (default TRUE).

## Value

A tibble with rows per contrast (default) or a named list if
`output = "list"`.
