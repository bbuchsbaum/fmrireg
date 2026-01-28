# Compute linear model contrast statistics (public API)

A polished wrapper around the internal fast contrast engine. Accepts
name-based contrast specifications and returns a consistent tibble by
default.

## Usage

``` r
compute_lm_contrasts(
  B,
  XtXinv,
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

- B:

  Numeric matrix (p x V) of coefficients.

- XtXinv:

  Numeric matrix (p x p), inverse cross-product of the design.

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

## Details

Inputs mirror a standard GLM contrast computation: provide `B` (p x V),
`XtXinv` (p x p), residual degrees of freedom `df`, and either `sigma`
or `sigma2` (per-voxel noise). Contrasts can be specified as a single
mixed `contrasts` list or separately as `t_contrasts` (vectors) and
`f_contrasts` (matrices). When contrast vectors/matrices are named (for
vectors) or have column names (for matrices), names are matched to
`columns` to determine the appropriate design indices; otherwise lengths
must match the number of parameters `p`.
