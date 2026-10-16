# Fit Contrasts for Linear Model (Default Method)

This function calculates contrasts for a fitted linear model.

## Usage

``` r
# Default S3 method
fit_contrasts(object, conmat, colind, se = TRUE, cov.unscaled = NULL, ...)
```

## Arguments

- object:

  The fitted linear model object.

- conmat:

  The contrast matrix or contrast vector.

- colind:

  The subset column indices in the design associated with the contrast.

- se:

  Whether to compute standard errors, t-statistics, and p-values
  (default: TRUE).

- cov.unscaled:

  Optional unscaled coefficient covariance \\(X'X)^{-1}\\ of `object`,
  to reuse across several contrasts on the same fit (default: computed
  from `object`).

- ...:

  Additional arguments (unused)

## Value

A list containing the following elements:

- `conmat`: Contrast matrix.

- `sigma`: Residual standard error.

- `df.residual`: Degrees of freedom for residuals.

- `estimate`: Estimated contrasts.

- `se`: Standard errors of the contrasts (if `se = TRUE`).

- `stat`: t-statistics for the contrasts (if `se = TRUE`).

- `prob`: Probabilities associated with the t-statistics (if
  `se = TRUE`).

- `stat_type`: Type of the statistics calculated.
