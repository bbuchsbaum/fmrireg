# Fit Linear Model Contrasts

This function computes contrasts and beta statistics for a fitted linear
model.

## Usage

``` r
fit_lm_contrasts(fit, conlist, fcon, vnames, se = TRUE)
```

## Arguments

- fit:

  A fitted linear model object.

- conlist:

  A list of contrast matrices.

- fcon:

  A list of F-contrasts.

- vnames:

  Variable names corresponding to the model coefficients.

- se:

  Logical. Whether to compute standard errors. Default is `TRUE`.

## Value

A list containing contrasts, beta statistics, and the fitted model.
