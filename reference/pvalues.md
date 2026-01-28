# Compute P-values from Meta-Analysis

Compute P-values from Meta-Analysis

## Usage

``` r
pvalues(object, two_tailed = TRUE)
```

## Arguments

- object:

  An fmri_meta object

- two_tailed:

  Logical. Use two-tailed test (default: TRUE)

## Value

Matrix of p-values

## Examples

``` r
toy_meta <- structure(
  list(
    coefficients = matrix(c(0.2, -0.1), nrow = 1),
    se = matrix(c(0.05, 0.08), nrow = 1)
  ),
  class = "fmri_meta"
)
pvalues(toy_meta)
#>              [,1]      [,2]
#> [1,] 6.334248e-05 0.2112995
```
