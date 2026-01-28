# Compute Z-scores from Meta-Analysis

Compute Z-scores from Meta-Analysis

## Usage

``` r
zscores(object)
```

## Arguments

- object:

  An fmri_meta object

## Value

Matrix of z-scores

## Examples

``` r
toy_meta <- structure(
  list(
    coefficients = matrix(c(0.2, -0.1), nrow = 1),
    se = matrix(c(0.05, 0.08), nrow = 1)
  ),
  class = "fmri_meta"
)
zscores(toy_meta)
#>      [,1]  [,2]
#> [1,]    4 -1.25
```
