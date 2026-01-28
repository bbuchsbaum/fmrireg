# Extract Standard Errors from Meta-Analysis

Extract Standard Errors from Meta-Analysis

## Usage

``` r
se(object)
```

## Arguments

- object:

  An fmri_meta object

## Value

Matrix of standard errors

## Examples

``` r
toy_meta <- structure(
  list(se = matrix(c(0.05, 0.06), nrow = 1)),
  class = "fmri_meta"
)
se(toy_meta)
#>      [,1] [,2]
#> [1,] 0.05 0.06
```
