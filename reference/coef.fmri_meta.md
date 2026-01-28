# Extract Coefficients from Meta-Analysis

Extract Coefficients from Meta-Analysis

## Usage

``` r
# S3 method for class 'fmri_meta'
coef(object, ...)
```

## Arguments

- object:

  An fmri_meta object

- ...:

  Additional arguments

## Value

Matrix of coefficients

## Examples

``` r
toy_meta <- structure(
  list(coefficients = matrix(c(0.2, -0.1), nrow = 1)),
  class = "fmri_meta"
)
coef(toy_meta)
#>      [,1] [,2]
#> [1,]  0.2 -0.1
```
