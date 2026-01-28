# Extract Image/Volume for Coefficient

Extract Image/Volume for Coefficient

## Usage

``` r
coef_image(object, coef = 1, statistic = c("estimate", "se", "z", "p"))
```

## Arguments

- object:

  An fmri_meta object

- coef:

  Coefficient name or index

- statistic:

  Type of statistic to extract ("estimate", "se", "z", "p")

## Value

NeuroVol object or matrix

## Examples

``` r
toy_meta <- structure(
  list(
    coefficients = matrix(c(0.3, 0.1), nrow = 1,
      dimnames = list(NULL, c("A", "B"))),
    se = matrix(c(0.05, 0.06), nrow = 1)
  ),
  class = "fmri_meta"
)
coef_image(toy_meta, coef = "A")
#>   A 
#> 0.3 
```
