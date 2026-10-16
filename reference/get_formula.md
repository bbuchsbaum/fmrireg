# Retrieve the formula underlying a model object

Generic used to expose the original modelling formula for fitted
objects.

Functions for creating and manipulating fmri_model objects

## Usage

``` r
get_formula(x, ...)

# S3 method for class 'fmri_model'
get_formula(x, ...)

# S3 method for class 'fmri_model'
get_formula(x, ...)
```

## Arguments

- x:

  The object to extract a formula from.

- ...:

  Additional arguments passed to methods.

## Value

A formula.

## Examples

``` r
fm <- fmrireg:::.demo_fmri_model()
get_formula(fm)
#> .y ~ condition + drift + block - 1
#> <environment: 0x563e0b4634d0>
```
