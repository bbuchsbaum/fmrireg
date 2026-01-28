# Design Matrix Method for fmri_lm Objects

Extract the design matrix from an fmri_lm object by delegating to its
model component.

## Usage

``` r
# S3 method for class 'fmri_lm'
design_matrix(x, ...)
```

## Arguments

- x:

  An fmri_lm object

- ...:

  Additional arguments passed to the design_matrix method for the model

## Value

The design matrix from the fmri_lm object's model
