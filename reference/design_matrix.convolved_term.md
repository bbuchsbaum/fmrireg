# Design Matrix for Convolved Terms

Extract the design matrix from a convolved term object, optionally
filtered by block ID.

## Usage

``` r
# S3 method for class 'convolved_term'
design_matrix(x, blockid = NULL, ...)
```

## Arguments

- x:

  A convolved_term object

- blockid:

  Optional numeric vector specifying which blocks/runs to include

- ...:

  Additional arguments (not used)

## Value

A matrix containing the convolved design matrix
