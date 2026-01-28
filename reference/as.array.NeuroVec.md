# Coerce NeuroVec to base array

This method allows NeuroVec objects from neuroim2 to be converted to
base R arrays using the standard as.array() function. This is
particularly useful in testing and data manipulation contexts.

## Usage

``` r
# S3 method for class 'NeuroVec'
as.array(x, ...)
```

## Arguments

- x:

  A NeuroVec object from neuroim2

- ...:

  Additional arguments (currently unused)

## Value

A base R array containing the data from the NeuroVec object
