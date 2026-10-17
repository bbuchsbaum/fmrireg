# Print fitted fMRI model objects

Summarises a fitted model: estimation scope, design size, inference
degrees of freedom, contrasts and noise/robust configuration.

## Usage

``` r
# S3 method for class 'fmri_lm'
print(x, ...)

# S3 method for class 'fmri_model'
print(x, ...)

# S3 method for class 'fmri_rlm'
print(x, ...)
```

## Arguments

- x:

  An fmri_rlm object

- ...:

  Additional arguments passed to print.fmri_lm

## Value

Invisibly returns the input object x
