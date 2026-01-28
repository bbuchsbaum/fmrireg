# Write Meta-Analysis Results

Write Meta-Analysis Results

## Usage

``` r
# S3 method for class 'fmri_meta'
write_results(x, path = ".", prefix = "meta", format = c("nifti", "h5"), ...)
```

## Arguments

- x:

  An fmri_meta object

- path:

  Output directory

- prefix:

  File name prefix

- format:

  Output format ("nifti" or "h5")

- ...:

  Additional arguments passed to write_results.fmri_lm

## Value

Invisible list of created files
