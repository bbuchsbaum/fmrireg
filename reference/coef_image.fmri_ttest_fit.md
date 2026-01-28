# Extract Coefficient Image from fmri_ttest_fit

Creates a NeuroVol image from coefficients in an fmri_ttest_fit object.

## Usage

``` r
# S3 method for class 'fmri_ttest_fit'
coef_image(object, coef = 1, statistic = c("estimate", "se", "z", "p"), ...)
```

## Arguments

- object:

  An fmri_ttest_fit object

- coef:

  Character or integer; coefficient to extract

- statistic:

  Character string; type of statistic to extract:

  - "estimate": Coefficient estimates (beta values)

  - "se": Standard errors (if available)

  - "z": Z-scores

  - "p": P-values

- ...:

  Additional arguments (e.g., mask to apply)

## Value

NeuroVol object or numeric vector
