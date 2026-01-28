# Flip Sign of Coefficients

Reverses the sign of coefficient-like outputs in a fit object. Useful
for switching between A-B and B-A conventions.

## Usage

``` r
flip_sign(fit, coef = NULL)
```

## Arguments

- fit:

  An fmri_ttest_fit or similar object

- coef:

  Character vector of coefficient names to flip (default: all)

## Value

Modified fit object with flipped signs
