# Parse Contrast Formula

Accepts a one-sided formula whose right-hand side is a signed sum of
coefficient names, each optionally multiplied by a number, e.g.
`~ groupold - groupyoung` or `~ 0.5 * a + 0.5 * b - c`. Unknown
coefficient names and any other expression are errors.

## Usage

``` r
parse_contrast_formula(formula, object)
```

## Arguments

- formula:

  Contrast formula

- object:

  fmri_meta object

## Value

Numeric vector of contrast weights
