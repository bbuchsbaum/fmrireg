# Extract Statistical Measures from an fmri_lm Object

This function extracts statistical measures (e.g., estimates, standard
errors) from an `fmri_lm` object.

## Usage

``` r
pull_stat_revised(x, type, element)
```

## Arguments

- x:

  An `fmri_lm` object.

- type:

  The type of statistic to extract: `"betas"`, `"contrasts"`, or `"F"`.

- element:

  The specific element to extract, such as `"estimate"`, `"se"`,
  `"stat"`, or `"prob"`.

## Value

A tibble containing the requested statistical measures.
