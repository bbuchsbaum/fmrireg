# Reducer: return the fitted model unchanged

The escape hatch when you want the entire `fmri_lm` object per subject.
Note this is the largest possible output.

## Usage

``` r
reduce_identity()
```

## Value

A reducer function suitable for
[`fmri_template()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_template.md).

## See also

Other reducers:
[`reduce_betas()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_betas.md),
[`reduce_contrasts()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_contrasts.md),
[`reduce_write_results()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_write_results.md)
