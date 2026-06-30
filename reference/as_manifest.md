# Coerce bindings to an fmri_manifest

Wraps a list of bindings (or a manifest `data.frame`) into an
`fmri_manifest` for
[`instantiate()`](https://bbuchsbaum.github.io/fmrireg/reference/instantiate.md).
This is the generic, BIDS-free way to describe subjects; see
[`from_bids()`](https://bbuchsbaum.github.io/fmrireg/reference/from_bids.md)
for the BIDS populator.

## Usage

``` r
as_manifest(x)
```

## Arguments

- x:

  A list of bindings (named lists, each with at least `id`, `scans`,
  `TR`, `run_length`) or a manifest `data.frame`.

## Value

An object of class `fmri_manifest`.

## See also

[`from_bids()`](https://bbuchsbaum.github.io/fmrireg/reference/from_bids.md),
[`instantiate()`](https://bbuchsbaum.github.io/fmrireg/reference/instantiate.md)

## Examples

``` r
m <- as_manifest(list(
  list(id = "sub-01", scans = matrix(rnorm(80 * 2), 80, 2), TR = 2,
       run_length = c(40, 40),
       events = data.frame(onset = c(5, 45),
                           condition = factor(c("A", "B")), run = c(1, 2)))))
length(m)
#> [1] 1
```
