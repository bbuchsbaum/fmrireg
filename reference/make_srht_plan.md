# Build SRHT plan

The plan's `scale` is `1 / sqrt(m)`, which makes the sketch normalised
(\\E\\S r\\^2 = \\r\\^2\\) because the Walsh-Hadamard transform is
applied unnormalised.

## Usage

``` r
make_srht_plan(Tlen, m)
```

## Arguments

- Tlen:

  Integer series length.

- m:

  Integer number of sketch rows.

## Value

A list with the sign, permutation, row-selection and scale components
used by
[`srht_apply()`](https://bbuchsbaum.github.io/fmrireg/reference/srht_apply.md).
