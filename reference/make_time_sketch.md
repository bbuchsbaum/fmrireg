# Build a time sketch matrix S (m x T)

All sketches are normalised so that \\E\\S r\\^2 = \\r\\^2\\.

## Usage

``` r
make_time_sketch(Tlen, ctrl)
```

## Arguments

- Tlen:

  Integer time length

- ctrl:

  List(method = "gaussian"\|"countsketch"\|"srht", m)

## Value

A dense or sparse sketch matrix for Gaussian/CountSketch methods; `NULL`
for `"srht"`, which is applied through a plan
([`srht_apply()`](https://bbuchsbaum.github.io/fmrireg/reference/srht_apply.md)).
