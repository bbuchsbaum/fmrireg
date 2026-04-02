# Build a time sketch matrix S (m x T)

Build a time sketch matrix S (m x T)

## Usage

``` r
make_time_sketch(Tlen, ctrl)
```

## Arguments

- Tlen:

  Integer time length

- ctrl:

  List(method = "gaussian"\|"countsketch"\|"srht"\|"ihs", m, iters)

## Value

A dense or sparse sketch matrix for Gaussian/CountSketch methods; `NULL`
for `"srht"` and `"ihs"` because those methods are applied via
plans/operators.
