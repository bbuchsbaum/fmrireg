# Apply SRHT to a T x k matrix using a plan

Apply SRHT to a T x k matrix using a plan

## Usage

``` r
srht_apply(M, plan)
```

## Arguments

- M:

  Numeric matrix with `plan$T` rows.

- plan:

  A plan from
  [`make_srht_plan()`](https://bbuchsbaum.github.io/fmrireg/reference/make_srht_plan.md).

## Value

The `plan$m` x `ncol(M)` sketched matrix.
