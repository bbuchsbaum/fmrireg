# Deprecated OLS and Welch t-test kernels

`ols_t_cpp()` and `welch_t_cpp()` are deprecated and will be removed in
the next minor release.
[`fmri_ttest()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_ttest.md)
(`engine = "classic"` or `"welch"`) and
[`fmri_ols_fit()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_ols_fit.md)
now use the fmrigds reducers `"ols:voxelwise"` and `"ttest:welch"`,
which also drop non-finite subjects per feature.

## Usage

``` r
ols_t_cpp(Y, X)

welch_t_cpp(Y, g_in)
```

## Arguments

- Y:

  S x P matrix (subjects x features).

- X:

  S x K design matrix.

- g_in:

  Length-S vector of group indicators (1/2 or 0/1).

## Value

`ols_t_cpp()`: list with beta, se, t (K x P), df, n_obs and ok.
`welch_t_cpp()`: list with muA, muB, t, df (Welch), nA and nB.
