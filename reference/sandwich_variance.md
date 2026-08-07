# Sandwich Variance Estimation in fmrireg

This documentation describes the sandwich variance estimation
capabilities in fmrireg, which provide robust standard errors for
regression coefficients when model assumptions are violated.

## Background

The sandwich variance estimator (also known as the Huber-White
estimator) provides valid standard errors even when the residuals
exhibit heteroscedasticity or other violations of the classical linear
model assumptions.

## Mathematical Details

The sandwich estimator is computed as: \$\$V\_{sandwich} = (X'X)^{-1} X'
\Omega X (X'X)^{-1}\$\$

where \\\Omega\\ is a diagonal matrix with squared residuals on the
diagonal. For robust regression with weights \\w_i\\, the weighted
version is: \$\$V\_{sandwich} = (X'WX)^{-1} X'W \Omega WX
(X'WX)^{-1}\$\$

## Usage in fmrireg

Sandwich variance estimation is **not** currently applied by
[`fmri_lm`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md).
The helpers documented here are available for direct use, but no fitting
path calls them: the robust path reports a model-based
weighted-least-squares variance, \\\hat\sigma^2\_{robust} \\
\mathrm{diag}((X'WX)^{-1})\\, not a sandwich.

Applying a sandwich to the robust path correctly would require the
M-estimator influence-function form, whose weights and influence matrix
are voxel-dependent, rather than the ordinary-least-squares expression
above.

## Effective Degrees of Freedom

Robust down-weighting reduces the residual degrees of freedom
multiplicatively: \$\$df\_{effective} = (n - p) \times \frac{\sum_i
w_i}{n}.\$\$

AR order does **not** enter this formula. fmrireg whitens the data and
then fits by ordinary least squares, so when the filter is adequate the
whitened residuals already carry \\n - p\\ degrees of freedom and
reducing them again would double-count. A filter that leaves residual
correlation behind does cost degrees of freedom, and that shortfall is
measured by the Satterthwaite approximation \$\$\nu = \mathrm{tr}(RV)^2
/ \mathrm{tr}(RVRV)\$\$ applied to the post-whitening covariance \\V\\,
which reduces exactly to \\n - p\\ when \\V = I\\.

## Implementation Notes

- Small sample corrections are applied (n/(n-p) factor)

- For multi-voxel data, computation is vectorized for efficiency

- Compatible with all contrast types (t, F, custom)

## References

Huber, P. J. (1967). The behavior of maximum likelihood estimates under
nonstandard conditions. Proceedings of the Fifth Berkeley Symposium on
Mathematical Statistics and Probability.

White, H. (1980). A heteroskedasticity-consistent covariance matrix
estimator and a direct test for heteroskedasticity. Econometrica, 48(4),
817-838.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit model with robust regression
cfg <- fmri_lm_control(
  robust = list(
    type = "bisquare",
    c_tukey = 4.685
  )
)

fit <- fmri_lm(model, dataset, config = cfg)

# Standard errors in fit$betas will use sandwich variance
# P-values will use effective degrees of freedom
} # }
```
