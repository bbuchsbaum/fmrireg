# Fit Group-Level Meta-Analysis

Performs voxelwise or ROI-based meta-analysis on group fMRI data using
fixed-effects, random-effects, or robust methods. Supports
meta-regression with covariates for group comparisons and other
moderator analyses.

## Usage

``` r
fmri_meta(
  data,
  formula = ~1,
  method = c("pm", "fe", "dl", "reml"),
  robust = c("none", "huber", "t"),
  weights = c("ivw", "equal", "custom"),
  weights_custom = NULL,
  combine = NULL,
  contrasts = NULL,
  return_cov = NULL,
  chunk_size = 10000,
  n_threads = .fmrireg_num_threads(),
  verbose = TRUE
)
```

## Arguments

- data:

  A group_data object created by
  [`group_data`](https://bbuchsbaum.github.io/fmrireg/reference/group_data.md)

- formula:

  Formula specifying the meta-regression model. Default is ~ 1
  (intercept only). Use ~ 1 + group for group comparisons, or include
  continuous covariates.

- method:

  Character string specifying the meta-analysis method:

  - "fe": Fixed-effects (inverse variance weighted)

  - "pm": Paule-Mandel random-effects (default, good for whole-brain)

  - "dl": DerSimonian-Laird random-effects

  - "reml": Not implemented. Accepted for backward compatibility; it
    runs the Paule-Mandel estimator, says so with a message, and the fit
    records `method = "pm"`.

- robust:

  Character string specifying robust estimation:

  - "none": No robust estimation (default)

  - "huber": Huber M-estimator with IRLS

  - "t": Student-t mixture model (for heavy-tailed distributions)

- weights:

  Character string specifying weighting scheme:

  - "ivw": Inverse variance weighting (default). The subject-level
    variances are taken as known and `method` chooses the
    between-subject variance estimator.

  - "equal": Equal weights for all subjects

  - "custom": User-provided weights (must supply weights argument)

  "equal" and "custom" weights are *relative*: the model is weighted
  least squares, \\Var(y_i) = \sigma^2 / w_i\\, with the residual scale
  \\\sigma^2\\ estimated per feature (as `lm(weights = w)` does; equal
  weights give the one-sample / OLS group t-statistic). The
  subject-level standard errors and `method` are then not used, and
  `tau2`, `Q` and `I2` are `NA`.

- weights_custom:

  Numeric vector or matrix of positive custom weights (required if
  `weights = "custom"`). If a vector, length must equal the number of
  subjects. If a matrix, it must be subjects x features (analysed voxels
  or ROIs).

- combine:

  For t-statistic-only data, combination method ("stouffer", "fisher",
  or "lancaster"). Stouffer combines z-scores and supports equal or
  custom weighting (via `weights`); inverse-variance weights are
  unavailable without effect standard errors. Fisher uses equal weights.
  Lancaster implements a weighted Fisher method by mapping weights to
  per-subject degrees of freedom.

- contrasts:

  Optional numeric vector or matrix specifying fit-time exact contrasts.
  If a vector is provided, its names must match the column names of the
  design matrix X. A matrix should have columns corresponding to
  predictors and rows corresponding to contrasts.

- return_cov:

  Optional; kept for backward compatibility. The packed upper-triangular
  Var(beta) per feature is now always returned under `$cov` (for
  subject-level effect data) so that
  [`contrast()`](https://bbuchsbaum.github.io/fmrireg/reference/contrast.md)
  gives exact standard errors.

- chunk_size:

  Number of voxels to process at once (default: 10000)

- n_threads:

  Number of OpenMP threads for the compiled fitter. Defaults to
  `options(fmrireg.num_threads)`, else the environment variable
  `FMRIREG_NUM_THREADS`, else 0 (the OpenMP default).

- verbose:

  Logical. Print progress messages (default: TRUE)

## Value

An fmri_meta object containing:

- coefficients: Meta-regression coefficients

- se: Standard errors

- tau2: Between-study variance (for random-effects)

- I2: I-squared heterogeneity statistic

- Q: Cochran's Q statistic

- model: Model specification

- data: Input group_data object

For fmrigds-backed group data (from
[`group_data()`](https://bbuchsbaum.github.io/fmrireg/reference/group_data.md))
the fit is the fmrigds reducer `"fmrireg:meta_reg"`. When such data hold
several subject-level contrasts, each is fitted and the result is an
`fmri_meta_list`: a list of `fmri_meta` objects named by contrast.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simple fixed-effects meta-analysis
fit <- fmri_meta(gd, method = "fe")

# Random-effects with group comparison
fit <- fmri_meta(gd, formula = ~ 1 + group, method = "pm")

# Robust meta-regression with continuous covariate
fit <- fmri_meta(gd, formula = ~ 1 + age + sex, method = "reml", robust = "huber")

# Stouffer's Z for t-statistics only
fit <- fmri_meta(gd_tstat, combine = "stouffer")

# Exact post-hoc contrasts by storing covariance
fit_cov <- fmri_meta(gd, formula = ~ 1 + group, method = "pm", return_cov = "tri")
con <- contrast(fit_cov, c("(Intercept)" = 0, groupB = 1))

# Exact fit-time contrast without storing covariance
fit_con <- fmri_meta(gd, formula = ~ 1 + group, method = "pm",
                     contrasts = c("(Intercept)" = 0, groupB = 1))
} # }
```
