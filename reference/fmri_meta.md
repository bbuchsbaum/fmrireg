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
  n_threads = getOption("fmrireg.num_threads", 0),
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

  - "reml": Restricted maximum likelihood random-effects

- robust:

  Character string specifying robust estimation:

  - "none": No robust estimation (default)

  - "huber": Huber M-estimator with IRLS

  - "t": Student-t mixture model (for heavy-tailed distributions)

- weights:

  Character string specifying weighting scheme:

  - "ivw": Inverse variance weighting (default)

  - "equal": Equal weights for all subjects

  - "custom": User-provided weights (must supply weights argument)

- weights_custom:

  Numeric vector or matrix of custom weights (required if
  `weights = "custom"`). If a vector, length must equal the number of
  subjects. If a matrix, it must be subjects x features.

- combine:

  For t-statistic-only data, combination method ("stouffer", "fisher",
  or "lancaster"). Stouffer combines z-scores and supports equal,
  inverse-variance, or custom weighting (via `weights`). Fisher uses
  equal weights. Lancaster implements a weighted Fisher method by
  mapping weights to per-subject degrees of freedom.

- contrasts:

  Optional numeric vector or matrix specifying fit-time exact contrasts.
  If a vector is provided, its names must match the column names of the
  design matrix X. A matrix should have columns corresponding to
  predictors and rows corresponding to contrasts.

- return_cov:

  Optional. If set to "tri", returns the packed upper-triangular
  Var(beta) per feature under `$cov` to enable exact post-hoc contrasts
  via
  [`contrast()`](https://bbuchsbaum.github.io/fmrireg/reference/contrast.md).

- chunk_size:

  Number of voxels to process at once (default: 10000)

- n_threads:

  Number of parallel threads to use. Defaults to fmrireg.num_threads
  option.

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
con <- contrast(fit_cov, c("(Intercept)" = 0, group = 1))

# Exact fit-time contrast without storing covariance
fit_con <- fmri_meta(gd, formula = ~ 1 + group, method = "pm",
                     contrasts = c("(Intercept)" = 0, group = 1))
} # }
```
