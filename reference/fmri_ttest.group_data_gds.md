# fmri_ttest method for gds-backed group_data

Delegates to fmrigds reducers for meta or classical (OLS) engines and
wraps outputs into an fmri_ttest_fit-compatible object.

## Usage

``` r
fmri_ttest.group_data_gds(
  gd,
  formula = ~1,
  engine = c("auto", "meta", "classic", "welch"),
  paired = FALSE,
  mu0 = 0,
  contrast = NULL,
  mc = NULL,
  alpha = 0.05,
  sign = c("AminusB", "BminusA"),
  mask = NULL,
  voxelwise_cov = NULL,
  center_voxelwise = TRUE,
  voxel_name = "voxel_cov",
  weights = c("ivw", "equal", "custom"),
  weights_custom = NULL,
  combine = NULL
)
```

## Arguments

- gd:

  A group_data object or data frame with subject data

- formula:

  R formula for between-subjects design. Examples:

  - ~ 1: One-sample t-test

  - ~ 1 + group: Two-sample t-test

  - ~ 1 + group + age: ANCOVA with covariate

- engine:

  Character string; analysis engine to use:

  - "auto": Automatically select based on data (default)

  - "meta": Use inverse-variance meta-analysis (requires SE)

  - "classic": Use OLS/Student t-test

  - "welch": Use Welch t-test for unequal variances

- paired:

  Logical; if TRUE, perform paired t-test on differences (default:
  FALSE)

- mu0:

  Numeric scalar; constant to test against for one-sample tests
  (default: 0)

- contrast:

  Optional named numeric vector for linear combination of coefficients

- mc:

  Character string or NULL; multiple comparisons correction:

  - NULL: No correction (default)

  - "bh": Benjamini-Hochberg FDR

  - "by": Benjamini-Yekutieli FDR

  - "spatial_fdr": Spatially-aware FDR

- alpha:

  Numeric scalar; FDR level if mc is not NULL (default: 0.05)

- sign:

  Character string; sign convention for group differences:

  - "AminusB": A - B (default)

  - "BminusA": B - A

- mask:

  Optional mask object or path

- voxelwise_cov:

  Optional S x P matrix of voxelwise covariates

- center_voxelwise:

  Logical; center voxelwise covariate per feature (default: TRUE)

- voxel_name:

  Character string; name for voxelwise coefficient (default:
  "voxel_cov")

- weights:

  Character string for meta-engine weighting: "ivw" (default, uses
  provided SE), "equal" (equal weights), or "custom" (supply
  `weights_custom`). Ignored for classic/welch engines.

- weights_custom:

  Numeric vector (length S) or matrix (S x P) of custom weights when
  `weights = "custom"`.

- combine:

  Optional. When using the meta engine with t-only inputs (i.e.,
  per-subject t-statistics and df), specify the t-combination method:
  "stouffer", "fisher", or "lancaster". Passed through to
  [`fmri_meta()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_meta.md)
  when delegating to the meta engine on group_data\_\*.
