# Create Group Dataset from NIfTI Files

Creates a group dataset from NIfTI files containing effect sizes and
their standard errors or variances. Supports various input
configurations including beta/SE pairs, beta/variance pairs, or
t-statistics with degrees of freedom.

## Usage

``` r
group_data_from_nifti(
  beta_paths = NULL,
  se_paths = NULL,
  var_paths = NULL,
  t_paths = NULL,
  df = NULL,
  subjects = NULL,
  covariates = NULL,
  mask = NULL,
  target_space = NULL,
  validate = TRUE
)
```

## Arguments

- beta_paths:

  Character vector of paths to beta/effect size NIfTI files

- se_paths:

  Character vector of paths to standard error NIfTI files

- var_paths:

  Character vector of paths to variance NIfTI files (alternative to
  se_paths)

- t_paths:

  Character vector of paths to t-statistic NIfTI files

- df:

  Degrees of freedom (scalar or vector). Required if using t_paths.

- subjects:

  Character vector of subject identifiers. If NULL, extracted from file
  paths.

- covariates:

  Data frame of subject-level covariates (optional)

- mask:

  Path to mask NIfTI file or mask object (optional but recommended)

- target_space:

  Path to template NIfTI for spatial alignment checking

- validate:

  Logical. Validate that all files exist and have matching dimensions
  (default: TRUE)

## Value

A group_data_nifti object

## Examples

``` r
if (FALSE) { # \dontrun{
# From FSL FEAT output (COPE and VARCOPE files)
gd <- group_data_from_nifti(
  beta_paths = Sys.glob("feat_output/sub-*/cope1.nii.gz"),
  var_paths = Sys.glob("feat_output/sub-*/varcope1.nii.gz"),
  subjects = sprintf("sub-%02d", 1:20),
  mask = "group_mask.nii.gz"
)

# From SPM contrast images
gd <- group_data_from_nifti(
  beta_paths = Sys.glob("SPM/sub*/con_0001.nii"),
  se_paths = Sys.glob("SPM/sub*/se_0001.nii"),
  mask = "SPM/mask.nii"
)

# From t-statistics only (for Stouffer's Z combination)
gd <- group_data_from_nifti(
  t_paths = Sys.glob("stats/sub-*/tstat1.nii.gz"),
  df = 100,  # Or vector of per-subject df
  mask = "group_mask.nii.gz"
)
} # }
```
