# Create Group Dataset for Meta-Analysis

Generic constructor that creates a group dataset from various input
formats for use in group-level meta-analysis with
[`fmri_meta`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_meta.md).

## Usage

``` r
group_data(data, format = c("auto", "h5", "nifti", "csv", "fmrilm"), ...)
```

## Arguments

- data:

  Input data. Format depends on the `format` argument:

  - For "h5": Character vector of HDF5 file paths

  - For "nifti": List or data frame with beta/SE/variance paths

  - For "csv": Path to CSV file or data frame

  - For "fmrilm": not supported; see below

- format:

  Character string specifying the input format. One of "auto" (default),
  "h5", "nifti", "csv", or "fmrilm". If "auto", attempts to detect
  format. `format = "fmrilm"` (and a list of `fmri_lm` fits) is not
  supported and fails immediately: write each subject's fit with
  [`write_results`](https://bbuchsbaum.github.io/fmrireg/reference/write_results.md)
  and load the resulting beta/SE maps with `format = "nifti"` (or HDF5
  files with `format = "h5"`).

- ...:

  Additional arguments passed to format-specific constructors

## Value

A group_data object suitable for meta-analysis

## Examples

``` r
if (FALSE) { # \dontrun{
# From HDF5 files created by write_results.fmri_lm
gd <- group_data(
  h5_paths,
  format = "h5",
  subjects = subject_ids,
  covariates = covariates_df,
  contrast = "FaceVsPlace"
)

# From NIfTI files
gd <- group_data(
  list(beta = beta_paths, se = se_paths),
  format = "nifti",
  subjects = subject_ids,
  mask = "group_mask.nii.gz"
)
} # }
```
