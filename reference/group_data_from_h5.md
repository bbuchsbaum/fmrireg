# Create Group Dataset from HDF5 Files

Creates a group dataset from HDF5 files produced by
[`write_results.fmri_lm`](https://bbuchsbaum.github.io/fmrireg/reference/write_results.fmri_lm.md).
These files use the fmristore LabeledVolumeSet format for efficient
storage of multiple statistical maps.

## Usage

``` r
group_data_from_h5(
  paths,
  subjects = NULL,
  covariates = NULL,
  mask = NULL,
  contrast = NULL,
  stat = c("beta", "se"),
  validate = TRUE
)
```

## Arguments

- paths:

  Character vector of HDF5 file paths, one per subject

- subjects:

  Character vector of subject identifiers. If NULL, extracted from file
  paths.

- covariates:

  Data frame of subject-level covariates (optional)

- mask:

  Path to mask file or mask object (optional)

- contrast:

  Character string specifying which contrast to extract (for
  multi-contrast files)

- stat:

  Character vector of statistics to extract (e.g., c("beta", "se",
  "tstat"))

- validate:

  Logical. Validate that all files exist and contain expected data
  (default: TRUE)

## Value

A group_data_h5 object

## Examples

``` r
if (FALSE) { # \dontrun{
# Read HDF5 files from write_results.fmri_lm
subjects <- data.frame(
  subject = sprintf("sub-%02d", 1:20),
  group = rep(c("young", "old"), each = 10),
  age = c(rnorm(10, 25, 3), rnorm(10, 70, 5))
)

h5_paths <- sprintf("derivatives/sub-%02d_task-nback_desc-GLMstatmap_bold.h5", 1:20)

gd <- group_data_from_h5(
  h5_paths,
  subjects = subjects$subject,
  covariates = subjects[c("group", "age")],
  contrast = "FaceVsPlace",
  stat = c("beta", "se")
)
} # }
```
