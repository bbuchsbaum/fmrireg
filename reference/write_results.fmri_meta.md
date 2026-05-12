# Write Meta-Analysis Results

Exports spatial meta-analysis maps with the same core export contract
used by
[`write_results.fmri_lm`](https://bbuchsbaum.github.io/fmrireg/reference/write_results.fmri_lm.md):
BIDS-style filenames, HDF5 and NIfTI image outputs, JSON metadata
sidecars, vectorized `format`, explicit overwrite handling, and atomic
finalization.

## Usage

``` r
# S3 method for class 'fmri_meta'
write_results(
  x,
  path = NULL,
  subject = NULL,
  task = NULL,
  space = NULL,
  desc = "Meta",
  format = c("h5"),
  strategy = c("by_stat", "by_coefficient"),
  coefficients = NULL,
  coefficient_match = c("auto", "exact", "regex"),
  coefficient_stats = c("beta", "se", "z", "pval"),
  heterogeneity = TRUE,
  overwrite = FALSE,
  validate_inputs = TRUE,
  ...
)
```

## Arguments

- x:

  An fmri_meta object

- path:

  Output directory path. If NULL, uses current working directory.

- subject:

  Optional subject or group identifier. Unlike subject-level GLM
  exports, meta-analysis outputs may omit this for group-level maps.

- task:

  Optional task identifier.

- space:

  Optional spatial reference label.

- desc:

  Description of the analysis (default: "Meta").

- format:

  Output format(s). Use `"h5"` for HDF5 LabeledVolumeSet outputs,
  `"nifti"` for NIfTI outputs, or a character vector to write both.

- strategy:

  Storage strategy. `"by_stat"` writes one file per statistic with
  coefficients along the 4th dimension; `"by_coefficient"` writes one
  file per coefficient with statistics along the 4th dimension.

- coefficients:

  Character vector of coefficient names to save. NULL saves all
  coefficients.

- coefficient_match:

  How to match `coefficients`: `"auto"` first matches exact names and
  treats unmatched selectors as regular expressions; `"exact"` uses
  literal names only; `"regex"` treats all selectors as regular
  expressions.

- coefficient_stats:

  Character vector of coefficient statistics to save. Supported values
  are `"beta"`, `"se"`, `"z"`, and `"pval"`.

- heterogeneity:

  Logical. Save heterogeneity maps (`tau2`, `I2`, `Q`, `Q_df`) when
  available.

- overwrite:

  Logical. Overwrite existing files (default: FALSE).

- validate_inputs:

  Logical. Validate fmri_meta object structure (default: TRUE).

- ...:

  Additional arguments passed to internal functions.

## Value

Invisible list of created files
