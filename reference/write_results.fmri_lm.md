# Write Results from fMRI Linear Model

Exports statistical maps from an fmri_lm object with BIDS-compliant
naming and JSON metadata sidecars. Supported image outputs are HDF5
LabeledVolumeSet files and NIfTI files; an fmrigds representation can
also be requested.

## Usage

``` r
# S3 method for class 'fmri_lm'
write_results(
  x,
  path = NULL,
  subject = NULL,
  task = NULL,
  space = NULL,
  desc = "GLM",
  format = c("h5"),
  strategy = c("by_stat", "by_contrast"),
  save_betas = TRUE,
  contrasts = NULL,
  contrast_match = c("auto", "exact", "regex"),
  contrast_stats = c("beta", "tstat", "pval", "se"),
  overwrite = FALSE,
  validate_inputs = TRUE,
  ...
)
```

## Arguments

- x:

  An fmri_lm object containing fitted model results

- path:

  Output directory path. If NULL, uses current working directory

- subject:

  Subject identifier (e.g., "01", "1001"). Required.

- task:

  Task identifier (e.g., "nback", "rest"). Required for BIDS compliance.

- space:

  Spatial reference (e.g., "MNI152NLin2009cAsym"). Optional but
  recommended.

- desc:

  Description of the analysis (default: "GLM")

- format:

  Output format(s). Use `"h5"` for BIDS-style HDF5 outputs, `"nifti"`
  for BIDS-style NIfTI outputs, `"gds"` for fmrigds-compatible assays
  plus an `.rds` plan, or a character vector to write multiple formats.

- strategy:

  Storage strategy: "by_stat" (group contrasts by statistic) or
  "by_contrast" (separate files)

- save_betas:

  Logical. Save raw regressor betas (default: TRUE)

- contrasts:

  Character vector of contrast names to save. NULL saves all contrasts

- contrast_match:

  How to match `contrasts`: `"auto"` first matches exact contrast names
  and treats unmatched selectors as regular expressions; `"exact"` uses
  literal names only; `"regex"` treats all selectors as regular
  expressions.

- contrast_stats:

  Character vector of contrast statistics to save (default: c("beta",
  "tstat", "pval", "se"))

- overwrite:

  Logical. Overwrite existing files (default: FALSE)

- validate_inputs:

  Logical. Validate fmrilm object structure (default: TRUE)

- ...:

  Additional arguments passed to internal functions

## Value

Invisible list of file paths created

## Examples

``` r
if (FALSE) { # \dontrun{
# Save all results using default settings
write_results(fitted_model, subject = "01", task = "nback")

# Save only specific contrasts and statistics  
write_results(fitted_model, 
              subject = "01", task = "nback", space = "MNI152NLin2009cAsym",
              contrasts = c("FacesVsPlaces", "GoVsNoGo"),
              contrast_stats = c("beta", "tstat"))
} # }
```
