# Write Results from fMRI Linear Model

Exports statistical maps from an fmri_lm object to HDF5 files with
BIDS-compliant naming and JSON metadata sidecars using the fmristore
LabeledVolumeSet infrastructure.

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
  strategy = c("by_stat", "by_contrast"),
  save_betas = TRUE,
  contrasts = NULL,
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

- strategy:

  Storage strategy: "by_stat" (group contrasts by statistic) or
  "by_contrast" (separate files)

- save_betas:

  Logical. Save raw regressor betas (default: TRUE)

- contrasts:

  Character vector of contrast names to save. NULL saves all contrasts

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
