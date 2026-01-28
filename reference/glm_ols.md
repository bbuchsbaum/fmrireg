# GLM OLS Estimation Convenience Function

A convenience wrapper around `estimate_betas` for ordinary least squares
(OLS) estimation. This function provides a simplified interface for
fitting GLMs using OLS on matrix datasets.

## Usage

``` r
glm_ols(
  dataset,
  model_obj,
  basis_obj,
  basemod = NULL,
  block = ~1,
  progress = TRUE,
  ...
)
```

## Arguments

- dataset:

  A `matrix_dataset` object containing the fMRI time series data

- model_obj:

  An `event_model` object specifying the experimental design

- basis_obj:

  An HRF basis object (e.g., from
  [`fmrihrf::HRF_SPMG1`](https://bbuchsbaum.github.io/fmrihrf/reference/HRF_objects.html),
  `HRF_FIR`, etc.)

- basemod:

  A `baseline_model` instance to regress out of data before beta
  estimation (default: NULL)

- block:

  A formula specifying the block factor (default: ~ 1 for single block)

- progress:

  Logical; show progress bar (default: TRUE)

- ...:

  Additional arguments passed to `estimate_betas`

## Value

A list of class "fmri_betas" containing the estimated coefficients

## Details

**Use Cases:**

- **Condition-level estimation**: Estimates average responses for each
  experimental condition

- **General linear modeling**: Standard GLM approach for group-level or
  condition-level effects

- **Multi-trial averaging**: Combines trials of the same condition to
  estimate mean responses

For single-trial estimation where each trial gets its own beta estimate,
use
[`glm_lss()`](https://bbuchsbaum.github.io/fmrireg/reference/glm_lss.md)
instead.

## See also

[`estimate_betas`](https://bbuchsbaum.github.io/fmrireg/reference/estimate_betas.md)
for the underlying estimation function,
[`glm_lss`](https://bbuchsbaum.github.io/fmrireg/reference/glm_lss.md)
for single trial estimation

## Examples

``` r
if (FALSE) { # \dontrun{
# Create event model and data
event_data <- data.frame(
  onset = c(10, 30, 50, 70),
  condition = factor(c("A", "B", "A", "B")),
  run = rep(1, 4)
)
sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
model_obj <- event_model(onset ~ hrf(condition), 
                        data = event_data, 
                        block = ~ run, 
                        sampling_frame = sframe)

# Create data matrix (100 timepoints, 10 voxels)
Y <- matrix(rnorm(1000), 100, 10)

# Create matrix_dataset with event table
dset <- matrix_dataset(Y, TR = 2, run_length = 100, event_table = event_data)

# Fit with OLS - estimates average response for each condition
fit <- glm_ols(dset, model_obj, fmrihrf::HRF_SPMG1)
dim(fit$betas_ran)  # 2 conditions x 10 voxels
} # }
```
