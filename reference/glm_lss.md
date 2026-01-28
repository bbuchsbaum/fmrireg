# GLM LSS Estimation Convenience Function (Single Trial Estimation)

A convenience wrapper around `estimate_betas` for least squares separate
(LSS) estimation. **This is primarily designed for single trial
estimation**, where each individual trial/event gets its own separate
beta estimate rather than averaging across trials of the same condition.

## Usage

``` r
glm_lss(
  dataset,
  model_obj,
  basis_obj,
  basemod = NULL,
  block = ~1,
  use_cpp = FALSE,
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

- use_cpp:

  Deprecated. The C++ implementation has been retired. This parameter is
  ignored; fmrilss is always used.

- progress:

  Logical; show progress bar (default: TRUE)

- ...:

  Additional arguments passed to `estimate_betas`

## Value

A list of class "fmri_betas" containing the estimated trial-wise
coefficients

## Details

**Primary Use Case - Single Trial Estimation:**

- **Trial-wise beta estimation**: Each trial gets its own beta
  coefficient

- **Single trial analysis**: Useful for decoding, representational
  similarity analysis (RSA)

- **Trial-by-trial variability**: Captures individual trial responses
  rather than condition averages

- **Avoiding trial averaging**: Preserves trial-specific information
  that would be lost in standard GLM

**Method Details:** LSS (Least Squares Separate) fits a separate model
for each trial, where the trial of interest gets its own regressor while
all other trials of the same condition are modeled together. This
approach avoids the collinearity issues that would arise from including
separate regressors for every trial simultaneously.

For standard condition-level estimation (averaging trials within
conditions), use
[`glm_ols()`](https://bbuchsbaum.github.io/fmrireg/reference/glm_ols.md)
instead.

## See also

[`estimate_betas`](https://bbuchsbaum.github.io/fmrireg/reference/estimate_betas.md)
for the underlying estimation function,
[`glm_ols`](https://bbuchsbaum.github.io/fmrireg/reference/glm_ols.md)
for condition-level estimation

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

# Fit with LSS - estimates separate beta for each individual trial
fit <- glm_lss(dset, model_obj, fmrihrf::HRF_SPMG1)
dim(fit$betas_ran)  # 4 trials x 10 voxels (NOT averaged by condition)

# This is useful for:
# - Decoding analysis (predicting condition from single trial patterns)
# - RSA (representational similarity analysis)
# - Studying trial-by-trial variability
} # }
```
