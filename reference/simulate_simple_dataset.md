# Simulate Complete fMRI Dataset

This function simulates a complete fMRI dataset by combining
task-related signals with a declared structured-noise model. It returns
both the clean signals and the noisy data.

## Usage

``` r
simulate_simple_dataset(
  ncond,
  nreps = 12,
  TR = 1.5,
  snr = 0.5,
  hrf = fmrihrf::HRF_SPMG1,
  seed = NULL
)
```

## Arguments

- ncond:

  Number of conditions to simulate

- nreps:

  Number of repetitions per condition (default is 12)

- TR:

  Repetition time in seconds (default is 1.5)

- snr:

  Finite positive signal-to-noise ratio, defined as the empirical
  standard deviation of the clean condition signals divided by the
  empirical standard deviation of the complete structured-noise matrix
  (default is 0.5)

- hrf:

  Hemodynamic response function to use (default is fmrihrf::HRF_SPMG1)

- seed:

  Optional seed for reproducibility (default is NULL). Reusing a seed
  with a different `snr` preserves the clean signal and standardized
  noise realization while changing only the noise scale.

## Value

A list containing:

- clean: The simulated signals without noise (from simulate_bold_signal)

- noisy: The signals with added noise

- noise: The simulated noise component

- onsets: Trial onset times

- conditions: Condition labels for each trial

## Examples

``` r
# Simulate a dataset with 3 conditions
data <- simulate_simple_dataset(ncond = 3, TR = 2, snr = 0.5)

# Plot clean and noisy data
old_par <- par(no.readonly = TRUE)
par(mfrow = c(2,1))
matplot(data$clean$mat[,1], data$clean$mat[,-1], type = "l",
        main = "Clean Signal", xlab = "Time (s)", ylab = "BOLD")
matplot(data$noisy[,1], data$noisy[,-1], type = "l",
        main = "Noisy Signal", xlab = "Time (s)", ylab = "BOLD")

par(old_par)
```
