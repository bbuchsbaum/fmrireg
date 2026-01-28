# Simulate fMRI Time Series

This function simulates an fMRI time series for multiple experimental
conditions with specified parameters. It generates a realistic
event-related design with randomized inter-stimulus intervals and
condition orders.

## Usage

``` r
simulate_bold_signal(
  ncond,
  hrf = fmrihrf::HRF_SPMG1,
  nreps = 12,
  amps = rep(1, ncond),
  isi = c(3, 6),
  ampsd = 0,
  TR = 1.5
)
```

## Arguments

- ncond:

  The number of conditions to simulate.

- hrf:

  The hemodynamic response function to use (default is
  fmrihrf::HRF_SPMG1).

- nreps:

  The number of repetitions per condition (default is 12).

- amps:

  A vector of amplitudes for each condition (default is a vector of 1s
  with length ncond).

- isi:

  A vector of length 2 specifying the range of inter-stimulus intervals
  to sample from (default is c(3, 6) seconds).

- ampsd:

  The standard deviation of the amplitudes (default is 0).

- TR:

  The repetition time of the fMRI acquisition (default is 1.5 seconds).

## Value

A list with the following components:

- onset: A vector of the onset times for each trial

- condition: A vector of condition labels for each trial

- mat: A matrix containing the simulated fMRI time series:

  - Column 1: Time points (in seconds)

  - Columns 2:(ncond+1): Simulated BOLD responses for each condition

## Examples

``` r
# Simulate 3 conditions with different amplitudes
sim <- simulate_bold_signal(ncond = 3, amps = c(1, 1.5, 2), TR = 2)

# Plot the simulated time series
matplot(sim$mat[,1], sim$mat[,-1], type = "l", 
        xlab = "Time (s)", ylab = "BOLD Response")

```
