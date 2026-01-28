# Compute an HRF smoothing kernel

This function computes a temporal similarity matrix from a series of
hemodynamic response functions.

## Usage

``` r
hrf_smoothing_kernel(
  len,
  TR = 2,
  form = onset ~ trialwise(),
  buffer_scans = 3L,
  normalise = TRUE,
  method = c("gram", "cosine")
)
```

## Arguments

- len:

  The number of scans.

- TR:

  The repetition time (default is 2 seconds).

- form:

  the `trialwise` formula expression, see examples.

- buffer_scans:

  The number of scans to buffer before and after the event.

- normalise:

  Whether to normalise the kernel.

- method:

  The method to use for computing the kernel.

## Value

a smoothing matrix

## Examples

``` r
form <- onset ~ trialwise(basis="gaussian")
sk <- hrf_smoothing_kernel(100, TR=1.5, form)
```
