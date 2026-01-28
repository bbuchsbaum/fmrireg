# Compute DVARS (Temporal Derivative of Timecourses)

DVARS measures the root mean square of the temporal derivative across
all voxels, providing a single quality metric per volume. High DVARS
indicates rapid signal changes, often associated with motion artifacts.

## Usage

``` r
compute_dvars(Y, normalize = TRUE)
```

## Arguments

- Y:

  Numeric matrix of fMRI data (time x voxels).

- normalize:

  Logical. If TRUE, normalize by median DVARS. Default TRUE.

## Value

Numeric vector of length nrow(Y) with DVARS values. The first timepoint
is set to the median of subsequent values.

## Examples

``` r
# Simulate fMRI data with a motion spike
set.seed(123)
Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
Y[50, ] <- Y[50, ] + 5  # Add spike at volume 50
dvars <- compute_dvars(Y)
#> Error in compute_dvars(Y): could not find function "compute_dvars"
# dvars[50] will be elevated compared to other volumes

# Identify suspect volumes (normalized DVARS > 1.5)
suspect <- which(dvars > 1.5)
#> Error: object 'dvars' not found
cat("Suspect volumes:", suspect, "\n")
#> Error: object 'suspect' not found
```
