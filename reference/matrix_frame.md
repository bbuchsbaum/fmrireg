# Build an fmri_frame from an in-memory matrix

fmrireg fits models to
[`fmridataset::fmri_frame()`](https://bbuchsbaum.github.io/fmridataset/reference/fmri_frame.html)
objects. This helper builds one from a time-by-feature matrix plus the
run structure and event table that the design machinery needs.
Observation IDs are derived deterministically as `run-<r>-vol-<index>`;
feature IDs default to the matrix column names when they are unique and
to `feature-<j>` otherwise.

## Usage

``` r
matrix_frame(
  datamat,
  TR,
  run_length,
  event_table = NULL,
  censor = NULL,
  feature_ids = NULL,
  assay = "signal"
)
```

## Arguments

- datamat:

  A numeric matrix with one row per acquired volume and one column per
  feature (voxel, vertex, parcel, or component).

- TR:

  Repetition time in seconds; one value, or one per run.

- run_length:

  Integer vector giving the number of volumes in each run. Must sum to
  `nrow(datamat)`.

- event_table:

  Optional data frame of events (onsets, conditions, and the block
  variable). It is stored on the frame as a keyed
  [`fmridataset::event_table()`](https://bbuchsbaum.github.io/fmridataset/reference/event_table.html);
  an `event_id` column is added when absent.

- censor:

  Optional censoring indicator with one entry per volume: a logical
  vector, a 0/1 vector, or a vector of 1-based volume indices. Stored as
  the logical `censor` observation column that
  [`fmridataset::temporal_schema()`](https://bbuchsbaum.github.io/fmridataset/reference/temporal-schema.html)
  recognises.

- feature_ids:

  Optional stable feature IDs (one per column).

- assay:

  Name of the assay holding `datamat`.

## Value

An `fmri_frame` whose feature space is an
[`fmridataset::index_space()`](https://bbuchsbaum.github.io/fmridataset/reference/index_space.html).

## See also

[`neurovec_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/neurovec_frame.md),
[`nifti_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/nifti_frame.md),
[`latent_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/latent_frame.md),
[`fmridataset::fmri_frame()`](https://bbuchsbaum.github.io/fmridataset/reference/fmri_frame.html)

## Examples

``` r
Y <- matrix(rnorm(80 * 3), 80, 3)
events <- data.frame(onset = c(5, 25, 45, 65),
                     condition = factor(c("A", "B", "A", "B")),
                     run = c(1, 1, 2, 2))
frame <- matrix_frame(Y, TR = 2, run_length = c(40, 40), event_table = events)
fmridataset::temporal_schema(frame)$run_lengths
#> run-1 run-2 
#>    40    40 
```
