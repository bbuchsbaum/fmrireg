# Data Structures and Sampling Frames

## Introduction: Linking Data and Design

Effective fMRI analysis requires associating the measured brain activity
(the imaging data) with crucial metadata, including:

- **Temporal Structure:** When each scan was acquired (TR) and how scans
  are grouped into runs.
- **Spatial Structure:** Which brain locations (voxels) are included in
  the analysis (mask).
- **Experimental Design:** Timing and properties of experimental events
  or conditions.

`fmrireg` represents all of this with one container, the `fmri_frame`
from the companion `fmridataset` package. A frame is an *observations by
features* matrix (scans by voxels, vertices, parcels, or components)
whose acquisition timing lives in the observation metadata, whose
spatial identity lives in a typed *feature space*, and whose
experimental events live in a keyed event table. Every modeling function
in `fmrireg` (`fmri_lm`, `estimate_betas`, `estimate_hrf`,
`fmri_latent_lm`, …) takes an `fmri_frame` as its dataset.

This vignette describes how to build a frame from the data you already
have and how to read its timing, events, and spatial layout back. For a
first complete fit, read
[`vignette("fmrireg", package = "fmrireg")`](https://bbuchsbaum.github.io/fmrireg/articles/fmrireg.md);
return here when you need to change how scans are stored or divided into
runs. For the frame model itself (views, laziness, persistence) see
[`vignette("fmridataset", package = "fmridataset")`](https://bbuchsbaum.github.io/fmridataset/articles/fmridataset.html).

| data_you_have | constructor | feature_space | choose_it_when |
|:---|:---|:---|:---|
| NeuroVec runs already in memory | neurovec_frame() | volume_space | The full volumes comfortably fit in memory |
| NIfTI runs and a mask on disk | nifti_frame() | volume_space | Runs are large and should be read lazily |
| Time-by-feature matrix | matrix_frame() | index_space | Data are ROI, surface, or already tabular |
| Latent components plus spatial loadings | latent_frame() | basis_space | A fmristore LatentNeuroVec already exists |
| fMRIPrep derivatives (BIDS) | fmridataset::read_bids_bold() | volume_space | One subject’s preprocessed BOLD runs are on disk |

Choose a frame constructor from the storage you already have {.table}

## The `sampling_frame`

Recall the `sampling_frame` object (introduced in the Overview and
detailed in other vignettes). It defines the temporal structure that the
design machinery expects:

- `blocklens`: A vector specifying the number of scans (time points) in
  each run.
- `TR`: The repetition time (time between scans) in seconds.

``` r

sframe_example <- sampling_frame(blocklens = c(150, 160), TR = 2.0)
print(sframe_example)
#> Sampling frame
#> - Blocks: 2 
#> - Scans: 310 (per block: 150, 160 )
#> - TR: 2 s
#> - Duration: 619 s
```

A frame does not store a `sampling_frame`. It stores the same facts as
two observation columns, `run_id` and `TR`, and
[`fmridataset::as_sampling_frame()`](https://bbuchsbaum.github.io/fmridataset/reference/temporal-schema.html)
rebuilds the `sampling_frame` from them whenever the design code needs
one. Because the columns are the truth, the timing follows subsetting
and reordering of the frame for free.

## In-Memory Volumetric Data (`neurovec_frame`)

Use this when your fMRI runs are loaded as
[`neuroim2::NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.html)
objects in your R session.

**Key Arguments:**

- `scans`: A `NeuroVec`, or a *list* of them, one for each run.
- `mask`: A
  [`neuroim2::LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.html)
  (or a `NeuroVol` whose non-zero voxels form the mask).
- `TR`: Repetition time (seconds), one value or one per run.
- `run_length` (Optional): Vector of run lengths; if omitted, inferred
  from the fourth dimension of each scan.
- `event_table` (Optional): A `data.frame` containing experimental
  design information.
- `censor` (Optional): One entry per scan marking volumes to exclude
  (logical, 0/1, or 1-based indices).

``` r

# Create minimal example data (2 runs)
d <- c(5, 5, 5, 20) # Small dimensions for example
mask_vol <- neuroim2::LogicalNeuroVol(array(TRUE, d[1:3]), neuroim2::NeuroSpace(d[1:3]))

scan1 <- neuroim2::NeuroVec(array(rnorm(prod(d)), d), neuroim2::NeuroSpace(d))
scan2 <- neuroim2::NeuroVec(array(rnorm(prod(d)), d), neuroim2::NeuroSpace(d))

# Example event table
events_df <- data.frame(
  onset = c(5, 15, 5, 15),
  condition = factor(c("A", "B", "A", "B")),
  run = c(1, 1, 2, 2)
)

# Create the frame
mem_frame <- neurovec_frame(scans = list(scan1, scan2),
                            mask = mask_vol,
                            TR = 2.0,
                            # run_length automatically inferred as c(20, 20)
                            event_table = events_df)

dim(mem_frame)              # scans by masked voxels
#> [1]  40 125
class(fmridataset::space(mem_frame))     # a volume_space on the mask's grid
#> [1] "volume_space"  "feature_space"
as_sampling_frame(mem_frame)
#> Sampling frame
#> - Blocks: 2 
#> - Scans: 40 (per block: 20, 20 )
#> - TR: 2 s
#> - Duration: 79 s
```

The masked voxel series are extracted once, with
[`neuroim2::series()`](https://bbuchsbaum.github.io/neuroim2/reference/series-methods.html),
so a `SparseNeuroVec` does not need to be densified.

## File-Based Volumetric Data (`nifti_frame`)

This is often the most practical option for typical fMRI analyses where
data resides in files.

**Key Arguments:**

- `scans`: A character vector of file paths to the 4D fMRI image files
  (e.g., `.nii.gz`), one path per run.
- `mask`: The path to a 3D mask image (or a `LogicalNeuroVol`, or a
  [`fmridataset::volume_space`](https://bbuchsbaum.github.io/fmridataset/reference/volume_space.html)).
- `TR`: Repetition time (seconds).
- `run_length` (Optional): Volumes per run; if omitted, read from each
  file’s header.
- `event_table` (Optional): A `data.frame` with experimental design
  info.
- `base_path` (Optional): A path to prepend to relative file paths in
  `scans` and `mask`.

The hidden fixture above creates tiny valid NIfTI files so the
constructor below genuinely executes. In a study, start with your
existing run and mask paths; the important teaching code is the
constructor itself.

``` r

# Pass only filenames to 'scans' and 'mask', and specify the directory in 'base_path'
file_frame <- nifti_frame(scans = c(scan1_filename, scan2_filename),
                          mask = mask_filename,
                          TR = 1.5,
                          run_length = c(20, 25), # Must match time dim of files
                          event_table = events_df,
                          base_path = tmp_dir)    # Set base_path to the temp directory

dim(file_frame)
#> [1]  45 125
temporal_schema(file_frame)$run_lengths
#> run-1 run-2 
#>    20    25
```

The frame wraps a
[`fmridataset::nifti_array_source()`](https://bbuchsbaum.github.io/fmridataset/reference/nifti_array_source.html):
headers and the mask are read at construction, but no BOLD volumes are
read until a fit needs them, and reads are pushed down per file and per
masked voxel. If the files change on disk after construction, the next
read fails with a structured stale-source error rather than silently
returning different values.

## Matrix Data (`matrix_frame`)

Use this if your fMRI data is already represented as a 2D matrix where
rows are time points and columns are voxels or components (e.g., after
surface projection or ROI averaging).

**Key Arguments:**

- `datamat`: The numeric matrix (time x features).
- `TR`: Repetition time (seconds).
- `run_length`: Vector specifying the number of rows (time points)
  belonging to each run.
- `event_table` (Optional): One row per experimental event, with onset,
  condition, run, and any modulators used by the model. Its row count
  usually differs from the number of scans; onsets and run labels must
  instead be valid for the declared run lengths.
- `feature_ids` (Optional): Stable feature IDs; unique column names are
  used when present.

``` r

# Example matrix (100 time points, 50 features/voxels)
# Two runs of 50 time points each
time_points <- 100
features <- 50
run_len <- c(50, 50)
example_matrix <- matrix(rnorm(time_points * features), time_points, features)

# Example event table for matrix data
events_mat_df <- data.frame(
  onset = c(seq(5, 45, by=10), seq(5, 45, by=10)),
  condition = factor(rep(c("C", "D"), 10)),
  run = rep(1:2, each = 5)
)

mat_frame <- matrix_frame(datamat = example_matrix,
                          TR = 2.5,
                          run_length = run_len,
                          event_table = events_mat_df)

dim(mat_frame)
#> [1] 100  50
class(fmridataset::space(mat_frame))   # an index_space: no spatial geometry
#> [1] "index_space"   "feature_space"
as_sampling_frame(mat_frame)
#> Sampling frame
#> - Blocks: 2 
#> - Scans: 100 (per block: 50, 50 )
#> - TR: 2.5 s
#> - Duration: 248.75 s
```

For `matrix_frame`, the concept of a spatial mask is implicit; all
columns provided in `datamat` are included, and the feature space is a
plain `index_space`.

## Latent Data (`latent_frame`)

This constructor is for data that has undergone dimensionality reduction
(e.g., PCA, ICA). It takes a
[`fmristore::LatentNeuroVec`](https://rdrr.io/pkg/fmrilatent/man/LatentNeuroVec.html),
which stores the basis (latent components over time) and loadings
(spatial maps of components), and builds a frame whose assay holds the
component scores and whose feature space is a
[`fmridataset::basis_space`](https://bbuchsbaum.github.io/fmridataset/reference/basis_space.html)
carrying the loadings as its synthesis operator.

**Key Arguments:**

- `x`: A `LatentNeuroVec` object from the `fmristore` package.
- `TR`: Repetition time (seconds).
- `run_length`: Vector specifying run lengths (must sum to the time
  dimension of `x`).
- `event_table` (Optional): Experimental design `data.frame`.

``` r

# Conceptual example (requires fmristore package and a LatentNeuroVec)
# Assuming 'my_latent_neuro_vec' is a LatentNeuroVec object representing
# 20 components over 300 time points (2 runs of 150)

latent_frame_obj <- latent_frame(
  my_latent_neuro_vec,
  TR = 2.0,
  run_length = c(150, 150),
  event_table = some_event_df
)
collect_assay(latent_frame_obj)              # scores: time x components
basis_synthesis(fmridataset::space(latent_frame_obj))     # loadings: voxels x components
```

Model fitting happens in component space;
[`fmri_latent_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_latent_lm.md)
and the `latent_sketch` engine reconstruct voxel-wise coefficients
through the loadings when you ask for them.

## Reading a frame back

Once created, a frame is the primary data input for `fmrireg`’s modeling
functions, and the same public accessors work whatever the storage:

``` r

# Timing: derived from the run_id / TR observation columns
schema <- temporal_schema(mem_frame)
schema$run_lengths
#> run-1 run-2 
#>    20    20
schema$TR
#> run-1 run-2 
#>     2     2
blocklens(as_sampling_frame(mem_frame))
#> [1] 20 20

# Events: the keyed event table (an event_id column is added when absent)
head(event_data(mem_frame$tables$events))
#> # A tibble: 4 × 4
#>   event_id onset condition   run
#>   <chr>    <dbl> <fct>     <dbl>
#> 1 event-1      5 A             1
#> 2 event-2     15 B             1
#> 3 event-3      5 A             2
#> 4 event-4     15 B             2

# Data: the dense scans-by-features matrix, read under a memory budget
dim(collect_assay(mem_frame))
#> [1]  40 125

# Space: the typed feature space and, for volumes, a map back to the grid
n_features(fmridataset::space(mem_frame))
#> [1] 125
first_volume <- spatial_map(mem_frame, 1)
class(first_volume)
#> [1] "DenseNeuroVol"
#> attr(,"package")
#> [1] "neuroim2"
```

Modeling functions take the frame directly:

- `event_model(..., sampling_frame = as_sampling_frame(frame))`
- `baseline_model(..., sframe = as_sampling_frame(frame))`
- `fmri_lm(formula, block, dataset = frame)`
- `estimate_betas(frame, fixed = ..., ran = ..., block = ...)`

Frames are lazy where they can be.
`filter_obs(frame, run_id == "run-2")` and `frame[, 1:100]` return views
that share the source and read nothing until `collect_assay()` or a fit
asks for values.

## Next

- [`vignette("a_08_simulation", package = "fmrireg")`](https://bbuchsbaum.github.io/fmrireg/articles/a_08_simulation.md)
  — Simulating fMRI data
- [`vignette("a_09_linear_model", package = "fmrireg")`](https://bbuchsbaum.github.io/fmrireg/articles/a_09_linear_model.md)
  — fMRI Linear Model (GLM)
- [`vignette("fmridataset", package = "fmridataset")`](https://bbuchsbaum.github.io/fmridataset/articles/fmridataset.html)
  — Frames, views, and the temporal contract
