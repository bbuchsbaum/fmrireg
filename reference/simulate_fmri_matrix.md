# Simulate fMRI Time Courses, Return Shared Onsets + Column-Specific Amplitudes/Durations

Generates \\n\\ time-series (columns) with a single set of onsets, but
*resampled* amplitudes/durations for each column if `amplitude_sd>0` or
`duration_sd>0`. Each column also gets independent noise. The result is
a list containing:

- `time_series`: a `matrix_dataset` with \\T \times n\\. The
  `event_table` uses the first column's amplitude/duration draws.

- `ampmat`: an \\n\\events \times n\\ matrix of per-column amplitudes.

- `durmat`: an \\n\\events \times n\\ matrix of per-column durations.

- `hrf_info`: info about the HRF.

- `noise_params`: info about noise generation (type + AR coefficients +
  SD).

## Usage

``` r
simulate_fmri_matrix(
  n = 1,
  total_time = 240,
  TR = 2,
  hrf = fmrihrf::HRF_SPMG1,
  n_events = 10,
  onsets = NULL,
  isi_dist = c("even", "uniform", "exponential"),
  isi_min = 2,
  isi_max = 6,
  isi_rate = 0.25,
  durations = 0,
  duration_sd = 0,
  duration_dist = c("lognormal", "gamma"),
  amplitudes = 1,
  amplitude_sd = 0,
  amplitude_dist = c("lognormal", "gamma", "gaussian"),
  single_trial = FALSE,
  noise_type = c("none", "white", "ar1", "ar2"),
  noise_ar = NULL,
  noise_sd = 1,
  random_seed = NULL,
  verbose = FALSE,
  buffer = 16
)
```

## Arguments

- n:

  Number of time-series (columns).

- total_time:

  Numeric. Total scan length (seconds).

- TR:

  Numeric. Repetition time (seconds).

- hrf:

  Hemodynamic response function, e.g.
  [`fmrihrf::HRF_SPMG1`](https://bbuchsbaum.github.io/fmrihrf/reference/HRF_objects.html).

- n_events:

  Number of events (ignored if `onsets` is provided).

- onsets:

  Optional numeric vector of event onsets. If `NULL`, will be generated.

- isi_dist:

  One of `"even"`, `"uniform"`, or `"exponential"`. Default is `"even"`
  so events are evenly spaced from 0..total_time.

- isi_min, isi_max:

  For `isi_dist="uniform"`.

- isi_rate:

  For `isi_dist="exponential"`.

- durations:

  Numeric, scalar or length-`n_events`. If `duration_sd>0`, random
  sampling is done per column.

- duration_sd:

  Numeric. If \>0, random variation in durations.

- duration_dist:

  `"lognormal"` or `"gamma"` (strictly positive).

- amplitudes:

  Numeric, scalar or length-`n_events`. If `amplitude_sd>0`, random
  sampling is done per column.

- amplitude_sd:

  Numeric. If \>0, random variation in amplitudes.

- amplitude_dist:

  `"lognormal"`, `"gamma"`, or `"gaussian"` (can be negative).

- single_trial:

  If TRUE, each event is a separate single-trial regressor that gets
  summed.

- noise_type:

  `"none"`, `"white"`, `"ar1"`, or `"ar2"`.

- noise_ar:

  Numeric vector for AR(1) or AR(2). If missing or insufficient,
  defaults are used (0.3 for AR(1); c(0.3,0.2) for AR(2)).

- noise_sd:

  Std dev of the noise.

- random_seed:

  Optional integer for reproducibility.

- verbose:

  If TRUE, prints messages.

- buffer:

  Numeric seconds appended to the end of the time grid to avoid edge
  truncation (default: 16).

## Value

A list containing:

- `time_series`:

  A `matrix_dataset` with \\T \times n\\ data and `event_table` for the
  *first* column's random draws.

- `ampmat`:

  An \\n\\events \times n\\ numeric matrix of amplitudes.

- `durmat`:

  An \\n\\events \times n\\ numeric matrix of durations.

- `hrf_info`:

  A list with HRF metadata.

- `noise_params`:

  A list describing noise generation.

## Details

- If `noise_type="ar1"` and you do not provide `noise_ar`, we default to
  `c(0.3)`.

- If `noise_type="ar2"` and you do not provide a 2-element `noise_ar`,
  we default to `c(0.3, 0.2)`.

- Onsets are either provided or generated once for all columns.

- **Amplitudes/durations** are re-sampled *inside the loop* so each
  column can differ randomly. The final arrays `ampmat` and `durmat`
  each have one column per time-series.

- The `matrix_dataset`'s `event_table` records the first column's
  amplitudes/durations. If you need each column's, see `ampmat` and
  `durmat`.
