# Fit a Linear Regression Model for fMRI Data Analysis

`fmri_lm` is a generic for fitting fMRI regression models. The default
interface accepts a model formula and dataset. An alternative method can
be used with a preconstructed `fmri_model` object that already contains
the design and data.

## Usage

``` r
fmri_lm(formula, ...)

# S3 method for class 'formula'
fmri_lm(
  formula,
  block,
  baseline_model = NULL,
  dataset,
  durations = 0,
  drop_empty = TRUE,
  robust = FALSE,
  robust_options = NULL,
  ar_options = NULL,
  volume_weights_options = NULL,
  soft_subspace_options = NULL,
  strategy = c("runwise", "chunkwise"),
  nchunks = 10,
  use_fast_path = TRUE,
  progress = FALSE,
  ar_voxelwise = FALSE,
  parallel_voxels = FALSE,
  cor_struct = NULL,
  cor_iter = NULL,
  cor_global = NULL,
  ar1_exact_first = NULL,
  ar_p = NULL,
  robust_psi = NULL,
  robust_max_iter = NULL,
  robust_scale_scope = NULL,
  volume_weights = NULL,
  nuisance_projection = NULL,
  parallel_chunks = FALSE,
  ...
)

# S3 method for class 'fmri_model'
fmri_lm(
  formula,
  dataset = NULL,
  robust = FALSE,
  robust_options = NULL,
  ar_options = NULL,
  volume_weights_options = NULL,
  soft_subspace_options = NULL,
  strategy = c("runwise", "chunkwise"),
  nchunks = 10,
  use_fast_path = TRUE,
  progress = FALSE,
  ar_voxelwise = FALSE,
  parallel_voxels = FALSE,
  cor_struct = NULL,
  cor_iter = NULL,
  cor_global = NULL,
  ar1_exact_first = NULL,
  ar_p = NULL,
  robust_psi = NULL,
  robust_max_iter = NULL,
  robust_scale_scope = NULL,
  volume_weights = NULL,
  nuisance_projection = NULL,
  parallel_chunks = FALSE,
  ...
)
```

## Arguments

- formula:

  A model formula describing the event structure or an `fmri_model`
  object.

- ...:

  Additional method arguments. Recognized engine arguments include
  `engine`, `engine_args`, and `lowrank`; see Details.

- block:

  The model formula for block structure.

- baseline_model:

  (Optional) A `baseline_model` object. Default is `NULL`.

- dataset:

  An `fmri_dataset` object containing the time-series data.

- durations:

  A vector of event durations. Default is `0`.

- drop_empty:

  Logical. Whether to remove factor levels with zero size. Default is
  `TRUE`.

- robust:

  Logical or character. Either `FALSE` (no robust fitting), `TRUE` (use
  Huber), or one of `"huber"` or `"bisquare"`. Default is `FALSE`.

- robust_options:

  List of robust fitting options. See Details.

- ar_options:

  List of autoregressive modeling options. See Details.

- volume_weights_options:

  List of volume weighting options. See
  [`fmri_lm_control`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm_control.md).

- soft_subspace_options:

  List of soft subspace projection options. See
  [`fmri_lm_control`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm_control.md).

- strategy:

  The data splitting strategy, either `"runwise"` or `"chunkwise"`.
  Default is `"runwise"`.

- nchunks:

  Number of data chunks when strategy is `"chunkwise"`. This controls
  memory partitioning; chunks are processed sequentially unless
  `parallel_chunks = TRUE`. Default is `10`.

- use_fast_path:

  Logical. If `TRUE` (the default), use the fast matrix engine, which
  supports OLS, AR whitening, robust, and preprocessing. `FALSE` selects
  the formula/lm reference engine (no robust or full preprocessing
  support); it is retained mainly as a parity oracle.

- progress:

  Logical. Whether to display a progress bar during model fitting.
  Default is `FALSE`.

- ar_voxelwise:

  Logical. Estimate AR parameters voxel-wise (overrides
  `ar_options$voxelwise`).

- parallel_voxels:

  Logical. Parallelize across voxels where supported; this does not
  control chunkwise execution.

- cor_struct:

  Character. Shorthand for `ar_options$struct` (e.g., "ar1", "ar2",
  "arp").

- cor_iter:

  Integer. Shorthand for `ar_options$iter_gls`.

- cor_global:

  Logical. Shorthand for `ar_options$global`.

- ar1_exact_first:

  Logical. Shorthand for `ar_options$exact_first`.

- ar_p:

  Integer. Shorthand for `ar_options$p`.

- robust_psi:

  Character. Shorthand for `robust_options$type` (e.g., "huber",
  "bisquare").

- robust_max_iter:

  Integer. Shorthand for `robust_options$max_iter`.

- robust_scale_scope:

  Character. Shorthand for `robust_options$scale_scope` ("run",
  "global", or "voxel").

- volume_weights:

  Logical or character. Simple interface for volume weighting:

  - `TRUE`: Enable with default method ("inverse_squared")

  - `"inverse_squared"`, `"soft_threshold"`, `"tukey"`: Enable with
    specific method

  - `FALSE` or `NULL`: Disable (default)

  For fine-grained control, use `volume_weights_options` instead.

- nuisance_projection:

  Matrix, character path, or NULL. Simple interface for soft subspace
  projection:

  - Matrix: Use as nuisance timeseries directly

  - Character: Path to NIfTI mask for WM/CSF voxels

  - `NULL`: Disable (default)

  For fine-grained control (lambda selection, warnings), use
  `soft_subspace_options`.

- parallel_chunks:

  Logical. For `strategy = "chunkwise"`, process chunks with
  [`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html)
  using the active `future` plan. Default is `FALSE`.

## Value

An object of class `fmri_lm`.

A fitted linear regression model for fMRI data analysis.

## Details

`robust_options` may contain:

- `type`: Character or logical. Type of robust fitting (`FALSE`,
  `"huber"`, `"bisquare"`)

- `k_huber`: Numeric. Tuning constant for Huber's psi (default: 1.345)

- `c_tukey`: Numeric. Tuning constant for Tukey's bisquare psi (default:
  4.685)

- `max_iter`: Integer. Maximum IRLS iterations (default: 2)

- `scale_scope`: Character. Scope for scale estimation (`"run"` or
  `"global"`)

- `reestimate_phi`: Logical. Whether to re-estimate AR parameters after
  robust fitting

`ar_options` may contain:

- `struct`: Character. Correlation structure (`"iid"`, `"ar1"`, `"ar2"`,
  `"arp"`)

- `p`: Integer. AR order when `struct = "arp"`

- `iter_gls`: Integer. Number of GLS iterations (default: 1)

- `global`: Logical. Use global AR coefficients (default: FALSE)

- `voxelwise`: Logical. Estimate AR parameters voxel-wise (default:
  FALSE)

- `exact_first`: Logical. Apply exact AR(1) scaling to first sample
  (default: FALSE)

Built-in fast engines are selected through `...`:

- `engine = "latent_sketch"` uses the sketched GLM path. Configure it
  with `lowrank = lowrank_control(...)`. The alias `engine = "sketch"`
  is normalized to `"latent_sketch"`.

- `engine = "rrr_gls"` uses the reduced-rank-regression GLS path.
  Configure it with `engine_args = list(...)`. This engine supports
  shared AR whitening from `ar_options` and inference for event/task
  parameters only; nuisance or baseline terms may be estimated, but
  standard-error and contrast inference is restricted to event/task
  coefficients.

For `engine = "rrr_gls"`, `engine_args` may contain:

- `rank`: Positive integer rank. If `NULL` and `rank_mode = "fixed"`,
  the engine uses the full available task rank.

- `rank_mode`: One of `"fixed"`, `"energy"`, or `"rss_budget"` (default:
  `"fixed"`).

- `energy_keep`: Fraction of task singular-value energy to retain when
  `rank_mode = "energy"` (default: `0.99`).

- `rss_budget`: Non-negative residual-sum-of-squares budget used when
  `rank_mode = "rss_budget"`.

- `se_mode`: Standard-error mode, either `"conditional"` or
  `"bootstrap"` (default: `"conditional"`). Bootstrap standard errors
  resample task-subspace residuals in blocks.

- `bootstrap_n`: Number of bootstrap resamples when
  `se_mode = "bootstrap"` (default: `200`).

- `bootstrap_block_size`: Positive integer block size for bootstrap
  resampling. If `NULL`, a block size of 1 is used.

- `bootstrap_seed`: Optional integer seed for reproducible bootstrap
  resampling.

- `contrast_policy`: How to handle post-hoc contrasts that include
  non-event coefficients: `"warn_drop"`, `"drop"`, or `"error"`
  (default: `"warn_drop"`).

The aliases `energy` and `nboot` are accepted as legacy shorthand for
`energy_keep` and `bootstrap_n`, respectively. Prefer the canonical
names in new code.

## See also

[`fmri_dataset`](https://bbuchsbaum.github.io/fmridataset/reference/fmri_dataset.html),
[`fmri_lm_fit`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm_fit.md),
[`fmri_lm_control`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm_control.md)

## Examples

``` r

facedes <- subset(read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), 
header=TRUE), face_gen != "n/a")
facedes$face_gen <- droplevels(factor(facedes$face_gen))
sframe <- sampling_frame(rep(430/2,6), TR=2)
ev <- event_model(onset ~ hrf(face_gen, basis="gaussian"), data=facedes, 
block= ~ run, sampling_frame=sframe)
globonsets <- fmrihrf::global_onsets(sframe, facedes$onset, facedes$run)
reg1_signal <- regressor(globonsets[facedes$face_gen == "male"], hrf=fmrihrf::HRF_GAUSSIAN)
reg2_signal <- regressor(globonsets[facedes$face_gen == "female"], hrf=fmrihrf::HRF_GAUSSIAN)
time <- samples(sframe, global=TRUE)
y1 <- fmrihrf::evaluate(reg1_signal, time)*1.5
y2 <- fmrihrf::evaluate(reg2_signal, time)*3.0
y <- y1+y2
ys1 <- y + rnorm(length(y), sd=.02)
ys2 <- y + rnorm(length(y), sd=.02)

h <<- gen_hrf(fmrihrf::hrf_bspline, N=7, span=25)
dset <- matrix_dataset(cbind(ys1,ys2), TR=2, 
                       run_length=fmrihrf::blocklens(sframe), 
                       event_table=facedes)
flm <- fmri_lm(onset ~ hrf(face_gen, 
                           basis=gen_hrf(fmrihrf::hrf_bspline, N=7, span=25)), 
               block = ~ run, 
               strategy="chunkwise", nchunks=1, dataset=dset)

if (FALSE) { # \dontrun{
fit_rrr <- fmri_lm(
  onset ~ hrf(condition), block = ~ run, dataset = dset,
  ar_options = list(struct = "ar1"),
  engine = "rrr_gls",
  engine_args = list(
    rank_mode = "energy",
    energy_keep = 0.99,
    se_mode = "bootstrap",
    bootstrap_n = 200L,
    bootstrap_block_size = 24L,
    bootstrap_seed = 1L
  )
)
standard_error(fit_rrr, type = "estimates")
} # }
```
