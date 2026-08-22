# Sketched GLM: Temporal Sketches, Pooled AR, and Landmarks

``` r

library(fmrireg)
library(neuroim2)
library(fmrihrf)
set.seed(1)
```

This vignette introduces the `latent_sketch` GLM path: shared
autoregressive (AR) whitening, temporal sketching, and an optional
landmark-based spatial approximation. These are approximation controls,
not universal accelerators. Their accuracy depends on the design, noise,
and—for landmarks—the spatial smoothness of the coefficient maps.

What you’ll see below

- A quick overview of the approach and when to use it
- Global vs parcel-pooled AR whitening
- Time sketches that preserve GLM geometry at much smaller m ~ 6-8\*p
  - SRHT (Subsampled Randomized Hadamard Transform): a structured
    Johnson-Lindenstrauss embedding built from random signs, a
    Walsh-Hadamard transform, and row sampling. It approximately
    preserves least-squares geometry on the design column space; the
    example below checks the resulting coefficients directly.
  - IHS (Iterative Hessian Sketch): an iterative sketch-and-solve
    method. Each iteration solves a small sketched normal-equation
    system to refine the coefficients.
- Optional Nyström-style spatial extension (solve on L landmarks, extend
  to all voxels)

At a glance

- Evidence shown here: a like-for-like coefficient comparison using both
  correlation and normalized root-mean-square error (NRMSE).
- Public controls: `noise_spec(..., pooling = "global" | "parcel")` for
  AR pooling and
  [`lowrank_control()`](https://bbuchsbaum.github.io/fmrireg/reference/lowrank_control.md)
  for sketch rows, iterations, landmarks, and neighbors.

Approach in one paragraph

- We first prewhiten in time, either globally or with coefficients
  pooled by parcel. A temporal sketch then replaces the time dimension
  `T` with `m` rows. The optional landmark path solves only at `L` voxel
  locations and interpolates coefficients to the remaining voxels with
  sparse heat-kernel weights. This vignette measures approximation
  error; it does not benchmark runtime, because its deliberately small
  dataset is not representative of a production performance regime.

### Set up a synthetic dataset

We start with a small 8 x 8 x 3 volume and 120 time-points so the
example runs quickly. A two-condition event design is convolved with the
canonical HRF to produce the design matrix.

``` r

TR <- 2; Tlen <- 120; dim3 <- c(8L, 8L, 3L)
space4d <- NeuroSpace(c(dim3, Tlen))
maskVol <- LogicalNeuroVol(array(TRUE, dim3), NeuroSpace(dim3))

onsets  <- seq(20, 200, by = 30)
conds   <- factor(rep(c("A", "B"), length.out = length(onsets)))
events_df <- data.frame(onset = onsets, condition = conds, run = 1L)

sframe <- sampling_frame(blocklens = Tlen, TR = TR)
em <- event_model(onset ~ hrf(condition), data = events_df,
                  block = ~ run, sampling_frame = sframe)
X <- as.matrix(design_matrix(em))
p <- ncol(X); V <- prod(dim3)
```

### Simulate data with AR(1) noise

The observed data is `Y = X B + noise` where the noise follows an AR(1)
process with `rho = 0.3`. This gives the sketched solvers declared
temporal dependence to estimate and whiten.

``` r

task_cols <- which(grepl("condition|hrf", colnames(X), ignore.case = TRUE))
B_true <- matrix(0, p, V)

# Landmark interpolation assumes spatial structure. Construct smooth task
# maps in the same voxel order used by the array-packing loop below.
coords_sim <- expand.grid(
  z = seq_len(dim3[3]), y = seq_len(dim3[2]), x = seq_len(dim3[1])
)[, c("x", "y", "z")]
xyz01 <- sweep(as.matrix(coords_sim), 2, 1, "-")
xyz01 <- sweep(xyz01, 2, pmax(dim3 - 1, 1), "/")
smooth_maps <- rbind(
  sin(pi * xyz01[, 1]) * cos(pi * xyz01[, 2]) * (0.5 + xyz01[, 3]),
  cos(pi * xyz01[, 1]) * sin(pi * xyz01[, 2]) * (1.0 - 0.5 * xyz01[, 3])
)
for (j in seq_along(task_cols)) {
  B_true[task_cols[j], ] <- smooth_maps[1L + (j - 1L) %% nrow(smooth_maps), ]
}

Y <- X %*% B_true + ar1_noise(Tlen, V, sd = 0.15)
```

### Pack into an fmri dataset

We reshape the simulated matrix into a 4D NeuroVec for
[`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md).

``` r

vec  <- NeuroVec(arr, space4d)
dset <- fmri_mem_dataset(scans = list(vec), mask = maskVol,
                         TR = TR, event_table = events_df)
```

### Choose the approximation engine

`fmrireg` has two built-in approximation paths with different controls.
`engine = "latent_sketch"` is the sketched GLM engine used in most of
this vignette; it is configured by
[`lowrank_control()`](https://bbuchsbaum.github.io/fmrireg/reference/lowrank_control.md).
`engine = "rrr_gls"` is the reduced-rank-regression GLS engine; it is
configured through `engine_args` and is the path that exposes
block-bootstrap task standard errors. The `rrr_gls` inference surface is
deliberately narrow: standard errors and post-hoc contrasts apply to
event/task coefficients, not baseline or nuisance terms.

### Build parcels

Parcels are used for parcel-pooled AR whitening. The landmark extension
below is coordinate-based and does **not** enforce parcel boundaries.
Here we k-means cluster voxel coordinates solely to demonstrate AR
pooling.

``` r

coords  <- expand.grid(x = seq_len(dim3[1]),
                       y = seq_len(dim3[2]),
                       z = seq_len(dim3[3]))
parcel_fit <- kmeans(
  coords, centers = 40, iter.max = 1000L, nstart = 10L,
  algorithm = "Lloyd"
)
parcels <- ClusteredNeuroVol(maskVol, parcel_fit$cluster)
```

### Fit with Global AR + SRHT sketch

`engine = "latent_sketch"` activates the sketched path.
[`lowrank_control()`](https://bbuchsbaum.github.io/fmrireg/reference/lowrank_control.md)
sets the sketch method and size; here we use SRHT with `m = 8p` rows.
`pooling = "global"` estimates one shared AR(1) coefficient across the
mask.

``` r

low_srht <- lowrank_control(
  parcels     = parcels,
  time_sketch = list(method = "srht", m = min(8L * p, Tlen))
)

fit_srht_global <- fmri_lm(
  onset ~ hrf(condition), block = ~ run,
  dataset    = dset,
  engine     = "latent_sketch",
  lowrank    = low_srht,
  control    = fmri_lm_control(noise = noise_spec("ar1", pooling = "global"))
)
```

### Fit with Parcel AR + SRHT sketch

Switching to `pooling = "parcel"` estimates a separate AR coefficient
per parcel, then shrinks them toward the global mean
(`shrink_c0 = 100`). This stabilizes whitening in small or noisy
parcels.

``` r

fit_srht_group <- fmri_lm(
  onset ~ hrf(condition), block = ~ run,
  dataset    = dset,
  engine     = "latent_sketch",
  lowrank    = low_srht,
  control    = fmri_lm_control(
    noise = noise_spec("ar1", pooling = "parcel", parcels = parcels,
                       shrink_c0 = 100L)
  )
)
```

### Fit with IHS sketch

IHS (Iterative Hessian Sketch) is an alternative to SRHT. It refines its
solution over multiple iterations. Treat `m` and `iters` as accuracy
controls and validate them against an unsketched fit for the target
design.

``` r

low_ihs <- lowrank_control(
  parcels     = parcels,
  time_sketch = list(method = "ihs", m = max(6L * p, p + 4L), iters = 3L)
)

fit_ihs <- fmri_lm(
  onset ~ hrf(condition), block = ~ run,
  dataset    = dset,
  engine     = "latent_sketch",
  lowrank    = low_ihs,
  control    = fmri_lm_control(noise = noise_spec("ar1", pooling = "global"))
)
```

### Compare sketching to an exact OLS baseline

The AR examples above illustrate whitening behavior, but they are not a
like-for-like comparison to exact OLS because the estimators differ. For
a direct sketch-fidelity check, we disable AR, fit the same model with
and without sketching, and compare the task-effect coefficients only. In
this moderate-SNR setting, correlations above 0.9 indicate that the
sketch is preserving the task part of the GLM solution well.

``` r

fit_exact_iid <- fmri_lm(
  onset ~ hrf(condition), block = ~ run,
  dataset = dset
)

fit_srht_iid <- fmri_lm(
  onset ~ hrf(condition), block = ~ run,
  dataset = dset,
  engine = "latent_sketch",
  lowrank = low_srht
)

B_exact_iid <- t(fit_exact_iid$result$betas$data[[1]]$estimate[[1]])

cor_srht_iid_task <- cor(
  as.numeric(B_exact_iid[task_cols, , drop = FALSE]),
  as.numeric(fit_srht_iid$betas_fixed[task_cols, , drop = FALSE])
)
cor_srht_global_task <- cor(
  as.numeric(B_exact_iid[task_cols, , drop = FALSE]),
  as.numeric(fit_srht_global$betas_fixed[task_cols, , drop = FALSE])
)
cor_srht_group_task <- cor(
  as.numeric(B_exact_iid[task_cols, , drop = FALSE]),
  as.numeric(fit_srht_group$betas_fixed[task_cols, , drop = FALSE])
)
cor_ihs_task <- cor(
  as.numeric(B_exact_iid[task_cols, , drop = FALSE]),
  as.numeric(fit_ihs$betas_fixed[task_cols, , drop = FALSE])
)

stopifnot(cor_srht_iid_task > 0.9)

list(
  cor_srht_iid_task = cor_srht_iid_task,
  cor_srht_global_task = cor_srht_global_task,
  cor_srht_group_task = cor_srht_group_task,
  cor_ihs_task = cor_ihs_task
)
#> $cor_srht_iid_task
#> [1] 0.9661459
#> 
#> $cor_srht_global_task
#> [1] 0.9739018
#> 
#> $cor_srht_group_task
#> [1] 0.9645555
#> 
#> $cor_ihs_task
#> [1] 0.9446414
```

The exact iid-versus-sketch comparison is the meaningful fidelity check
here, and it should stay high. The AR-whitened fits usually track the
same task effects reasonably well, but they are solving a different
estimation problem, so their correlations to exact OLS are not expected
to match the iid sketch baseline.

## Reduced-rank GLS with bootstrap standard errors

Use `rrr_gls` when you want a low-rank task fit with task-level standard
errors from either the conditional analytic approximation or a residual
block bootstrap. The bootstrap is useful when short-TR data leave
temporal autocorrelation that ordinary iid standard errors do not
represent well.

``` r

fit_rrr_boot <- fmri_lm(
  onset ~ hrf(condition), block = ~ run,
  dataset = dset,
  control = fmri_lm_control(noise = noise_spec("ar1")),
  engine = "rrr_gls",
  engine_args = list(
    rank_mode = "energy",
    energy_keep = 0.99,
    se_mode = "bootstrap",
    bootstrap_n = 20L,
    bootstrap_block_size = 8L,
    bootstrap_seed = 1L
  )
)
```

For a real analysis, increase `bootstrap_n` substantially; `200` is a
reasonable starting point. The vignette uses a smaller value so the
example remains quick.

``` r

rrr_se <- standard_error(fit_rrr_boot, type = "estimates")

data.frame(
  engine = attr(fit_rrr_boot, "engine"),
  se_mode = fit_rrr_boot$rrr$se_mode,
  rank_used = fit_rrr_boot$rrr$rank_used,
  finite_standard_errors = all(is.finite(as.matrix(rrr_se)))
)
#>    engine   se_mode rank_used finite_standard_errors
#> 1 rrr_gls bootstrap         2                   TRUE
```

## Landmark extension (full-voxel, global AR)

Here we isolate the spatial approximation: both reference and landmark
fits use the same SRHT plan size and global AR specification. The only
difference is that the landmark fit solves at 48 of 192 voxels and
interpolates the rest. The simulated task maps are spatially smooth
because interpolating independent voxel effects would be an ill-posed
target.

``` r

# Choose landmarks (L) and neighbors (k)
L <- 48L; k_nn <- 8L
low_lm <- lowrank_control(
  parcels = parcels,
  landmarks = L,
  k_neighbors = k_nn,
  time_sketch = list(method = "srht", m = min(8L*p, Tlen))
)

fit_lm_srht <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset,
                       engine = "latent_sketch",
                       lowrank = low_lm,
                       control = fmri_lm_control(
                         noise = noise_spec("ar1", pooling = "global")
                       ))

# Compare task coefficients to the like-for-like full-voxel sketch.
task_reference <- as.numeric(
  fit_srht_global$betas_fixed[task_cols, , drop = FALSE]
)
task_landmarks <- as.numeric(
  fit_lm_srht$betas_fixed[task_cols, , drop = FALSE]
)
cor_landmarks <- cor(task_reference, task_landmarks)
rmse_landmarks <- sqrt(mean((task_reference - task_landmarks)^2))
nrmse_landmarks <- rmse_landmarks / stats::sd(task_reference)

stopifnot(cor_landmarks > 0.85, nrmse_landmarks < 0.60)

c(
  correlation = cor_landmarks,
  RMSE = rmse_landmarks,
  NRMSE = nrmse_landmarks
)
#> correlation        RMSE       NRMSE 
#>   0.8955085   0.2233649   0.4624006
```

The correlation checks spatial ordering, while NRMSE detects scale
errors that correlation alone would miss. The thresholds above are an
explicit contract for this smooth synthetic field, not a general
guarantee for real images. Validate both metrics on representative pilot
data before adopting landmarks in an analysis.

``` r

# Compare one task-effect slice on a common scale, then show signed error.
task_index <- task_cols[1]
slice_index <- ceiling(dim3[3] / 2)
reference_volume <- array(
  fit_srht_global$betas_fixed[task_index, ], dim = dim3
)
landmark_volume <- array(
  fit_lm_srht$betas_fixed[task_index, ], dim = dim3
)
reference_slice <- reference_volume[, , slice_index]
landmark_slice <- landmark_volume[, , slice_index]
error_slice <- landmark_slice - reference_slice

effect_limit <- max(abs(c(reference_slice, landmark_slice)))
error_limit <- max(abs(error_slice))
effect_palette <- hcl.colors(64, "Blue-Red 3")
error_palette <- hcl.colors(64, "Blue-Red 3")

old_par <- par(mfrow = c(1, 3), mar = c(3, 3, 3, 1))
image(
  reference_slice, zlim = c(-effect_limit, effect_limit),
  col = effect_palette, asp = 1, axes = FALSE,
  main = "Full-voxel reference"
)
image(
  landmark_slice, zlim = c(-effect_limit, effect_limit),
  col = effect_palette, asp = 1, axes = FALSE,
  main = "Landmark extension"
)
image(
  error_slice, zlim = c(-error_limit, error_limit),
  col = error_palette, asp = 1, axes = FALSE,
  main = "Error: landmark - reference"
)
```

![Three matched spatial slices showing the full-voxel task estimate, its
landmark extension, and signed landmark-minus-reference
error.](sketched-ar_files/figure-html/plot-landmark-error-1.png)

``` r

par(old_par)
```

The common color scale in the first two panels makes attenuation
visible; the third panel localizes signed interpolation error. Smooth
regions are recovered most closely, while gradients and boundaries show
the largest departures. This is why a whole-map correlation alone is
insufficient.

Notes:

- Landmarks are selected by k-means over voxel coordinates and mapped to
  nearest mask voxels.
- The extension uses coordinate k-NN heat-kernel weights with
  row-normalization. It can smooth across anatomical or parcel
  boundaries.
- Residual variance is propagated from landmark residuals through
  squared interpolation weights; this is an approximation, not a
  calibrated voxelwise-inference guarantee.

### Choosing landmarks: practical guidance

- How many (L)?
  - There is no universal ratio. Start conservatively, validate
    correlation and NRMSE on representative regions, then reduce `L`
    only while those predeclared tolerances continue to hold.
- Where to place them?
  - k-means centers on voxel coordinates (as above) are a reproducible
    default when the RNG is seeded.
  - Farthest-point (a.k.a. max-min) sampling gives more uniform coverage
    than k-means; prefer it if you have irregular masks.
  - Boundaries: include landmarks near anatomical borders and inspect
    those regions explicitly; the current interpolator does not honor
    boundaries.
- Neighborhood size (k_neighbors) and bandwidth (h)
  - Use k_neighbors in \[8, 32\]; 16 is a good default. Larger k smooths
    more but risks crossing boundaries.
  - The default bandwidth uses the median distance to the k-th neighbor;
    adjust up to smooth more, down to sharpen.
- Diagnostics
  - Correlate landmark-extended betas with a smaller exact or sketched
    full-voxel solve on an ROI.
  - Visualize a few beta maps near sharp boundaries (sulci, ROI borders)
    to check that edges are preserved.
  - Do not use interpolated standard errors or test statistics for
    strict voxelwise inference without a separate calibration study.

## What the comparisons establish

- The iid SRHT comparison isolates temporal-sketch error from AR-model
  changes.
- The landmark comparison isolates spatial interpolation while holding
  the temporal sketch and AR specification fixed.
- Results on smooth synthetic fields do not establish accuracy at sharp
  boundaries or for spatially irregular effects.
- No runtime claim is made here. Benchmark the exact workload, including
  data access and preprocessing, before choosing an approximation for
  speed.

## When to use which knob

- Reduce `T` with a temporal sketch only after validating
  task-coefficient error against the unsketched estimator.
- Reduce `V` with landmarks only when coefficient fields are spatially
  smooth enough to satisfy a declared map-error tolerance.
- Omit landmarks when voxelwise detail or sharp boundaries are central
  to the scientific question.

## Pros and cons

Pros

- Temporal and spatial approximation controls are explicit in
  [`lowrank_control()`](https://bbuchsbaum.github.io/fmrireg/reference/lowrank_control.md).
- Global and parcel-pooled AR models are expressed through the same
  public
  [`noise_spec()`](https://bbuchsbaum.github.io/fmrireg/reference/noise_spec.md)
  interface.
- Correlation and scale-sensitive error can be validated separately.

Cons

- Landmark extension introduces spatial shrinkage that is data-dependent
  and can cross boundaries; if you need strict voxelwise inference, omit
  landmarks.
- Sketching slightly perturbs degrees of freedom (uses m - p); choose m
  conservatively if p is large or X is ill-conditioned.

## Practical defaults

- Time sketch: `list(method = "srht", m = min(8*p, T))` or IHS with
  `m ~ 6*p`, `iters = 2-3`.
- AR order: 1 is a solid default; use 2 if residuals show lag-2
  structure.
- Parcel AR: use
  `noise_spec("ar1", pooling = "parcel", parcels = parcels)`; tune
  shrinkage for the target data rather than assuming one universal
  value.
- Landmarks: choose `L` and `k_neighbors` by a declared correlation and
  NRMSE contract on representative data.

## Reproducibility tips

- Set a seed before building SRHT plans or k-means landmarks
  (`set.seed(...)`).
- Keep m and iters fixed across runs to compare apples-to-apples.
- Validate on a small ROI first: check correlations vs exact, and
  histogram of t-maps in null tissue.

## Next

- [`vignette("a_09_linear_model", package = "fmrireg")`](https://bbuchsbaum.github.io/fmrireg/articles/a_09_linear_model.md)
  — fMRI Linear Model (GLM)
- [`vignette("group_analysis", package = "fmrireg")`](https://bbuchsbaum.github.io/fmrireg/articles/group_analysis.md)
  — Group Analysis
