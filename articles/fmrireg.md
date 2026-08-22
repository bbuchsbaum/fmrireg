# Overview and Workflow

## Introduction

If you work with fMRI data, you already know the core analysis is
multiple regression applied voxel-by-voxel — but specifying the design
matrix, convolving events with hemodynamic response functions, and
managing multi-run structures can be tedious and error-prone.

The `fmrireg` package gives you a concise, formula-based interface for
these tasks, inspired by [`lm()`](https://rdrr.io/r/stats/lm.html) and
[`glm()`](https://rdrr.io/r/stats/glm.html). You describe your
experimental design and HRF choices in a formula, and fmrireg handles
convolution, baseline modeling, contrast computation, and model fitting.

This vignette starts with a complete, small fit and then unpacks the
objects behind it. You begin with a time-by-voxel matrix and an event
table; you finish with one coefficient estimate per condition and voxel.

## A complete fit first

For a compact example, we create two voxels with different known
responses to faces and scenes. In your own analysis, `bold` would be the
matrix extracted from preprocessed images or regions of interest.

``` r

quick_events <- data.frame(
  onset = seq(10, 80, by = 10),
  condition = factor(rep(c("face", "scene"), 4)),
  run = 1L
)
quick_frame <- sampling_frame(blocklens = 80, TR = 2)
quick_event_model <- event_model(
  onset ~ hrf(condition), block = ~ run,
  data = quick_events, sampling_frame = quick_frame
)
quick_design <- as.matrix(design_matrix(quick_event_model))
quick_truth <- matrix(
  c(1.0, 0.6, 0.7, 1.1), nrow = ncol(quick_design),
  dimnames = list(colnames(quick_design), c("voxel1", "voxel2"))
)
set.seed(2026)
bold <- quick_design %*% quick_truth + matrix(rnorm(80 * 2, sd = 0.2), 80, 2)
quick_dataset <- matrix_dataset(
  bold, TR = 2, run_length = 80, event_table = quick_events
)
```

``` r

face_vs_scene <- pair_contrast(
  ~ condition == "face", ~ condition == "scene",
  name = "face_vs_scene"
)
fit <- fmri_lm(
  onset ~ hrf(condition, contrasts = face_vs_scene), block = ~ run,
  dataset = quick_dataset
)

coefficient_table <- data.frame(
  condition = levels(quick_events$condition),
  truth_voxel1 = quick_truth[, 1],
  estimate_voxel1 = as.numeric(coef(fit)[1, ]),
  truth_voxel2 = quick_truth[, 2],
  estimate_voxel2 = as.numeric(coef(fit)[2, ])
)
knitr::kable(coefficient_table, digits = 2,
             row.names = FALSE,
             caption = "Known and recovered condition effects")
```

| condition | truth_voxel1 | estimate_voxel1 | truth_voxel2 | estimate_voxel2 |
|:----------|-------------:|----------------:|-------------:|----------------:|
| face      |          1.0 |            0.95 |          0.7 |            0.51 |
| scene     |          0.6 |            0.76 |          1.1 |            1.14 |

Known and recovered condition effects {.table}

``` r


quick_contrast_numeric <- data.frame(
  voxel = colnames(quick_truth),
  estimate = as.numeric(coef(fit, type = "contrasts")$face_vs_scene),
  std_error = as.numeric(standard_error(fit, type = "contrasts")$face_vs_scene),
  statistic = as.numeric(stats(fit, type = "contrasts")$face_vs_scene),
  p_value = as.numeric(p_values(fit, type = "contrasts")$face_vs_scene)
)
quick_contrast <- transform(
  quick_contrast_numeric,
  p_value = format.pval(p_value, digits = 3, eps = 1e-4)
)
knitr::kable(quick_contrast, digits = 3,
             row.names = FALSE,
             caption = "Face-minus-scene contrast with uncertainty")
```

| voxel  | estimate | std_error | statistic | p_value |
|:-------|---------:|----------:|----------:|:--------|
| voxel1 |    0.440 |     0.053 |     8.319 | \<1e-04 |
| voxel2 |   -0.379 |     0.062 |    -6.116 | \<1e-04 |

Face-minus-scene contrast with uncertainty {.table}

The first table separates the known effect sizes from their estimates.
The second table asks a scientific question—face minus scene—and reports
its estimate, standard error, test statistic, and p-value for each
voxel. Here the first voxel favors faces while the second favors scenes,
recovering the pattern used to generate the data. The next sections
unpack the dataset, event, baseline, and contrast objects used for
larger analyses.

## Key Components and Workflow

A typical fMRI analysis workflow with `fmrireg` involves these main
stages:

1.  **Data Representation:** Defining how your fMRI data and
    experimental design are organized.
2.  **Model Building:** Specifying the regression model using formulas,
    including:
    - Defining Hemodynamic Response Functions (HRFs).
    - Creating event-related regressors (`event_model`).
    - Modeling baseline components like drift and nuisance factors
      (`baseline_model`).
    - Defining contrasts of interest (`contrast`).
3.  **Model Estimation:** Fitting the specified model to the data (e.g.,
    using `fmri_lm`).
4.  **Advanced Topics:** Exploring single-trial estimation, simulation,
    group analysis, or accelerated engines.

### 1. Data Representation

Before modeling, you need to represent your fMRI data and its structure.
`fmrireg` uses objects like:

- `fmri_dataset` / `fmri_mem_dataset`: Encapsulates the 4D fMRI data
  (potentially file-backed or in-memory), mask information, TR, run
  structure, and the experimental design table.
- `sampling_frame`: Defines the temporal structure – the number of scans
  per run (`blocklens`) and the repetition time (`TR`).

``` r

# Example: one run, 120 scans, TR=2s
sframe <- sampling_frame(blocklens = 120, TR = 2)
sframe
#> Sampling frame
#> - Blocks: 1 
#> - Scans: 120 (per block: 120 )
#> - TR: 2 s
#> - Duration: 239 s
```

### 2. Model Building

#### a) Hemodynamic Response Function (HRF)

The HRF models the BOLD signal response to a brief neural event.
`fmrireg` provides various pre-defined HRFs (e.g., `HRF_SPMG1`,
`HRF_GAUSSIAN`) and functions to create or modify them (`gen_hrf`,
`hrf_blocked`, `hrf_lagged`).

``` r

# Plot the SPM canonical HRF
time <- seq(0, 24, by = 0.2)
plot(time, HRF_SPMG1(time), type = 'l', xlab = "Time (s)", ylab = "BOLD Response",
     main = "SPM Canonical HRF (HRF_SPMG1)")
```

![Line plot of the SPM Canonical HRF (BOLD response over time in
seconds)](fmrireg_files/figure-html/hrf_example-1.png)

- **Related Topic:** Hemodynamic Response Functions

#### b) Event Model (`event_model`)

This defines the task-related part of your GLM. You specify how
experimental events (onsets, conditions, durations) map to regressors
using a formula syntax. The
[`hrf()`](https://bbuchsbaum.github.io/fmridesign/reference/hrf.html) or
`afni_hrf()` functions within the formula link event variables to
specific HRF shapes.

``` r

# A simple design table with two stimulus types
design_table <- data.frame(
  onset = c(10, 40, 70, 100, 130, 160, 190, 220),
  stim_type = factor(rep(c("face", "scene"), 4)),
  run_id = rep(1, 8)
)
```

``` r

# Model the effect of stim_type, convolved with the SPM canonical HRF
event_mod <- event_model(onset ~ hrf(stim_type),
                         data = design_table,
                         block = ~ run_id,
                         sampling_frame = sframe)
event_mod
```

- **Related Topic:** Event Models

#### c) Baseline Model (`baseline_model`)

This models non-task-related components like scanner drift, run
intercepts, and other nuisance factors (e.g., motion parameters).

``` r

# Example: Cubic B-spline drift model + run intercepts
# (Requires 'sframe' defined earlier)
base_mod <- baseline_model(basis = "bs", degree = 3, sframe = sframe)
print(base_mod)
#> ================================================
#>            Baseline Model                       
#> ================================================
#>   Drift Components                           
#>     * Name: baseline_bs_3                
#>     * Basis type: bs                     
#>     * Degree: 3                          
#>     * Drift columns: 3                   
#> 
#>   Additional Components                    
#>     * Constant columns: 1                
#>     * Nuisance columns: 0                
#> 
#>   Model Summary                            
#>     * Total columns: 4                   
#> 
#>   Design Matrix Preview                    
#>      0.000  0.000  0.000  1.000          
#>      0.025  0.000  0.000  1.000          
#>      0.049  0.001  0.000  1.000          
#> ================================================
```

- **Related Topic:** Baseline Models

#### d) Contrasts (`pair_contrast`, `unit_contrast`, etc.)

Contrasts define specific comparisons between condition estimates.

``` r

# Example: Contrast comparing 'face' vs 'scene' levels of 'stim_type'
face_vs_scene <- pair_contrast(~ stim_type == "face", ~ stim_type == "scene", name = "FvS")
print(face_vs_scene)
#> contrast: FvS 
#>  A:  ~stim_type == "face" 
#>  B:  ~stim_type == "scene"
```

- **See Vignette:**
  [`vignette("a_05_contrasts", package = "fmrireg")`](https://bbuchsbaum.github.io/fmrireg/articles/a_05_contrasts.md)

### 3. Model Estimation (`fmri_model`, `fmri_lm`)

You combine the event and baseline models into a full `fmri_model` and
then fit it to your `fmri_dataset` using estimation functions like
`fmri_lm` (for standard GLM) or `estimate_betas` (for single-trial
betas).

``` r

# Create three example voxel series over the same one-run design
set.seed(42)
overview_data <- matrix_dataset(
  matrix(rnorm(120 * 3), 120, 3),
  TR = 2, run_length = 120, event_table = design_table
)

# Combine the design with its dataset, then fit it
full_model <- fmri_model(event_mod, base_mod, dataset = overview_data)
overview_fit <- fmri_lm(full_model)
overview_coef_matrix <- coef(overview_fit)
overview_coef_table <- data.frame(
  condition = levels(design_table$stim_type),
  voxel_1 = overview_coef_matrix[, 1],
  voxel_2 = overview_coef_matrix[, 2],
  voxel_3 = overview_coef_matrix[, 3],
  check.names = FALSE
)
knitr::kable(
  overview_coef_table, digits = 2, row.names = FALSE,
  caption = "Condition coefficients from the object-building route"
)
```

| condition | voxel_1 | voxel_2 | voxel_3 |
|:----------|--------:|--------:|--------:|
| face      |   -0.09 |   -0.14 |   -0.11 |
| scene     |    0.05 |    0.06 |    0.10 |

Condition coefficients from the object-building route {.table}

This second fit uses the object-building route rather than the compact
formula route shown first. See the linear-model vignette for AR
prewhitening, robust row weighting, uncertainty estimates, and model
diagnostics.

- **See Vignettes:**
  - [`vignette("a_09_linear_model", package = "fmrireg")`](https://bbuchsbaum.github.io/fmrireg/articles/a_09_linear_model.md)
  - Linear Models and Trialwise Analysis

### 4. Advanced Topics

`fmrireg` also supports:

- **Simulation:** Creating artificial fMRI data
  (`simulate_fmri_matrix`).
  - [`vignette("a_08_simulation", package = "fmrireg")`](https://bbuchsbaum.github.io/fmrireg/articles/a_08_simulation.md)
- **Low-level Regressors:** Manually creating and manipulating
  regressors (`regressor`).

#### Acceleration: Sketched / Nyström Analysis

For large problems, you can accelerate voxelwise GLM estimation with: -
SRHT/IHS temporal sketching to reduce timepoints while preserving LS
geometry. - Nyström/landmark spatial extension to solve on a subset of
voxels and extend with parcel‑aware kernels.

See
[`vignette("sketched-ar", package = "fmrireg")`](https://bbuchsbaum.github.io/fmrireg/articles/sketched-ar.md)
for guidance, tradeoffs, and practical defaults.

#### Group-Level Meta-Analysis

Combine subject‑level maps or ROI summaries using fixed/random‑effects
meta‑regression: - Inputs: HDF5/NIfTI (beta+SE or t+df) or ROI CSV. -
Random‑effects: Paule–Mandel/DerSimonian–Laird/REML alias; optional
robust Huber. - t‑only combine: Stouffer/Fisher/Lancaster (supports
equal/IVW/custom weights where applicable). - Exact contrasts: at
fit‑time (contrasts=…) or post‑hoc via `return_cov = "tri"` +
[`contrast()`](https://bbuchsbaum.github.io/fmrireg/reference/contrast.md). -
Multiple comparisons: BH FDR or spatial FDR (block/parcel aware).

See
[`vignette("group_analysis", package = "fmrireg")`](https://bbuchsbaum.github.io/fmrireg/articles/group_analysis.md)
for end‑to‑end examples.

## Getting Started

This overview provides a map of the `fmrireg` package. For detailed
explanations and examples, please refer to the specific vignettes linked
above.

## Next Steps

- [`vignette("a_10_dataset", package = "fmrireg")`](https://bbuchsbaum.github.io/fmrireg/articles/a_10_dataset.md)
  — Dataset management
- [`vignette("a_08_simulation", package = "fmrireg")`](https://bbuchsbaum.github.io/fmrireg/articles/a_08_simulation.md)
  — Build controlled examples
- [`vignette("a_09_linear_model", package = "fmrireg")`](https://bbuchsbaum.github.io/fmrireg/articles/a_09_linear_model.md)
  — Fit and diagnose a GLM
- [`vignette("a_05_contrasts", package = "fmrireg")`](https://bbuchsbaum.github.io/fmrireg/articles/a_05_contrasts.md)
  — Test planned hypotheses
