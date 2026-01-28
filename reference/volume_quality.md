# Volume Quality Metrics and Temporal Weighting

Functions for computing volume-level quality metrics (DVARS) and
converting them to weights for weighted least squares fitting. This
implements "soft scrubbing" - downweighting bad volumes rather than hard
censoring them.

## Conceptual Overview

Traditional fMRI artifact removal uses hard censoring ("scrubbing"),
where volumes exceeding a threshold are completely removed. This loses
temporal information and can create discontinuities.

Soft scrubbing instead assigns each volume a weight between 0 and 1
based on its quality. High-quality volumes receive full weight;
artifacts receive reduced weight. This preserves temporal continuity
while still downweighting problematic data.

## DVARS as Quality Metric

DVARS (Derivative of VARiance across voxels) measures the root mean
square of the temporal derivative across all voxels: \$\$DVARS_t =
\sqrt{\frac{1}{V} \sum_v (Y\_{t,v} - Y\_{t-1,v})^2}\$\$

Interpretation:

- Low DVARS: signal changed smoothly from previous volume (good)

- High DVARS: signal changed rapidly (possible artifact or motion)

- Normalized DVARS: values ~1 are typical; \>1.5 suggests artifacts

## Weighting Methods

Three methods convert DVARS to weights, offering different trade-offs:

- inverse_squared:

  Formula: `w = 1 / (1 + dvars^2)`

      Properties: Smooth, continuous decay. A volume with DVARS=1 gets
      weight 0.5. Most conservative choice - provides gentle downweighting
      even for moderately elevated DVARS.

      Use when: You want smooth, gradual artifact handling without
      sharp transitions.

- soft_threshold:

  Formula: Sigmoid decay above threshold

      Properties: Volumes below threshold get full weight; above
      threshold, weights decay smoothly. Steepness parameter controls
      how rapidly weights drop.

      Use when: You have a clear idea of what "acceptable" DVARS looks
      like (e.g., 1.5x median) and want to preserve good volumes fully.

- tukey:

  Formula: Tukey bisquare `(1 - u^2)^2` for `|u| <= 1`

      Properties: Complete downweighting for extreme values. Volumes
      beyond 2x threshold get zero weight. Most aggressive choice.

      Use when: You have clear artifacts that should be fully excluded,
      similar to hard scrubbing but with smooth transitions.

## When to Use Volume Weighting vs Robust Fitting

Both approaches handle artifacts, but through different mechanisms:

**Volume weighting** downweights entire timepoints uniformly across all
voxels. Best when artifacts affect the whole brain (head motion, scanner
spikes).

**Robust fitting** (Huber/Tukey bisquare) downweights outlier residuals
voxel-by-voxel. Best when artifacts are spatially localized or when
temporal structure matters.

**Combined approach**: Use volume weighting for global artifacts plus
robust fitting for residual voxel-level outliers. Set
`robust = TRUE, volume_weights = TRUE` in
[`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md).

## Typical Usage

For use within
[`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md),
see the convenience parameter `volume_weights` or the more detailed
`volume_weights_options`.
