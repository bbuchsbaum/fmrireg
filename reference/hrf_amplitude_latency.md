# Response Amplitude and Latency from a Multi-Basis Fit

For every condition modelled with a multi-basis HRF (canonical plus
temporal derivative, `basis = "spmg2"`, or the quadrature pair
`fmrihrf::HRF_SPMGQ`), combine the fitted basis coefficients into a
response amplitude and a latency per voxel, using
`fmrihrf::amplitude_latency()`.

The canonical coefficient alone underestimates a response that peaks
earlier or later than the canonical HRF, so condition differences in it
mix amplitude with timing. The coefficients of a time-shifted response
trace a known curve in coefficient space; each voxel's coefficients are
matched to the nearest direction on that curve, which gives its latency,
and their size along that direction gives the amplitude.

## Usage

``` r
hrf_amplitude_latency(
  object,
  term = NULL,
  shifts = seq(-4, 4, by = 0.05),
  method = c("curve", "boost"),
  min_fit = 0
)
```

## Arguments

- object:

  A fitted `fmri_lm` object.

- term:

  Optional character vector of term tags (as in
  [`fmridesign::design_colmap()`](https://bbuchsbaum.github.io/fmridesign/reference/design_colmap.html))
  to restrict the result to. By default every event term with at least
  two basis functions is used.

- shifts:

  Candidate latencies in seconds (positive = later than the canonical
  response). Estimates are limited to this range.

- method:

  `"curve"` (default) matches coefficients to the shift curve and
  returns amplitude and latency. `"boost"` returns the
  derivative-boosted amplitude of Calhoun et al. (2004) and no latency.

- min_fit:

  Voxels whose coefficients match the shift curve with a cosine below
  this value get `NA` (default 0).

## Value

An object of class `hrf_amplitude_latency`: a list with numeric matrices
`amplitude`, `latency` (seconds) and `fit` (cosine between the
coefficients and the matched curve point), each with one row per voxel
and one column per condition, plus `conditions`, a data frame giving
each column's `term`, `condition` and `basis`. `latency` and `fit` are
`NA` for `method = "boost"`.

## Details

Coefficients are matched in the metric of the basis functions' Gram
matrix, which is appropriate for well-separated events; with heavily
overlapping responses the estimates remain consistent but are not the
most efficient, and a response that is shifted by several seconds is
only approximately spanned by the basis, so part of it can be absorbed
by neighbouring regressors. Latency is relative to the model's own
timing: it includes any offset between the event onsets and the
acquisition times assumed by the sampling frame. Latency is poorly
determined where the amplitude is near zero, so mask the latency map
(for example with a test on the amplitude) before interpreting it. The
summaries are non-linear in the coefficients and have no closed-form
standard error here: take the per-subject maps to a group analysis.

This function needs an fmrihrf that provides `amplitude_latency()`.

## References

Calhoun, V. D., Stevens, M. C., Pearlson, G. D., & Kiehl, K. A. (2004).
fMRI analysis with the general linear model: removal of latency-induced
amplitude bias by incorporation of hemodynamic derivative terms.
*NeuroImage*, 22(1), 252-257.
[doi:10.1016/j.neuroimage.2003.12.029](https://doi.org/10.1016/j.neuroimage.2003.12.029)

Henson, R. N. A., Price, C. J., Rugg, M. D., Turner, R., & Friston, K.
J. (2002). Detecting latency differences in event-related BOLD
responses: application to words versus nonwords and initial versus
repeated face presentations. *NeuroImage*, 15(1), 83-97.
[doi:10.1006/nimg.2001.0940](https://doi.org/10.1006/nimg.2001.0940)

Liao, C. H., Worsley, K. J., Poline, J.-B., Aston, J. A. D., Duncan, G.
H., & Evans, A. C. (2002). Estimating the delay of the fMRI response.
*NeuroImage*, 16(3), 593-606.
[doi:10.1006/nimg.2002.1096](https://doi.org/10.1006/nimg.2002.1096)

## Examples

``` r
# \donttest{
if (exists("amplitude_latency", asNamespace("fmrihrf"))) {
  set.seed(1)
  events <- data.frame(
    onset = seq(4, 280, by = 16),
    condition = factor(rep(c("A", "B"), 9)),
    run = 1
  )
  # Simulate one voxel in which B peaks 1.5 s later than A
  time <- fmrihrf::samples(fmrihrf::sampling_frame(150, TR = 2))
  resp <- function(on) as.numeric(fmrihrf::evaluate(fmrihrf::regressor(on), time))
  y <- 100 + resp(events$onset[events$condition == "A"]) +
    resp(events$onset[events$condition == "B"] + 1.5) + rnorm(150, sd = 0.02)
  dset <- matrix_frame(cbind(y), TR = 2, run_length = 150, event_table = events)
  fit <- fmri_lm(onset ~ hrf(condition, basis = "spmg2"), block = ~ run,
                 dataset = dset)
  al <- hrf_amplitude_latency(fit)
  al$amplitude
  al$latency
}
# }
```
