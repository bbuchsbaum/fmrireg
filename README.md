# fmrireg

[![R-CMD-check](https://github.com/bbuchsbaum/fmrireg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/fmrireg/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/bbuchsbaum/fmrireg/branch/main/graph/badge.svg)](https://app.codecov.io/gh/bbuchsbaum/fmrireg?branch=main)
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

experimental regressor construction, and univariate fMRI regression using a
formula-based interface inspired by `lm()` and `glm()`.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("bbuchsbaum/fmrireg")
```

## Quick start

```r
library(fmrireg)

# Define the temporal structure: 2 runs, 100 scans each, TR = 2s
sframe <- sampling_frame(blocklens = c(100, 100), TR = 2)

# Build an event model from an experimental design table
emod <- event_model(onset ~ hrf(condition, basis = "spmg1"),
                    data = design_table,
                    block = ~ run,
                    sampling_frame = sframe)

# Add a baseline model with polynomial drift
bmod <- baseline_model(basis = "bs", degree = 5, sframe = sframe)

# Combine into a full fMRI model and fit
fmod <- fmri_model(emod, bmod)
fit  <- fmri_lm(fmod, dataset = dset)
```

## Key features

- **HRF library** -- SPM canonical, gamma, Gaussian, B-spline, and custom
  basis sets, with decorators for lag, block, and normalization.
- **Formula interface** -- Specify event models with `onset ~ hrf(...)` syntax;
  supports categorical events, continuous modulators, and multi-basis expansions.
- **Contrast system** -- Flexible contrasts via formulas, including pairwise,
  polynomial, and F-contrasts.
- **Robust estimation** -- OLS, iteratively reweighted least squares (IWLS),
  and sandwich variance estimators.
- **AR correction** -- Autoregressive noise modeling via the
  [fmriAR](https://github.com/bbuchsbaum/fmriAR) package.
- **Performance** -- C++ solvers (Rcpp/RcppArmadillo) with optional
  multithreading via RcppParallel.

## Documentation

Full documentation and tutorials are available at
<https://bbuchsbaum.github.io/fmrireg/>. Vignettes include:

- [Package Overview](https://bbuchsbaum.github.io/fmrireg/articles/Overview.html)
- [Statistical Contrasts](https://bbuchsbaum.github.io/fmrireg/articles/a_05_contrasts.html)
- [Simulation](https://bbuchsbaum.github.io/fmrireg/articles/a_08_simulation.html)
- [Linear Modeling](https://bbuchsbaum.github.io/fmrireg/articles/a_09_linear_model.html)
- [Dataset Management](https://bbuchsbaum.github.io/fmrireg/articles/a_10_dataset.html)
- [Group Analysis](https://bbuchsbaum.github.io/fmrireg/articles/group_analysis.html)
- [Functional Connectivity](https://bbuchsbaum.github.io/fmrireg/articles/functional_connectivity.html)
- [Sketched GLM](https://bbuchsbaum.github.io/fmrireg/articles/sketched-ar.html)
- [Benchmark Datasets](https://bbuchsbaum.github.io/fmrireg/articles/benchmark_datasets.html)

## Reporting

PDF report generation for `fmri_lm` fits is provided by the separate
[fmrireport](https://github.com/bbuchsbaum/fmrireport) package:

```r
fmrireport::report(fit, output_dir = "results")
```

## Performance configuration

The internal C++ routines use
[RcppParallel](https://rcppcore.github.io/RcppParallel/). Control the thread
count with:

```r
options(fmrireg.num_threads = 4)
```

or set the environment variable `FMRIREG_NUM_THREADS` before loading the
package.

## Citation

If you use fmrireg in your research, please cite:

```
Buchsbaum, B. R. (2025). fmrireg: Regression Analysis of Functional
Magnetic Resonance Imaging Data. R package version 0.1.2.
https://github.com/bbuchsbaum/fmrireg
```

## License

GPL (>= 2)
