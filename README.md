# fmrireg

[![R-CMD-check](https://github.com/bbuchsbaum/fmrireg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/fmrireg/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/bbuchsbaum/fmrireg/branch/main/graph/badge.svg)](https://app.codecov.io/gh/bbuchsbaum/fmrireg?branch=main)
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)


# Fit the model
fit <- fmri_lm(emod, dataset = bm$core_data)
```

## Documentation

See comprehensive examples and tutorials in the [vignettes](https://bbuchsbaum.github.io/fmrireg/articles/index.html):

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

PDF report generation for `fmri_lm` fits now lives in the separate
`fmrireport` package. Use `fmrireport::report(fit, ...)` to render analysis
reports.

## Performance Configuration

The internal C++ routines use [RcppParallel](https://rcppcore.github.io/RcppParallel/). You can control the number of threads by setting the R option `fmrireg.num_threads` or the environment variable `FMRIREG_NUM_THREADS` before loading the package. If either is set, `fmrireg` calls `RcppParallel::setThreadOptions()` when it loads.

## Development Status

`fmrireg` is currently in active development. While the core functionality is stable, the API may change as we continue to improve the package. Please [file issues](https://github.com/bbuchsbaum/fmrireg/issues) for bugs or feature requests.

## Citation

If you use `fmrireg` in your research, please cite:

```
Buchsbaum, B. R. (2024). fmrireg: fMRI Analysis in R.
R package version 0.1.0. https://github.com/bbuchsbaum/fmrireg
```
