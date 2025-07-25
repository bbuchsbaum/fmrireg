# fmrireg

[![R-CMD-check](https://github.com/bbuchsbaum/fmrireg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/fmrireg/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/bbuchsbaum/fmrireg/branch/main/graph/badge.svg)](https://app.codecov.io/gh/bbuchsbaum/fmrireg?branch=main)
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

- **Integration**: Works with AFNI and other neuroimaging software
- **Visualization**: Rich plotting capabilities for designs and results

## Installation

You can install the development version of fmrireg from [GitHub](https://github.com/) with:

```r
# install.packages("devtools")
devtools::install_github("bbuchsbaum/fmrireg")
```

## Quick Start

```r
library(fmrireg)

# Load a benchmark dataset
bm <- load_benchmark_dataset("BM_Canonical_HighSNR")

# Create an event model
events <- bm$core_data$event_table
sframe <- sampling_frame(blocklens = bm$core_data$run_length, TR = bm$core_data$TR)
emod <- event_model(onset ~ hrf(condition), data = events, 
                   block = ~run, sampling_frame = sframe)

# Fit the model
fit <- fmri_lm(emod, dataset = bm$core_data)
```

## Documentation

See comprehensive examples and tutorials in the [vignettes](https://bbuchsbaum.github.io/fmrireg/articles/index.html):

- [Hemodynamic Response Functions](https://bbuchsbaum.github.io/fmrireg/articles/a_01_hemodynamic_response.html)
- [Building Regressors](https://bbuchsbaum.github.io/fmrireg/articles/a_02_regressor.html)
- [Event Models](https://bbuchsbaum.github.io/fmrireg/articles/a_04_event_models.html)
- [Statistical Contrasts](https://bbuchsbaum.github.io/fmrireg/articles/a_05_contrasts.html)
- [Linear Modeling](https://bbuchsbaum.github.io/fmrireg/articles/a_09_linear_model.html)

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

