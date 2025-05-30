fmrireg
================

# In Development

[![codecov](https://codecov.io/github/bbuchsbaum/fmrireg/branch/master/graphs/badge.svg)](https://codecov.io/github/bbuchsbaum/fmrireg)

[![Build
Status](https://img.shields.io/travis/bbuchsbaum/fmrireg/master.svg)](https://travis-ci.org/bbuchsbaum/fmrireg)

The `fmrireg` package is currently in development. The goal of the
project is to provide functions for carrying out regression anaysis with
*f*MRI data. It’s Not Ready for Primetime, but still have a look at
vignettes:

## Vignettes

See examples of use of `fmrireg` in the
[vignettes](https://bbuchsbaum.github.io/fmrireg/articles/index.html).

## Thread configuration

The internal C++ routines use [RcppParallel](https://rcppcore.github.io/RcppParallel/). You can control the number of threads by setting the R option `fmrireg.num_threads` or the environment variable `FMRIREG_NUM_THREADS` before loading the package. If either is set, `fmrireg` calls `RcppParallel::setThreadOptions()` when it loads.

