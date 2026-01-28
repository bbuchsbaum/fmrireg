# Multiresponse bootstrap linear model

Performs block bootstrap resampling for multiresponse linear models,
particularly useful for fMRI time series data where temporal
dependencies exist. The function implements a block bootstrap approach
to maintain the temporal correlation structure within the data.

## Usage

``` r
multiresponse_bootstrap_lm(
  form,
  data_env,
  conlist,
  vnames,
  fcon,
  modmat = NULL,
  block_size = 30,
  boot_rows = FALSE,
  nboot = 100,
  event_indices
)
```

## Arguments

- form:

  Formula for the linear model. Required if modmat is NULL.

- data_env:

  Environment containing the data for the linear model.

- conlist:

  List of contrasts to be computed for each bootstrap sample.

- vnames:

  Vector of variable names.

- fcon:

  Contrasts for fixed effects.

- modmat:

  Optional pre-computed model matrix. If provided, form is ignored.

- block_size:

  Size of the blocks for the bootstrap (default: 30). Should be large
  enough to capture temporal dependencies but small enough to allow
  sufficient randomization.

- boot_rows:

  Logical flag indicating whether to bootstrap rows (default: FALSE).

- nboot:

  Number of bootstrap iterations (default: 100).

- event_indices:

  Indices of events for computing beta covariances.

## Value

A list containing:

- original: The fitted original model with contrasts

- con_cov: Covariance matrices for contrasts (if contrasts provided)

- beta_cov: Covariance matrices for beta estimates

- nboot: Number of bootstrap iterations performed

- bootstrap: Logical indicating this is a bootstrap result

## Details

The function performs the following steps:

1.  Fits the original linear model

2.  Implements block bootstrap resampling of residuals

3.  Reconstructs response variables using fitted values and resampled
    residuals

4.  Computes contrasts for each bootstrap sample

The block bootstrap approach helps preserve temporal dependencies in the
data by resampling blocks of consecutive observations rather than
individual observations.
