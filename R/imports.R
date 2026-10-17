# Package Imports
# This file consolidates all package imports for the fmrireg package
# to ensure proper namespace management and CRAN compliance

#' @importFrom stats arima.sim as.formula coef coefficients cor cov df.residual
#' @importFrom stats fitted lm lm.fit lsfit mad median model.matrix pchisq pf
#' @importFrom stats pnorm predict pt qnorm quantile resid rexp rgamma rlnorm
#' @importFrom stats rnorm runif sd setNames terms toeplitz var
#' @importFrom utils data head read.csv setTxtProgressBar str txtProgressBar
#' @importFrom methods as
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom neuroim2 space
NULL

# Non-standard-evaluation names that R CMD check cannot resolve (roxyglobals
# generates R/globals.R for the rest).
utils::globalVariables(c(
  # design_plot(): ggplot2 aesthetics
  ".time", "value", "condition",
  # tidy.fmri_meta(): dplyr column
  "n_significant",
  # packaged dataset loaded by .ensure_benchmark_data_loaded()
  "fmri_benchmark_datasets"
))
