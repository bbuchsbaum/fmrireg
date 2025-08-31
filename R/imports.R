# Package Imports
# This file consolidates all package imports for the fmrireg package
# to ensure proper namespace management and CRAN compliance

#' @importFrom stats arima.sim as.formula coef coefficients cor cov df.residual
#' @importFrom stats fitted lm lm.fit lsfit mad median model.matrix pchisq pf
#' @importFrom stats pnorm predict pt qnorm quantile resid rexp rgamma rlnorm
#' @importFrom stats rnorm runif sd terms toeplitz var
#' @importFrom utils data head read.csv setTxtProgressBar str txtProgressBar
#' @importFrom methods as
#' @importFrom magrittr %>%
#' @importFrom ggplot2 geom_tile scale_y_reverse scale_fill_gradientn
#' @importFrom ggplot2 scale_fill_gradient2 annotate theme element_blank element_text
#' @importFrom dplyr starts_with
#' @importFrom rlang .data
#' @importFrom neuroim2 space
NULL

# Global variables used in ggplot2 aesthetics and data manipulation
# This prevents R CMD check NOTEs about undefined global variables
utils::globalVariables(c(
  # Variables used in .correlation_map_common
  "Var1", "Var2", "Correlation",
  # Variables used in design_plot
  ".time", ".block", "value", "condition",
  # Variables used in design_map.fmri_model
  "scan_number", "Regressor", "Value",
  # Variables used in autoplot
  ".data",
  # Variables used in meta-analysis functions
  "n_significant",
  # Variables from chunkwise processing
  "data",
  # Variables from robust fitting
  "proj_global_robustly_weighted", "multiresponse_arma", "auto_whiten",
  # Package data
  "fmri_benchmark_datasets"
))