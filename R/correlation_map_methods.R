#' Correlation Map Methods for fmrireg
#'
#' These methods provide correlation heatmap visualizations for various model objects.
#' They are thin wrappers around methods from fmridesign when appropriate.

#' @rdname correlation_map
#' @method correlation_map baseline_model
#' @param method Correlation method: "pearson" (default) or "spearman"
#' @param half_matrix Logical; if TRUE, show only lower triangle (default: FALSE)
#' @param absolute_limits Logical; if TRUE, set color limits to \[-1,1\] (default: TRUE)
#' @export
correlation_map.baseline_model <- function(x,
                                          method          = c("pearson", "spearman"),
                                          half_matrix     = FALSE,
                                          absolute_limits = TRUE,
                                          ...) {
  # Extract the design matrix
  DM <- as.matrix(design_matrix(x))
  
  # Call the internal helper function that's defined in autoplot-methods.R
  .correlation_map_common(DM, method = method, half_matrix = half_matrix,
                          absolute_limits = absolute_limits, ...)
}