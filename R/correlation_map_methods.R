#' Correlation Map Methods for fmrireg
#'
#' These methods provide correlation heatmap visualizations for various model objects.

#' @rdname correlation_map
#' @method correlation_map baseline_model
#' @param method Correlation method: "pearson" (default) or "spearman"
#' @param half_matrix Logical; if TRUE (default), show only the lower triangle
#' @param absolute_limits Logical; if TRUE, set color limits to \[-1,1\] (default: TRUE)
#' @param within_run Logical; centre columns within runs and drop run
#'   intercepts before correlating (default: TRUE)
#' @param label_values Logical or NULL; print correlations in the cells
#'   (NULL: all cells when there are at most 12 columns)
#' @return A ggplot2 object containing the correlation heatmap visualization
#' @export
correlation_map.baseline_model <- function(x,
                                          method          = c("pearson", "spearman"),
                                          half_matrix     = TRUE,
                                          absolute_limits = TRUE,
                                          within_run      = TRUE,
                                          label_values    = NULL,
                                          ...) {
  DM <- as.matrix(design_matrix(x))
  pseudo <- list(baseline_model = x)
  .correlation_map_common(DM, method = method, half_matrix = half_matrix,
                          absolute_limits = absolute_limits,
                          info = .column_info(pseudo, DM),
                          runs = .run_lengths(pseudo),
                          within_run = within_run,
                          label_values = label_values, ...)
}
