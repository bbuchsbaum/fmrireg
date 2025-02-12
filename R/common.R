#' Internal helper to build a correlation heatmap from a numeric matrix
#'
#' @param DM A numeric matrix of regressors (columns).
#' @param method Correlation method, passed to `stats::cor()` (e.g. "pearson" or "spearman").
#' @param half_matrix Logical; if TRUE, show only the lower-triangle (including diagonal).
#' @param absolute_limits Logical; if TRUE, force the fill color scale to -1..+1.
#' @param ... Additional arguments passed to `geom_tile()`.
#'
#' @return A ggplot2 object.
#' @keywords internal
.correlation_map_common <- function(DM,
                                    method         = c("pearson", "spearman"),
                                    half_matrix    = FALSE,
                                    absolute_limits = TRUE,
                                    ...) {
  method <- match.arg(method)
  stopifnot(is.matrix(DM), ncol(DM) >= 2)  # must have at least 2 columns to correlate
  
  # 1) Compute correlation matrix
  cormat <- stats::cor(DM, use = "pairwise.complete.obs", method = method)
  
  # 2) Optionally mask to show only the lower triangle.
  #    (Below we keep the diagonal, so i >= j.)
  if (half_matrix) {
    cormat[upper.tri(cormat)] <- NA  # set upper triangle to NA
  }
  
  # 3) Convert to a long data frame:
  #    row var becomes `Var1`, column var -> `Var2`, correlation -> `value`
  #    We remove any rows with NA if half_matrix=TRUE
  df_long <- as.data.frame(as.table(cormat), stringsAsFactors = FALSE)
  names(df_long) <- c("Var1", "Var2", "Correlation")
  if (half_matrix) {
    df_long <- df_long[!is.na(df_long$Correlation), ]
  }
  
  # Ensure factor ordering matches the original column order
  varnames <- colnames(DM)
  df_long$Var1 <- factor(df_long$Var1, levels = varnames)
  df_long$Var2 <- factor(df_long$Var2, levels = varnames)
  
  # 4) Base ggplot
  plt <- ggplot2::ggplot(df_long, 
                         ggplot2::aes(x = Var2, y = Var1, fill = Correlation)) +
    ggplot2::geom_tile(...)
  
  # Flip the y-axis so the first row is at the top (like many correlation-lower-triangle plots)
  plt <- plt + ggplot2::scale_y_discrete(limits = rev(levels(df_long$Var1)))
  
  # 5) Decide color-scale limits
  fillmin <- if (absolute_limits) -1 else min(df_long$Correlation, na.rm=TRUE)
  fillmax <- if (absolute_limits)  1 else max(df_long$Correlation, na.rm=TRUE)
  
  plt <- plt + 
    ggplot2::scale_fill_gradient2(
      limits   = c(fillmin, fillmax),
      midpoint = 0,
      low  = "blue",
      mid  = "white",
      high = "red"
    )
  
  # 6) Theming
  plt <- plt + 
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      axis.title  = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      panel.grid  = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  plt
}  