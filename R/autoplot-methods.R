#' Autoplot method for Reg objects
#' 
#' Creates a ggplot visualization of an fMRI regressor object.
#' 
#' @param object A `Reg` object (or one inheriting from it, like `regressor`).
#' @param grid Optional numeric vector specifying time points (seconds) for evaluation. 
#'   If NULL, a default grid is generated based on the object's onsets and span.
#' @param precision Numeric precision for HRF evaluation if `grid` needs generation or 
#'   if internal evaluation requires it (passed to `evaluate`).
#' @param method Evaluation method passed to `evaluate`.
#' @param ... Additional arguments (currently unused).
#' 
#' @return A ggplot object.
#' 
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap labs theme_minimal autoplot
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate across where
#' @export
#' @method autoplot Reg
#' @rdname autoplot
autoplot.Reg <- function(object, grid = NULL, precision = 0.1, method = "conv", ...) {
  
  # Determine evaluation grid if not provided
  if (is.null(grid)) {
    if (length(onsets(object)) > 0) {
        min_onset <- min(onsets(object), na.rm = TRUE)
        max_onset <- max(onsets(object), na.rm = TRUE)
        # Extend grid by HRF span, ensure start isn't negative if possible
        grid_start <- max(0, min_onset - 5) 
        grid_end <- max_onset + attr(object$hrf, "span") + 5
        # Use a reasonable step based on span or default to precision
        grid_step <- min(precision, attr(object$hrf, "span") / 100, 0.5) 
        grid <- seq(grid_start, grid_end, by = grid_step)
    } else {
        # Default grid if no onsets (e.g., null regressor)
        grid <- seq(0, (attr(object$hrf, "span") %||% 24), by = precision)
    }
  }
  
  # Evaluate the regressor
  eval_data <- fmrihrf::evaluate(object, grid, precision = precision, method = method)
  nb <- fmrihrf::nbasis(object)
  
  # Prepare data frame for plotting
  if (nb == 1) {
    plot_df <- data.frame(Time = grid, Response = as.vector(eval_data))
    p <- ggplot(plot_df, aes(x = .data$Time, y = .data$Response)) +
         geom_line() +
         labs(title = paste("Regressor:", attr(object$hrf, "name") %||% "Custom"), 
              x = "Time (seconds)", y = "Amplitude")
  } else {
    colnames(eval_data) <- paste0("Basis_", 1:nb)
    plot_df <- as.data.frame(eval_data)
    plot_df$Time <- grid
    
    plot_df_long <- plot_df |>
      tidyr::pivot_longer(cols = starts_with("Basis_"), 
                          names_to = "Basis", 
                          values_to = "Response")
                          
    p <- ggplot(plot_df_long, aes(x = .data$Time, y = .data$Response)) +
         geom_line() +
         facet_wrap(~ .data$Basis) +
         labs(title = paste("Regressor:", attr(object$hrf, "name") %||% "Custom", "(Multi-basis)"),
              x = "Time (seconds)", y = "Amplitude")
  }
  
  p + theme_minimal()
} 




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
#' @noRd
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