###############################################################################
## fmrimodel.R
##
## This file creates the overall fMRI model by combining an event model 
## (describing experimental events) and a baseline model (modeling drift, 
## nuisance, and block effects).
##
## The file provides two main functions:
##   - create_fmri_model(): Creates an fMRI model from a formula, block formula,
##     (optionally) a baseline_model, and an fmri_frame.
##   - fmri_model(): Combines an event_model and a baseline_model into an fmri_model.
##
## Additional functions build design matrices, compute contrasts, print, and plot
## the overall fMRI model.
###############################################################################

## ============================================================================
## Section 1: fMRI Model Construction Functions
## ============================================================================

#' Create an fMRI Model
#'
#' This function creates an fMRI model consisting of an event model and a baseline model.
#'
#' @param formula The model formula for experimental events.
#' @param block The model formula for block structure.
#' @param baseline_model (Optional) A \code{baseline_model} object. Default is \code{NULL}.
#' @param dataset An \code{fmri_frame} object containing the time-series data.
#' @param drop_empty Logical. Whether to remove factor levels with zero size. Default is \code{TRUE}.
#' @param durations A vector of event durations. Default is \code{0}.
#' @return An \code{fmri_model} object.
#' @export
#' @examples
#' \dontrun{
#' # Assuming you have an fmri_frame object named ds and a formula for events:
#' fmri_mod <- create_fmri_model(formula = onset ~ hrf(x) + hrf(y),
#'                               block = ~ run,
#'                               dataset = ds,
#'                               drop_empty = TRUE,
#'                               durations = rep(0, nrow(ds$event_table)))
#' }
create_fmri_model <- function(formula, block, baseline_model = NULL, dataset, drop_empty = TRUE, durations = 0) {
  assert_that(is.formula(formula), msg = "'formula' must be a formula")
  assert_that(is.formula(block), msg = "'block' must be a formula")
  assert_that(inherits(dataset, "fmri_frame"), msg = "'dataset' must be an 'fmri_frame'")
  assert_that(is.numeric(durations), msg = "'durations' must be numeric")
  formula <- .fmrireg_inject_registered_bases(formula)
  
  # Replicate durations if a single value is provided.
  if (length(durations) == 1) {
    durations <- rep(durations, nrow(.dset_event_table(dataset)))
  }
  
  # Resolve conflict: use a temporary variable to hold the baseline model.
  if (is.null(baseline_model)) {
    base_model_obj <- baseline_model(
      basis = "bs",
      degree = max(ceiling(median(.dset_run_lengths(dataset)) / 100), 3),
      sframe = .dset_sampling_frame(dataset)
    )
  } else {
    assert_that(inherits(baseline_model, "baseline_model"),
                msg = "'baseline_model' must have class 'baseline_model'")
    base_model_obj <- baseline_model
  }
  
  ev_model <- event_model(
    x = formula,
    block = block,
    data = .dset_event_table(dataset),
    sampling_frame = .dset_sampling_frame(dataset),
    drop_empty = drop_empty,
    durations = durations
  )
  
  fmri_model(ev_model, base_model_obj, dataset)
}


#' Construct an fMRI Regression Model
#'
#' This function constructs an fMRI regression model consisting of an event model
#' and a baseline model. The resulting model can be used for the analysis of fMRI data.
#'
#' @param event_model An object of class "event_model" representing the event-related part of the fMRI regression model.
#' @param baseline_model An object of class "baseline_model" representing the baseline-related part of the fMRI regression model.
#' @param dataset An \code{fmri_frame} used to build the model.
#' @return An object of class \code{fmri_model} containing the event and baseline models along with the dataset.
#' @export
#' @seealso event_model, baseline_model
fmri_model <- function(event_model, baseline_model, dataset) {
  assert_that(inherits(event_model, "event_model"))
  assert_that(inherits(baseline_model, "baseline_model"))
  assert_that(inherits(dataset, "fmri_frame"))

  fmodel <- list(event_model = event_model,
                 baseline_model = baseline_model,
                 dataset = dataset)
  class(fmodel) <- "fmri_model"
  fmodel
}


#' (Internal) Prediction Matrix
#'
#' This function is intended to compute a prediction matrix for the model.
#' (Currently a stub.)
#'
#' @param x An fmri_model object.
#' @return (Not implemented)
#' @keywords internal
#' @noRd
prediction_matrix <- function(x) {
  stop("not implemented")
}


## ============================================================================
## Section 2: Design Matrix and Environment for fMRI Models
## ============================================================================

#' Design Matrix for fMRI Models
#' 
#' Extract the combined design matrix from an fMRI model containing both event and baseline terms.
#' 
#' @param x An fmri_model object
#' @param blockid Optional numeric vector specifying which blocks/runs to include
#' @param ... Additional arguments (not used)
#' @return A tibble containing the combined design matrix with event and baseline terms
#' @examples
#' fm <- fmrireg:::.demo_fmri_model()
#' head(design_matrix(fm))
#' @method design_matrix fmri_model
#' @export
#' @importFrom tibble as_tibble
design_matrix.fmri_model <- function(x, blockid = NULL, ...) {
  suppressMessages(
    tibble::as_tibble(
      cbind(
        design_matrix(x$event_model, blockid),
        design_matrix(x$baseline_model, blockid)
      ),
      .name_repair = "check_unique"
    )
  )
}


#' @importFrom tibble as_tibble
#' @keywords internal
#' @noRd
design_env.fmri_model <- function(x, blockid = NULL, ...) {
  stop("Not implemented")
}


## ============================================================================
## Section 3: Accessor Functions for fMRI Models
## ============================================================================

#' @export
terms.fmri_model <- function(x, ...) {
  c(terms(x$event_model), terms(x$baseline_model))
}

#' @export
#' @autoglobal
cells.event_model <- function(x, ...) {
  eterms <- terms(x)
  if (length(eterms) == 0L) {
    return(tibble::tibble())
  }
  parts <- lapply(eterms, function(term) tibble::as_tibble(cells(term, ...)))
  dplyr::bind_rows(parts)
}

#' @export
#' @autoglobal
cells.fmri_model <- function(x, ...) {
  c1 <- tibble::as_tibble(cells(x$event_model, ...))
  if (nrow(c1) > 0L) c1$type <- "event"
  c2 <- tibble::as_tibble(cells(x$baseline_model, ...))
  if (nrow(c2) > 0L) c2$type <- "baseline"
  out <- dplyr::bind_rows(c1, c2)
  cols <- intersect(c("index", "type"), names(out))
  if (length(cols) > 0L) {
    out <- dplyr::relocate(out, dplyr::all_of(cols))
  }
  out
}

#' @export
blocklens.fmri_model <- function(x, ...) {
  fmrihrf::blocklens(x$event_model)
}

#' @export
event_terms.fmri_model <- function(x, ...) {
  terms(x$event_model)
}

#' @export
baseline_terms.fmri_model <- function(x, ...) {
  terms(x$baseline_model)
}

#' @export
contrast_weights.fmri_model <- function(x, ...) {
  contrast_weights(x$event_model, ...)
}

#' @export
conditions.fmri_model <- function(x, ...) {
  unlist(lapply(terms(x), function(t) conditions(t)), use.names = FALSE)
}

#' @export
conditions.baseline_model <- function(x, ...) {
  unlist(lapply(terms(x), function(t) conditions(t)), use.names = FALSE)
}


## ============================================================================
## Section 4: Plot and Print Methods for fMRI Models
## ============================================================================

#' Plot the regressors of an fmri_model over time
#'
#' Draws the design-matrix columns as time series in stacked small multiples
#' that share one time axis. Each event regressor gets its own row, coloured
#' by condition (designs with more than 12 event columns get one row per event
#' term instead). Baseline regressors are drawn in grey, one row per baseline
#' term (drift, nuisance). Run-specific columns are drawn only within their own
#' run, and no line is joined across a run boundary. Run intercepts are
#' omitted, since a constant carries no temporal information.
#'
#' @param x An \code{fmri_model}.
#' @param baseline Logical; include baseline rows (default \code{TRUE}).
#' @param ... Unused.
#' @return A ggplot2 object.
#' @export
plot.fmri_model <- function(x, baseline = TRUE, ...) {
  DM <- as.matrix(design_matrix(x))
  info <- .column_info(x, DM)
  n <- nrow(DM)
  runs <- .run_lengths(x) %||% n
  run_id <- rep(seq_along(runs), runs)
  sf <- .model_sframe(x)
  time <- if (!is.null(sf)) fmrihrf::samples(sf, global = TRUE) else seq_len(n)

  is_event <- info$kind == "event"
  keep <- (is_event | (baseline & info$kind != "intercept"))
  n_ev <- sum(is_event)
  per_term <- n_ev > 12L
  info$panel <- ifelse(is_event, if (per_term) info$group else info$label, info$group)
  cond_cols <- .condition_colours(info)
  info$colour <- ifelse(is_event, cond_cols[info$condition], .fmrireg_ink[["muted"]])
  if (per_term) {
    # a collapsed row with many conditions (e.g. trialwise) cannot use colour
    # for identity; draw it in one colour
    many <- tapply(info$condition[is_event], info$group[is_event],
                   function(x) length(unique(x)) > 8L)
    crowded <- names(many)[many]
    info$colour[is_event & info$group %in% crowded] <- .fmrireg_categorical(1)
  }
  info$colour[is.na(info$colour)] <- .fmrireg_ink[["muted"]]

  sel <- which(keep)
  vals <- DM[, sel, drop = FALSE]
  # run-specific columns exist only inside their run
  for (k in seq_along(sel)) {
    r <- info$run[sel[k]]
    if (!is.na(r)) vals[run_id != r, k] <- NA
  }
  df <- data.frame(
    time   = rep(time, times = length(sel)),
    value  = as.vector(vals),
    grp    = paste(rep(info$column[sel], each = n), rep(run_id, times = length(sel))),
    panel  = factor(rep(info$panel[sel], each = n), levels = unique(info$panel[sel])),
    colour = rep(info$colour[sel], each = n)
  )
  df <- df[!is.na(df$value), ]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$value,
                                        group = .data$grp, colour = .data$colour)) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::scale_colour_identity() +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$panel), scales = "free_y",
                        switch = "y") +
    ggplot2::scale_y_continuous(breaks = .row_breaks, labels = .axis_num)

  if (length(runs) > 1L) {
    starts <- c(1L, cumsum(runs)[-length(runs)] + 1L)
    dt <- if (length(time) > 1L) diff(time[1:2]) else 1
    bounds <- time[starts[-1]] - dt / 2
    mids <- (time[starts] + time[cumsum(runs)]) / 2
    p <- p + ggplot2::geom_vline(xintercept = bounds, colour = .fmrireg_ink[["muted"]],
                                 linewidth = 0.3) +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  sec.axis = ggplot2::dup_axis(breaks = mids,
                                                               labels = paste("Run", seq_along(runs)),
                                                               name = NULL))
  } else {
    p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
  }

  omitted <- sum(info$kind == "intercept")
  p +
    ggplot2::labs(title = "Model regressors over time",
                  caption = if (!is.null(sf)) sprintf("Regressors drawn as sampled at each scan (TR = %s s); each row has its own vertical scale.",
                                                      paste(unique(sf$TR), collapse = "/")) else NULL,
                  subtitle = sprintf("%d event columns%s%s", n_ev,
                                     if (baseline) sprintf(", %d baseline columns (grey)",
                                                           sum(keep & !is_event)) else "",
                                     if (baseline && omitted) sprintf("; %d run intercepts not shown", omitted) else ""),
                  x = "Time (s)", y = NULL) +
    .plot_theme(grid = "none") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.9)),
      axis.text.x.top = ggplot2::element_text(colour = .fmrireg_ink[["primary"]]),
      panel.spacing.y = ggplot2::unit(0.4, "lines")
    )
}

#' @export
#' @rdname print
#' @return Invisibly returns the input object x
print.fmri_model <- function(x, ...) {
  # Header with fancy border
  cat("\n=============================================")
  cat("\n             fMRI Model                     ")
  cat("\n=============================================")
  
  # Event Model Section
  cat("\n Event Model                                ")
  cat("\n---------------------------------------------")
  cat("\n| Formula:", crayon::cyan(Reduce(paste, deparse(x$event_model$model_spec$formula))))
  
  # Event Model Summary
  cat("\n| Summary:")
  cat("\n|   - Terms:", crayon::yellow(length(terms(x$event_model))))
  cat("\n|   - Events:", crayon::yellow(nrow(x$event_model$model_spec$event_table)))
  cat("\n|   - Design Columns:", crayon::yellow(length(conditions(x$event_model))))
  cat("\n|   - Blocks:", crayon::yellow(length(unique(x$event_model$blockids))))
  
  # Baseline Model Section (if present)
  if (!is.null(x$baseline_model)) {
    cat("\n---------------------------------------------")
    cat("\n| Baseline Model                           |")
    cat("\n| Components:")
    
    # Drift term info
    if (!is.null(x$baseline_model$drift_term)) {
      drift_name <- x$baseline_model$drift_term$varname
      basis_type <- x$baseline_model$drift_spec$basis
      degree <- x$baseline_model$drift_spec$degree
      drift_cols <- ncol(design_matrix(x$baseline_model$drift_term))
      cat("\n|   - Drift:", crayon::magenta(drift_name))
      cat("\n|     - Type:", crayon::blue(basis_type))
      cat("\n|     - Degree:", crayon::blue(degree))
      cat("\n|     - Columns:", crayon::yellow(drift_cols))
    }
    
    # Block term info
    if (!is.null(x$baseline_model$block_term)) {
      const_cols <- ncol(design_matrix(x$baseline_model$block_term))
      cat("\n|   - Block Terms:", crayon::yellow(const_cols), "columns")
    }
    
    # Nuisance term info
    if (!is.null(x$baseline_model$nuisance_term)) {
      nuis_cols <- ncol(design_matrix(x$baseline_model$nuisance_term))
      cat("\n|   - Nuisance Terms:", crayon::yellow(nuis_cols), "columns")
    }
  }
  
  # Total Model Summary
  cat("\n---------------------------------------------")
  cat("\n| Total Model                              |")
  total_cols <- ncol(design_matrix(x))
  cat("\n|   - Total Design Columns:", crayon::yellow(total_cols))
  
  # Footer
  cat("\n=============================================\n")
}

#' correlation_map.fmri_model
#'
#' @description
#' Heatmap of pairwise correlations between the columns of an \code{fmri_model}'s
#' combined event + baseline design matrix, grouped into panels by model term
#' (event terms, drift, nuisance).
#'
#' By default (\code{within_run = TRUE}) each column is centred within each
#' run and the run-intercept columns are dropped before correlating. This
#' shows the collinearity that matters for estimation once run means are
#' modelled, and removes spurious correlations that zero-padding creates
#' between columns belonging to different runs. Cells for pairs of columns
#' that live in different runs are left blank.
#'
#' @param x An \code{fmri_model}.
#' @param method Correlation method (e.g., "pearson", "spearman").
#' @param half_matrix Logical; if TRUE (default), display only the lower triangle.
#' @param absolute_limits Logical; if TRUE, the colour scale spans -1 to 1;
#'   otherwise it spans plus or minus the largest absolute correlation shown.
#' @param within_run Logical; centre columns within runs and drop run
#'   intercepts before correlating (default \code{TRUE}). \code{FALSE}
#'   correlates the raw columns.
#' @param label_values Logical or \code{NULL}; print correlations in the cells.
#'   \code{NULL} (default) labels every cell when there are at most 12 columns,
#'   and otherwise only cells with |r| >= 0.3 that involve an event regressor.
#' @param ... Additional arguments passed to \code{\link[ggplot2]{geom_tile}}.
#' @return A ggplot2 object.
#' @export
correlation_map.fmri_model <- function(x,
                                       method          = c("pearson", "spearman"),
                                       half_matrix     = TRUE,
                                       absolute_limits = TRUE,
                                       within_run      = TRUE,
                                       label_values    = NULL,
                                       ...) {
  DM <- as.matrix(design_matrix(x))
  .correlation_map_common(DM, method = method, half_matrix = half_matrix,
                          absolute_limits = absolute_limits,
                          info = .column_info(x, DM), runs = .run_lengths(x),
                          within_run = within_run,
                          label_values = label_values, ...)
}


#' Heatmap visualization of the combined fmri_model design matrix
#'
#' @description
#' Produces a single heatmap of all columns in the design matrix of an
#' \code{fmri_model}, which merges the event_model and baseline_model
#' regressors. Rows are scans (time runs downward); columns are regressors,
#' kept in model order and grouped into panels by term (event terms, drift,
#' intercept, nuisance). Run boundaries are marked with horizontal rules, and
#' the subtitle reports the rank and condition number of the design.
#'
#' By default each column is divided by its maximum absolute value, so
#' regressors on very different scales (an intercept, a motion trace, a
#' convolved event train) are all legible on one colour scale. Set
#' \code{scale_columns = FALSE} to show raw values. Run-intercept columns
#' are drawn in grey, since they carry no information beyond run membership.
#'
#' The condition number is computed after scaling every column to unit
#' length, so it reflects collinearity rather than the units of the columns.
#'
#' @param x An \code{fmri_model} object.
#' @param block_separators Logical; if \code{TRUE}, draw rules between runs.
#' @param rotate_x_text Logical; if \code{TRUE}, rotate x-axis labels to vertical.
#' @param fill_midpoint Numeric or \code{NULL}; centre of the diverging colour
#'   scale. Defaults to 0.
#' @param fill_limits Numeric vector of length 2 or \code{NULL}; passed to the fill scale
#'   \code{limits=} argument. This can clip or expand the color range.
#' @param scale_columns Logical; if \code{TRUE} (default), scale each column to
#'   a maximum absolute value of 1.
#' @param ... Additional arguments passed to \code{\link[ggplot2]{geom_raster}}.
#'
#' @return A ggplot2 plot object.
#' @export
design_map.fmri_model <- function(x,
                                  block_separators = TRUE,
                                  rotate_x_text    = TRUE,
                                  fill_midpoint    = NULL,
                                  fill_limits      = NULL,
                                  scale_columns    = TRUE,
                                  ...) {
  DM <- as.matrix(design_matrix(x))
  info <- .column_info(x, DM)
  n_scans <- nrow(DM)
  p <- ncol(DM)

  # rank and condition number on unit-length columns
  norms <- sqrt(colSums(DM^2, na.rm = TRUE))
  norms[norms == 0] <- 1
  sv <- svd(sweep(DM, 2, norms, "/"), nu = 0, nv = 0)$d
  tol <- max(dim(DM)) * max(sv) * .Machine$double.eps
  rk <- sum(sv > tol)
  kappa_val <- if (rk == p) max(sv) / min(sv) else Inf

  if (scale_columns) {
    mx <- apply(abs(DM), 2, max, na.rm = TRUE)
    mx[!is.finite(mx) | mx == 0] <- 1
    DM <- sweep(DM, 2, mx, "/")
  }
  vals <- DM
  icpt <- info$kind == "intercept"
  vals[, icpt][vals[, icpt] != 0] <- NA

  df <- data.frame(
    scan      = rep(seq_len(n_scans), times = p),
    regressor = factor(rep(info$label, each = n_scans), levels = info$label),
    group     = factor(rep(info$group, each = n_scans), levels = unique(info$group)),
    value     = as.vector(vals)
  )

  lim <- fill_limits %||% if (scale_columns) c(-1, 1) else {
    m <- max(abs(DM), na.rm = TRUE); c(-m, m)
  }
  plt <- ggplot2::ggplot(df, ggplot2::aes(x = .data$regressor, y = .data$scan,
                                          fill = .data$value)) +
    ggplot2::geom_raster(...) +
    ggplot2::facet_grid(cols = ggplot2::vars(.data$group),
                        scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradient2(
      low = .fmrireg_diverging[["low"]], mid = .fmrireg_diverging[["mid"]],
      high = .fmrireg_diverging[["high"]], midpoint = fill_midpoint %||% 0,
      limits = lim, oob = .squish, na.value = "#BDBDBD",
      name = if (scale_columns) paste0("Value ", intToUtf8(247L), "\ncolumn's\nmax |value|") else "Value"
    ) +
    # the run of a run-specific column is visible from its block, so the
    # axis shows the short name
    ggplot2::scale_x_discrete(expand = c(0, 0), labels = function(x) {
      sub(paste0(" ", .glyph("dot"), " run [0-9]+$"), "", x)
    })

  runs <- .run_lengths(x)
  if (!is.null(runs) && length(runs) > 1L && sum(runs) == n_scans) {
    starts <- c(1L, cumsum(runs)[-length(runs)] + 1L)
    mids <- starts + (runs - 1) / 2
    plt <- plt + ggplot2::scale_y_reverse(
      breaks = mids, labels = paste("Run", seq_along(runs)),
      sec.axis = ggplot2::dup_axis(breaks = c(starts, n_scans),
                                   labels = c(starts, n_scans), name = "Scan"),
      expand = c(0, 0)
    )
    if (block_separators) {
      plt <- plt + ggplot2::geom_hline(yintercept = starts[-1] - 0.5,
                                       colour = "white", linewidth = 0.9)
    }
  } else {
    plt <- plt + ggplot2::scale_y_reverse(expand = c(0, 0), name = "Scan")
  }

  cond_txt <- if (rk < p) {
    sprintf("rank %d of %d: RANK DEFICIENT", rk, p)
  } else {
    sprintf("full rank, condition number %s", format(signif(kappa_val, 3), big.mark = ","))
  }
  plt +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Design matrix",
                  subtitle = sprintf("%d scans %s %d regressors; %s", n_scans, .glyph("times"), p, cond_txt),
                  caption = if (any(icpt)) "Grey: run intercepts." else NULL) +
    .plot_theme(grid = "none") +
    ggplot2::theme(
      axis.text.x = if (rotate_x_text) {
        ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5,
                              size = ggplot2::rel(if (p > 60) 0.55 else 0.8))
      } else {
        ggplot2::element_text()
      },
      axis.text.y.left = ggplot2::element_text(colour = .fmrireg_ink[["primary"]]),
      axis.title.y.right = ggplot2::element_text(size = ggplot2::rel(0.8)),
      panel.spacing.x = ggplot2::unit(0.3, "lines"),
      strip.text = ggplot2::element_text(hjust = 0, face = "bold", size = ggplot2::rel(0.75)),
      strip.clip = "off",
      legend.key.height = ggplot2::unit(1.6, "lines"),
      legend.key.width = ggplot2::unit(0.6, "lines")
    )
}
