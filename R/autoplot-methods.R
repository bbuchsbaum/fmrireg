#' Autoplot method for Reg objects
#'
#' Plots the predicted BOLD time course of an fMRI regressor: the event train
#' convolved with its HRF. Event onsets are marked as ticks beneath the curve,
#' so the hemodynamic lag is visible. When responses to nearby events overlap,
#' each event's own response is drawn faintly under the summed curve.
#' Multi-basis HRFs are drawn as stacked panels that share the time axis, each
#' with its own vertical scale.
#'
#' @param object A `Reg` object (or one inheriting from it, like `regressor`).
#' @param grid Optional numeric vector specifying time points (seconds) for evaluation.
#'   If NULL, a default grid spans the first onset to the end of the last
#'   event's response (including its duration).
#' @param precision Numeric precision for HRF evaluation if `grid` needs generation or
#'   if internal evaluation requires it (passed to `evaluate`).
#' @param method Evaluation method passed to `evaluate`.
#' @param title Optional plot title, e.g. the condition this regressor models.
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 autoplot
#' @export
#' @method autoplot Reg
#' @rdname autoplot
autoplot.Reg <- function(object, grid = NULL, precision = 0.1, method = "conv",
                         title = NULL, ...) {
  # Use fmrihrf::onsets - Reg methods live on the fmrihrf generic, while the
  # package-level `onsets` import is from fmridesign and does not dispatch for Reg.
  ons <- fmrihrf::onsets(object)
  durs <- rep_len(object$duration %||% 0, length(ons))
  amps <- rep_len(object$amplitude %||% 1, length(ons))
  span <- attr(object$hrf, "span") %||% object$span %||% 24
  if (is.null(grid)) {
    if (length(ons) > 0) {
      grid_start <- max(0, min(ons, na.rm = TRUE) - 5)
      # start at 0 when the first event is near the start of the run
      if (grid_start < 0.1 * (max(ons + durs, na.rm = TRUE) + span)) grid_start <- 0
      grid_end <- max(ons + durs, na.rm = TRUE) + span + 5
      grid_step <- min(precision, span / 100, 0.5)
      grid <- seq(grid_start, grid_end, by = grid_step)
    } else {
      grid <- seq(0, span, by = precision)
    }
  }

  eval_data <- as.matrix(fmrihrf::evaluate(object, grid, precision = precision,
                                           method = method))
  nb <- fmrihrf::nbasis(object)
  hrf_name <- attr(object$hrf, "name") %||% "custom HRF"
  basis_labels <- .basis_labels(hrf_name, nb)
  lev <- basis_labels

  df <- data.frame(
    time  = rep(grid, times = nb),
    value = as.vector(eval_data),
    basis = factor(rep(basis_labels, each = length(grid)), levels = lev)
  )
  ons_in <- ons[ons >= min(grid) & ons <= max(grid)]
  # onset ticks once, under the bottom panel
  ticks <- data.frame(onset = ons_in,
                      basis = factor(rep(basis_labels[nb], length(ons_in)), levels = lev))

  # individual event responses, drawn only when they overlap
  overlap <- length(ons) > 1L && any(diff(sort(ons)) < span)
  singles <- NULL
  if (overlap && length(ons) <= 60L) {
    singles <- do.call(rbind, lapply(seq_along(ons), function(i) {
      ri <- fmrihrf::regressor(onsets = ons[i], hrf = object$hrf, duration = durs[i],
                               amplitude = amps[i], span = span)
      vi <- as.matrix(fmrihrf::evaluate(ri, grid, precision = precision, method = method))
      data.frame(time = rep(grid, times = nb), value = as.vector(vi),
                 basis = factor(rep(basis_labels, each = length(grid)), levels = lev),
                 event = i)
    }))
    singles <- singles[abs(singles$value) > 1e-6 * max(abs(eval_data)), ]
  }

  modulated <- length(unique(round(amps, 8))) > 1L
  blocks <- any(durs > 0)
  subtitle <- paste0(
    sprintf("%d event%s", length(ons), if (length(ons) == 1L) "" else "s"),
    if (blocks) sprintf(", duration %s s", paste(unique(signif(durs, 3)), collapse = "/")) else "",
    if (modulated) ", amplitude-modulated" else "",
    "; ticks mark onsets"
  )
  # peak of a single unit event, so readers can relate betas to response size
  unit_peak <- tryCatch({
    h1 <- as.matrix(object$hrf(seq(0, span, by = 0.1)))[, 1]
    h1[which.max(abs(h1))]
  }, error = function(e) NA_real_)
  # are the basis functions orthogonal? If not, the first coefficient alone
  # is not the response amplitude
  ortho_note <- NULL
  if (nb > 1) {
    # correlation of the plotted regressor columns, as a reader would compute it
    rr <- suppressWarnings(stats::cor(eval_data)[1, -1])
    if (all(is.finite(rr))) {
      j <- which.max(abs(rr))
      ortho_note <- if (max(abs(rr)) < 0.05) {
        "The basis regressors are uncorrelated."
      } else {
        sprintf("The %s and %s regressors correlate (r = %.2f), so the first coefficient alone is not the response amplitude.",
                basis_labels[1], basis_labels[j + 1], rr[j])
      }
    }
  }
  notes <- c(
    if (!is.null(singles)) "Thin lines: each event alone; thick line: their sum.",
    if (nb > 1) paste("Each basis function has its own vertical scale.", ortho_note %||% ""),
    if (nb == 1 && is.finite(unit_peak)) sprintf("A single unit-amplitude event peaks at %s.",
                                                 formatC(unit_peak, format = "fg", digits = 3))
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$value)) +
    ggplot2::geom_hline(yintercept = 0, colour = .fmrireg_ink[["rule"]], linewidth = 0.4)
  if (!is.null(singles)) {
    p <- p + ggplot2::geom_line(data = singles, ggplot2::aes(group = .data$event),
                                colour = "#9E9E9E", linewidth = 0.3)
  }
  xb <- pretty(range(grid), n = 7)
  xb <- sort(unique(c(xb[xb >= min(grid) & xb <= max(grid)])))
  p <- p +
    ggplot2::geom_rug(data = ticks, ggplot2::aes(x = .data$onset),
                      inherit.aes = FALSE, sides = "b",
                      colour = .fmrireg_ink[["secondary"]],
                      length = ggplot2::unit(0.05, "npc"), linewidth = 0.5) +
    ggplot2::geom_line(colour = .fmrireg_categorical(1), linewidth = 0.7) +
    ggplot2::scale_x_continuous(breaks = xb, limits = range(grid), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(labels = .axis_num) +
    ggplot2::labs(title = title %||% sprintf("Regressor (%s HRF)", hrf_name),
                  subtitle = subtitle, caption = .wrap_notes(notes),
                  x = "Time (s)", y = if (nb == 1) "Predicted BOLD (a.u.)" else NULL) +
    .plot_theme() +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 16, 8, 8))

  if (nb > 1) {
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(.data$basis),
                                 scales = "free_y", switch = "y")
  }
  p
}

#' Readable names for HRF basis functions
#' @keywords internal
#' @noRd
.basis_labels <- function(hrf_name, nb) {
  if (nb == 1L) {
    return("response")
  }
  nm <- toupper(hrf_name %||% "")
  if (nb == 2L && grepl("SPMG2", nm)) {
    return(c("canonical", "temporal deriv."))
  }
  if (nb == 3L && grepl("SPMG3", nm)) {
    return(c("canonical", "temporal deriv.", "dispersion deriv."))
  }
  paste("basis", seq_len(nb))
}


#' Internal helper to build a correlation heatmap from a numeric matrix
#'
#' @param DM A numeric matrix of regressors (columns).
#' @param method Correlation method, passed to `stats::cor()`.
#' @param half_matrix Logical; if TRUE, show only the lower triangle.
#' @param absolute_limits Logical; if TRUE, force the fill scale to -1..+1.
#' @param info Optional column metadata from `.column_info()`.
#' @param runs Optional run lengths (scans), used by `within_run`.
#' @param within_run Logical; centre within runs and drop run intercepts.
#' @param label_values Logical or NULL; print r in cells.
#' @param ... Additional arguments passed to `geom_tile()`.
#'
#' @return A ggplot2 object.
#' @keywords internal
#' @noRd
.correlation_map_common <- function(DM,
                                    method          = c("pearson", "spearman"),
                                    half_matrix     = TRUE,
                                    absolute_limits = TRUE,
                                    info            = NULL,
                                    runs            = NULL,
                                    within_run      = TRUE,
                                    label_values    = NULL,
                                    ...) {
  method <- match.arg(method)
  stopifnot(is.matrix(DM), ncol(DM) >= 2)
  if (is.null(info)) {
    info <- data.frame(column = colnames(DM),
                       label = make.unique(.pretty_regressor_names(colnames(DM))),
                       group = "", kind = "other", run = NA_integer_,
                       stringsAsFactors = FALSE)
  }

  use_runs <- within_run && !is.null(runs) && length(runs) > 1L && sum(runs) == nrow(DM)
  run_id <- if (use_runs) rep(seq_along(runs), runs) else rep(1L, nrow(DM))
  # VIFs always describe the design as fit: run intercepts in, columns centred
  # within runs, whatever the display settings
  has_runs <- !is.null(runs) && sum(runs) == nrow(DM)
  info_vif <- info[info$kind != "intercept", , drop = FALSE]
  DM_vif <- .centre_within(DM[, info$kind != "intercept", drop = FALSE],
                           if (has_runs) rep(seq_along(runs), runs) else rep(1L, nrow(DM)))
  if (within_run) {
    keep <- info$kind != "intercept"
    DM <- DM[, keep, drop = FALSE]
    info <- info[keep, , drop = FALSE]
    for (r in unique(run_id)) {
      rows <- run_id == r
      DM[rows, ] <- sweep(DM[rows, , drop = FALSE], 2,
                          colMeans(DM[rows, , drop = FALSE]), "-")
    }
  }
  p <- ncol(DM)
  if (p < 2) {
    stop("Fewer than two columns remain to correlate after removing run intercepts.",
         call. = FALSE)
  }
  const <- apply(DM, 2, function(v) {
    s <- stats::sd(v)
    !is.finite(s) || s <= 1e-12 * max(1, max(abs(v)))
  })
  if (all(const)) {
    stop("Every design column is constant (after centring within runs); there is nothing to correlate.",
         call. = FALSE)
  }

  cor_q <- function(M) suppressWarnings(stats::cor(M, use = "pairwise.complete.obs", method = method))
  cormat <- cor_q(DM)
  ri <- info$run
  if (use_runs && any(!is.na(ri))) {
    # a pair involving a run-specific column is correlated on that run's rows,
    # so zero rows from the other runs do not dilute it
    for (r in sort(unique(ri[!is.na(ri)]))) {
      rows <- run_id == r
      Cr <- cor_q(DM[rows, , drop = FALSE])
      in_r <- which(ri == r)
      cormat[in_r, ] <- Cr[in_r, ]
      cormat[, in_r] <- Cr[, in_r]
    }
    diff_run <- outer(ri, ri, function(a, b) !is.na(a) & !is.na(b) & a != b)
    cormat[diff_run] <- NA
  } else {
    diff_run <- matrix(FALSE, p, p)
  }
  diag(cormat) <- NA
  const_cell <- outer(const, const, `|`) & !diff_run
  diag(const_cell) <- FALSE
  cormat[const_cell] <- NA
  shown <- !is.na(cormat) | const_cell
  if (half_matrix) {
    shown[upper.tri(shown)] <- FALSE
  }

  vif <- .event_vif(DM_vif, info_vif)

  labels <- info$label
  groups <- info$group
  glev <- unique(groups)
  is_ev <- info$kind == "event"

  df <- data.frame(
    row   = factor(rep(labels, times = p), levels = rev(labels)),
    col   = factor(rep(labels, each = p), levels = labels),
    rgrp  = factor(rep(groups, times = p), levels = glev),
    cgrp  = factor(rep(groups, each = p), levels = glev),
    ev    = rep(is_ev, times = p) | rep(is_ev, each = p),
    r     = as.vector(cormat),
    const = as.vector(const_cell)
  )
  df <- df[as.vector(shown), ]

  lim <- if (absolute_limits) c(-1, 1) else {
    m <- max(abs(df$r), na.rm = TRUE); c(-m, m)
  }

  plt <- ggplot2::ggplot(df, ggplot2::aes(x = .data$col, y = .data$row)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$r), colour = .fmrireg_ink[["rule"]],
                       linewidth = if (p <= 60) 0.2 else 0, ...) +
    ggplot2::scale_fill_gradient2(
      low = .fmrireg_diverging[["low"]], mid = "white",
      high = .fmrireg_diverging[["high"]], midpoint = 0, limits = lim,
      name = "r", na.value = "#BDBDBD",
      breaks = if (absolute_limits) c(-1, -0.5, 0, 0.5, 1) else ggplot2::waiver()
    ) +
    ggplot2::scale_x_discrete(expand = c(0, 0), drop = TRUE, labels = function(x) {
      # each event column carries its VIF where the reader looks for it
      v <- vif[x]
      ifelse(is.na(v), x, sprintf("%s (VIF %s%s)", x,
                                  ifelse(is.infinite(v), "Inf", formatC(v, format = "f", digits = 1)),
                                  ifelse(v > 5, "!", "")))
    }) +
    ggplot2::scale_y_discrete(expand = c(0, 0), drop = TRUE)

  if (is.null(label_values)) label_values <- NA
  ok <- !is.na(df$r)
  lab <- if (isTRUE(label_values) || (is.na(label_values) && p <= 12)) {
    df[ok, ]
  } else if (is.na(label_values)) {
    # label the substantial correlations; above 40 columns only those that
    # involve an event regressor, to keep the map readable
    df[ok & abs(df$r) >= 0.3 & (df$ev | p <= 40), ]
  } else {
    df[0, ]
  }
  if (nrow(lab)) {
    lab$txt <- sub("^(-?)0\\.", "\\1.", sprintf(if (p <= 12) "%.2f" else "%.1f", lab$r))
    lab$ink <- ifelse(abs(lab$r) > 0.6, "white", .fmrireg_ink[["primary"]])
    plt <- plt + ggplot2::geom_text(
      data = lab, ggplot2::aes(label = .data$txt, colour = .data$ink),
      size = if (p <= 12) 3 else if (p <= 30) 2.6 else 2
    ) + ggplot2::scale_colour_identity()
  }

  if (length(glev) > 1L) {
    # names of very narrow groups would overprint their neighbours
    width <- table(factor(groups, levels = glev))
    narrow <- names(width)[width / p < 0.04]
    ev_groups <- unique(groups[is_ev])
    # narrow event groups keep an abbreviated name (they are what readers look
    # for); narrow baseline groups are left unlabelled
    lab_fun <- function(x) {
      x <- as.character(x)
      ifelse(x %in% narrow,
             ifelse(x %in% ev_groups, abbreviate(x, 4, named = FALSE), ""),
             x)
    }
    plt <- plt + ggplot2::facet_grid(rows = ggplot2::vars(.data$rgrp),
                                     cols = ggplot2::vars(.data$cgrp),
                                     scales = "free", space = "free",
                                     switch = "y", drop = TRUE,
                                     labeller = ggplot2::as_labeller(lab_fun))
  }

  subtitle <- sprintf("%s correlation, %s; %d columns",
                      if (method == "pearson") "Pearson" else "Spearman",
                      if (within_run) "within run, intercepts removed" else "raw columns", p)
  notes <- .vif_note(vif)
  if (nrow(lab) && !isTRUE(label_values)) {
    notes <- c(notes, sprintf("Numbers: |r| >= 0.3%s.", if (p > 40) ", event regressors only" else ""))
  }
  if (use_runs && any(!is.na(ri)) && any(diff_run)) {
    notes <- c(notes, "Blank: columns from different runs. Pairs with a run-specific column use that run's scans.")
  }
  if (any(const)) {
    notes <- c(notes, paste0("Grey: constant column (", paste(labels[const], collapse = ", "), ")."))
  }
  txt_size <- if (p <= 30) 0.8 else if (p <= 80) 0.6 else 0.45
  plt +
    ggplot2::labs(x = NULL, y = NULL, title = "Regressor correlations", subtitle = subtitle,
                  caption = .wrap_notes(notes)) +
    .plot_theme(grid = "none") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5,
                                          size = ggplot2::rel(txt_size)),
      axis.text.y = ggplot2::element_text(size = ggplot2::rel(txt_size)),
      panel.spacing = ggplot2::unit(0.3, "lines"),
      panel.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(hjust = 0, size = ggplot2::rel(0.75)),
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, face = "bold",
                                                size = ggplot2::rel(0.75)),
      strip.clip = "off",
      legend.key.height = ggplot2::unit(1.6, "lines"),
      legend.key.width = ggplot2::unit(0.6, "lines")
    )
}

#' Centre each column within each run
#' @keywords internal
#' @noRd
.centre_within <- function(M, run_id) {
  for (r in unique(run_id)) {
    rows <- run_id == r
    M[rows, ] <- sweep(M[rows, , drop = FALSE], 2, colMeans(M[rows, , drop = FALSE]), "-")
  }
  M
}

#' Variance inflation factors of the event columns
#'
#' Computed by regressing each event column on every other non-constant
#' column of the (run-centred) design. Columns aliased with others get `Inf`.
#' @return Named numeric vector (labels as names), possibly empty.
#' @keywords internal
#' @noRd
.event_vif <- function(Z, info) {
  ev <- which(info$kind == "event")
  if (!length(ev)) return(stats::setNames(numeric(0), character(0)))
  sds <- apply(Z, 2, stats::sd)
  usable <- which(is.finite(sds) & sds > 1e-12 * pmax(1, apply(abs(Z), 2, max)))
  Zu <- Z[, usable, drop = FALSE]
  q <- qr(Zu)
  aliased <- usable[q$pivot[seq_len(ncol(Zu)) > q$rank]]
  basis <- setdiff(usable, aliased)
  out <- vapply(ev, function(j) {
    if (!(j %in% usable) || j %in% aliased) return(Inf)
    others <- setdiff(basis, j)
    if (!length(others)) return(1)
    res <- stats::lm.fit(Z[, others, drop = FALSE], Z[, j])$residuals
    tss <- sum(Z[, j]^2)
    1 / max(sum(res^2) / tss, .Machine$double.eps)
  }, numeric(1))
  stats::setNames(out, info$label[ev])
}

#' One caption line summarising event VIFs
#' @keywords internal
#' @noRd
.vif_note <- function(vif) {
  if (!length(vif)) return(character(0))
  fmt <- function(v) ifelse(is.infinite(v), "Inf (aliased)", formatC(v, format = "f", digits = 1))
  hi <- vif[vif > 5]
  paste0("VIF: variance inflation of each event regressor in the joint design (shown with its column",
         if (length(vif) > 8L) {
           fin <- vif[is.finite(vif)]
           sprintf("; median %s, max %s", fmt(stats::median(fin)), fmt(max(vif)))
         } else "",
         "). ",
         if (any(is.infinite(vif))) paste0(paste(names(vif)[is.infinite(vif)], collapse = ", "), ": Inf (aliased). ") else "",
         if (length(hi)) sprintf("%d above 5 (marked !) suggest collinearity.", length(hi))
         else "None above 5.")
}
