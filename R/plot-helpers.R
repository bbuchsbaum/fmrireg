## ============================================================================
## Shared plotting infrastructure
##
## Every fmrireg plot draws its theme, palettes, and regressor labels from
## here, so the design map, correlation map, regressor plots, and fit plots
## read as one family.
## ============================================================================

#' Okabe-Ito categorical palette (colour-blind safe), yellow dropped for
#' legibility on white. Conditions keep their hue by position in the design,
#' never by rank. Past eight conditions, hues are reused and identity must
#' come from labels.
#' @keywords internal
#' @noRd
.fmrireg_categorical <- function(n) {
  base <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
            "#E69F00", "#56B4E9", "#882255", "#555555")
  rep_len(base, n)
}

#' Diverging endpoints (ColorBrewer RdBu) with a neutral light-grey midpoint.
#' @keywords internal
#' @noRd
.fmrireg_diverging <- c(low = "#2166AC", mid = "#F4F4F4", high = "#B2182B")

#' Ink colours for non-data elements.
#' @keywords internal
#' @noRd
.fmrireg_ink <- c(primary = "#1F1F1F", secondary = "#595959",
                  muted = "#8C8C8C", rule = "#D9D9D9", faint = "#EDEDED")

#' The fmrireg ggplot2 theme
#'
#' A restrained theme used by all fmrireg plots: hairline major grid only,
#' left-aligned titles, and muted axis ink so the data carry the contrast.
#' Add it to any ggplot, or override it by adding your own theme afterwards.
#'
#' fmrireg plots use this theme by default. To use a different base theme for
#' every fmrireg plot, set `options(fmrireg.plot_theme = <theme>)`, for example
#' `options(fmrireg.plot_theme = ggplot2::theme_bw())`.
#'
#' @param base_size Base font size in points.
#' @param base_family Base font family.
#' @param grid Which major gridlines to keep: `"y"`, `"x"`, `"xy"`, or `"none"`.
#' @return A [ggplot2::theme()] object.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_fmrireg()
#' @export
theme_fmrireg <- function(base_size = 11, base_family = "",
                          grid = c("y", "x", "xy", "none")) {
  grid <- match.arg(grid)
  ink <- .fmrireg_ink
  hair <- ggplot2::element_line(colour = ink[["faint"]], linewidth = 0.35)
  blank <- ggplot2::element_blank()
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = ink[["primary"]]),
      plot.title = ggplot2::element_text(size = ggplot2::rel(1.15), face = "bold",
                                         margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(colour = ink[["secondary"]],
                                            size = ggplot2::rel(0.9),
                                            margin = ggplot2::margin(b = 8)),
      plot.caption = ggplot2::element_text(colour = ink[["muted"]],
                                           size = ggplot2::rel(0.8), hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      axis.title = ggplot2::element_text(colour = ink[["secondary"]],
                                         size = ggplot2::rel(0.9)),
      axis.text = ggplot2::element_text(colour = ink[["secondary"]],
                                        size = ggplot2::rel(0.85)),
      axis.ticks = blank,
      panel.grid.minor = blank,
      panel.grid.major.x = if (grid %in% c("x", "xy")) hair else blank,
      panel.grid.major.y = if (grid %in% c("y", "xy")) hair else blank,
      # one strip style for the family: panel titles bold and left-aligned,
      # row labels (strips switched to the left) plain, horizontal, right-aligned
      strip.text = ggplot2::element_text(colour = ink[["primary"]], face = "bold",
                                         size = ggplot2::rel(0.85), hjust = 0),
      strip.text.y.left = ggplot2::element_text(colour = ink[["primary"]],
                                                face = "plain", angle = 0, hjust = 1,
                                                size = ggplot2::rel(0.85)),
      strip.placement = "outside",
      legend.title = ggplot2::element_text(colour = ink[["secondary"]],
                                           size = ggplot2::rel(0.85)),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.85)),
      legend.key.size = ggplot2::unit(0.9, "lines"),
      panel.spacing = ggplot2::unit(0.9, "lines"),
      plot.margin = ggplot2::margin(8, 10, 8, 8)
    )
}

#' Base theme for a package plot, honouring `options(fmrireg.plot_theme)`
#' @keywords internal
#' @noRd
.plot_theme <- function(grid = "y") {
  user <- getOption("fmrireg.plot_theme")
  if (inherits(user, "theme")) user else theme_fmrireg(grid = grid)
}

#' Human-readable names for raw design-matrix column names (fallback path)
#'
#' Used when term metadata is unavailable. `base_poly2_block_1` ->
#' `poly2 - run 1`, `constant_1` -> `run 1`, `nuis#01_3` -> `nuis 3 - run 1`,
#' `rt_c_rt_c` -> `rt_c`, `condition_condition.faces` -> `faces`.
#' @keywords internal
#' @noRd
.pretty_regressor_names <- function(x) {
  dot <- paste0(" ", .glyph("dot"), " ")
  out <- x
  is_baseline <- grepl("^(base_|constant|nuis#)", out)
  for (i in which(!is_baseline)) {
    us <- gregexpr("_", out[i], fixed = TRUE)[[1]]
    if (us[1] < 0) next
    for (k in us) {
      term <- substr(out[i], 1, k - 1)
      rest <- substr(out[i], k + 1, nchar(out[i]))
      if (identical(rest, term)) {
        out[i] <- term
        break
      }
      if (grepl("^[^_]+\\.", rest) && startsWith(term, sub("\\..*$", "", rest))) {
        out[i] <- .strip_level_prefixes(rest)
        break
      }
    }
  }
  out <- sub("^base_([A-Za-z]+)([0-9]+)_block_([0-9]+)$", paste0("\\1\\2", dot, "run \\3"), out)
  out <- sub("^base_(.*)_block_([0-9]+)$", paste0("\\1", dot, "run \\2"), out)
  out <- sub("^constant_([0-9]+)$", "run \\1", out)
  out <- sub("^constant$", "intercept", out)
  out <- sub("^nuis#0*([0-9]+)_([0-9]+)$", paste0("nuis \\2", dot, "run \\1"), out)
  out
}

#' `condition.faces_task.a` -> `faces x a`; `rt_c` -> `rt_c`
#' @keywords internal
#' @noRd
.strip_level_prefixes <- function(x) {
  vapply(x, function(s) {
    parts <- strsplit(s, "_(?=[^_]+\\.)", perl = TRUE)[[1]]
    parts <- sub("^[^.]+\\.", "", parts)
    paste(parts, collapse = paste0(" ", .glyph("times"), " "))
  }, character(1), USE.NAMES = FALSE)
}

#' One row of metadata per design-matrix column
#'
#' The single source of truth for how plots name, group, and colour design
#' columns. Built from term metadata (`col_indices`, `conditions()`), not by
#' parsing column names.
#'
#' Multi-basis event terms are laid out condition-major in the design matrix
#' (all basis functions of the first condition, then the next condition),
#' which is how `fitted_hrf()` indexes them; labels follow that layout.
#'
#' @return data.frame with columns `column` (raw name), `label` (unique display
#'   label), `group` (term or baseline group), `kind` (`"event"`, `"drift"`,
#'   `"intercept"`, `"nuisance"`, `"other"`), `condition` (event condition
#'   label, without basis), `basis` (integer, 1 for single-basis), `nbasis`,
#'   and `run` (run index for run-specific columns, `NA` otherwise).
#' @keywords internal
#' @noRd
.column_info <- function(fmrimod, DM = as.matrix(design_matrix(fmrimod))) {
  cols <- colnames(DM)
  p <- length(cols)
  info <- data.frame(column = cols, label = .pretty_regressor_names(cols),
                     group = "other", kind = "other", condition = NA_character_,
                     basis = 1L, nbasis = 1L, run = NA_integer_,
                     stringsAsFactors = FALSE)

  # event terms
  em <- fmrimod$event_model
  ci <- NULL
  if (!is.null(em)) {
    ev_dm <- design_matrix(em)
    ci <- attr(ev_dm, "col_indices")
    ev_cols <- colnames(ev_dm)
    eterms <- terms(em)
  }
  if (!is.null(ci)) {
    for (tn in names(ci)) {
      idx <- match(ev_cols[ci[[tn]]], cols)
      idx <- idx[!is.na(idx)]
      if (!length(idx)) next
      info$group[idx] <- tn
      info$kind[idx] <- "event"
      conds <- tryCatch(fmridesign::conditions(eterms[[tn]]), error = function(e) NULL)
      nb <- if (length(conds)) length(idx) / length(conds) else NA
      if (length(conds) && isTRUE(nb == round(nb))) {
        lab <- .strip_level_prefixes(conds)
        info$condition[idx] <- rep(lab, each = nb)
        info$basis[idx] <- rep(seq_len(nb), times = length(lab))
        info$nbasis[idx] <- as.integer(nb)
      } else {
        info$condition[idx] <- info$label[idx]
      }
    }
  }
  # a condition label used by two terms is qualified with its term
  ev <- info$kind == "event"
  clash <- character(0)
  if (any(ev)) {
    cond_terms <- tapply(info$group[ev], info$condition[ev], function(g) length(unique(g)))
    clash <- names(cond_terms)[cond_terms > 1]
  }
  qual <- ev & info$condition %in% clash
  info$condition[qual] <- paste0(info$group[qual], ": ", info$condition[qual])
  info$label[ev] <- ifelse(info$nbasis[ev] > 1,
                           paste0(info$condition[ev], " [b", info$basis[ev], "]"),
                           info$condition[ev])

  # baseline terms
  bl_terms <- if (!is.null(fmrimod$baseline_model)) terms(fmrimod$baseline_model) else list()
  bl_kind <- c(drift = "drift", block = "intercept", nuisance = "nuisance")
  for (tn in names(bl_terms)) {
    tc <- colnames(design_matrix(bl_terms[[tn]]))
    idx <- match(tc, cols)
    idx <- idx[!is.na(idx)]
    k <- if (tn %in% names(bl_kind)) bl_kind[[tn]] else "other"
    info$group[idx] <- if (k == "other") tn else k
    info$kind[idx] <- k
  }

  # run membership: a column is run-specific when it is zero outside one run
  runs <- .run_lengths(fmrimod)
  if (!is.null(runs) && length(runs) > 1L && sum(runs) == nrow(DM)) {
    run_id <- rep(seq_along(runs), runs)
    nz <- DM != 0
    for (j in seq_len(p)) {
      r <- unique(run_id[nz[, j]])
      if (length(r) == 1L) info$run[j] <- r
    }
  }
  info$label <- make.unique(info$label, sep = " #")
  info
}

#' Named colours for event conditions, in design order
#'
#' Every fit and design plot takes condition colours from here, so "faces" is
#' the same hue in the regressor, HRF, and time-course plots.
#' @keywords internal
#' @noRd
.condition_colours <- function(info) {
  conds <- unique(info$condition[info$kind == "event"])
  stats::setNames(.fmrireg_categorical(length(conds)), conds)
}

#' Run boundaries (in scans) for a model's sampling frame
#' @keywords internal
#' @noRd
.run_lengths <- function(fmrimod) {
  sf <- .model_sframe(fmrimod)
  if (is.null(sf)) {
    return(NULL)
  }
  as.integer(fmrihrf::blocklens(sf))
}

#' @keywords internal
#' @noRd
.model_sframe <- function(x) {
  x$event_model$sampling_frame %||% x$baseline_model$sampling_frame %||%
    x$sampling_frame
}

#' Clamp out-of-range values to the scale limits (scales::squish without the dep)
#' @keywords internal
#' @noRd
.squish <- function(x, range = c(0, 1), ...) {
  pmin(pmax(x, range[1]), range[2])
}

#' Format numbers for axis labels with a consistent precision
#' @keywords internal
#' @noRd
.axis_num <- function(x) {
  ok <- is.finite(x)
  out <- rep("", length(x))
  if (!any(ok)) return(out)
  # one number of decimals for the whole axis (no "0.5" next to "0.25")
  dec <- vapply(x[ok], function(v) {
    for (d in 0:4) if (abs(round(v, d) - v) < 1e-8 * max(1, abs(v))) return(d)
    4L
  }, numeric(1))
  out[ok] <- formatC(x[ok], format = "f", digits = max(dec))
  out
}

#' Typographic glyphs, built from code points so R sources stay ASCII.
#' Restricted to Latin-1, which every graphics device can draw.
#' @keywords internal
#' @noRd
.glyph <- function(name) {
  # Latin-1 only: these render on every graphics device, including pdf()
  cp <- c(dot = 183L, times = 215L, sq = 178L, pm = 177L)
  intToUtf8(cp[[name]])
}

#' Wrap caption notes, one note per paragraph; NULL when there are none
#' @keywords internal
#' @noRd
.wrap_notes <- function(notes, width = 100) {
  notes <- trimws(notes[nzchar(notes)])
  if (!length(notes)) return(NULL)
  # balanced wrapping: split a long note into lines of similar length, so no
  # line is left holding a single orphaned word
  wrap1 <- function(x) {
    n_lines <- ceiling(nchar(x) / width)
    if (n_lines <= 1L) return(x)
    strwrap(x, width = ceiling(nchar(x) / n_lines) + 8L)
  }
  paste(unlist(lapply(notes, wrap1)), collapse = "\n")
}

#' Two or three readable y breaks for a small-multiple row: 0 plus one
#' rounded value per sign that the data reach
#' @keywords internal
#' @noRd
.row_breaks <- function(lims) {
  lims <- lims[is.finite(lims)]
  if (length(lims) < 2L) return(numeric(0))
  m <- max(abs(lims))
  if (m == 0) return(0)
  b <- signif(0.8 * m, 1)
  out <- c(if (lims[1] < -0.25 * b) -b, 0, if (lims[2] > 0.25 * b) b)
  out[out >= lims[1] & out <= lims[2]]
}
