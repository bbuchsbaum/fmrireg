## ============================================================================
## Plots for fitted fmri_lm models
## ============================================================================

#' Plot a fitted fMRI linear model
#'
#' Visual summaries of an \code{fmri_lm} fit. Five views are available:
#'
#' \describe{
#'   \item{\code{"betas"}}{With no \code{voxel}, the distribution of
#'     t-statistics across voxels for each event regressor, drawn as a sina
#'     plot (points spread in proportion to their density), with reference
#'     lines at the two-sided p < .001 (uncorrected) critical value and the
#'     share of voxels beyond it in each tail. With \code{voxel}, a coefficient
#'     plot of estimates with confidence intervals, one panel per event term
#'     so that regressors in different units (a condition and a parametric
#'     modulator) do not share an axis.}
#'   \item{\code{"contrasts"}}{As \code{"betas"}, for the fitted contrasts.}
#'   \item{\code{"hrf"}}{The estimated hemodynamic response for each condition:
#'     the term's HRF basis weighted by the fitted coefficients. For one voxel,
#'     the band is a pointwise confidence interval from the coefficient
#'     covariance, and a curve whose interval includes zero at its peak is
#'     drawn faded and dashed and labelled "n.s.". For several voxels, the line
#'     is the mean across voxels, the band a t-based confidence interval of
#'     that mean (voxels treated as independent), and up to 200 voxels are
#'     drawn faintly underneath. Peak latency is reported only for terms with
#'     more than one basis function; with a single basis the shape, and so the
#'     latency, is fixed by the HRF rather than estimated.}
#'   \item{\code{"timecourse"}}{Observed signal and model fit for one voxel,
#'     with event onsets in a lane per condition and residuals in a lower panel
#'     drawn on the same vertical scale.}
#'   \item{\code{"residuals"}}{Autocorrelation of the residuals at lags 1 to
#'     10 scans. By default this summarises up to 200 voxels spread evenly
#'     through the data (median and 10th-90th percentile per lag); with
#'     \code{voxel}, it shows those voxels. Least-squares residuals are
#'     negatively autocorrelated even when the noise is white, so the 95\%
#'     white-noise range is centred on the autocorrelation this design
#'     induces, not on zero. For fits with an AR noise model the residuals
#'     are also shown after prewhitening with the fitted AR coefficients,
#'     which is the check that the noise model removed the serial
#'     correlation. The subtitle reports a portmanteau test of lags 1-10 that
#'     accounts for both the design-induced mean and the correlation between
#'     lags that regression introduces (a generalised Ljung-Box test).}
#' }
#'
#' When \code{voxel} is \code{NULL} for the \code{"hrf"} and
#' \code{"timecourse"} views, the voxel with the largest absolute event
#' t-statistic is shown and named in the subtitle.
#'
#' Critical values and intervals use the degrees of freedom the fit used for
#' inference (\code{fit$result$df$inference}), which differ from the nominal
#' residual df for AR, robust, or effective-df corrected fits.
#'
#' In the \code{"timecourse"} and \code{"residuals"} views the fitted values
#' are an ordinary least squares fit of the full design (events plus baseline)
#' to the voxel: run by run for runwise fits, jointly otherwise. For OLS fits
#' these are the model's own fitted values; for AR or robust fits they are the
#' unweighted counterpart, on the scale of the data. For runwise fits the
#' fitted values use run-specific coefficients, whereas the estimates view
#' shows coefficients pooled across runs. The reported partial R-squared is
#' \eqn{1 - RSS_{full} / RSS_{baseline}}: the share of the variance left by
#' the baseline model that the event regressors explain.
#'
#' @param object,x An \code{fmri_lm} object.
#' @param type One of \code{"betas"} (default), \code{"contrasts"},
#'   \code{"hrf"}, \code{"timecourse"}, \code{"residuals"}. The first two
#'   name the parameter family, as in \code{\link{stats}()} and
#'   \code{\link[=coef.fmri_lm]{coef()}}. \code{"estimates"} is accepted
#'   as a synonym of \code{"betas"}.
#' @param voxel Optional integer voxel index (or indices, for
#'   \code{"betas"}, \code{"contrasts"}, \code{"hrf"} and
#'   \code{"residuals"}).
#' @param level Confidence level for intervals and bands.
#' @param sample_at Time points (s) at which to evaluate the HRF, for
#'   \code{type = "hrf"}.
#' @param direct_labels Logical; for \code{type = "hrf"}, label each curve
#'   next to it (default) or, if \code{FALSE}, use a legend. A legend reads
#'   better in small or many-panel figures.
#' @param ... Passed on from \code{plot()} to \code{autoplot()}.
#' @return A ggplot2 object. \code{plot()} prints it and returns it invisibly.
#' @examples
#' \donttest{
#' set.seed(1)
#' ev <- data.frame(onset = seq(10, 180, by = 12), run = 1,
#'                  condition = factor(rep(c("a", "b"), length.out = 15)))
#' Y <- matrix(rnorm(100 * 5), 100, 5)
#' dset <- matrix_frame(Y, TR = 2, run_length = 100, event_table = ev)
#' fit <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset)
#' ggplot2::autoplot(fit)
#' ggplot2::autoplot(fit, type = "timecourse", voxel = 1)
#' }
#' @export
#' @method autoplot fmri_lm
#' @rdname autoplot.fmri_lm
autoplot.fmri_lm <- function(object,
                             type = c("betas", "contrasts", "hrf", "timecourse",
                                      "residuals"),
                             voxel = NULL, level = 0.95,
                             sample_at = seq(0, 20, by = 0.25),
                             direct_labels = TRUE, ...) {
  # "estimates" was this view's name before the accessor family was renamed
  # (#217); it stays a silent synonym, since a plot is labelled and cannot be
  # mistaken for the estimates themselves.
  type <- match.arg(type[1], c("betas", "contrasts", "hrf", "timecourse",
                               "residuals", "estimates"))
  if (identical(type, "estimates")) type <- "betas"
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be a single number between 0 and 1.", call. = FALSE)
  }
  if (!is.logical(direct_labels) || length(direct_labels) != 1L || is.na(direct_labels)) {
    stop("`direct_labels` must be TRUE or FALSE.", call. = FALSE)
  }
  switch(type,
    betas      = .plot_fit_stats(object, "betas", voxel, level),
    contrasts  = .plot_fit_stats(object, "contrasts", voxel, level),
    hrf        = .plot_fit_hrf(object, voxel, sample_at, level, direct_labels),
    timecourse = .plot_fit_timecourse(object, voxel),
    residuals  = .plot_fit_residual_acf(object, voxel)
  )
}

#' @export
#' @rdname autoplot.fmri_lm
plot.fmri_lm <- function(x, type = c("betas", "contrasts", "hrf", "timecourse",
                                     "residuals"),
                         voxel = NULL, ...) {
  p <- autoplot.fmri_lm(x, type = type, voxel = voxel, ...)
  print(p)
  invisible(p)
}

## ---------------------------------------------------------------------------

#' Voxel with the largest absolute event t-statistic
#' @keywords internal
#' @noRd
.peak_voxel <- function(fit) {
  tt <- as.matrix(stats(fit, type = "betas"))
  score <- suppressWarnings(apply(abs(tt), 1, max, na.rm = TRUE))
  score[!is.finite(score)] <- -Inf
  which.max(score)
}

#' @keywords internal
#' @noRd
.check_voxels <- function(fit, voxel, single = FALSE) {
  nvox <- nrow(coef(fit))
  ok <- is.numeric(voxel) && length(voxel) >= 1L && !anyNA(voxel) &&
    all(voxel >= 1) && all(voxel <= nvox) && all(voxel == round(voxel))
  if (!ok) {
    stop(sprintf("`voxel` must contain integer indices between 1 and %d.", nvox),
         call. = FALSE)
  }
  voxel <- unique(as.integer(voxel))
  if (single && length(voxel) != 1L) {
    stop("`voxel` must be a single index for this view.", call. = FALSE)
  }
  voxel
}

#' Degrees of freedom used for inference, per voxel
#' @keywords internal
#' @noRd
.fit_inference_df <- function(fit, nvox) {
  d <- fit$result$df$inference %||% fit$result$rdf %||% Inf
  rep_len(as.numeric(d), nvox)
}

#' @keywords internal
#' @noRd
.df_text <- function(d) {
  d <- d[is.finite(d)]
  if (!length(d)) return("Inf")
  r <- range(round(d, 1))
  if (diff(r) < 0.5) format(round(r[1])) else paste0(r[1], "-", r[2])
}

#' One-line description of how the fit was estimated
#' @keywords internal
#' @noRd
.fit_strategy_text <- function(fit) {
  nr <- length(.run_lengths(fit$model) %||% 1L)
  if (identical(fit$strategy, "runwise") && nr > 1L) {
    sprintf("Runwise fit: coefficients estimated per run and pooled across %d runs.", nr)
  } else {
    "Fit jointly across all scans."
  }
}

#' Sina offsets: spread points vertically in proportion to their local density
#' @keywords internal
#' @noRd
.sina_offset <- function(x, width = 0.38) {
  n <- length(x)
  if (n < 3L || stats::sd(x) == 0) return(rep(0, n))
  d <- stats::density(x, n = 512)
  dens <- stats::approx(d$x, d$y, xout = x)$y
  # deterministic low-discrepancy sequence instead of an RNG
  u <- ((seq_len(n) * 0.6180339887) %% 1)[order(order(x))]
  (u - 0.5) * 2 * width * dens / max(dens)
}

#' @keywords internal
#' @noRd
.pct <- function(x) {
  ifelse(x == 0, "0%", ifelse(x < 0.01, "<1%", sprintf("%.0f%%", 100 * x)))
}

#' Push label positions apart so neighbours are at least `gap` apart
#' @keywords internal
#' @noRd
.dodge_labels <- function(y, gap) {
  if (length(y) < 2L) return(y)
  o <- order(y, decreasing = TRUE)
  z <- y[o]
  for (k in seq_along(z)[-1]) {
    if (z[k - 1] - z[k] < gap) z[k] <- z[k - 1] - gap
  }
  out <- y
  out[o] <- z
  out
}

#' @keywords internal
#' @noRd
.plot_fit_stats <- function(fit, type, voxel, level) {
  tt <- tryCatch(as.matrix(stats(fit, type = type)), error = function(e) NULL)
  if (is.null(tt) || !length(tt) || ncol(tt) == 0L) {
    stop(sprintf("This fit has no %s to plot.", type), call. = FALSE)
  }
  se <- as.matrix(standard_error(fit, type = type))
  # Take the estimates from coef(), which shares the voxels x terms layout of
  # stats() and standard_error(), rather than rebuilding them as t * se (t is
  # set to 0 where se is numerically zero).
  est <- as.matrix(coef(fit, type = type))
  stopifnot(identical(dim(est), dim(tt)), identical(colnames(est), colnames(tt)))
  nvox <- nrow(tt)
  dfs <- .fit_inference_df(fit, nvox)
  what <- if (type == "betas") "Regressor" else "Contrast"

  if (type == "betas") {
    info <- .column_info(fit$model)
    info <- info[match(colnames(tt), info$column), ]
    labels <- ifelse(is.na(info$label), .pretty_regressor_names(colnames(tt)), info$label)
    groups <- ifelse(is.na(info$group), "", info$group)
  } else {
    labels <- colnames(tt)
    groups <- rep("contrasts", length(labels))
  }
  labels <- make.unique(labels, sep = " #")
  ink <- .fmrireg_ink

  if (is.null(voxel)) {
    crit <- stats::qt(1 - 0.001 / 2, stats::median(dfs))
    df <- data.frame(
      term = factor(rep(labels, each = nvox), levels = rev(labels)),
      t    = as.vector(tt)
    )
    df <- df[is.finite(df$t), ]
    estimable <- colSums(is.finite(tt)) > 0
    pos <- colMeans(tt > crit, na.rm = TRUE)
    neg <- colMeans(tt < -crit, na.rm = TRUE)
    ann <- data.frame(term = factor(labels, levels = rev(labels)),
                      pos = ifelse(estimable, .pct(pos), "aliased"),
                      neg = ifelse(estimable, .pct(neg), ""))
    xr <- range(c(df$t, -crit, crit), na.rm = TRUE)
    show_neg <- any(neg > 0, na.rm = TRUE)
    x_pos <- xr[2] + 0.11 * diff(xr)
    x_neg <- xr[2] + 0.21 * diff(xr)
    x_end <- if (show_neg) x_neg else x_pos

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$t, y = .data$term)) +
      ggplot2::geom_vline(xintercept = 0, colour = ink[["rule"]], linewidth = 0.5) +
      ggplot2::geom_vline(xintercept = c(-crit, crit), colour = ink[["muted"]],
                          linewidth = 0.3, linetype = "22")
    if (nvox <= 5000) {
      df$y <- as.numeric(df$term)
      for (lv in unique(df$y)) {
        sel <- df$y == lv
        df$y[sel] <- lv + .sina_offset(df$t[sel])
      }
      alpha <- max(0.15, min(0.8, 12 / sqrt(nvox)))
      p <- p + ggplot2::geom_point(data = df, ggplot2::aes(y = .data$y),
                                   size = 1.3, alpha = alpha, stroke = 0,
                                   colour = ink[["primary"]])
    } else {
      p <- p + ggplot2::geom_bin_2d(ggplot2::aes(alpha = ggplot2::after_stat(.data$ncount)),
                                    fill = ink[["primary"]],
                                    bins = c(120, length(labels)), drop = TRUE) +
        ggplot2::scale_alpha_continuous(range = c(0.05, 1), guide = "none")
    }
    head_y <- length(labels) + 0.55
    p +
      ggplot2::geom_text(data = ann, ggplot2::aes(x = x_pos, label = .data$pos),
                         hjust = 1, size = 3, colour = ink[["secondary"]]) +
      (if (show_neg) ggplot2::geom_text(data = ann, ggplot2::aes(x = x_neg, label = .data$neg),
                                        hjust = 1, size = 3, colour = ink[["secondary"]])) +
      ggplot2::annotate("text", x = c(x_pos, x_neg)[seq_len(1 + show_neg)], y = head_y,
                        label = sprintf(c("t > %.2f", "t < %.2f"), c(crit, -crit))[seq_len(1 + show_neg)],
                        hjust = 1, vjust = 0,
                        size = 2.9, colour = ink[["secondary"]]) +
      ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.6, 0.9))) +
      ggplot2::coord_cartesian(xlim = c(xr[1], x_end + 0.01 * diff(xr)), clip = "off") +
      ggplot2::labs(
        title = sprintf("%s t-statistics across %s voxels", what,
                        format(nvox, big.mark = ",")),
        subtitle = sprintf("Dashed: |t| = %.2f (p < .001, two-sided). Right: %% of voxels beyond.",
                           crit),
        caption = .wrap_notes(c(
          paste0(sprintf("df = %s. Under the null, 0.05%% of voxels fall beyond each line.",
                         .df_text(dfs)),
                 if (!show_neg) sprintf(" No voxel has t below %.2f.", -crit) else ""),
          .fit_strategy_text(fit))),
        x = "t-statistic", y = NULL) +
      .plot_theme(grid = "x") +
      ggplot2::theme(plot.margin = ggplot2::margin(8, 16, 8, 8))
  } else {
    voxel <- .check_voxels(fit, voxel)
    df <- data.frame(
      voxel = factor(rep(paste("voxel", voxel), times = length(labels)),
                     levels = paste("voxel", voxel)),
      term  = factor(rep(labels, each = length(voxel)), levels = rev(labels)),
      group = factor(rep(groups, each = length(voxel)), levels = unique(groups)),
      est   = as.vector(est[voxel, , drop = FALSE]),
      se    = as.vector(se[voxel, , drop = FALSE]),
      t     = as.vector(tt[voxel, , drop = FALSE]),
      dfi   = rep(dfs[voxel], times = length(labels))
    )
    q <- stats::qt(1 - (1 - level) / 2, df$dfi)
    df$lo <- df$est - q * df$se
    df$hi <- df$est + q * df$se
    pv <- 2 * stats::pt(-abs(df$t), df$dfi)
    df$txt <- ifelse(is.finite(df$t),
                     sprintf("t = %.1f, %s", df$t,
                             ifelse(pv < 0.001, "p < .001",
                                    sprintf("p = %s", sub("^0", "", sprintf("%.3f", pv))))),
                     "not estimable (aliased)")
    df$est[!is.finite(df$est)] <- NA_real_
    df <- df[!(is.na(df$est) & df$txt != "not estimable (aliased)"), ]
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$est, y = .data$term)) +
      ggplot2::geom_vline(xintercept = 0, colour = ink[["rule"]], linewidth = 0.5) +
      ggplot2::geom_linerange(data = df[is.finite(df$est), ],
                              ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
                              linewidth = 0.8, colour = ink[["secondary"]]) +
      ggplot2::geom_point(data = df[is.finite(df$est), ], size = 2.4,
                          colour = ink[["primary"]]) +
      ggplot2::geom_text(ggplot2::aes(x = Inf, label = .data$txt), hjust = -0.08,
                         size = 2.9, colour = ink[["secondary"]]) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(
        title = if (length(voxel) == 1L) sprintf("%s estimates, voxel %d", what, voxel)
                else sprintf("%s estimates", what),
        subtitle = sprintf("Point: estimate. Line: %g%% CI. Each panel has its own units.",
                           100 * level),
        caption = .wrap_notes(c(sprintf("df = %s.", .df_text(dfs[voxel])),
                                .fit_strategy_text(fit))),
        x = expression(paste("Estimate (", beta, ")")), y = NULL) +
      .plot_theme(grid = "x") +
      ggplot2::theme(plot.margin = ggplot2::margin(8, 96, 8, 8))
    if (length(voxel) > 1L) {
      p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$group, .data$voxel),
                                   ncol = length(voxel), scales = "free")
    } else if (length(unique(groups)) > 1L) {
      # one panel per term, each with its own x axis because terms differ in units
      p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$group), ncol = 1L, scales = "free")
    }
    p
  }
}

#' HRF basis matrix and coefficient blocks for each event term
#'
#' Mirrors `fitted_hrf()`: the term's HRF evaluated at `sample_at` (one column
#' per basis function) and, per condition, the positions of that condition's
#' coefficients among the event coefficients (condition-major layout).
#' @keywords internal
#' @noRd
.hrf_terms <- function(fit, sample_at) {
  em <- fit$model$event_model
  eterms <- terms(em)
  ci <- attr(design_matrix(em), "col_indices")
  info <- .column_info(fit$model)
  ev_names <- colnames(design_matrix(em))
  out <- list()
  for (tn in names(eterms)) {
    et <- eterms[[tn]]
    if (inherits(et, "feature_term") || is.null(ci[[tn]])) next
    spec <- attr(et, "hrfspec") %||% et$hrfspec
    hrf_fun <- if (!is.null(spec$hrf)) spec$hrf else fmrihrf::HRF_SPMG1
    G <- as.matrix(hrf_fun(sample_at))
    ind <- as.integer(ci[[tn]])
    nb <- ncol(G)
    ncond <- length(ind) / nb
    if (ncond != round(ncond)) next
    conds <- info$condition[match(ev_names[ind], info$column)]
    blocks <- lapply(seq_len(ncond), function(k) ind[(k - 1) * nb + seq_len(nb)])
    names(blocks) <- conds[(seq_len(ncond) - 1) * nb + 1]
    out[[tn]] <- list(G = G, blocks = blocks,
                      hrf_name = attr(hrf_fun, "name") %||% "HRF")
  }
  out
}

#' Unscaled coefficient covariance that applies to one voxel
#'
#' A shared `cov.unscaled` when the fit has one; otherwise the covariance of
#' the voxel's cluster for fits that carry one per cluster
#' (`covariance_by_cluster` indexed by `cluster_voxels`, e.g. by_cluster
#' sketch fits). `NULL` when neither is available.
#' @keywords internal
#' @noRd
.voxel_cov_unscaled <- function(fit, voxel) {
  cu <- fit$result$cov.unscaled
  if (is.matrix(cu)) return(cu)
  covs <- fit$result$covariance_by_cluster
  groups <- fit$result$cluster_voxels
  if (is.null(covs) || is.null(groups)) return(NULL)
  hit <- which(vapply(groups, function(g) voxel %in% g, logical(1)))
  if (length(hit) != 1L) return(NULL)
  as.matrix(covs[[hit]])
}

#' @keywords internal
#' @noRd
.plot_fit_hrf <- function(fit, voxel, sample_at, level, direct_labels = TRUE) {
  auto <- is.null(voxel)
  voxel <- if (auto) .peak_voxel(fit) else .check_voxels(fit, voxel)
  terms_h <- .hrf_terms(fit, sample_at)
  if (!length(terms_h)) {
    stop("This fit has no event terms with an HRF to plot.", call. = FALSE)
  }
  # coef() and standard_error() are both voxels x terms.
  beta <- coef(fit)
  se <- as.matrix(standard_error(fit, type = "betas"))
  stopifnot(identical(dim(beta), dim(se)))
  nvox <- nrow(beta)
  dfs <- .fit_inference_df(fit, nvox)
  single <- length(voxel) == 1L
  cu <- if (single) .voxel_cov_unscaled(fit, voxel) else NULL
  ink <- .fmrireg_ink
  max_lines <- 200L
  line_vox <- if (length(voxel) > max_lines) {
    voxel[unique(round(seq(1, length(voxel), length.out = max_lines)))]
  } else {
    voxel
  }

  curves <- list()
  voxlines <- list()
  aliased <- character(0)
  for (tn in names(terms_h)) {
    G <- terms_h[[tn]]$G
    for (cond in names(terms_h[[tn]]$blocks)) {
      idx <- terms_h[[tn]]$blocks[[cond]]
      B <- beta[voxel, idx, drop = FALSE]
      # a condition with no estimable coefficient is not drawn at all: a flat
      # zero line would look like an estimated null response
      if (!any(is.finite(B))) {
        aliased <- c(aliased, cond)
        next
      }
      B[!is.finite(B)] <- 0
      H <- G %*% t(B)                           # time x voxels
      if (single) {
        m <- drop(H)
        half <- NA_real_
        s_j <- se[voxel, idx]
        if (is.matrix(cu) && all(idx <= nrow(cu)) && all(is.finite(s_j))) {
          # V = D cor(cu) D reproduces the reported standard errors exactly
          Cj <- suppressWarnings(stats::cov2cor(cu[idx, idx, drop = FALSE]))
          if (all(is.finite(Cj))) {
            V <- (s_j %o% s_j) * Cj
            half <- stats::qt(1 - (1 - level) / 2, dfs[voxel]) *
              sqrt(pmax(rowSums((G %*% V) * G), 0))
          }
        }
      } else {
        m <- rowMeans(H)
        half <- stats::qt(1 - (1 - level) / 2, ncol(H) - 1) *
          apply(H, 1, stats::sd) / sqrt(ncol(H))
        keep <- match(line_vox, voxel)
        voxlines[[length(voxlines) + 1L]] <- data.frame(
          term = tn, condition = cond, time = rep(sample_at, length(keep)),
          value = as.vector(H[, keep, drop = FALSE]), vox = rep(line_vox, each = nrow(H)))
      }
      pk <- which.max(abs(m))
      lo <- m - half
      hi <- m + half
      sig <- if (all(is.na(half))) TRUE else isTRUE(lo[pk] > 0 || hi[pk] < 0)
      curves[[length(curves) + 1L]] <- data.frame(
        term = tn, condition = cond, time = sample_at, mean = m, lo = lo, hi = hi,
        sig = sig, nbasis = ncol(G), peak_t = sample_at[pk], peak_v = m[pk])
    }
  }
  if (!length(curves)) {
    stop("No condition in this fit has an estimable response for the selected voxel(s).",
         call. = FALSE)
  }
  df <- do.call(rbind, curves)
  conds <- unique(df$condition)
  df$condition <- factor(df$condition, levels = conds)
  # terms whose every condition is aliased have no curves and no panel
  df$term <- factor(df$term, levels = intersect(names(terms_h), unique(df$term)))
  df$sig <- factor(df$sig, levels = c(TRUE, FALSE))
  info <- .column_info(fit$model)
  cols <- .condition_colours(info)
  cols <- cols[names(cols) %in% conds]
  missing <- setdiff(conds, names(cols))
  cols <- c(cols, stats::setNames(rep(ink[["secondary"]], length(missing)), missing))

  # direct labels just right of each curve's peak, dodged vertically, on a
  # translucent white ground so they stay legible where they cross a line
  labs <- unique(df[, c("term", "condition", "peak_t", "peak_v", "sig", "nbasis")])
  labs$txt <- ifelse(labs$sig == "FALSE", paste(labs$condition, "(n.s.)"),
                     ifelse(labs$nbasis > 1, sprintf("%s  %.1f s", labs$condition, labs$peak_t),
                            as.character(labs$condition)))
  labs$x <- labs$peak_t + 0.04 * diff(range(sample_at))
  labs$y <- labs$peak_v
  for (tn in levels(df$term)) {
    sel <- which(labs$term == tn)
    placed <- .place_curve_labels(df[df$term == tn, ], labs[sel, ])
    labs$x[sel] <- placed$x
    labs$y[sel] <- placed$y
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$mean,
                                        colour = .data$condition)) +
    ggplot2::geom_hline(yintercept = 0, colour = ink[["rule"]], linewidth = 0.5)
  if (!single && length(voxlines)) {
    vl <- do.call(rbind, voxlines)
    vl$condition <- factor(vl$condition, levels = conds)
    vl$term <- factor(vl$term, levels = levels(df$term))
    p <- p + ggplot2::geom_line(data = vl,
                                ggplot2::aes(y = .data$value,
                                             group = interaction(.data$condition, .data$vox)),
                                alpha = 0.12, linewidth = 0.3)
  }
  if (any(is.finite(df$lo))) {
    # where bands of several conditions overlap, fills blend into colours that
    # are in no key, so keep the fill faint and trace each band's edges
    crowded <- tapply(df$condition, df$term, function(x) length(unique(x)) > 1L)
    band_alpha <- if (any(crowded)) 0.07 else 0.15
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi,
                                               fill = .data$condition),
                                  colour = NA, alpha = band_alpha) +
      ggplot2::scale_fill_manual(values = cols, guide = "none")
    # trace band edges only for a single voxel; with several voxels the faint
    # per-voxel lines already show the spread
    if (any(crowded) && single) {
      edges <- rbind(transform(df, edge = "lo", v = df$lo), transform(df, edge = "hi", v = df$hi))
      p <- p + ggplot2::geom_line(data = edges,
                                  ggplot2::aes(y = .data$v,
                                               group = interaction(.data$condition, .data$edge)),
                                  linewidth = 0.3, linetype = "22", alpha = 0.7)
    }
  }
  p <- p +
    ggplot2::geom_line(ggplot2::aes(linetype = .data$sig, alpha = .data$sig),
                       linewidth = 0.8) +
    ggplot2::scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"),
                                   guide = "none", drop = FALSE) +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.55),
                                guide = "none", drop = FALSE) +
    ggplot2::geom_point(data = labs, ggplot2::aes(x = .data$peak_t, y = .data$peak_v),
                        size = 1.6, show.legend = FALSE) +
    (if (direct_labels) .halo_text(labs, ggplot2::aes(x = .data$x, y = .data$y, label = .data$txt))) +
    ggplot2::scale_colour_manual(values = cols, guide = if (direct_labels) "none" else "legend",
                                 name = NULL,
                                 labels = if (direct_labels) ggplot2::waiver() else
                                   stats::setNames(labs$txt[match(names(cols), labs$condition)], names(cols))) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(labels = .axis_num,
                                expand = ggplot2::expansion(mult = c(0.08, 0.12)))

  if (length(unique(df$term)) > 1L) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$term), scales = "free_y")
  }

  basis_txt <- paste(unique(vapply(terms_h, `[[`, "", "hrf_name")), collapse = "/")
  multi <- any(df$nbasis > 1)
  sub <- if (single) {
    sprintf("Voxel %d%s, %s basis. Band: %g%% CI.", voxel,
            if (auto) " (largest |t|)" else "", basis_txt, 100 * level)
  } else {
    sprintf("Mean of %d voxels, %s basis. Band: %g%% CI of the mean.",
            length(voxel), basis_txt, 100 * level)
  }
  notes <- c(
    if (single) "Dashed, faded curve: its interval includes 0 at the peak (n.s.).",
    if (!single) sprintf("Faint lines: %s. Voxels are treated as independent.",
                         if (length(voxel) > max_lines) sprintf("%d of the voxels", max_lines)
                         else "each voxel"),
    if (multi) "Dots mark each peak; labels give its estimated latency."
    else "One basis function fixes the response shape, so peak latency is not estimated.",
    if (length(aliased)) sprintf("Not estimable (aliased), not drawn: %s.",
                                 paste(unique(aliased), collapse = ", "))
  )
  p + ggplot2::labs(title = "Estimated hemodynamic response", subtitle = sub,
                    caption = .wrap_notes(notes, width = 120),
                    x = "Time since onset (s)",
                    y = "Response (a.u. per unit regressor)") +
    .plot_theme() +
    ggplot2::theme(legend.position = if (direct_labels) "none" else "top",
                   legend.justification = "left")
}

#' Place one direct label per curve where no other curve or label passes
#'
#' Puts each label beside its curve's peak dot (right, else left), nudging it
#' up or down only as far as needed for its text box to clear the other
#' conditions' mean curves and the labels already placed. Widths
#' are estimated from the character count in data units, which is adequate
#' at typical figure sizes.
#' @keywords internal
#' @noRd
.place_curve_labels <- function(curves, labs, char_frac = 0.0145, clear_frac = 0.045) {
  xr <- range(curves$time)
  yr <- range(c(curves$mean, curves$lo, curves$hi), na.rm = TRUE)
  dx <- diff(xr)
  dy <- diff(yr)
  if (dx <= 0 || dy <= 0) return(labs[, c("x", "y")])
  h <- clear_frac * dy
  gap <- 0.02 * dx
  boxes <- list()
  out <- data.frame(x = labs$x, y = labs$y)
  ord <- order(-abs(labs$peak_v))
  for (i in ord) {
    others <- curves[curves$condition != labs$condition[i], ]
    own <- curves[curves$condition == labs$condition[i], ]
    w <- nchar(labs$txt[i]) * char_frac * dx
    # beside the peak dot first (right, then left), then nudged up or down;
    # only other conditions' mean curves and other labels are obstacles
    cand <- expand.grid(side = c(1, -1), lift = c(0, 1, -1, 2, -2, 3, -3),
                        push = c(0, 0.05, 0.1))
    placed <- FALSE
    for (k in seq_len(nrow(cand))) {
      x0 <- if (cand$side[k] > 0) {
        labs$peak_t[i] + gap + cand$push[k] * dx
      } else {
        labs$peak_t[i] - gap - cand$push[k] * dx - w
      }
      if (x0 < xr[1] || x0 + w > xr[2]) next
      # text sits just above the dot so its own falling curve passes below it
      yc <- labs$peak_v[i] + (cand$lift[k] + 0.9) * h
      if (yc > yr[2] + 2 * h || yc < yr[1] - h) next
      span <- others$time >= x0 & others$time <= x0 + w
      own_span <- own$time >= x0 & own$time <= x0 + w
      hit_curve <- any(abs(others$mean[span] - yc) < h, na.rm = TRUE) ||
        any(abs(own$mean[own_span] - yc) < 0.7 * h, na.rm = TRUE)
      hit_label <- any(vapply(boxes, function(b) {
        b$x0 < x0 + w && x0 < b$x1 && abs(b$y - yc) < 2 * h
      }, logical(1)))
      if (!hit_curve && !hit_label) {
        out$x[i] <- x0
        out$y[i] <- yc
        boxes[[length(boxes) + 1L]] <- list(x0 = x0, x1 = x0 + w, y = yc)
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      out$x[i] <- labs$peak_t[i] + gap
      out$y[i] <- labs$peak_v[i] + 0.9 * h
      boxes[[length(boxes) + 1L]] <- list(x0 = out$x[i], x1 = out$x[i] + w, y = out$y[i])
    }
  }
  out
}

#' Coloured text with a white halo
#'
#' A white copy, fractionally larger, is drawn under the text, so the halo
#' covers only the glyphs and never a rectangle of data. Labels are placed
#' clear of other curves by `.place_curve_labels()`; the halo only guards
#' against gridlines and the label's own curve.
#' @keywords internal
#' @noRd
.halo_text <- function(data, mapping, size = 3.1, hjust = 0) {
  list(
    ggplot2::geom_text(data = data, mapping = mapping, colour = "white",
                       size = size + 0.12, hjust = hjust,
                       vjust = 0.5, show.legend = FALSE),
    ggplot2::geom_text(data = data, mapping = mapping, size = size, hjust = hjust,
                       vjust = 0.5, show.legend = FALSE)
  )
}

#' Fitted values of the full design for one voxel, run by run or jointly
#' @keywords internal
#' @noRd
.voxel_refit <- function(y, Xe, Xb, run_id, runwise) {
  n <- length(y)
  fitted <- base_only <- rep(NA_real_, n)
  pieces <- if (runwise) split(seq_len(n), run_id) else list(seq_len(n))
  for (rows in pieces) {
    keep_b <- colSums(Xb[rows, , drop = FALSE] != 0) > 0
    Xbr <- Xb[rows, keep_b, drop = FALSE]
    Xfull <- cbind(Xe[rows, , drop = FALSE], Xbr)
    fitted[rows] <- stats::lm.fit(Xfull, y[rows])$fitted.values
    base_only[rows] <- if (ncol(Xbr)) stats::lm.fit(Xbr, y[rows])$fitted.values else mean(y[rows])
  }
  rss_full <- sum((y - fitted)^2)
  rss_base <- sum((y - base_only)^2)
  r2 <- if (rss_base <= 1e-10 * max(sum(y^2), .Machine$double.xmin)) {
    NA_real_
  } else {
    1 - rss_full / rss_base
  }
  list(fitted = fitted, resid = y - fitted, partial_r2 = r2)
}

#' @keywords internal
#' @noRd
.voxel_series <- function(dataset, voxel) {
  one <- tryCatch(dataset[, voxel], error = function(e) NULL)
  if (!is.null(one)) {
    y <- tryCatch(as.matrix(.dset_data_matrix(one)), error = function(e) NULL)
    if (!is.null(y) && ncol(y) == 1L) return(as.numeric(y[, 1]))
  }
  as.numeric(as.matrix(.dset_data_matrix(dataset))[, voxel])
}

#' Observed series, refit, and time axis for one voxel
#' @keywords internal
#' @noRd
.voxel_fit_parts <- function(fit, voxel) {
  if (is.null(fit$dataset)) {
    stop("This fit does not carry its dataset; the voxel's data cannot be drawn.",
         call. = FALSE)
  }
  y <- .voxel_series(fit$dataset, voxel)
  model <- fit$model
  Xe <- as.matrix(design_matrix(model$event_model))
  Xb <- as.matrix(design_matrix(model$baseline_model))
  sf <- .model_sframe(model)
  time <- fmrihrf::samples(sf, global = TRUE)
  runs <- as.integer(fmrihrf::blocklens(sf))
  run_id <- rep(seq_along(runs), runs)
  runwise <- identical(fit$strategy, "runwise") && length(runs) > 1L
  rf <- .voxel_refit(y, Xe, Xb, run_id, runwise)
  list(y = y, rf = rf, time = time, runs = runs, run_id = run_id, runwise = runwise)
}

#' @keywords internal
#' @noRd
.plot_fit_timecourse <- function(fit, voxel) {
  auto <- is.null(voxel)
  voxel <- if (auto) .peak_voxel(fit) else .check_voxels(fit, voxel, single = TRUE)
  vp <- .voxel_fit_parts(fit, voxel)
  y <- vp$y
  rf <- vp$rf
  time <- vp$time
  runs <- vp$runs
  run_id <- vp$run_id
  model <- fit$model
  ink <- .fmrireg_ink

  # Onset lanes live in their own panel, placed far from the data values so a
  # breaks function can give that panel no axis ticks or gridlines.
  lane_off <- 1e9
  ev <- .event_onsets(model)
  panel_lv <- c("Signal", if (!is.null(ev) && nrow(ev)) "Onsets", "Residual")
  yr <- range(c(y, rf$fitted), na.rm = TRUE)
  lane_h <- 0.1 * diff(yr)

  df <- data.frame(
    time = rep(time, 3),
    value = c(y, rf$fitted, rf$resid),
    series = rep(c("observed", "fitted", "residual"), each = length(y)),
    run = rep(run_id, 3),
    panel = factor(rep(c("Signal", "Signal", "Residual"), each = length(y)),
                   levels = panel_lv)
  )
  dt <- if (length(time) > 1L) stats::median(diff(time)) else 1
  x_lab <- max(time) + dt

  breaks_fun <- function(lims) {
    if (min(lims) > lane_off / 2) return(numeric(0))
    b <- pretty(lims, n = 4)
    b[b >= lims[1] & b <= lims[2]]
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$value,
                                        group = .data$run)) +
    ggplot2::geom_hline(data = data.frame(panel = factor("Residual", levels = panel_lv), y = 0),
                        ggplot2::aes(yintercept = .data$y),
                        colour = ink[["rule"]], linewidth = 0.5)
  if (length(runs) > 1L) {
    starts <- c(1L, cumsum(runs)[-length(runs)] + 1L)
    bounds <- time[starts[-1]] - dt / 2
    mids <- (time[starts] + time[cumsum(runs)]) / 2
    p <- p + ggplot2::geom_vline(xintercept = bounds, colour = ink[["muted"]],
                                 linewidth = 0.3) +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  sec.axis = ggplot2::dup_axis(breaks = mids,
                                                               labels = paste("Run", seq_along(runs)),
                                                               name = NULL))
  } else {
    p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
  }
  obs <- df[df$series == "observed", ]
  fitd <- df[df$series == "fitted", ]
  p <- p +
    ggplot2::geom_line(data = obs, colour = ink[["muted"]], linewidth = 0.35) +
    ggplot2::geom_line(data = df[df$series == "residual", ],
                       colour = ink[["muted"]], linewidth = 0.35) +
    ggplot2::geom_line(data = fitd, colour = ink[["primary"]], linewidth = 0.6)

  # direct labels for the two signal series, at their final values
  n <- length(y)
  end_y <- .dodge_labels(c(y[n], rf$fitted[n]), 0.07 * diff(yr))
  end_lab <- data.frame(x = x_lab, y = end_y, txt = c("observed", "fitted"),
                        ink = c(ink[["secondary"]], ink[["primary"]]),
                        panel = factor("Signal", levels = panel_lv))
  p <- p + ggplot2::geom_text(data = end_lab,
                              ggplot2::aes(x = .data$x, y = .data$y, label = .data$txt,
                                           colour = .data$ink),
                              inherit.aes = FALSE, hjust = 0, size = 2.8) +
    ggplot2::scale_colour_identity()

  if (!is.null(ev) && nrow(ev)) {
    info <- .column_info(model)
    cols <- .condition_colours(info)
    lv <- levels(ev$condition)
    miss <- setdiff(lv, names(cols))
    cols <- c(cols, stats::setNames(rep(ink[["secondary"]], length(miss)), miss))
    lane_y <- lane_off + lane_h * rev(seq_along(lv))
    ev$y <- lane_y[as.integer(ev$condition)]
    ev$col <- unname(cols[as.character(ev$condition)])
    ev$panel <- factor("Onsets", levels = panel_lv)
    lane_lab <- data.frame(x = x_lab, y = lane_y, txt = lv,
                           col = unname(cols[lv]),
                           panel = factor("Onsets", levels = panel_lv))
    p <- p +
      ggplot2::geom_segment(data = ev,
                            ggplot2::aes(x = .data$onset, xend = .data$onset,
                                         y = .data$y - 0.38 * lane_h, yend = .data$y + 0.38 * lane_h,
                                         colour = .data$col),
                            inherit.aes = FALSE, linewidth = 0.6) +
      ggplot2::geom_text(data = lane_lab,
                         ggplot2::aes(x = .data$x, y = .data$y, label = .data$txt,
                                      colour = .data$col),
                         inherit.aes = FALSE, hjust = 0, size = 2.8) +
      # pad the lane panel so it keeps its height under space = "free_y"
      ggplot2::geom_blank(data = data.frame(time = min(time),
                                            value = lane_off + lane_h * c(0.3, length(lv) + 0.7),
                                            run = 1L,
                                            panel = factor("Onsets", levels = panel_lv)))
  }

  r2_txt <- if (is.na(rf$partial_r2)) "n/a" else sprintf("%.2f", rf$partial_r2)
  p + ggplot2::facet_grid(rows = ggplot2::vars(.data$panel), scales = "free_y",
                          space = "free_y", switch = "y") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::scale_y_continuous(breaks = breaks_fun, labels = .axis_num,
                                expand = ggplot2::expansion(mult = 0.04)) +
    ggplot2::labs(
      title = sprintf("Observed and fitted signal, voxel %d%s", voxel,
                      if (auto) " (largest |t|)" else ""),
      subtitle = sprintf("Partial R%s of events given baseline = %s", .glyph("sq"), r2_txt),
      caption = .wrap_notes(c(
        paste0(sprintf("Fitted: %s OLS of the full design", if (vp$runwise) "per-run" else "joint"),
               if (vp$runwise) " (run-specific coefficients; the estimates view pools runs)" else "",
               ". Signal and residual share one scale, in data units."))),
      x = "Time (s)", y = NULL) +
    .plot_theme() +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 70, 8, 8),
                   axis.text.x.top = ggplot2::element_text(colour = ink[["primary"]]),
                   panel.spacing.y = ggplot2::unit(0.5, "lines"))
}

#' Data for a set of voxels (time x voxels), reading only those columns
#' @keywords internal
#' @noRd
.voxel_matrix <- function(dataset, voxels) {
  sub <- tryCatch(dataset[, voxels], error = function(e) NULL)
  if (!is.null(sub)) {
    y <- tryCatch(as.matrix(.dset_data_matrix(sub)), error = function(e) NULL)
    if (!is.null(y) && ncol(y) == length(voxels)) return(y)
  }
  as.matrix(.dset_data_matrix(dataset))[, voxels, drop = FALSE]
}

#' Residual autocorrelation of several voxels, with its white-noise expectation
#'
#' Fits the full design by least squares (run by run for runwise fits, jointly
#' otherwise) and returns the autocorrelation of the residuals at lags
#' 1..lag_max, pooled across runs by run length. With `phi` (AR coefficients
#' per run) the data and design are first prewhitened run by run, so the
#' residuals are GLS residuals.
#'
#' Least-squares residuals e = M y (M = I - H) are autocorrelated even when the
#' noise is white; for white noise E[r_k] ~= sum_t M[t, t + k] / sum_t M[t, t]
#' within a run. `expect` returns that reference, pooled the same way.
#' @keywords internal
#' @noRd
.residual_acf <- function(Ymat, X, run_id, runwise, lag_max, phi = NULL) {
  runs <- split(seq_len(nrow(Ymat)), run_id)
  Yw <- vector("list", length(runs))
  Xw <- vector("list", length(runs))
  for (r in seq_along(runs)) {
    rows <- runs[[r]]
    y <- Ymat[rows, , drop = FALSE]
    x <- X[rows, , drop = FALSE]
    if (!is.null(phi)) {
      y <- .ar_whiten(y, phi[[r]])
      x <- .ar_whiten(x, phi[[r]])
    }
    Yw[[r]] <- y
    Xw[[r]] <- x
  }
  groups <- if (runwise) as.list(seq_along(runs)) else list(seq_along(runs))
  acf_sum <- 0
  exp_sum <- 0
  cov_sum <- 0
  n_tot <- 0
  for (g in groups) {
    y <- do.call(rbind, Yw[g])
    x <- do.call(rbind, Xw[g])
    x <- x[, colSums(x != 0) > 0, drop = FALSE]
    seg <- rep(seq_along(g), vapply(Yw[g], nrow, integer(1)))
    q <- qr(x)
    Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
    E <- y - Q %*% crossprod(Q, y)
    for (s in unique(seg)) {
      rows <- which(seg == s)
      n <- length(rows)
      # the run's block of M = I - QQ', from Q alone (never the full n x n M)
      Qs <- Q[rows, , drop = FALSE]
      trM <- n - sum(Qs^2)
      e_k <- vapply(seq_len(lag_max), function(k) {
        if (k >= n) return(0)
        -sum(Qs[seq_len(n - k), , drop = FALSE] * Qs[(k + 1):n, , drop = FALSE]) / trM
      }, numeric(1))
      acf_sum <- acf_sum + n * .acf_cols(E[rows, , drop = FALSE], lag_max)
      exp_sum <- exp_sum + n * e_k
      cov_sum <- cov_sum + n^2 * .acf_null_cov(Qs, lag_max)
      n_tot <- n_tot + n
    }
  }
  list(acf = acf_sum / n_tot, expect = exp_sum / n_tot, cov = cov_sum / n_tot^2,
       n = n_tot)
}

#' Null covariance of residual autocorrelations at lags 1..L
#'
#' For residuals e = R eps of white noise, with Ms = R R' = I - Qs Qs' the
#' residual-maker restricted to one run (Qs: that run's rows of Q), r_k ~ e' A_k e / e' e with A_k the symmetrised lag-k
#' shift. Treating the denominator as its expectation tr(Ms),
#' Cov(r_j, r_k) = 2 tr(A_j Ms A_k Ms) / tr(Ms)^2. Regression correlates the
#' lags, which an i.i.d. 1/n approximation ignores.
#' @keywords internal
#' @noRd
.acf_null_cov <- function(Qs, lag_max) {
  # Expanding M = I - Q Q' gives, for lags j, k >= 1,
  #   tr(A_j M A_k M) = [j == k] (n - j) / 2 - 2 <A_j Q, A_k Q> + tr(G_j G_k),
  # with G_k = Q' A_k Q, so only n x p and p x p objects are needed.
  n <- nrow(Qs)
  trM <- n - sum(Qs^2)
  AQ <- lapply(seq_len(lag_max), function(k) {
    if (k >= n) return(matrix(0, n, ncol(Qs)))
    up <- rbind(Qs[(k + 1):n, , drop = FALSE], matrix(0, k, ncol(Qs)))
    down <- rbind(matrix(0, k, ncol(Qs)), Qs[seq_len(n - k), , drop = FALSE])
    (up + down) / 2
  })
  G <- lapply(AQ, function(a) crossprod(Qs, a))
  S <- matrix(0, lag_max, lag_max)
  for (j in seq_len(lag_max)) {
    for (k in j:lag_max) {
      tr <- (if (j == k) max(n - j, 0) / 2 else 0) - 2 * sum(AQ[[j]] * AQ[[k]]) +
        sum(G[[j]] * t(G[[k]]))
      S[j, k] <- S[k, j] <- 2 * tr / trM^2
    }
  }
  S
}

#' Autocorrelation at lags 1..lag_max of each column of E
#' @keywords internal
#' @noRd
.acf_cols <- function(E, lag_max) {
  E <- sweep(E, 2, colMeans(E))
  den <- colSums(E^2)
  den[den == 0] <- NA_real_
  n <- nrow(E)
  # lags x voxels, also for a single voxel
  do.call(rbind, lapply(seq_len(lag_max), function(k) {
    if (k >= n) return(rep(NA_real_, ncol(E)))
    colSums(E[seq_len(n - k), , drop = FALSE] * E[(k + 1):n, , drop = FALSE]) / den
  }))
}

#' Apply an AR(p) prewhitening filter to each column
#' @keywords internal
#' @noRd
.ar_whiten <- function(E, phi) {
  p <- length(phi)
  if (!p || nrow(E) <= p + 2L) return(E)
  out <- E
  for (k in seq_len(p)) {
    out[(p + 1):nrow(E), ] <- out[(p + 1):nrow(E), ] - phi[k] * E[(p + 1 - k):(nrow(E) - k), ]
  }
  out[-seq_len(p), , drop = FALSE]
}

#' AR coefficients the fit used, one vector per run (or NULL)
#' @keywords internal
#' @noRd
.fit_ar_phi <- function(fit, nrun) {
  a <- fit$ar_coef
  if (is.null(a)) return(NULL)
  a <- tryCatch(lapply(if (is.list(a)) a else list(a), function(v) as.numeric(v)),
                error = function(e) NULL)
  if (!length(a) || !all(lengths(a) >= 1L) || anyNA(unlist(a))) return(NULL)
  rep_len(a, nrun)
}

#' @keywords internal
#' @noRd
.plot_fit_residual_acf <- function(fit, voxel, max_lag = 10L) {
  nvox <- nrow(coef(fit))
  chosen <- if (is.null(voxel)) "default" else "user"
  voxel <- if (is.null(voxel)) {
    unique(round(seq(1, nvox, length.out = min(nvox, 200L))))
  } else {
    .check_voxels(fit, voxel)
  }
  if (is.null(fit$dataset)) {
    stop("This fit does not carry its dataset; residuals cannot be computed.", call. = FALSE)
  }
  model <- fit$model
  X <- cbind(as.matrix(design_matrix(model$event_model)),
             as.matrix(design_matrix(model$baseline_model)))
  runs <- .run_lengths(model)
  run_id <- rep(seq_along(runs), runs)
  runwise <- identical(fit$strategy, "runwise") && length(runs) > 1L
  lag_max <- min(max_lag, min(runs) - 3L)
  if (lag_max < 1L) {
    stop("Runs are too short to estimate residual autocorrelation.", call. = FALSE)
  }
  Ymat <- .voxel_matrix(fit$dataset, voxel)
  phi <- .fit_ar_phi(fit, length(runs))

  series <- list(OLS = .residual_acf(Ymat, X, run_id, runwise, lag_max))
  if (!is.null(phi)) {
    series$prewhitened <- .residual_acf(Ymat, X, run_id, runwise, lag_max, phi = phi)
  }
  single <- length(voxel) == 1L
  ns <- length(series)
  dodge <- if (ns > 1L) c(-0.17, 0.17) else 0
  box_w <- if (ns > 1L) 0.15 else 0.36

  lb_count <- function(sr) {
    # portmanteau test on the deviations from the design-induced expectation,
    # using their exact null covariance (lags are correlated after regression).
    # AR coefficients are pooled across voxels, not fit to each series, so no
    # degrees of freedom are subtracted.
    d <- sr$acf - sr$expect
    Si <- tryCatch(solve(sr$cov), error = function(e) NULL)
    q <- if (is.null(Si)) {
      colSums(d^2 / diag(sr$cov))
    } else {
      colSums(d * (Si %*% d))
    }
    pv <- stats::pchisq(q, df = lag_max, lower.tail = FALSE)
    c(k = sum(pv < 0.05, na.rm = TRUE), p = if (single) pv[1] else NA_real_)
  }
  pts <- bands <- list()
  verdict <- character(0)
  for (i in seq_len(ns)) {
    nm <- names(series)[i]
    sr <- series[[i]]
    A <- sr$acf
    x <- seq_len(lag_max) + dodge[i]
    pts[[i]] <- data.frame(
      series = nm, x = x,
      r = if (single) A[, 1] else apply(A, 1, stats::median, na.rm = TRUE),
      lo = if (single) NA_real_ else apply(A, 1, stats::quantile, 0.1, na.rm = TRUE),
      hi = if (single) NA_real_ else apply(A, 1, stats::quantile, 0.9, na.rm = TRUE))
    half <- stats::qnorm(0.975) * sqrt(pmax(diag(sr$cov), 0))
    bands[[i]] <- data.frame(series = nm, x = x, c = sr$expect, lo = sr$expect - half,
                             hi = sr$expect + half)
    lb <- lb_count(sr)
    verdict <- c(verdict, if (single) {
      sprintf("%s p = %s", nm, formatC(lb[["p"]], format = "fg", digits = 2))
    } else {
      sprintf("%s %d/%d", nm, as.integer(lb[["k"]]), length(voxel))
    })
  }
  pts <- do.call(rbind, pts)
  bands <- do.call(rbind, bands)
  lv <- names(series)
  pts$series <- factor(pts$series, levels = lv)
  bands$series <- factor(bands$series, levels = lv)
  ink <- .fmrireg_ink
  inks <- stats::setNames(c(ink[["muted"]], ink[["primary"]])[seq_len(ns)], lv)
  shapes <- stats::setNames(c(21, 16)[seq_len(ns)], lv)

  p <- ggplot2::ggplot(pts, ggplot2::aes(x = .data$x, y = .data$r)) +
    ggplot2::geom_rect(data = bands,
                       ggplot2::aes(xmin = .data$x - box_w, xmax = .data$x + box_w,
                                    ymin = .data$lo, ymax = .data$hi),
                       inherit.aes = FALSE, fill = ink[["faint"]]) +
    ggplot2::geom_segment(data = bands,
                          ggplot2::aes(x = .data$x - box_w, xend = .data$x + box_w,
                                       y = .data$c, yend = .data$c),
                          inherit.aes = FALSE, colour = ink[["muted"]], linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0, colour = ink[["rule"]], linewidth = 0.4)
  if (!single) {
    p <- p + ggplot2::geom_linerange(ggplot2::aes(ymin = .data$lo, ymax = .data$hi,
                                                  colour = .data$series), linewidth = 0.7)
  }
  p +
    ggplot2::geom_point(ggplot2::aes(shape = .data$series, colour = .data$series),
                        size = 2.4, fill = "white", stroke = 0.9) +
    ggplot2::scale_shape_manual(values = shapes, name = NULL) +
    ggplot2::scale_colour_manual(values = inks, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq_len(lag_max), expand = c(0, 0),
                                limits = c(0.5, lag_max + 0.5)) +
    ggplot2::scale_y_continuous(labels = .axis_num) +
    ggplot2::labs(
      title = sprintf("Residual autocorrelation, %s",
                      if (single) sprintf("voxel %d", voxel)
                      else sprintf("%d %svoxels", length(voxel), if (chosen == "user") "selected " else "")),
      subtitle = if (single) {
        paste0("Portmanteau test, lags 1-", lag_max, ": ", paste(verdict, collapse = ", "))
      } else {
        sprintf("Voxels failing the portmanteau test (p < .05): %s (%s expected by chance)",
                paste(verdict, collapse = ", "), format(0.05 * length(voxel), digits = 2))
      },
      caption = .wrap_notes(c(
        paste0("Grey: 95% white-noise range for each series under this design, centred on the dash; the test uses the exact null covariance of the lags.",
               if (!single) " Points: median across voxels; lines: 10th-90th percentile." else ""),
        if (!is.null(phi)) sprintf("Prewhitened with the model's AR coefficients (pooled across voxels; lag 1 by run: %s).",
                                  paste(formatC(vapply(phi, `[`, numeric(1), 1), format = "f", digits = 2),
                                        collapse = ", "))),
        width = 120),
      x = "Lag (scans)", y = "Autocorrelation") +
    .plot_theme() +
    ggplot2::theme(legend.position = if (ns > 1L) "top" else "none",
                   legend.justification = "left",
                   legend.margin = ggplot2::margin(0, 0, 0, 0))
}

#' Global onset times and condition labels for all event terms
#'
#' Condition labels match `.column_info()`, so onsets share the condition
#' colours used by the other plots.
#' @keywords internal
#' @noRd
.event_onsets <- function(model) {
  em <- model$event_model
  sf <- em$sampling_frame
  level_order <- character(0)
  out <- lapply(terms(em), function(t) {
    if (is.null(t$onsets) || !length(t$onsets)) return(NULL)
    lev <- t$condition_levels
    ids <- t$condition_ids
    if (!is.null(lev) && !is.null(ids) && length(ids) == length(t$onsets)) {
      lev <- .strip_level_prefixes(as.character(lev))
      cond <- lev[ids]
    } else {
      lev <- t$varname %||% "events"
      cond <- rep(lev, length(t$onsets))
    }
    level_order <<- c(level_order, lev)
    data.frame(onset = fmrihrf::global_onsets(sf, t$onsets, t$blockids),
               condition = cond)
  })
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  if (is.null(out)) return(NULL)
  # parametric terms share onsets with the factor term; keep one tick per onset
  out <- out[!duplicated(out$onset), ]
  out$condition <- factor(out$condition,
                          levels = intersect(unique(level_order), unique(out$condition)))
  out
}
