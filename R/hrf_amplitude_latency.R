# Amplitude and latency of multi-basis responses from a fitted model

#' Response Amplitude and Latency from a Multi-Basis Fit
#'
#' @description
#' For every condition modelled with a multi-basis HRF (canonical plus
#' temporal derivative, `basis = "spmg2"`, or the quadrature pair
#' `fmrihrf::HRF_SPMGQ`), combine the fitted basis coefficients into a
#' response amplitude and a latency per voxel, using
#' `fmrihrf::amplitude_latency()`.
#'
#' The canonical coefficient alone underestimates a response that peaks
#' earlier or later than the canonical HRF, so condition differences in it mix
#' amplitude with timing. The coefficients of a time-shifted response trace a
#' known curve in coefficient space; each voxel's coefficients are matched to
#' the nearest direction on that curve, which gives its latency, and their
#' size along that direction gives the amplitude.
#'
#' @param object A fitted `fmri_lm` object.
#' @param term Optional character vector of term tags (as in
#'   `fmridesign::design_colmap()`) to restrict the result to. By default
#'   every event term with at least two basis functions is used.
#' @param shifts Candidate latencies in seconds (positive = later than the
#'   canonical response). Estimates are limited to this range.
#' @param method `"curve"` (default) matches coefficients to the shift curve
#'   and returns amplitude and latency. `"boost"` returns the
#'   derivative-boosted amplitude of Calhoun et al. (2004) and no latency.
#' @param min_fit Voxels whose coefficients match the shift curve with a
#'   cosine below this value get `NA` (default 0).
#'
#' @return An object of class `hrf_amplitude_latency`: a list with numeric
#'   matrices `amplitude`, `latency` (seconds) and `fit` (cosine between the
#'   coefficients and the matched curve point), each with one row per voxel
#'   and one column per condition, plus `conditions`, a data frame giving
#'   each column's `term`, `condition` and `basis`. `latency` and `fit` are
#'   `NA` for `method = "boost"`.
#'
#' @details
#' Coefficients are matched in the metric of the basis functions' Gram
#' matrix, which is appropriate for well-separated events; with heavily
#' overlapping responses the estimates remain consistent but are not the
#' most efficient, and a response that is shifted by several seconds is only
#' approximately spanned by the basis, so part of it can be absorbed by
#' neighbouring regressors. Latency is relative to the model's own timing:
#' it includes any offset between the event onsets and the acquisition times
#' assumed by the sampling frame. Latency is poorly determined where the amplitude is near
#' zero, so mask the latency map (for example with a test on the amplitude)
#' before interpreting it. The summaries are non-linear in the coefficients
#' and have no closed-form standard error here: take the per-subject maps to
#' a group analysis.
#'
#' This function needs an fmrihrf that provides `amplitude_latency()`.
#'
#' @references
#' Calhoun, V. D., Stevens, M. C., Pearlson, G. D., & Kiehl, K. A. (2004).
#' fMRI analysis with the general linear model: removal of latency-induced
#' amplitude bias by incorporation of hemodynamic derivative terms.
#' *NeuroImage*, 22(1), 252-257. \doi{10.1016/j.neuroimage.2003.12.029}
#'
#' Henson, R. N. A., Price, C. J., Rugg, M. D., Turner, R., & Friston, K. J.
#' (2002). Detecting latency differences in event-related BOLD responses:
#' application to words versus nonwords and initial versus repeated face
#' presentations. *NeuroImage*, 15(1), 83-97. \doi{10.1006/nimg.2001.0940}
#'
#' Liao, C. H., Worsley, K. J., Poline, J.-B., Aston, J. A. D., Duncan, G. H.,
#' & Evans, A. C. (2002). Estimating the delay of the fMRI response.
#' *NeuroImage*, 16(3), 593-606. \doi{10.1006/nimg.2002.1096}
#'
#' @examples
#' \donttest{
#' if (exists("amplitude_latency", asNamespace("fmrihrf"))) {
#'   set.seed(1)
#'   events <- data.frame(
#'     onset = seq(4, 280, by = 16),
#'     condition = factor(rep(c("A", "B"), 9)),
#'     run = 1
#'   )
#'   # Simulate one voxel in which B peaks 1.5 s later than A
#'   time <- fmrihrf::samples(fmrihrf::sampling_frame(150, TR = 2))
#'   resp <- function(on) as.numeric(fmrihrf::evaluate(fmrihrf::regressor(on), time))
#'   y <- 100 + resp(events$onset[events$condition == "A"]) +
#'     resp(events$onset[events$condition == "B"] + 1.5) + rnorm(150, sd = 0.02)
#'   dset <- matrix_frame(cbind(y), TR = 2, run_length = 150, event_table = events)
#'   fit <- fmri_lm(onset ~ hrf(condition, basis = "spmg2"), block = ~ run,
#'                  dataset = dset)
#'   al <- hrf_amplitude_latency(fit)
#'   al$amplitude
#'   al$latency
#' }
#' }
#' @export
hrf_amplitude_latency <- function(object, term = NULL,
                                  shifts = seq(-4, 4, by = 0.05),
                                  method = c("curve", "boost"), min_fit = 0) {
  method <- match.arg(method)
  if (!inherits(object, "fmri_lm")) {
    stop("`object` must be a fitted fmri_lm object.", call. = FALSE)
  }
  if (!exists("amplitude_latency", envir = asNamespace("fmrihrf"), inherits = FALSE)) {
    stop("hrf_amplitude_latency() needs fmrihrf::amplitude_latency(); the installed fmrihrf (",
         as.character(utils::packageVersion("fmrihrf")), ") does not have it. Update fmrihrf.",
         call. = FALSE)
  }
  estimator <- get("amplitude_latency", envir = asNamespace("fmrihrf"))

  event_model <- object$model$event_model
  betas <- stats::coef(object)
  colmap <- as.data.frame(fmridesign::design_colmap(event_model))
  need <- c("name", "term_tag", "condition", "basis_ix", "basis_total")
  if (!all(need %in% names(colmap))) {
    stop("The event model's column map lacks basis information.", call. = FALSE)
  }
  colmap <- colmap[colmap$name %in% colnames(betas) &
                     !is.na(colmap$basis_total) & colmap$basis_total >= 2L, , drop = FALSE]
  if (!is.null(term)) {
    unknown <- setdiff(term, unique(colmap$term_tag))
    if (length(unknown) > 0L) {
      stop("No multi-basis event term named: ", paste(unknown, collapse = ", "),
           ". Available: ", paste(unique(colmap$term_tag), collapse = ", "), call. = FALSE)
    }
    colmap <- colmap[colmap$term_tag %in% term, , drop = FALSE]
  }
  if (nrow(colmap) == 0L) {
    stop("The model has no event term with a multi-basis HRF ",
         "(e.g. hrf(cond, basis = \"spmg2\")).", call. = FALSE)
  }

  event_terms <- fmridesign::event_terms(event_model)
  term_tags <- vapply(event_terms, function(tm) {
    tag <- attr(tm, "term_tag")
    if (is.null(tag)) NA_character_ else as.character(tag)[1L]
  }, character(1))

  curve_id <- paste(colmap$term_tag, colmap$condition, sep = "\r")
  curves <- unique(curve_id)
  amp <- lat <- fit <- matrix(NA_real_, nrow = nrow(betas), ncol = length(curves))
  info <- data.frame(term = character(length(curves)),
                     condition = character(length(curves)),
                     basis = character(length(curves)), stringsAsFactors = FALSE)

  for (j in seq_along(curves)) {
    rows <- colmap[curve_id == curves[[j]], , drop = FALSE]
    rows <- rows[order(rows$basis_ix), , drop = FALSE]
    tag <- rows$term_tag[1L]
    k <- rows$basis_total[1L]
    if (!identical(as.integer(rows$basis_ix), seq_len(k))) {
      stop(sprintf("Condition '%s' of term '%s' does not have basis columns 1 to %d.",
                   rows$condition[1L], tag, k), call. = FALSE)
    }
    tm <- event_terms[[match(tag, term_tags)]]
    hrf <- attr(tm, "hrfspec")$hrf
    if (!inherits(hrf, "HRF") || fmrihrf::nbasis(hrf) != k) {
      stop(sprintf("Could not recover the %d-function HRF basis of term '%s'.", k, tag),
           call. = FALSE)
    }
    est <- estimator(betas[, rows$name, drop = FALSE], hrf, shifts = shifts,
                     method = method, min_fit = min_fit)
    amp[, j] <- est$amplitude
    lat[, j] <- est$latency
    fit[, j] <- est$fit
    info$term[j] <- tag
    info$condition[j] <- as.character(rows$condition[1L])
    info$basis[j] <- if ("basis_name" %in% names(rows)) as.character(rows$basis_name[1L]) else NA_character_
  }
  labels <- paste(info$term, info$condition, sep = "_")
  colnames(amp) <- colnames(lat) <- colnames(fit) <- labels
  structure(list(amplitude = amp, latency = lat, fit = fit, conditions = info,
                 method = method),
            class = "hrf_amplitude_latency")
}

#' @export
print.hrf_amplitude_latency <- function(x, ...) {
  cat(sprintf("Response amplitude and latency (%s): %d voxel%s x %d condition%s\n",
              x$method, nrow(x$amplitude), if (nrow(x$amplitude) == 1L) "" else "s",
              ncol(x$amplitude), if (ncol(x$amplitude) == 1L) "" else "s"))
  summ <- data.frame(
    x$conditions,
    median_amplitude = apply(x$amplitude, 2L, stats::median, na.rm = TRUE),
    median_latency = apply(x$latency, 2L, stats::median, na.rm = TRUE),
    row.names = NULL
  )
  print(summ, digits = 3)
  invisible(x)
}
