###############################################################################
## reducers.R
##
## The reducer protocol for the multi-subject fan-out workflow. A reducer is a
## function `function(fit, job)` that turns a fitted `fmri_lm` into the
## per-subject output: it may write results to disk (and return the paths) or
## return a small in-memory summary (e.g. a tidy contrast table) to the driver.
##
## Builtin reducers are tagged with class "fmri_reducer" so they are exempt from
## the closure-capture warning in `fmri_template()` (we trust that they capture
## only small, portable configuration).
##
## See `docs/plans/model-templates-multisubject-fanout.md`.
###############################################################################

#' @keywords internal
#' @noRd
new_reducer <- function(fn) {
  stopifnot(is.function(fn))
  class(fn) <- c("fmri_reducer", "function")
  fn
}

#' Apply a reducer to a fitted model
#'
#' Internal helper used by [run()]. A \code{NULL} reducer returns the fitted
#' object unchanged (equivalent to [reduce_identity()]).
#'
#' @param reducer A reducer function or \code{NULL}.
#' @param fit A fitted \code{fmri_lm} object.
#' @param job The [fmri_job] that produced \code{fit}.
#' @return The reducer's output.
#' @keywords internal
#' @noRd
apply_reducer <- function(reducer, fit, job) {
  if (is.null(reducer)) return(fit)
  reducer(fit, job)
}

#' Reducer: return the fitted model unchanged
#'
#' The escape hatch when you want the entire \code{fmri_lm} object per subject.
#' Note this is the largest possible output.
#'
#' @return A reducer function suitable for [fmri_template()].
#' @family reducers
#' @export
reduce_identity <- function() {
  new_reducer(function(fit, job) fit)
}

#' @keywords internal
#' @noRd
# Reshape a (terms x voxels) coefficient matrix to a long data frame with
# columns (term, voxel, <value>). Order is column-major, matching as.vector().
.coef_to_long <- function(m, value = "estimate") {
  m <- as.matrix(m)
  terms <- rownames(m)
  if (is.null(terms)) terms <- paste0("term", seq_len(nrow(m)))
  out <- data.frame(
    term  = rep(terms, times = ncol(m)),
    voxel = rep(seq_len(ncol(m)), each = nrow(m)),
    stringsAsFactors = FALSE
  )
  out[[value]] <- as.vector(m)
  out
}

#' Reducer: tidy contrast table
#'
#' Extracts the fitted contrasts as a long data frame (one row per
#' contrast x voxel), optionally prefixed with the job id so results stack
#' cleanly across subjects. Surfaces template-level contrasts (those passed to
#' [fmri_template()], applied by [run_job()]) when present, otherwise the
#' model's formula-embedded contrasts via [coef()]. Returns a zero-row frame
#' when no contrasts are defined.
#'
#' @param add_id Logical; prepend a \code{job_id} column. Default \code{TRUE}.
#' @return A reducer function suitable for [fmri_template()].
#' @family reducers
#' @export
reduce_contrasts <- function(add_id = TRUE) {
  force(add_id)
  new_reducer(function(fit, job) {
    tc <- attr(fit, "template_contrasts")
    if (!is.null(tc)) {
      out <- .tidy_template_contrasts(tc)
    } else {
      con <- tryCatch(coef(fit, type = "contrasts"), error = function(e) NULL)
      if (is.null(con) || prod(dim(as.matrix(con))) == 0) {
        out <- data.frame(term = character(0), voxel = integer(0),
                          estimate = numeric(0), stringsAsFactors = FALSE)
      } else {
        # coef(type="contrasts") is voxels x contrasts; transpose to contrasts x voxels.
        out <- .coef_to_long(t(as.matrix(con)), "estimate")
      }
    }
    if (isTRUE(add_id)) out <- cbind(job_id = job$id, out, stringsAsFactors = FALSE)
    out
  })
}

#' @keywords internal
#' @noRd
# Turn the list returned by fit_contrasts() (one entry per contrast, each with
# per-voxel estimate/se/stat) into a long (term, voxel, estimate, se, stat) frame.
.tidy_template_contrasts <- function(tc) {
  parts <- lapply(names(tc), function(nm) {
    e <- tc[[nm]]
    est <- as.numeric(e$estimate)
    df <- data.frame(term = e$name %||% nm, voxel = seq_along(est),
                     estimate = est, stringsAsFactors = FALSE)
    if (!is.null(e$se)) df$se <- as.numeric(e$se)
    if (!is.null(e$stat)) df$stat <- as.numeric(e$stat)
    df
  })
  do.call(rbind, parts)
}

#' Reducer: tidy beta-estimate table
#'
#' Extracts the fitted condition estimates as a long data frame (one row per
#' term x voxel), with standard error and a Wald statistic when available.
#' Built on [coef()] / [standard_error()] for robustness.
#'
#' @param add_id Logical; prepend a \code{job_id} column. Default \code{TRUE}.
#' @param include_baseline Logical; include baseline/nuisance terms. Default
#'   \code{FALSE}.
#' @return A reducer function suitable for [fmri_template()].
#' @family reducers
#' @export
reduce_betas <- function(add_id = TRUE, include_baseline = FALSE) {
  force(add_id); force(include_baseline)
  new_reducer(function(fit, job) {
    est <- as.matrix(coef(fit, type = "betas", include_baseline = include_baseline))
    # Normalize to terms x voxels: terms are the *named* dimension. coef() returns
    # terms x voxels (rownames = terms) without baseline, but voxels x terms
    # (colnames = terms) with baseline -- transpose the latter.
    if (is.null(rownames(est)) && !is.null(colnames(est))) est <- t(est)
    out <- .coef_to_long(est, "estimate")
    se <- tryCatch(standard_error(fit, type = "estimates"), error = function(e) NULL)
    # standard_error() is voxels x terms; attach only when every estimate term is
    # present (so the baseline-included case, whose SE is unavailable, is skipped
    # rather than mis-aligned).
    if (!is.null(se) && !is.null(rownames(est)) &&
        length(rownames(est)) > 0 && all(rownames(est) %in% colnames(se))) {
      se_mat <- t(as.matrix(se[, rownames(est), drop = FALSE]))
      out$se <- as.vector(se_mat)
      out$stat <- out$estimate / out$se
    }
    if (isTRUE(add_id)) out <- cbind(job_id = job$id, out, stringsAsFactors = FALSE)
    out
  })
}

#' Reducer: write per-subject results to disk (BIDS-keyed)
#'
#' Calls [write_results()] on the fitted model, keying the output filenames by
#' the job's \code{meta} (\code{subject}, \code{task}, \code{space}). Returns the
#' vector of created file paths. This is the default reducer for large studies:
#' workers write compact statistical maps that [group_data()] / [fmri_meta()]
#' can pick up for the group level.
#'
#' @param format Output format passed to [write_results()].
#' @param stats Which contrast statistics to write (mapped to
#'   \code{contrast_stats}).
#' @param contrasts Optional restriction to specific contrasts.
#' @param path Output directory (passed to [write_results()]).
#' @param desc BIDS \code{desc-} label.
#' @param overwrite Overwrite existing files.
#' @return A reducer function suitable for [fmri_template()].
#' @note [write_results()] requires BIDS entities: each job's \code{meta} should
#'   carry \code{subject} (defaults to the job id), \code{task}, and
#'   \code{space}. A missing \code{task}/\code{space} surfaces as an error when
#'   the reducer runs.
#' @family reducers
#' @export
reduce_write_results <- function(format = c("h5", "nifti", "gds"),
                                 stats = c("beta", "tstat"),
                                 contrasts = NULL,
                                 path = NULL,
                                 desc = "GLM",
                                 overwrite = FALSE) {
  format <- match.arg(format)
  force(stats); force(contrasts); force(path); force(desc); force(overwrite)
  new_reducer(function(fit, job) {
    meta <- job$meta
    write_results(
      fit,
      path = path,
      subject = meta$subject %||% job$id,
      task = meta$task,
      space = meta$space,
      desc = desc,
      format = format,
      contrasts = contrasts,
      contrast_stats = stats,
      overwrite = overwrite
    )
  })
}
