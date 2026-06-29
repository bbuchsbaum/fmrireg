###############################################################################
## from_bids.R
##
## Discover BIDS-formatted data via the `bidser` package and build a per-subject
## manifest of bindings that `instantiate()` turns into jobs. This is the
## convenience populator for the common "BIDS in -> fan out -> BIDS out" case;
## `as_manifest()` is the generic, BIDS-free alternative.
##
## bidser is an optional dependency (Suggests). nvols / run_length is read from
## the BOLD header via RNifti pending bidser#73 (an nvols accessor).
##
## See `docs/plans/model-templates-multisubject-fanout.md`.
###############################################################################

#' Build a manifest from a BIDS project
#'
#' Uses \pkg{bidser} to discover preprocessed BOLD scans, events, confounds, TR,
#' and run lengths for each subject of a BIDS dataset, returning an
#' \code{fmri_manifest} (a list of per-subject bindings) suitable for
#' [instantiate()].
#'
#' The design formula's variables must match the columns of the BIDS
#' \code{events.tsv} (e.g. \code{trial_type}); a \code{run} column is added for
#' the block structure. Confound columns are selected by \code{confounds} (a
#' character vector or \code{bidser::confound_set()} result) passed straight to
#' \code{bidser::read_confounds()}.
#'
#' @param proj A \code{bidser::bids_project} (opened with \code{fmriprep = TRUE}
#'   for derivative discovery).
#' @param task BIDS task label.
#' @param space BIDS space for preprocessed scans (e.g.
#'   \code{"MNI152NLin2009cAsym"}); \code{NULL} matches any.
#' @param confounds Confound selection passed to \code{bidser::read_confounds()}
#'   as \code{cvars} (character vector or a \code{confound_set}); \code{NULL}
#'   reads no confounds.
#' @param mask A brain-mask path (length-1 character) applied to every subject's
#'   file-backed dataset. Required for file-backed BIDS data, which needs a mask
#'   (auto-resolution of the fmriprep brain mask is a planned enhancement).
#' @param desc BIDS \code{desc-} label of the preprocessed scans
#'   (default \code{"preproc"}).
#' @param subjects Optional character vector of subject labels (bare ids, e.g.
#'   \code{"01"}); default is all participants.
#' @param ... Forwarded to \code{bidser::read_confounds()} (e.g. \code{npcs},
#'   \code{clean}, \code{na_action}).
#' @return An \code{fmri_manifest}: a list of bindings, one per subject.
#' @seealso [as_manifest()], [instantiate()], [preflight()]
#' @export
#' @examples
#' \dontrun{
#' proj <- bidser::bids_project("study/", fmriprep = TRUE)
#' mani <- from_bids(proj, task = "stroop", space = "MNI152NLin2009cAsym",
#'                   confounds = bidser::confound_set("motion6"),
#'                   mask = "study/derivatives/.../space-MNI..._desc-brain_mask.nii.gz")
#' jobs <- instantiate(template, mani)
#' }
from_bids <- function(proj, task, space = NULL, confounds = NULL,
                      mask = NULL, desc = "preproc", subjects = NULL, ...) {
  if (!requireNamespace("bidser", quietly = TRUE)) {
    stop("from_bids() requires the 'bidser' package; install it or use as_manifest().",
         call. = FALSE)
  }
  if (!requireNamespace("RNifti", quietly = TRUE)) {
    stop("from_bids() requires 'RNifti' to read scan headers (run lengths).",
         call. = FALSE)
  }
  assert_that(inherits(proj, "bids_project"),
              msg = "'proj' must be a bidser::bids_project")
  assert_that(is.character(task), length(task) == 1, msg = "'task' must be a single string")
  if (!is.null(mask)) {
    assert_that(is.character(mask), length(mask) == 1,
                msg = "'mask' must be a single file path")
  }

  subs <- subjects %||% bidser::participants(proj)
  assert_that(length(subs) > 0, msg = "no participants found in project")

  bindings <- lapply(subs, function(s) {
    .bids_subject_binding(proj, subid = s, task = task, space = space,
                          confounds = confounds, desc = desc, mask = mask, ...)
  })
  bindings <- Filter(Negate(is.null), bindings)
  if (length(bindings) == 0) {
    stop(sprintf("no usable subjects for task '%s' (no preprocessed scans found)", task),
         call. = FALSE)
  }
  as_manifest(bindings)
}

#' @keywords internal
#' @noRd
.bids_run_num <- function(paths) {
  vapply(basename(paths), function(b) {
    m <- regmatches(b, regexpr("run-[0-9]+", b))
    if (length(m) == 0L) NA_integer_ else as.integer(sub("run-", "", m))
  }, integer(1), USE.NAMES = FALSE)
}

#' @keywords internal
#' @noRd
.bids_subject_binding <- function(proj, subid, task, space, confounds, desc,
                                   mask = NULL, ...) {
  scans <- bidser::preproc_scans(proj, subid = subid, task = task,
                                 space = space %||% ".*", desc = desc,
                                 full_path = TRUE)
  scans <- unname(scans)
  if (length(scans) == 0L) {
    warning(sprintf("subject '%s': no preprocessed scans for task '%s'; skipping",
                    subid, task), call. = FALSE)
    return(NULL)
  }

  # Order scans by run so scans / run_length / events / confounds all align.
  runs <- .bids_run_num(scans)
  ord <- if (all(!is.na(runs))) order(runs) else seq_along(scans)
  scans <- scans[ord]
  runs <- runs[ord]

  # run_length (nvols) from the BOLD header -- no voxel data loaded.
  run_length <- vapply(scans, function(p) as.integer(RNifti::niftiHeader(p)$dim[5]),
                       integer(1))

  tr <- bidser::get_repetition_time(proj, subid = subid, task = task)
  tr <- as.numeric(tr)[1]

  events <- .bids_events(proj, subid, task)
  conf <- if (!is.null(confounds)) {
    .bids_confounds(proj, subid, task, confounds, ...)
  } else NULL

  list(
    id = paste0("sub-", subid),
    scans = scans,
    TR = tr,
    run_length = unname(run_length),
    events = events,
    confounds = conf,
    mask = mask,
    meta = list(subject = subid, task = task, space = space)
  )
}

#' @keywords internal
#' @noRd
.bids_events <- function(proj, subid, task) {
  ev <- bidser::read_events(proj, subid = subid, task = task)
  ev <- as.data.frame(ev)               # drop grouping
  if (nrow(ev) == 0L) return(data.frame())
  ord <- order(suppressWarnings(as.integer(ev$.run)))
  ev <- ev[ord, , drop = FALSE]
  parts <- lapply(seq_len(nrow(ev)), function(i) {
    d <- as.data.frame(ev$data[[i]])
    d$.file <- NULL                     # bidser internal column
    d$run <- suppressWarnings(as.integer(ev$.run[i]))
    d
  })
  do.call(rbind, parts)
}

#' @keywords internal
#' @noRd
.bids_confounds <- function(proj, subid, task, confounds, ...) {
  cf <- bidser::read_confounds(proj, subid = subid, task = task,
                               cvars = confounds, ...)
  cf <- as.data.frame(cf)               # drop grouping
  if (nrow(cf) == 0L) return(NULL)
  ord <- order(suppressWarnings(as.integer(cf$run)))
  cf <- cf[ord, , drop = FALSE]
  lapply(seq_len(nrow(cf)), function(i) as.matrix(cf$data[[i]]))
}

#' Coerce bindings to an fmri_manifest
#'
#' Wraps a list of bindings (or a manifest \code{data.frame}) into an
#' \code{fmri_manifest} for [instantiate()]. This is the generic, BIDS-free way
#' to describe subjects; see [from_bids()] for the BIDS populator.
#'
#' @param x A list of bindings (named lists, each with at least \code{id},
#'   \code{scans}, \code{TR}, \code{run_length}) or a manifest \code{data.frame}.
#' @return An object of class \code{fmri_manifest}.
#' @seealso [from_bids()], [instantiate()]
#' @export
#' @examples
#' m <- as_manifest(list(
#'   list(id = "sub-01", scans = matrix(rnorm(80 * 2), 80, 2), TR = 2,
#'        run_length = c(40, 40),
#'        events = data.frame(onset = c(5, 45),
#'                            condition = factor(c("A", "B")), run = c(1, 2)))))
#' length(m)
as_manifest <- function(x) {
  bindings <- .as_binding_list(x)
  bad <- which(vapply(bindings, function(b) is.null(b[["id"]]), logical(1)))
  if (length(bad) > 0) {
    stop(sprintf("binding(s) %s have no 'id'", paste(bad, collapse = ", ")),
         call. = FALSE)
  }
  structure(bindings, class = c("fmri_manifest", "list"))
}

#' @export
print.fmri_manifest <- function(x, ...) {
  ids <- vapply(x, function(b) as.character(b[["id"]]), character(1))
  cat(sprintf("<fmri_manifest> %d unit(s)\n", length(x)))
  tasks <- unique(unlist(lapply(x, function(b) b$meta$task)))
  if (length(tasks)) cat("  task(s): ", paste(tasks, collapse = ", "), "\n", sep = "")
  cat("  ids: ", paste(utils::head(ids, 8), collapse = ", "),
      if (length(ids) > 8) sprintf(" ... (+%d)", length(ids) - 8) else "", "\n", sep = "")
  invisible(x)
}
