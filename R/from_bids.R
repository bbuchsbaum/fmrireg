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
#' Files are paired by their BIDS entities, parsed from the file names: a
#' preprocessed series, its events file and its confounds file belong to the
#' same run only when \code{sub}, \code{ses}, \code{task} and \code{run} are
#' equal, and labels are compared exactly (\code{subjects = "1"} never selects
#' \code{sub-10}; \code{task = "faces"} never selects \code{task-facesloc}).
#' Runs are ordered by session, then run number, and numbered \code{1, 2, ...}
#' in that order; each session's runs are blocks of their own. Two series for
#' one run (another resolution, acquisition or echo) are an error unless
#' \code{res} picks one of them.
#'
#' The design formula's variables must match the columns of the BIDS
#' \code{events.tsv} (e.g. \code{trial_type}); a \code{run} column holding the
#' block number is added for the block structure. Confound columns are selected
#' by \code{confounds} (a character vector or \code{bidser::confound_set()}
#' result) passed straight to \code{bidser::read_confounds()}.
#'
#' The TR is the \code{RepetitionTime} of the series' own JSON sidecar, else of
#' the raw series' sidecar (or an inherited \code{task-<label>_bold.json}),
#' else the NIfTI header's. Sidecars of other files under \code{derivatives/}
#' (for example maps written by [write_results()]) are never read.
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
#'   (auto-resolution of the fmriprep brain mask is a planned enhancement). The
#'   functional brain masks paired with each run are listed in each binding's
#'   \code{meta$runs$mask}.
#' @param desc BIDS \code{desc-} label of the preprocessed scans
#'   (default \code{"preproc"}).
#' @param subjects Optional character vector of subject labels (bare ids, e.g.
#'   \code{"01"}); default is all participants.
#' @param sessions Optional character vector of session labels (bare, e.g.
#'   \code{"1"}); default is every session.
#' @param res Optional BIDS \code{res-} label, for derivatives written at
#'   several resolutions; \code{NULL} accepts a single resolution per run.
#' @param ... Forwarded to \code{bidser::read_confounds()} (e.g. \code{npcs},
#'   \code{clean}, \code{na_action}).
#' @return An \code{fmri_manifest}: a list of bindings, one per subject. Each
#'   binding's \code{meta$runs} is a data frame with one row per run
#'   (\code{session}, \code{run}, \code{scan}, \code{events},
#'   \code{mask}).
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
                      mask = NULL, desc = "preproc", subjects = NULL,
                      sessions = NULL, res = NULL, ...) {
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
  for (arg in c("space", "res", "desc")) {
    v <- get(arg)
    if (!is.null(v)) {
      assert_that(is.character(v), length(v) == 1,
                  msg = sprintf("'%s' must be a single string or NULL", arg))
    }
  }
  if (!is.null(mask)) {
    assert_that(is.character(mask), length(mask) == 1,
                msg = "'mask' must be a single file path")
  }

  subs <- as.character(subjects %||% bidser::participants(proj))
  subs <- sub("^sub-", "", subs)
  assert_that(length(subs) > 0, msg = "no participants found in project")
  if (!is.null(sessions)) sessions <- sub("^ses-", "", as.character(sessions))

  runs <- .bids_run_table(proj, task = task, space = space, desc = desc,
                          subjects = subs, sessions = sessions, res = res)

  bindings <- lapply(subs, function(s) {
    r <- runs[runs$sub == s, , drop = FALSE]
    if (nrow(r) == 0L) {
      warning(sprintf("subject '%s': no preprocessed scans for task '%s'; skipping",
                      s, task), call. = FALSE)
      return(NULL)
    }
    .bids_subject_binding(proj, r, task = task, space = space,
                          confounds = confounds, mask = mask, ...)
  })
  bindings <- Filter(Negate(is.null), bindings)
  if (length(bindings) == 0) {
    stop(sprintf("no usable subjects for task '%s' (no preprocessed scans found)", task),
         call. = FALSE)
  }
  as_manifest(bindings)
}

#' BIDS entities of file names
#'
#' One row per file, one column per requested key; `NA` where a name lacks the
#' entity. Values are exact strings (`sub-01` gives `"01"`).
#' @keywords internal
#' @noRd
.bids_entities <- function(paths, keys = c("sub", "ses", "task", "acq", "ce",
                                           "rec", "dir", "run", "echo", "space",
                                           "res", "den", "desc")) {
  name <- basename(paths)
  out <- lapply(keys, function(key) {
    m <- regmatches(name, regexec(paste0("(^|_)", key, "-([^_.]+)"), name))
    vapply(m, function(x) if (length(x)) x[3] else NA_character_, character(1))
  })
  names(out) <- keys
  as.data.frame(out, stringsAsFactors = FALSE)
}

#' Run key shared by a run's BOLD, events and confounds files
#'
#' Sessions and runs that a name leaves out are `"1"`, as in bidser.
#' @keywords internal
#' @noRd
.bids_run_key <- function(ent) {
  ses <- ifelse(is.na(ent$ses), "1", ent$ses)
  run <- ifelse(is.na(ent$run), "1", as.character(suppressWarnings(as.integer(ent$run))))
  paste(ent$sub, ses, ent$task, run, sep = "|")
}

#' Every selected run, with the files that belong to it
#'
#' Discovers the preprocessed series once through bidser and filters them by
#' exact entity equality; events files come from the raw tree (never from
#' `derivatives/`) and are paired by `.bids_run_key()`.
#' @keywords internal
#' @noRd
.bids_run_table <- function(proj, task, space, desc, subjects, sessions = NULL,
                            res = NULL) {
  scans <- unname(bidser::preproc_scans(proj, full_path = TRUE))
  scans <- scans[grepl("_bold\\.nii(\\.gz)?$", basename(scans))]
  ent <- .bids_entities(scans)
  keep <- !is.na(ent$sub) & ent$sub %in% subjects &
    !is.na(ent$task) & ent$task == task
  if (!is.null(desc)) keep <- keep & !is.na(ent$desc) & ent$desc == desc
  if (!is.null(space)) keep <- keep & !is.na(ent$space) & ent$space == space
  if (!is.null(sessions)) keep <- keep & !is.na(ent$ses) & ent$ses %in% sessions
  if (!is.null(res)) keep <- keep & !is.na(ent$res) & ent$res == sub("^res-", "", res)
  scans <- scans[keep]
  ent <- ent[keep, , drop = FALSE]

  runs <- data.frame(
    sub = ent$sub,
    session = ent$ses,
    run = suppressWarnings(as.integer(ent$run)),
    space = ent$space,
    res = ent$res,
    scan = scans,
    key = .bids_run_key(ent),
    stringsAsFactors = FALSE
  )
  dup <- duplicated(runs$key) | duplicated(runs$key, fromLast = TRUE)
  if (any(dup)) {
    several_res <- length(unique(runs$res[dup])) > 1L
    stop(paste0(
      "several preprocessed series for one run: ",
      paste(basename(runs$scan[dup]), collapse = ", "), ". ",
      if (several_res) "Choose a resolution with `res = `." else
        "Narrow the selection with `space = ` or `res = `."
    ), call. = FALSE)
  }
  runs <- runs[order(runs$sub, ifelse(is.na(runs$session), "", runs$session),
                     ifelse(is.na(runs$run), 1L, runs$run)), , drop = FALSE]

  # Raw events files, by run key. derivatives/ is skipped: it can hold copies
  # (or outputs) that are not the experiment's events.
  root <- proj$path
  ev_files <- list.files(root, "_events\\.tsv$", recursive = TRUE, full.names = FALSE)
  ev_files <- ev_files[!grepl("^derivatives(/|$)", ev_files)]
  ev_ent <- .bids_entities(ev_files, keys = c("sub", "ses", "task", "run"))
  ev_ok <- !is.na(ev_ent$sub) & !is.na(ev_ent$task)
  ev_files <- file.path(root, ev_files[ev_ok])
  ev_keys <- .bids_run_key(ev_ent[ev_ok, , drop = FALSE])
  if (anyDuplicated(ev_keys)) {
    clash <- ev_keys %in% ev_keys[duplicated(ev_keys)]
    stop(sprintf("several events files for one run: %s",
                 paste(basename(ev_files[clash]), collapse = ", ")), call. = FALSE)
  }
  runs$events <- ev_files[match(runs$key, ev_keys)]

  # The functional brain mask fMRIPrep writes beside each series.
  mask_path <- sub("_desc-[^_]+_bold\\.nii(\\.gz)?$", "_desc-brain_mask.nii.gz", runs$scan)
  runs$mask <- ifelse(file.exists(mask_path), mask_path, NA_character_)
  rownames(runs) <- NULL
  runs
}

#' @keywords internal
#' @noRd
.bids_subject_binding <- function(proj, runs, task, space, confounds,
                                  mask = NULL, ...) {
  subid <- runs$sub[1]
  scans <- runs$scan

  # run_length (nvols) from the BOLD header -- no voxel data loaded.
  run_length <- vapply(scans, function(p) as.integer(RNifti::niftiHeader(p)$dim[5]),
                       integer(1), USE.NAMES = FALSE)

  tr <- vapply(scans, .bids_repetition_time, numeric(1), root = proj$path,
               USE.NAMES = FALSE)
  if (any(abs(tr - tr[1]) > 1e-6)) {
    stop(sprintf("subject '%s': runs have different repetition times (%s)",
                 subid, paste(unique(tr), collapse = ", ")), call. = FALSE)
  }

  events <- .bids_events(runs)
  conf <- if (!is.null(confounds)) {
    .bids_confounds(proj, runs, task, confounds, ...)
  } else NULL

  list(
    id = paste0("sub-", subid),
    scans = scans,
    TR = tr[1],
    run_length = run_length,
    events = events,
    confounds = conf,
    mask = mask,
    meta = list(subject = subid, task = task, space = space,
                runs = runs[c("session", "run", "scan", "events", "mask")])
  )
}

#' Repetition time of one preprocessed series
#'
#' Its own sidecar, else the raw series' sidecar or an inherited
#' `task-<label>_bold.json` (subject, then dataset level), else the header.
#' @keywords internal
#' @noRd
.bids_repetition_time <- function(scan, root) {
  read_tr <- function(f) {
    if (!file.exists(f)) return(NULL)
    tr <- tryCatch(jsonlite::read_json(f)[["RepetitionTime"]], error = function(e) NULL)
    if (is.numeric(tr) && length(tr) == 1L && is.finite(tr) && tr > 0) tr else NULL
  }
  name <- basename(scan)
  ent <- .bids_entities(scan, keys = c("sub", "ses", "task"))
  raw_stem <- sub("_(space|res|den|desc)-.*$", "", sub("_bold\\.nii(\\.gz)?$", "", name))
  sub_dir <- file.path(root, paste0("sub-", ent$sub))
  ses_dir <- if (is.na(ent$ses)) sub_dir else file.path(sub_dir, paste0("ses-", ent$ses))
  task_json <- paste0("task-", ent$task, "_bold.json")
  candidates <- c(
    sub("\\.nii(\\.gz)?$", ".json", scan),
    file.path(ses_dir, "func", paste0(raw_stem, "_bold.json")),
    file.path(ses_dir, paste0("sub-", ent$sub,
                              if (!is.na(ent$ses)) paste0("_ses-", ent$ses),
                              "_", task_json)),
    file.path(sub_dir, paste0("sub-", ent$sub, "_", task_json)),
    file.path(root, task_json)
  )
  for (f in candidates) {
    tr <- read_tr(f)
    if (!is.null(tr)) return(as.numeric(tr))
  }
  tr <- as.numeric(RNifti::niftiHeader(scan)$pixdim[5])
  if (!is.finite(tr) || tr <= 0) {
    stop(sprintf("no RepetitionTime for '%s': no sidecar gives one and the header has none",
                 name), call. = FALSE)
  }
  warning(sprintf("'%s': no sidecar RepetitionTime; using the NIfTI header's (%g)",
                  name, tr), call. = FALSE)
  tr
}

#' @keywords internal
#' @noRd
.bids_events <- function(runs) {
  missing <- is.na(runs$events)
  if (all(missing)) return(data.frame())
  if (any(missing)) {
    stop(sprintf("no events file for: %s",
                 paste(basename(runs$scan[missing]), collapse = ", ")), call. = FALSE)
  }
  parts <- lapply(seq_len(nrow(runs)), function(i) {
    d <- utils::read.delim(runs$events[i], na.strings = c("n/a", "NA"),
                           stringsAsFactors = FALSE, check.names = FALSE)
    d$run <- rep(i, nrow(d))
    d
  })
  .rbind_fill(parts)
}

#' Row-bind data frames whose columns differ, filling with NA
#' @keywords internal
#' @noRd
.rbind_fill <- function(parts) {
  cols <- unique(unlist(lapply(parts, names)))
  parts <- lapply(parts, function(d) {
    for (nm in setdiff(cols, names(d))) d[[nm]] <- rep(NA, nrow(d))
    d[cols]
  })
  do.call(rbind, parts)
}

#' @keywords internal
#' @noRd
.bids_confounds <- function(proj, runs, task, confounds, ...) {
  subid <- runs$sub[1]
  cf <- bidser::read_confounds(proj, subid = paste0("^", subid, "$"),
                               task = paste0("^", task, "$"),
                               cvars = confounds, ...)
  cf <- as.data.frame(cf)               # drop grouping
  if (nrow(cf) == 0L) return(NULL)
  cf_ent <- data.frame(sub = subid, ses = cf$session, task = task, run = cf$run,
                       stringsAsFactors = FALSE)
  idx <- match(runs$key, .bids_run_key(cf_ent))
  if (anyNA(idx)) {
    stop(sprintf("subject '%s': no confounds for: %s", subid,
                 paste(basename(runs$scan[is.na(idx)]), collapse = ", ")),
         call. = FALSE)
  }
  lapply(idx, function(i) as.matrix(cf$data[[i]]))
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
