###############################################################################
## preflight.R
##
## Validate jobs BEFORE fanning a model out over many subjects, so problems
## surface on the driver in seconds rather than on a compute node after the
## queue has drained. Checks the design contract (formula columns present),
## run-length / TR consistency, nuisance dimensions, and (optionally) that
## file-backed scans exist.
##
## See `docs/plans/model-templates-multisubject-fanout.md`.
###############################################################################

#' Preflight-check jobs before fan-out
#'
#' Validates one or more [fmri_job]s against the structural contract their
#' template implies, collecting issues per job. Intended to run on the driver
#' before [run_jobs()] / [export_jobs()].
#'
#' Checks performed per job: template validity; every variable referenced by the
#' design \code{formula} and \code{block} is a column of that job's event table;
#' \code{TR} is positive; run lengths are consistent with the data
#' (\code{matrix_dataset}: rows match \code{sum(run_length)};
#' \code{fmri_dataset}: one \code{run_length} per scan file); nuisance regressor
#' rows match the total number of scans; and, when \code{check_files = TRUE},
#' that file-backed scans exist on disk.
#'
#' The design-column check uses \code{all.vars()} on the formula, so it is
#' deliberately conservative: a formula with a variable-valued HRF argument
#' (e.g. \code{hrf(x, basis = my_basis)}) may flag \code{my_basis} as a missing
#' column. Event tables supplied as file paths are not yet validated here.
#'
#' @param x A single [fmri_job] or a list of them.
#' @param check_files Logical; for file-backed jobs, verify scan paths exist.
#' @param on_issue One of \code{"warn"} (default), \code{"error"} (stop if any
#'   issue), or \code{"collect"} (return silently).
#' @return Invisibly, an object of class \code{fmri_preflight} with an
#'   \code{$issues} data frame (\code{job_id}, \code{message}) and \code{$ok}.
#' @seealso [validate_template()], [run_jobs()], [export_jobs()]
#' @export
#' @examples
#' tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
#' job <- instantiate(tmpl, list(id = "sub-01",
#'                               scans = matrix(rnorm(80 * 2), 80, 2), TR = 2,
#'                               run_length = c(40, 40),
#'                               events = data.frame(onset = c(5, 45),
#'                                                   condition = factor(c("A", "B")),
#'                                                   run = c(1, 2))))
#' preflight(job)
preflight <- function(x, check_files = FALSE, on_issue = c("warn", "error", "collect")) {
  on_issue <- match.arg(on_issue)
  jobs <- if (inherits(x, "fmri_job")) list(x) else x
  assert_that(is.list(jobs) && length(jobs) > 0,
              msg = "'x' must be an fmri_job or a non-empty list of fmri_jobs")

  rows <- lapply(jobs, function(job) {
    msgs <- .preflight_job(job, check_files = check_files)
    if (length(msgs) == 0) return(NULL)
    data.frame(job_id = job$id, message = msgs, stringsAsFactors = FALSE)
  })
  issues <- do.call(rbind, rows)
  if (is.null(issues)) {
    issues <- data.frame(job_id = character(0), message = character(0),
                         stringsAsFactors = FALSE)
  }

  report <- structure(
    list(issues = issues, ok = nrow(issues) == 0,
         n_jobs = length(jobs)),
    class = "fmri_preflight"
  )

  if (!report$ok) {
    txt <- paste(sprintf("  [%s] %s", issues$job_id, issues$message), collapse = "\n")
    if (identical(on_issue, "error")) {
      stop(sprintf("preflight found %d issue(s):\n%s", nrow(issues), txt), call. = FALSE)
    } else if (identical(on_issue, "warn")) {
      warning(sprintf("preflight found %d issue(s):\n%s", nrow(issues), txt),
              call. = FALSE)
    }
  }
  invisible(report)
}

#' @keywords internal
#' @noRd
.preflight_job <- function(job, check_files = FALSE) {
  msgs <- character(0)
  add <- function(m) msgs <<- c(msgs, m)

  if (!inherits(job, "fmri_job")) return("not an 'fmri_job'")
  ok_tmpl <- tryCatch({ validate_template(job$template); TRUE },
                      error = function(e) { add(conditionMessage(e)); FALSE })
  if (!ok_tmpl) return(msgs)

  tmpl <- job$template
  spec <- job$dataset_spec
  args <- spec$args

  # Design contract: formula/block variables must be present as event columns.
  et <- args$event_table
  req <- unique(c(all.vars(tmpl$formula), all.vars(tmpl$block)))
  have <- if (is.data.frame(et)) names(et) else character(0)
  miss <- setdiff(req, have)
  if (length(miss) > 0) {
    add(sprintf("event table is missing design column(s): %s",
                paste(miss, collapse = ", ")))
  }

  # TR
  tr <- args$TR
  if (is.null(tr) || !is.numeric(tr) || any(tr <= 0)) {
    add("TR is missing or non-positive")
  }

  # Run-length consistency
  rl <- args$run_length
  if (is.null(rl) || !is.numeric(rl) || any(rl <= 0)) {
    add("run_length is missing or non-positive")
  } else if (identical(spec$constructor, "matrix_dataset") && !is.null(args$datamat)) {
    n <- nrow(args$datamat)
    if (sum(rl) != n) {
      add(sprintf("sum(run_length)=%d does not match data rows=%d", sum(rl), n))
    }
  } else if (identical(spec$constructor, "fmri_dataset") && !is.null(args$scans)) {
    if (length(args$scans) != length(rl)) {
      add(sprintf("number of scans (%d) does not match run_length entries (%d)",
                  length(args$scans), length(rl)))
    }
  }

  # Nuisance dimensions (matrix/data.frame: total rows; per-run list: one entry
  # per run, each with run_length rows).
  if (!is.null(job$nuisance) && is.numeric(rl)) {
    nz <- job$nuisance
    if (is.matrix(nz) || is.data.frame(nz)) {
      if (nrow(nz) != sum(rl)) {
        add(sprintf("nuisance rows (%d) do not match total scans (%d)",
                    nrow(nz), sum(rl)))
      }
    } else if (is.list(nz)) {
      if (length(nz) != length(rl)) {
        add(sprintf("nuisance list length (%d) does not match number of runs (%d)",
                    length(nz), length(rl)))
      } else {
        for (k in seq_along(nz)) {
          nrk <- nrow(as.matrix(nz[[k]]))
          if (!is.null(nrk) && nrk != rl[k]) {
            add(sprintf("nuisance run %d has %d rows, expected %d", k, nrk, rl[k]))
          }
        }
      }
    }
  }

  # File existence (opt-in)
  if (isTRUE(check_files) && identical(spec$constructor, "fmri_dataset") &&
      !is.null(args$scans)) {
    base <- args$base_path %||% "."
    paths <- ifelse(grepl("^(/|[A-Za-z]:)", args$scans),
                    args$scans, file.path(base, args$scans))
    gone <- args$scans[!file.exists(paths)]
    if (length(gone) > 0) {
      add(sprintf("scan file(s) not found: %s", paste(gone, collapse = ", ")))
    }
  }

  msgs
}

#' @export
print.fmri_preflight <- function(x, ...) {
  if (isTRUE(x$ok)) {
    cat(sprintf("<fmri_preflight> %d job(s): all clear\n", x$n_jobs))
  } else {
    cat(sprintf("<fmri_preflight> %d job(s): %d issue(s)\n",
                x$n_jobs, nrow(x$issues)))
    for (i in seq_len(nrow(x$issues))) {
      cat(sprintf("  ! %s: %s\n", x$issues$job_id[i], x$issues$message[i]))
    }
  }
  invisible(x)
}
