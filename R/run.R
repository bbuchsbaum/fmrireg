###############################################################################
## run.R
##
## Execution layer for the multi-subject fan-out workflow. `run_job()` realizes
## one job's dataset, assembles the model, fits it, and applies the template's
## reducer. `run_jobs()` runs a batch sequentially or in parallel via the active
## `future` plan, isolating per-job errors so one bad subject does not abort the
## whole study.
##
## No scheduler-specific code lives here: parallelism is expressed purely
## through the `future` API, so a cluster is reached by the user choosing a
## `future` backend (e.g. `future.batchtools`).
##
## See `docs/plans/model-templates-multisubject-fanout.md`.
###############################################################################

#' Run a single job
#'
#' Realizes the dataset, assembles the model, fits it with the template's
#' control options (and engine, if any), and applies the template's reducer.
#'
#' @param job An [fmri_job].
#' @param progress Logical; show a fitting progress bar.
#' @return The reducer's output, or the fitted \code{fmri_lm} if the template has
#'   no reducer.
#' @seealso [run_jobs()], [build_model()]
#' @export
run_job <- function(job, progress = FALSE) {
  assert_that(inherits(job, "fmri_job"), msg = "'job' must be an 'fmri_job'")
  tmpl <- job$template
  dataset <- realize_dataset(job)
  model <- build_model(job, dataset)
  cfg <- tmpl$control

  fit <- if (!is.null(tmpl$engine)) {
    .fmri_lm_dispatch_engine(
      model = model, dataset = dataset, engine = tmpl$engine,
      lowrank = NULL, cfg = cfg, engine_args = tmpl$engine_args
    )
  } else {
    fmri_lm_fit(model, dataset, strategy = tmpl$strategy, cfg = cfg,
                progress = progress)
  }

  # Apply template-level contrasts so they actually drive output (reduce_contrasts
  # and any downstream consumer read them from this attribute). Formula-embedded
  # contrasts continue to flow through coef()/write_results() as before.
  if (!is.null(tmpl$contrasts)) {
    attr(fit, "template_contrasts") <- fit_contrasts(fit, tmpl$contrasts)
  }

  apply_reducer(tmpl$reducer, fit, job)
}

#' Run a batch of jobs
#'
#' Executes a list of [fmri_job]s, returning an \code{fmri_batch_result}. By
#' default each job is isolated: a failure is captured rather than aborting the
#' batch. Set \code{parallel = TRUE} to dispatch through the active
#' \code{future::plan()} (e.g. \code{multisession}, \code{cluster}, or a
#' \code{future.batchtools} cluster backend).
#'
#' @param jobs A single [fmri_job] or a list of them.
#' @param parallel Logical shorthand; \code{TRUE} selects the \code{"future"}
#'   backend, \code{FALSE} the \code{"sequential"} one. Ignored if \code{backend}
#'   is given.
#' @param progress Logical; per-job fitting progress.
#' @param on_error Either \code{"isolate"} (default; capture per-job errors) or
#'   \code{"stop"} (fail the batch on the first error).
#' @param backend Optional execution backend name (see [run_backends()]).
#'   Overrides \code{parallel}. The \code{"future"} backend honours the active
#'   \code{future::plan()} (including \code{future.batchtools} cluster plans).
#' @return An object of class \code{fmri_batch_result}.
#' @seealso [run_job()], [batch_values()], [batch_errors()], [register_run_backend()]
#' @export
#' @examples
#' \dontrun{
#' jobs <- instantiate(tmpl, manifest)
#' future::plan(future::multisession, workers = 4)
#' res <- run_jobs(jobs, parallel = TRUE)
#' values <- batch_values(res)
#' }
run_jobs <- function(jobs, parallel = FALSE, progress = FALSE,
                     on_error = c("isolate", "stop"), backend = NULL) {
  if (inherits(jobs, "fmri_job")) jobs <- list(jobs)
  assert_that(is.list(jobs) && length(jobs) > 0,
              msg = "'jobs' must be an fmri_job or a non-empty list of fmri_jobs")
  assert_that(all(vapply(jobs, inherits, logical(1), "fmri_job")),
              msg = "every element of 'jobs' must be an 'fmri_job'")
  on_error <- match.arg(on_error)

  run_one <- function(job) {
    if (identical(on_error, "stop")) {
      list(id = job$id, ok = TRUE, value = run_job(job, progress = progress),
           error = NULL)
    } else {
      tryCatch(
        list(id = job$id, ok = TRUE, value = run_job(job, progress = progress),
             error = NULL),
        error = function(e) list(id = job$id, ok = FALSE, value = NULL,
                                 error = conditionMessage(e))
      )
    }
  }

  backend_name <- backend %||% (if (isTRUE(parallel)) "future" else "sequential")
  results <- .get_run_backend(backend_name)(jobs, run_one)

  structure(
    list(
      results = results,
      ids = vapply(results, `[[`, character(1), "id"),
      ok = vapply(results, `[[`, logical(1), "ok")
    ),
    class = "fmri_batch_result"
  )
}

#' Extract successful values from a batch result
#'
#' @param x An \code{fmri_batch_result}.
#' @return A named list (by job id) of reducer outputs for jobs that succeeded.
#' @seealso [run_jobs()], [batch_errors()]
#' @export
batch_values <- function(x) {
  assert_that(inherits(x, "fmri_batch_result"))
  ok <- Filter(function(r) isTRUE(r$ok), x$results)
  stats::setNames(lapply(ok, `[[`, "value"), vapply(ok, `[[`, character(1), "id"))
}

#' Extract per-job errors from a batch result
#'
#' @param x An \code{fmri_batch_result}.
#' @return A named character vector (by job id) of error messages for jobs that
#'   failed; empty if all succeeded.
#' @seealso [run_jobs()], [batch_values()]
#' @export
batch_errors <- function(x) {
  assert_that(inherits(x, "fmri_batch_result"))
  failed <- Filter(function(r) !isTRUE(r$ok), x$results)
  stats::setNames(vapply(failed, function(r) r$error %||% NA_character_, character(1)),
                  vapply(failed, `[[`, character(1), "id"))
}

#' @export
print.fmri_batch_result <- function(x, ...) {
  n <- length(x$results)
  nok <- sum(x$ok)
  cat(sprintf("<fmri_batch_result> %d job(s): %d ok, %d failed\n", n, nok, n - nok))
  if (nok < n) {
    errs <- batch_errors(x)
    for (id in names(errs)) cat(sprintf("  ! %s: %s\n", id, errs[[id]]))
  }
  invisible(x)
}
