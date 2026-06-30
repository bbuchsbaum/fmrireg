###############################################################################
## run_backends.R
##
## Execution-backend registry for the fan-out workflow. This is the
## extensibility seam: [run_jobs()] dispatches the actual mapping of jobs to a
## named backend. Two builtins ship -- "sequential" and "future" (the latter
## reaching any cluster via the active `future` plan). Third parties can
## register additional backends (e.g. a direct batchtools/clustermq driver)
## WITHOUT this package depending on any scheduler.
##
## A backend is a function `function(jobs, run_one, ...)` that returns a list of
## per-job result records (whatever `run_one` produces). It must preserve order.
##
## See `docs/plans/model-templates-multisubject-fanout.md`.
###############################################################################

.run_backend_registry <- new.env(parent = emptyenv())

#' @keywords internal
#' @noRd
.ensure_builtin_backends <- function() {
  if (!exists("sequential", envir = .run_backend_registry, inherits = FALSE)) {
    assign("sequential", function(jobs, run_one, ...) lapply(jobs, run_one),
           envir = .run_backend_registry)
  }
  if (!exists("future", envir = .run_backend_registry, inherits = FALSE)) {
    assign("future", function(jobs, run_one, ...) {
      future.apply::future_lapply(jobs, run_one, future.seed = TRUE)
    }, envir = .run_backend_registry)
  }
  invisible(NULL)
}

#' Register an execution backend for run_jobs()
#'
#' Adds a named backend that [run_jobs()] can dispatch to. This is how to plug in
#' a custom scheduler driver without modifying \pkg{fmrireg}. The builtins
#' \code{"sequential"} and \code{"future"} are always available.
#'
#' @param name Backend name (a string).
#' @param fn A function \code{function(jobs, run_one, ...)} returning an
#'   order-preserving list of per-job results.
#' @param overwrite Allow replacing an existing backend of the same name.
#' @return The backend name, invisibly.
#' @seealso [run_jobs()], [run_backends()]
#' @export
#' @examples
#' register_run_backend("first_only", function(jobs, run_one, ...) {
#'   list(run_one(jobs[[1]]))
#' })
register_run_backend <- function(name, fn, overwrite = FALSE) {
  assert_that(is.character(name), length(name) == 1, nzchar(name),
              msg = "'name' must be a single non-empty string")
  assert_that(is.function(fn), msg = "'fn' must be a function(jobs, run_one, ...)")
  .ensure_builtin_backends()
  if (exists(name, envir = .run_backend_registry, inherits = FALSE) && !isTRUE(overwrite)) {
    stop(sprintf("run backend '%s' already registered (use overwrite = TRUE)", name),
         call. = FALSE)
  }
  assign(name, fn, envir = .run_backend_registry)
  invisible(name)
}

#' List registered execution backends
#'
#' @return A character vector of backend names.
#' @seealso [register_run_backend()], [run_jobs()]
#' @export
run_backends <- function() {
  .ensure_builtin_backends()
  sort(ls(.run_backend_registry))
}

#' @keywords internal
#' @noRd
.get_run_backend <- function(name) {
  .ensure_builtin_backends()
  if (!exists(name, envir = .run_backend_registry, inherits = FALSE)) {
    stop(sprintf("unknown run backend '%s'; available: %s",
                 name, paste(run_backends(), collapse = ", ")), call. = FALSE)
  }
  get(name, envir = .run_backend_registry)
}
