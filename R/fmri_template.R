###############################################################################
## fmri_template.R
##
## Subject-invariant model "template" objects for the multi-subject fan-out
## workflow. A template captures everything about a model that is constant
## across subjects (the design formula, baseline specification, contrasts,
## fitting control, and an output reducer); the per-subject data is bound in
## later via `instantiate()` to produce serializable `fmri_job` recipes.
##
## See `docs/plans/model-templates-multisubject-fanout.md`.
###############################################################################

#' Baseline specification (subject-invariant)
#'
#' Describes how the nuisance / baseline model should be built for every
#' subject, without binding any particular subject's confound values. The
#' concrete \code{baseline_model} is assembled per subject at instantiation
#' time from this spec, the subject's \code{sampling_frame}, and the subject's
#' confound matrix.
#'
#' @param degree Integer drift degree passed to \code{baseline_model()}.
#' @param basis Drift basis: one of \code{"bs"}, \code{"poly"}, \code{"ns"},
#'   \code{"constant"}.
#' @param confounds Optional confound selection used to populate the per-subject
#'   nuisance regressors. Either \code{NULL} (no confounds), a character vector
#'   of confound column names / patterns, or a \code{bidser} confound-set object.
#'   Resolved to actual values by \code{from_bids()} / \code{instantiate()}.
#' @param intercept Intercept handling passed to \code{baseline_model()}.
#' @param nuisance_check How to handle problematic nuisance regressors, passed
#'   to \code{baseline_model()}.
#' @return An object of class \code{baseline_spec}.
#' @seealso [fmri_template()], [baseline_model()]
#' @export
#' @examples
#' baseline_spec(degree = 3, basis = "bs")
baseline_spec <- function(degree = 3,
                          basis = c("bs", "poly", "ns", "constant"),
                          confounds = NULL,
                          intercept = c("runwise", "global", "none"),
                          nuisance_check = c("warn", "error", "drop", "none")) {
  basis <- match.arg(basis)
  intercept <- match.arg(intercept)
  nuisance_check <- match.arg(nuisance_check)
  assert_that(is.numeric(degree), length(degree) == 1, degree >= 0,
              msg = "'degree' must be a single non-negative number")
  # Fail fast at template-definition time, not on a worker node.
  if (basis %in% c("bs", "ns")) {
    assert_that(degree >= 3,
                msg = "'bs' and 'ns' drift bases require degree >= 3")
  }
  if (!is.null(confounds)) {
    assert_that(is.character(confounds) || inherits(confounds, "confound_set"),
                msg = "'confounds' must be NULL, a character vector, or a bidser confound_set")
  }
  structure(
    list(degree = degree, basis = basis, confounds = confounds,
         intercept = intercept, nuisance_check = nuisance_check),
    class = "baseline_spec"
  )
}

#' Define a subject-invariant model template
#'
#' Captures a complete model specification that is constant across subjects.
#' Per-subject data (BOLD paths, event table, run lengths, confounds) is bound
#' later with [instantiate()] to produce serializable [fmri_job] recipes that
#' can be executed locally, in parallel via \pkg{future}, or on a cluster.
#'
#' The template is pure, serializable data: it holds no voxel data and (by
#' design) no captured execution environment. The \code{reducer} is the one
#' field that can accidentally capture state; see \strong{Reducer
#' serializability} below.
#'
#' @section Reducer serializability:
#' A \code{reducer} should be a top-level / package function (or a closure that
#' captures nothing), so it survives serialization to a worker node. A closure
#' that captures local variables will drag those bindings along when the job is
#' written with \code{saveRDS()}; \code{fmri_template()} warns in that case.
#'
#' @param formula Event model formula, e.g. \code{onset ~ hrf(condition)}.
#' @param block Block / run-structure formula, e.g. \code{~ run}.
#' @param baseline A [baseline_spec()] describing the nuisance model.
#' @param durations Event durations passed through to the event model.
#' @param contrasts Optional contrast specification (e.g. from
#'   \code{contrast_set()}).
#' @param control An \code{fmri_lm_config} from [fmri_lm_control()] holding the
#'   robust / AR / preprocessing options for fitting.
#' @param strategy Fitting strategy: \code{"runwise"} or \code{"chunkwise"}.
#' @param engine Optional fitting engine name (see [register_engine()]).
#' @param engine_args Optional list of engine arguments.
#' @param reducer Optional function \code{function(fit, job)} that turns a fitted
#'   \code{fmri_lm} into the per-subject output (written to disk and/or returned).
#'   \code{NULL} (the default) returns the fitted object unchanged. See builtins
#'   such as [reduce_write_results()].
#' @return An object of class \code{fmri_template}.
#' @seealso [baseline_spec()], [instantiate()], [fmri_lm_control()]
#' @export
#' @examples
#' tmpl <- fmri_template(onset ~ hrf(condition), ~ run,
#'                       baseline = baseline_spec(degree = 3))
#' tmpl
fmri_template <- function(formula, block,
                          baseline = baseline_spec(),
                          durations = 0,
                          contrasts = NULL,
                          control = fmri_lm_control(),
                          strategy = c("runwise", "chunkwise"),
                          engine = NULL,
                          engine_args = list(),
                          reducer = NULL) {
  strategy <- match.arg(strategy)
  assert_that(inherits(formula, "formula"), msg = "'formula' must be a formula")
  assert_that(inherits(block, "formula"), msg = "'block' must be a formula")
  assert_that(inherits(baseline, "baseline_spec"),
              msg = "'baseline' must come from baseline_spec()")
  assert_that(is.numeric(durations), msg = "'durations' must be numeric")
  assert_that(inherits(control, "fmri_lm_config"),
              msg = "'control' must come from fmri_lm_control()")
  assert_that(is.list(engine_args), msg = "'engine_args' must be a list")
  if (!is.null(engine)) {
    assert_that(is.character(engine), length(engine) == 1,
                msg = "'engine' must be a single string or NULL")
  }
  if (!is.null(reducer)) {
    assert_that(is.function(reducer), msg = "'reducer' must be a function or NULL")
    warn_if_unserializable_fn(reducer, "reducer")
  }

  # Reset the formula/block environments to the global environment so jobs stay
  # compact and portable: a formula defined in a local scope would otherwise drag
  # that scope's bindings into every serialized job. Design variables resolve
  # against each subject's event table (passed as data), and term functions such
  # as hrf() against the loaded fmrireg namespace -- neither needs the captured
  # environment.
  environment(formula) <- globalenv()
  environment(block) <- globalenv()

  structure(
    list(formula = formula, block = block, baseline = baseline,
         durations = durations, contrasts = contrasts, control = control,
         strategy = strategy, engine = engine, engine_args = engine_args,
         reducer = reducer),
    class = "fmri_template"
  )
}

#' @keywords internal
#' @noRd
warn_if_unserializable_fn <- function(f, what = "function") {
  # Builtin reducers (see reducers.R) are trusted: they capture only small,
  # portable configuration.
  if (inherits(f, "fmri_reducer")) return(invisible(FALSE))
  e <- environment(f)
  if (is.null(e)) return(invisible(FALSE))          # primitive
  if (identical(e, globalenv()) || identical(e, baseenv()) ||
      identical(e, emptyenv()) || isNamespace(e)) {
    return(invisible(FALSE))
  }
  # A local closure environment: its bindings will be serialized with the job.
  if (length(ls(e, all.names = TRUE)) > 0) {
    warning(sprintf(
      paste0("'%s' is a closure that captures local variables; these will be ",
             "serialized with every job and may not be portable to a worker. ",
             "Prefer a top-level / package function."),
      what), call. = FALSE)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

#' Validate an fmri_template
#'
#' Re-checks the structural invariants of a template. Useful in preflight before
#' fanning a model out over many subjects.
#'
#' @param x An [fmri_template].
#' @return \code{TRUE} invisibly if valid; otherwise an error is raised.
#' @export
validate_template <- function(x) {
  assert_that(inherits(x, "fmri_template"), msg = "not an 'fmri_template'")
  assert_that(inherits(x$formula, "formula"), msg = "template formula is not a formula")
  assert_that(inherits(x$block, "formula"), msg = "template block is not a formula")
  assert_that(inherits(x$baseline, "baseline_spec"), msg = "template baseline is not a baseline_spec")
  assert_that(inherits(x$control, "fmri_lm_config"), msg = "template control is not an fmri_lm_config")
  if (!is.null(x$reducer)) {
    assert_that(is.function(x$reducer), msg = "template reducer is not a function")
    warn_if_unserializable_fn(x$reducer, "reducer")
  }
  invisible(TRUE)
}

#' @export
print.baseline_spec <- function(x, ...) {
  cat("<baseline_spec>\n")
  cat(sprintf("  basis: %s  degree: %s  intercept: %s\n",
              x$basis, x$degree, x$intercept))
  if (is.null(x$confounds)) {
    cat("  confounds: none\n")
  } else if (is.character(x$confounds)) {
    cat(sprintf("  confounds: %d column(s)/pattern(s)\n", length(x$confounds)))
  } else {
    cat("  confounds: <confound_set>\n")
  }
  invisible(x)
}

#' @export
print.fmri_template <- function(x, ...) {
  cat("<fmri_template>\n")
  cat("  formula:  ", deparse(x$formula), "\n")
  cat("  block:    ", deparse(x$block), "\n")
  cat("  strategy: ", x$strategy,
      if (!is.null(x$engine)) sprintf(" (engine: %s)", x$engine) else "", "\n", sep = "")
  cat("  baseline: ", sprintf("%s(degree=%s)", x$baseline$basis, x$baseline$degree), "\n")
  cat("  contrasts:", if (is.null(x$contrasts)) "none" else "set", "\n")
  cat("  reducer:  ", if (is.null(x$reducer)) "none (returns fitted object)" else "set", "\n")
  invisible(x)
}
