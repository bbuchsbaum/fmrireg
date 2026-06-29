###############################################################################
## fmri_job.R
##
## A serializable per-subject "job" recipe: an `fmri_template` bound to one
## subject's data via a `dataset_spec`. A job holds no voxel data and no
## captured execution environment, so it can be written with saveRDS() and
## shipped to a worker (a future, or a node in an array job) that never saw the
## driver session. The worker reconstructs the dataset + model lazily and fits.
##
## See `docs/plans/model-templates-multisubject-fanout.md`.
###############################################################################

#' Describe how to (re)construct a subject's dataset
#'
#' A small, serializable recipe for a dataset: the name of a dataset constructor
#' plus the arguments to call it with. For file-backed data the arguments are
#' paths and run lengths (no voxel data), which keeps the enclosing [fmri_job]
#' tiny and portable. The dataset is realized lazily on the worker.
#'
#' @param constructor Name of a dataset constructor (a string), e.g.
#'   \code{"fmri_dataset"} or \code{"matrix_dataset"}. Resolved at run time, so
#'   the data is not loaded when the spec is built.
#' @param args A named list of arguments passed to \code{constructor} (for
#'   \code{"fmri_dataset"}: \code{scans}, \code{TR}, \code{run_length},
#'   \code{event_table}, \code{mask}, \code{base_path}, ...).
#' @param source Either \code{"file"} (paths; nothing loaded until run) or
#'   \code{"inline"} (data already in \code{args}, e.g. a \code{matrix_dataset}).
#' @return An object of class \code{dataset_spec}.
#' @seealso [fmri_job()], [fmri_dataset()], [matrix_dataset()]
#' @export
#' @examples
#' dataset_spec("fmri_dataset",
#'              args = list(scans = c("run-1_bold.nii.gz", "run-2_bold.nii.gz"),
#'                          TR = 2, run_length = c(200, 200)),
#'              source = "file")
dataset_spec <- function(constructor, args = list(), source = c("file", "inline")) {
  source <- match.arg(source)
  assert_that(is.character(constructor), length(constructor) == 1, nzchar(constructor),
              msg = "'constructor' must be a single non-empty function name")
  assert_that(is.list(args), msg = "'args' must be a (named) list")
  structure(
    list(constructor = constructor, args = args, source = source),
    class = "dataset_spec"
  )
}

#' Construct a per-subject job recipe
#'
#' Binds an [fmri_template] to one subject's [dataset_spec]. The result is a
#' fully serializable recipe (no voxel data, no captured environment) that
#' [run()] turns into a fitted model and reduced output. Build jobs with
#' [instantiate()] rather than by hand in the common case.
#'
#' @param id A unique job identifier (e.g. a subject label).
#' @param template An [fmri_template].
#' @param dataset_spec A [dataset_spec] describing this subject's data.
#' @param meta Optional named list of metadata used for output keying
#'   (e.g. \code{list(subject = "01", task = "stroop", space = "MNI152...")}).
#' @param nuisance Optional per-subject nuisance / confound regressors fed to the
#'   baseline model at run time: \code{NULL}, a numeric matrix with one row per
#'   scan (split across runs), or a list of per-run matrices. (A deferred
#'   file-read spec is also accepted by [run()].)
#' @return An object of class \code{fmri_job}.
#' @seealso [instantiate()], [run()], [export_jobs()]
#' @export
#' @examples
#' tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
#' ds <- dataset_spec("fmri_dataset",
#'                    args = list(scans = "run-1_bold.nii.gz", TR = 2,
#'                                run_length = 200), source = "file")
#' fmri_job("sub-01", tmpl, ds, meta = list(subject = "01"))
fmri_job <- function(id, template, dataset_spec, meta = list(), nuisance = NULL) {
  assert_that(is.character(id), length(id) == 1, nzchar(id),
              msg = "'id' must be a single non-empty string")
  # ids become filenames in export_jobs()/run_one.R, so keep them path-safe.
  assert_that(!grepl("[/\\\\]", id),
              msg = "'id' must not contain path separators ('/' or '\\\\')")
  assert_that(!id %in% c(".", ".."),
              msg = "'id' must not be '.' or '..'")
  assert_that(inherits(template, "fmri_template"),
              msg = "'template' must be an 'fmri_template'")
  assert_that(inherits(dataset_spec, "dataset_spec"),
              msg = "'dataset_spec' must be a 'dataset_spec'")
  assert_that(is.list(meta), msg = "'meta' must be a list")
  structure(
    list(id = id, template = template, dataset_spec = dataset_spec,
         meta = meta, nuisance = nuisance),
    class = "fmri_job"
  )
}

#' Validate an fmri_job
#'
#' @param x An [fmri_job].
#' @return \code{TRUE} invisibly if valid; otherwise an error is raised.
#' @export
validate_job <- function(x) {
  assert_that(inherits(x, "fmri_job"), msg = "not an 'fmri_job'")
  assert_that(is.character(x$id), length(x$id) == 1, nzchar(x$id),
              msg = "job 'id' must be a single non-empty string")
  validate_template(x$template)
  assert_that(inherits(x$dataset_spec, "dataset_spec"),
              msg = "job 'dataset_spec' is not a 'dataset_spec'")
  invisible(TRUE)
}

#' @export
print.dataset_spec <- function(x, ...) {
  cat("<dataset_spec>\n")
  cat(sprintf("  constructor: %s()  source: %s\n", x$constructor, x$source))
  cat(sprintf("  args: %s\n", paste(names(x$args), collapse = ", ")))
  invisible(x)
}

#' @export
print.fmri_job <- function(x, ...) {
  cat(sprintf("<fmri_job> %s\n", x$id))
  cat(sprintf("  dataset: %s() [%s]\n", x$dataset_spec$constructor, x$dataset_spec$source))
  if (length(x$meta)) {
    kv <- vapply(names(x$meta), function(k) sprintf("%s=%s", k, as.character(x$meta[[k]])[1]),
                 character(1))
    cat("  meta:   ", paste(kv, collapse = ", "), "\n")
  }
  cat("  formula:", deparse(x$template$formula), "\n")
  invisible(x)
}
