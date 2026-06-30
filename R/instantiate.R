###############################################################################
## instantiate.R
##
## Bind an `fmri_template` to per-subject data ("bindings") to produce
## serializable `fmri_job` recipes, plus the realize helpers that reconstruct a
## dataset and assemble the full `fmri_model` on the worker at run time.
##
## A "binding" is a named list describing one analysis unit:
##   id          (required) unique identifier / subject label
##   scans       (required) a numeric matrix (inline) OR character paths (file)
##   TR          (required) repetition time
##   run_length  (required) integer vector of timepoints per run
##   events      (optional) event-table data.frame (design + block columns)
##   confounds   (optional) nuisance regressors (matrix or per-run list)
##   mask        (optional) mask path/object (file scans only)
##   base_path   (optional) base path for file scans
##   meta        (optional) named list for output keying (subject/task/space)
##
## See `docs/plans/model-templates-multisubject-fanout.md`.
###############################################################################

#' Instantiate a template into per-subject jobs
#'
#' Binds an [fmri_template] to one or more per-subject data bindings, producing
#' serializable [fmri_job] recipes. No data is loaded: file-backed scans are kept
#' as paths and realized lazily by [run()].
#'
#' @param template An [fmri_template].
#' @param x A single binding (a named list with at least \code{id}), a list of
#'   bindings, or a manifest \code{data.frame} (one row per analysis unit, with
#'   list-columns for \code{scans}/\code{events}/\code{run_length}/etc.).
#' @param ... Unused.
#' @return A single [fmri_job] if \code{x} is one binding, otherwise a list of
#'   [fmri_job]s.
#' @seealso [fmri_job()], [run()], [from_bids()]
#' @export
#' @examples
#' tmpl <- fmri_template(onset ~ hrf(condition), ~ run)
#' b <- list(id = "sub-01", scans = matrix(rnorm(80 * 3), 80, 3),
#'           TR = 2, run_length = c(40, 40),
#'           events = data.frame(onset = c(5, 25, 45, 65),
#'                               condition = factor(c("A", "B", "A", "B")),
#'                               run = c(1, 1, 2, 2)))
#' job <- instantiate(tmpl, b)
instantiate <- function(template, x, ...) {
  assert_that(inherits(template, "fmri_template"),
              msg = "'template' must be an 'fmri_template'")
  single <- .is_single_binding(x)
  bindings <- .as_binding_list(x)
  jobs <- lapply(bindings, function(b) .binding_to_job(template, b))
  if (single) jobs[[1]] else jobs
}

#' @keywords internal
#' @noRd
.is_single_binding <- function(x) {
  is.list(x) && !is.data.frame(x) && !is.null(x[["id"]])
}

#' @keywords internal
#' @noRd
.as_binding_list <- function(x) {
  if (is.data.frame(x)) {
    # `col[[i]]` already unwraps a list-column element (matrix / data.frame /
    # per-run list) and yields the scalar for an atomic column.
    return(lapply(seq_len(nrow(x)), function(i) lapply(x, function(col) col[[i]])))
  }
  if (.is_single_binding(x)) return(list(x))
  assert_that(is.list(x), msg = "'x' must be a binding, a list of bindings, or a manifest data.frame")
  x
}

#' @keywords internal
#' @noRd
.binding_to_job <- function(template, b) {
  assert_that(!is.null(b[["id"]]), msg = "each binding needs an 'id'")
  assert_that(!is.null(b[["scans"]]), msg = "each binding needs 'scans'")
  assert_that(!is.null(b[["TR"]]), msg = "each binding needs 'TR'")
  assert_that(!is.null(b[["run_length"]]), msg = "each binding needs 'run_length'")

  events <- b[["events"]]
  if (is.null(events)) events <- data.frame()

  scans <- b[["scans"]]
  if (is.character(scans)) {
    args <- list(scans = scans, TR = b[["TR"]],
                 run_length = b[["run_length"]], event_table = events)
    if (!is.null(b[["mask"]])) args$mask <- b[["mask"]]
    if (!is.null(b[["base_path"]])) args$base_path <- b[["base_path"]]
    spec <- dataset_spec("fmri_dataset", args = args, source = "file")
  } else if (is.matrix(scans) || is.numeric(scans)) {
    spec <- dataset_spec(
      "matrix_dataset",
      args = list(datamat = as.matrix(scans), TR = b[["TR"]],
                  run_length = b[["run_length"]], event_table = events),
      source = "inline"
    )
  } else {
    stop("binding 'scans' must be a numeric matrix (inline) or character paths (file)",
         call. = FALSE)
  }

  meta <- b[["meta"]]
  if (is.null(meta)) meta <- list(subject = as.character(b[["id"]]))

  fmri_job(id = as.character(b[["id"]]), template = template,
           dataset_spec = spec, meta = meta, nuisance = b[["confounds"]])
}

#' Realize the dataset described by a job
#'
#' Reconstructs the \code{fmri_dataset} from the job's [dataset_spec]. For
#' file-backed specs this is where data first becomes addressable (still lazily,
#' per the dataset backend).
#'
#' @param job An [fmri_job].
#' @return An \code{fmri_dataset}.
#' @seealso [build_model()], [run()]
#' @export
realize_dataset <- function(job) {
  assert_that(inherits(job, "fmri_job"), msg = "'job' must be an 'fmri_job'")
  spec <- job$dataset_spec
  # Allowlist the constructor: a job is a recipe that may have travelled from
  # disk, so only known dataset constructors may be invoked by name.
  if (!spec$constructor %in% .allowed_dataset_constructors) {
    stop(sprintf("dataset constructor '%s' is not allowed; expected one of: %s",
                 spec$constructor,
                 paste(.allowed_dataset_constructors, collapse = ", ")),
         call. = FALSE)
  }
  fn <- match.fun(spec$constructor)
  do.call(fn, spec$args)
}

#' @keywords internal
#' @noRd
.allowed_dataset_constructors <- c(
  "matrix_dataset", "fmri_dataset", "fmri_mem_dataset", "latent_dataset"
)

#' @keywords internal
#' @noRd
# Apply the baseline_spec confound *selector* to the per-subject confound values.
# `selector` is a character vector of column names; when given, the named columns
# are kept (and missing ones are an error), so the model includes only the
# intended invariant subset rather than every supplied confound.
.select_confounds <- function(nuisance, selector) {
  if (is.null(nuisance) || is.null(selector)) return(nuisance)
  if (!is.character(selector)) return(nuisance)  # e.g. a bidser confound_set used by from_bids()
  pick <- function(m) {
    m <- as.matrix(m)
    cn <- colnames(m)
    if (is.null(cn)) {
      stop("baseline confound selector was supplied but the subject's confounds have no column names",
           call. = FALSE)
    }
    miss <- setdiff(selector, cn)
    if (length(miss) > 0) {
      stop(sprintf("confound column(s) not found in subject data: %s",
                   paste(miss, collapse = ", ")), call. = FALSE)
    }
    m[, selector, drop = FALSE]
  }
  if (is.list(nuisance) && !is.data.frame(nuisance)) lapply(nuisance, pick) else pick(nuisance)
}

#' @keywords internal
#' @noRd
.resolve_nuisance <- function(nuisance, run_length) {
  if (is.null(nuisance)) return(NULL)
  if (is.list(nuisance) && !is.data.frame(nuisance)) return(nuisance) # already per-run
  m <- as.matrix(nuisance)
  if (nrow(m) != sum(run_length)) {
    stop(sprintf("nuisance rows (%d) do not match total scans (%d)",
                 nrow(m), sum(run_length)), call. = FALSE)
  }
  idx <- rep(seq_along(run_length), run_length)
  lapply(split(seq_len(nrow(m)), idx), function(i) m[i, , drop = FALSE])
}

#' Assemble the full fMRI model for a job
#'
#' Builds the per-subject \code{baseline_model} (from the template's
#' [baseline_spec], the dataset's sampling frame, and any nuisance regressors)
#' and combines it with the event model into an \code{fmri_model}.
#'
#' @param job An [fmri_job].
#' @param dataset The realized dataset (defaults to \code{realize_dataset(job)}).
#' @return An \code{fmri_model}.
#' @seealso [realize_dataset()], [run()]
#' @export
build_model <- function(job, dataset = realize_dataset(job)) {
  assert_that(inherits(job, "fmri_job"), msg = "'job' must be an 'fmri_job'")
  tmpl <- job$template
  sf <- dataset$sampling_frame
  bs <- tmpl$baseline
  nz <- .select_confounds(job$nuisance, bs$confounds)
  nl <- .resolve_nuisance(nz, fmrihrf::blocklens(sf))
  bmodel <- baseline_model(
    basis = bs$basis, degree = bs$degree, sframe = sf,
    intercept = bs$intercept, nuisance_list = nl,
    nuisance_check = bs$nuisance_check
  )
  create_fmri_model(
    formula = tmpl$formula, block = tmpl$block,
    baseline_model = bmodel, dataset = dataset,
    durations = tmpl$durations
  )
}
