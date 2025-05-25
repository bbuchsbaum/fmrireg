#' Build fmri_config Object from Validated IOR List (Contextual Validation)
#'
#' This function performs context-dependent checks and constructs
#' an \code{fmri_config} object. DSL-201 covers loading the BIDS
#' project specified in the configuration and verifying that the
#' dataset path exists.
#'
#' @param validated_ior A list produced by [parse_and_validate_config()].
#' @return An object of class \code{fmri_config}.
#' @keywords internal
build_config_from_ior <- function(validated_ior) {
  if (!isTRUE(attr(validated_ior, "validated_schema"))) {
    stop("Input 'validated_ior' must come from parse_and_validate_config().")
  }

  errors <- ValidationErrors$new()

  bids_path <- fs::path_expand(validated_ior$dataset$path)

  if (!fs::dir_exists(bids_path)) {
    errors$add_error("dataset$path", paste0("BIDS directory does not exist: '", bids_path, "'"))
    errors$stop_if_invalid("Cannot proceed without a valid BIDS project")
  }

  project <- NULL
  tryCatch({
    project <- bidser::bids_project(bids_path)
  }, error = function(e) {
    errors$add_error("dataset$path", paste0("Failed to load BIDS project at '", bids_path, "': ", e$message))
  })

  errors$stop_if_invalid("Cannot proceed without a valid BIDS project")

  bids_check_level <- validated_ior$validation_settings$bids_content_checks %||% "Warn"

  ## DSL-202: Validate subject, task, and run selections against BIDS content
  avail_subjects <- tryCatch(bidser::participants(project), error = function(e) character())
  avail_tasks    <- tryCatch(bidser::tasks(project), error = function(e) character())
  avail_runs     <- tryCatch(bidser::runs(project),  error = function(e) character())

  ds <- validated_ior$dataset

  ## Subjects
  included_subs <- if (is.null(ds$subjects$include)) avail_subjects else sub("^sub-", "", ds$subjects$include)
  excluded_subs <- if (is.null(ds$subjects$exclude)) character() else sub("^sub-", "", ds$subjects$exclude)
  resolved_subs <- setdiff(included_subs, excluded_subs)
  missing_subs  <- setdiff(resolved_subs, avail_subjects)
  if (length(missing_subs) > 0) {
    msg <- paste0("Subjects not found in dataset: ", paste(missing_subs, collapse = ", "))
    if (identical(bids_check_level, "Error")) {
      errors$add_error("dataset$subjects", msg)
    } else if (identical(bids_check_level, "Warn")) {
      warning(msg, call. = FALSE)
    }
    resolved_subs <- intersect(resolved_subs, avail_subjects)
  }
  if (length(resolved_subs) == 0) {
    msg <- "No matching subjects found in BIDS dataset"
    if (identical(bids_check_level, "Error")) {
      errors$add_error("dataset$subjects", msg)
    } else if (identical(bids_check_level, "Warn")) {
      warning(msg, call. = FALSE)
    }
  }

  ## Tasks
  requested_tasks <- if (is.null(ds$tasks)) avail_tasks else sub("^task-", "", ds$tasks)
  missing_tasks   <- setdiff(requested_tasks, avail_tasks)
  if (length(missing_tasks) > 0) {
    msg <- paste0("Tasks not found in dataset: ", paste(missing_tasks, collapse = ", "))
    if (identical(bids_check_level, "Error")) {
      errors$add_error("dataset$tasks", msg)
    } else if (identical(bids_check_level, "Warn")) {
      warning(msg, call. = FALSE)
    }
    requested_tasks <- intersect(requested_tasks, avail_tasks)
  }

  ## Runs
  requested_runs <- if (is.null(ds$runs)) avail_runs else sub("^run-", "", ds$runs)
  if (length(avail_runs) > 0) {
    missing_runs <- setdiff(requested_runs, avail_runs)
    if (length(missing_runs) > 0) {
      msg <- paste0("Runs not found in dataset: ", paste(missing_runs, collapse = ", "))
      if (identical(bids_check_level, "Error")) {
        errors$add_error("dataset$runs", msg)
      } else if (identical(bids_check_level, "Warn")) {
        warning(msg, call. = FALSE)
      }
      requested_runs <- intersect(requested_runs, avail_runs)
    }
  }

  errors$stop_if_invalid("BIDS content validation failed")

  config <- list(
    spec      = validated_ior,
    project   = project,
    subjects  = resolved_subs,
    tasks     = requested_tasks,
    runs      = requested_runs,
    validated = TRUE
  )
  class(config) <- "fmri_config"
  config
}

#' Load, Validate, and Build fMRI Analysis Configuration
#'
#' This is the main user-facing helper for the DSL. It parses the YAML
#' file, applies defaults, validates against the schema, then calls
#' [build_config_from_ior()] to perform BIDS checks.
#'
#' @param yaml_file Path to the YAML configuration file.
#' @return An \code{fmri_config} object.
#' @export
load_fmri_config <- function(yaml_file) {
  validated_ior <- parse_and_validate_config(yaml_file)
  build_config_from_ior(validated_ior)
}
