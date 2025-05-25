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

  # Initialize info containers
  events_info <- list(columns = character(), mapping = validated_ior$events)
  confounds_info <- list(available_columns = character(), groups = list())

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

  ## DSL-203: Validate existence of essential event columns
  if (length(resolved_subs) > 0 && length(requested_tasks) > 0) {
    rep_sub  <- resolved_subs[1]
    rep_task <- requested_tasks[1]
    rep_run  <- if (length(requested_runs) > 0) requested_runs[1] else NULL

    rep_events_df <- NULL
    tryCatch({
      ev_res <- bidser::read_events(
        project,
        subid = rep_sub,
        task  = rep_task,
        run   = rep_run
      )
      if (nrow(ev_res) > 0 && length(ev_res$data) > 0 && !is.null(ev_res$data[[1]])) {
        rep_events_df <- ev_res$data[[1]]
        events_info$columns <- names(rep_events_df)
      }
    }, error = function(e) {
      msg <- paste0("Failed to load events for representative subject/task/run: ", e$message)
      if (identical(bids_check_level, "Error")) {
        errors$add_error("events", msg)
      } else if (identical(bids_check_level, "Warn")) {
        warning(msg, call. = FALSE)
      }
    })

    if (!is.null(rep_events_df)) {
      essential_cols <- c(
        validated_ior$events$onset_column,
        validated_ior$events$duration_column,
        validated_ior$events$block_column
      )
      missing_cols <- setdiff(essential_cols, names(rep_events_df))
      if (length(missing_cols) > 0) {
        msg <- paste0(
          "Required event columns missing: ",
          paste(missing_cols, collapse = ", ")
        )
        if (identical(bids_check_level, "Error")) {
          errors$add_error("events", msg)
        } else if (identical(bids_check_level, "Warn")) {
          warning(msg, call. = FALSE)
        }
      }
    } else {
      msg <- "No events data found for representative subject/task/run"
      if (identical(bids_check_level, "Error")) {
        errors$add_error("events", msg)
      } else if (identical(bids_check_level, "Warn")) {
        warning(msg, call. = FALSE)
      }
    }

    ## DSL-204: Validate variables -> BIDS column mapping
    rep_confounds_df <- NULL
    tryCatch({
      conf_res <- bidser::read_confounds(
        project,
        subid = rep_sub,
        task  = rep_task,
        run   = rep_run
      )
      if (nrow(conf_res) > 0 && length(conf_res$data) > 0 && !is.null(conf_res$data[[1]])) {
        rep_confounds_df <- conf_res$data[[1]]
        confounds_info$available_columns <- names(rep_confounds_df)
      }
    }, error = function(e) {
      msg <- paste0("Failed to load confounds for representative subject/task/run: ", e$message)
      if (identical(bids_check_level, "Error")) {
        errors$add_error("confounds", msg)
      } else if (identical(bids_check_level, "Warn")) {
        warning(msg, call. = FALSE)
      }
    })

    event_cols <- if (!is.null(rep_events_df)) names(rep_events_df) else character()
    conf_cols  <- if (!is.null(rep_confounds_df)) names(rep_confounds_df) else character()

    for (var_nm in names(validated_ior$variables)) {
      var <- validated_ior$variables[[var_nm]]
      bcol <- var$bids_column
      path <- paste0("variables$", var_nm, "$bids_column")
      if (identical(var$role, "NuisanceSource")) {
        if (!(bcol %in% conf_cols)) {
          msg <- paste0("Confound column '", bcol, "' not found in confounds file.")
          if (identical(bids_check_level, "Error")) {
            errors$add_error(path, msg)
          } else if (identical(bids_check_level, "Warn")) {
            warning(msg, call. = FALSE)
          }
        }
      } else {
        if (!(bcol %in% event_cols)) {
          msg <- paste0("Event column '", bcol, "' not found in events file.")
          if (identical(bids_check_level, "Error")) {
            errors$add_error(path, msg)
          } else if (identical(bids_check_level, "Warn")) {
            warning(msg, call. = FALSE)
          }
        }
      }
    }

    ## DSL-205: Confound group resolution and baseline checks
    resolved_groups <- list()
    if (!is.null(validated_ior$confound_groups)) {
      for (grp_nm in names(validated_ior$confound_groups)) {
        grp <- validated_ior$confound_groups[[grp_nm]]
        matches <- character()
        if (!is.null(grp$select_by_pattern)) {
          for (pat in grp$select_by_pattern) {
            matches <- union(matches, grep(pat, conf_cols, value = TRUE))
          }
        }
        if (!is.null(grp$select_by_bids_column)) {
          missing_cols <- setdiff(grp$select_by_bids_column, conf_cols)
          if (length(missing_cols) > 0) {
            msg <- paste0(
              "Confound column(s) not found for group '", grp_nm,
              "': ", paste(missing_cols, collapse = ", ")
            )
            if (identical(bids_check_level, "Error")) {
              errors$add_error(paste0("confound_groups$", grp_nm), msg)
            } else if (identical(bids_check_level, "Warn")) {
              warning(msg, call. = FALSE)
            }
          }
          matches <- union(matches, intersect(grp$select_by_bids_column, conf_cols))
        }
        if (length(matches) == 0) {
          msg <- paste0("No confound columns resolved for group '", grp_nm, "'")
          if (identical(bids_check_level, "Error")) {
            errors$add_error(paste0("confound_groups$", grp_nm), msg)
          } else if (identical(bids_check_level, "Warn")) {
            warning(msg, call. = FALSE)
          }
        }
        resolved_groups[[grp_nm]] <- matches
        confounds_info$groups[[grp_nm]] <- matches
      }
    }

    cr_level <- validated_ior$validation_settings$cross_references %||% "Error"
    for (i in seq_along(validated_ior$models)) {
      model <- validated_ior$models[[i]]
      groups <- model$baseline$include_confound_groups %||% list()
      for (g in groups) {
        path <- paste0("models[", i, "]$baseline$include_confound_groups")
        if (!(g %in% names(resolved_groups))) {
          msg <- paste0("Confound group '", g, "' not defined")
          if (identical(cr_level, "Error")) {
            errors$add_error(path, msg)
          } else if (identical(cr_level, "Warn")) {
            warning(msg, call. = FALSE)
          }
        } else if (length(resolved_groups[[g]]) == 0) {
          msg <- paste0("Confound group '", g, "' resolved to no columns")
          if (identical(bids_check_level, "Error")) {
            errors$add_error(path, msg)
          } else if (identical(bids_check_level, "Warn")) {
            warning(msg, call. = FALSE)
          }
        }
      }
    }
  }

  errors$stop_if_invalid("BIDS content validation failed")

  roles_vec <- vapply(validated_ior$variables, `[[`, character(1), "role")
  variable_roles <- split(names(roles_vec), roles_vec)
  variable_roles <- lapply(variable_roles, unname)

  config <- list(
    spec      = validated_ior,
    project   = project,
    subjects  = resolved_subs,
    tasks     = requested_tasks,
    runs      = requested_runs,
    events_info = events_info,
    confounds_info = confounds_info,
    variable_roles = variable_roles,
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
