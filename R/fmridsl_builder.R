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
#' Load Events and Confounds for a Subject
#'
#' Helper used by the DSL model builder. Given a validated
#' `fmri_config` and a subject identifier, this loads all matching
#' event files and any requested confound variables. Run lengths and
#' the effective TR are derived from the configuration scan parameters.
#'
#' @param config Validated \code{fmri_config} object.
#' @param subject_id Subject ID string.
#'
#' @return List with elements \code{events_df}, \code{confounds_df},
#'   \code{run_lengths}, and \code{TR}.
#' @noRd
#' @keywords internal
load_and_prepare_subject_data <- function(config, subject_id) {
  if (!inherits(config, "fmri_config") || !isTRUE(config$validated)) {
    stop("'config' must be a validated fmri_config")
  }
  if (!subject_id %in% config$subjects) {
    stop(sprintf("Subject '%s' not listed in configuration.", subject_id))
  }

  tasks <- config$tasks
  runs  <- config$runs

  events_res <- bidser::read_events(
    config$project,
    subid = subject_id,
    task  = if (length(tasks) > 0) paste(tasks, collapse = "|") else NULL,
    run   = if (length(runs) > 0) paste(runs, collapse = "|") else NULL
  )

  if (nrow(events_res) == 0) {
    warning("No events found for subject", call. = FALSE)
    events_df <- data.frame()
  } else {
    events_df <- tidyr::unnest(events_res, data)
    events_df <- dplyr::arrange(events_df, .data$.run)
  }

  sp <- config$spec$dataset$scan_params
  default_TR <- sp$TR %||% NA_real_
  TR_over <- sp$TR_overrides %||% list()

  effective_TR <- default_TR
  if (length(TR_over) > 0 && nrow(events_df) > 0) {
    run_tags <- paste0(
      "sub-", subject_id, "_task-", events_df$.task,
      "_run-", sprintf("%02d", as.numeric(events_df$.run))
    )
    for (pat in names(TR_over)) {
      if (any(grepl(pat, run_tags))) {
        effective_TR <- TR_over[[pat]]
        break
      }
    }
  }
  if (is.na(effective_TR)) {
    stop("TR could not be determined for subject")
  }

  default_rl <- sp$run_lengths %||% list()
  rl_over <- sp$run_length_overrides %||% list()
  run_info <- unique(events_df[, c(".run", ".task")])
  final_rl <- integer(nrow(run_info))

  for (i in seq_len(nrow(run_info))) {
    rstr <- sprintf("run-%02d", as.numeric(run_info$.run[i]))
    tstr <- run_info$.task[i]
    matched <- FALSE
    for (pat in names(rl_over)) {
      test <- paste0("sub-", subject_id, "_task-", tstr, "_", rstr)
      if (grepl(pat, test)) {
        final_rl[i] <- rl_over[[pat]]
        matched <- TRUE
        break
      }
    }
    if (!matched) {
      final_rl[i] <- default_rl[[tstr]] %||% NA_integer_
    }
  }
  if (any(is.na(final_rl))) {
    stop("Run lengths could not be determined for all runs")
  }

  conf_columns <- unique(unlist(config$confounds_info$groups, use.names = FALSE))
  if (length(conf_columns) > 0) {
    conf_res <- bidser::read_confounds(
      config$project,
      subid = subject_id,
      task  = if (length(tasks) > 0) paste(tasks, collapse = "|") else NULL,
      run   = if (length(runs) > 0) paste(runs, collapse = "|") else NULL
    )
    if (nrow(conf_res) == 0) {
      confounds_df <- NULL
    } else {
      confounds_df <- tidyr::unnest(conf_res, data)
      confounds_df <- dplyr::arrange(confounds_df, as.numeric(run))
      confounds_df <- dplyr::select(confounds_df, dplyr::all_of(conf_columns))
    }
  } else {
    confounds_df <- NULL
  }

  list(
    events_df = events_df,
    confounds_df = confounds_df,
    run_lengths = as.integer(final_rl),
    TR = as.numeric(effective_TR)
  )
}

#' Instantiate an HRF Object from DSL Definition
#'
#' Given the name of an HRF defined in an `fmri_config`, this helper
#' returns the corresponding `HRF` object with any specified global
#' decorators applied. If the `hrfs` block is omitted from the
#' configuration, only a default `canonical` HRF (SPM canonical) is
#' available.
#'
#' @param hrf_name Name of the HRF to retrieve.
#' @param fmri_config A validated `fmri_config` object.
#'
#' @return An object of class `HRF`.
#' @keywords internal
get_hrf_from_dsl <- function(hrf_name, fmri_config) {
  defs <- fmri_config$hrfs %||% list(canonical = list(type = "SPMCanonical"))

  if (!hrf_name %in% names(defs)) {
    stop(sprintf("HRF '%s' not defined in configuration", hrf_name), call. = FALSE)
  }

  def <- defs[[hrf_name]]
  type <- def$type
  params <- def$parameters %||% list()

  base <- switch(type,
    SPMCanonical = HRF_SPMG1,
    SPMCanonicalDerivs = {
      derivs <- def$derivatives %||% character()
      if ("Dispersion" %in% derivs) {
        HRF_SPMG3
      } else {
        HRF_SPMG2
      }
    },
    GammaFunction = as_hrf(
      function(t) hrf_gamma(t,
        shape = params$shape %||% 6,
        rate  = params$rate  %||% 1
      ),
      name = "gamma",
      params = list(shape = params$shape %||% 6, rate = params$rate %||% 1)
    ),
    Gaussian = as_hrf(
      function(t) hrf_gaussian(t,
        mean = params$mean %||% 6,
        sd   = params$sd   %||% 2
      ),
      name = "gaussian",
      params = list(mean = params$mean %||% 6, sd = params$sd %||% 2)
    ),
    BSplineBasisHRF = do.call(hrfspline_generator, params),
    TentBasisHRF    = do.call(hrf_tent_generator, params),
    FourierBasisHRF = do.call(hrf_fourier_generator, params),
    DaguerreBasisHRF = do.call(hrf_daguerre_generator, params),
    CustomR = {
      fun_name <- def$definition
      if (is.null(fun_name)) {
        stop("CustomR HRF requires a 'definition' field", call. = FALSE)
      }
      fun <- NULL
      if (is.function(fun_name)) {
        fun <- fun_name
      } else if (exists(fun_name, mode = "function")) {
        fun <- get(fun_name, mode = "function")
      } else if (file.exists(fun_name)) {
        env <- new.env(parent = baseenv())
        sys.source(fun_name, envir = env)
        if (exists(fun_name, envir = env, mode = "function")) {
          fun <- get(fun_name, envir = env)
        } else {
          stop(sprintf("Function '%s' not found after sourcing", fun_name), call. = FALSE)
        }
      } else {
        stop(sprintf("Cannot resolve custom HRF function '%s'", fun_name), call. = FALSE)
      }
      as_hrf(fun, name = fun_name)
    },
    stop(sprintf("Unknown HRF type '%s'", type))
  )

  lag      <- def$lag %||% 0
  width    <- def$width %||% 0
  summate  <- def$summate %||% TRUE
  normalize <- def$normalize %||% FALSE

  gen_hrf(base, lag = lag, width = width,
          summate = summate, normalize = normalize)
}

#' Apply Transformation Operations
#'
#' Internal helper used by the DSL builder to process variable
#' transformations defined in an `fmri_config` object.
#'
#' @param x Vector or factor to transform.
#' @param ops List of operations (strings or lists) specifying the
#'   transformations to apply in order.
#' @param env Environment containing additional variables referenced by
#'   transformation ops.
#' @param block Optional grouping vector used for group-based operations
#'   like `demean-by-group`.
#' @return Transformed object.
#' @keywords internal
apply_transform_ops <- function(x, ops, env = parent.frame(), block = NULL) {
  for (op in ops) {
    if (is.character(op)) {
      if (op == "center") {
        mu <- mean(x, na.rm = TRUE)
        x <- x - mu
      } else if (op == "scale-sd") {
        sdv <- stats::sd(x, na.rm = TRUE)
        if (is.na(sdv) || sdv == 0) sdv <- 1e-6
        x <- x / sdv
      } else if (op == "z-score") {
        mu <- mean(x, na.rm = TRUE)
        sdv <- stats::sd(x, na.rm = TRUE)
        if (is.na(sdv) || sdv == 0) sdv <- 1e-6
        x <- (x - mu) / sdv
      } else if (op == "log") {
        x <- log(x)
      } else if (op == "exp") {
        x <- exp(x)
      } else if (op == "factorize") {
        x <- as.factor(x)
      } else if (op == "demean-by-group") {
        if (is.null(block)) {
          stop("demean-by-group requires grouping variable", call. = FALSE)
        }
        mus <- tapply(x, block, mean, na.rm = TRUE)
        x <- x - mus[as.character(block)]
      }
    } else if (is.list(op)) {
      type <- op$type
      if (type == "scale-within-group") {
        gname <- op$group_by_variable
        g <- env[[gname]]
        if (is.null(g)) {
          stop(sprintf("group_by_variable '%s' not found", gname), call. = FALSE)
        }
        mus <- tapply(x, g, mean, na.rm = TRUE)
        sds <- tapply(x, g, sd, na.rm = TRUE)
        sds[is.na(sds) | sds == 0] <- 1e-6
        x <- mapply(function(val, grp) {
          mu <- mus[as.character(grp)]
          sdv <- sds[as.character(grp)]
          (val - mu) / sdv
        }, x, g)
      } else if (type == "clip") {
        if (!is.null(op$min)) x <- pmax(op$min, x)
        if (!is.null(op$max)) x <- pmin(op$max, x)
      } else if (type == "recode-levels") {
        lvlmap <- op$level_map %||% list()
        if (!is.factor(x)) x <- factor(x)
        levs <- levels(x)
        for (nm in names(lvlmap)) {
          levs[levs == nm] <- lvlmap[[nm]]
        }
        levels(x) <- levs
      }
    }
  }
  x
}

#' Create Parametric Basis from DSL Definition
#'
#' @param basis_def Basis definition list from the DSL.
#' @param x Numeric vector used to build the basis.
#' @return A `ParametricBasis` object.
#' @keywords internal
create_basis_from_dsl <- function(basis_def, x) {
  params <- basis_def$parameters %||% list()
  type <- basis_def$type
  switch(type,
    Polynomial = Poly(x, degree = params$degree %||% 1),
    BSpline    = BSpline(x, degree = params$degree %||% 3),
    Standardized = Standardized(x),
    Identity   = Ident(x),
    stop(sprintf("Unknown basis type '%s'", type))
  )
}

#' Process Variables and Transformations into Environment
#'
#' Implements DSL2-301. Creates a unified environment containing all
#' raw variables, derived variables from the `transformations` block,
#' and basis-expanded variables for any parametric modulators.
#'
#' @param fmri_config Validated `fmri_config` object.
#' @param subject_data List returned by `load_and_prepare_subject_data()`.
#' @return Environment with model variables available by name.
#' @keywords internal
process_variables_and_transformations <- function(fmri_config, subject_data) {
  env <- new.env(parent = baseenv())

  events_df   <- subject_data$events_df
  conf_df     <- subject_data$confounds_df
  var_defs    <- fmri_config$spec$variables
  block_var   <- fmri_config$spec$events$block_column

  # Raw variables
  for (nm in names(var_defs)) {
    def <- var_defs[[nm]]
    col <- def$bids_column
    val <- if (identical(def$role, "NuisanceSource")) {
      if (is.null(conf_df)) NULL else conf_df[[col]]
    } else {
      events_df[[col]]
    }
    if (identical(def$role, "Factor")) {
      val <- as.factor(val)
    } else if (identical(def$role, "Numeric")) {
      val <- as.numeric(val)
    }
    env[[nm]] <- val
  }

  # Derived variables
  trans_defs <- fmri_config$spec$transformations %||% list()
  if (length(trans_defs) > 0) {
    for (nm in names(trans_defs)) {
      tdef <- trans_defs[[nm]]
      src  <- env[[tdef$source_variable]]
      env[[nm]] <- apply_transform_ops(src, tdef$ops, env, events_df[[block_var]])
    }
  }

  # Basis-expanded modulators
  term_defs <- fmri_config$spec$terms %||% list()
  for (tnm in names(term_defs)) {
    term <- term_defs[[tnm]]
    if (!is.null(term$modulator_basis)) {
      mod_var <- term$mod_var
      val <- env[[mod_var]]
      if (is.null(val)) next
      basis_obj <- create_basis_from_dsl(term$modulator_basis, val)
      deg_part <- term$modulator_basis$parameters$degree %||% ""
      name <- paste0(mod_var, "_", tolower(term$modulator_basis$type),
                     if (!identical(deg_part, "")) paste0("_deg", deg_part) else "")
      env[[name]] <- basis_obj$y
    }
  }

  env
}
#' Convert DSL Terms to hrfspec/covariatespec Objects
#'
#' Implements DSL2-401. For a given model definition and processed
#' variable environment, this helper creates a named list of term
#' specification objects suitable for `event_model()`.
#'
#' @param fmri_config Validated `fmri_config` object.
#' @param model A single model definition from `fmri_config$spec$models`.
#' @param var_env Environment returned by `process_variables_and_transformations()`.
#'
#' @return Named list of `hrfspec` and `covariatespec` objects.
#' @keywords internal
convert_terms_to_specs <- function(fmri_config, model, var_env) {
  stopifnot(inherits(fmri_config, "fmri_config"))
  term_defs <- fmri_config$spec$terms %||% list()
  term_names <- model$terms %||% character()

  res <- list()

  for (tnm in term_names) {
    tdef <- term_defs[[tnm]]
    if (is.null(tdef)) {
      stop(sprintf("Term '%s' not defined in configuration", tnm), call. = FALSE)
    }

    subset_expr <- if (!is.null(tdef$subset)) rlang::parse_expr(tdef$subset) else NULL

    if (identical(tdef$type, "NuisanceRegressors")) {
      vars <- tdef$nuisance_source_variables
      data_df <- as.data.frame(mget(vars, envir = var_env))
      spec <- covariate(!!!rlang::syms(vars), data = data_df, id = tnm, subset = subset_expr)
      res[[tnm]] <- spec
      next
    }

    hrf_name <- tdef$hrf %||% "canonical"
    hobj <- get_hrf_from_dsl(hrf_name, fmri_config$spec)
    lag_val <- tdef$lag %||% 0
    if (!is.null(lag_val) && lag_val != 0) {
      hobj <- gen_hrf(hobj, lag = lag_val)
    }

    if (identical(tdef$type, "ParametricModulation")) {
      mod_name <- tdef$mod_var
      if (!is.null(tdef$modulator_basis)) {
        deg_part <- tdef$modulator_basis$parameters$degree %||% ""
        mod_name <- paste0(mod_name, "_", tolower(tdef$modulator_basis$type),
                            if (!identical(deg_part, "")) paste0("_deg", deg_part) else "")
      }
      vars <- c(tdef$selector_vars, mod_name)
    } else if (identical(tdef$type, "EventRelated")) {
      vars <- tdef$event_variables
    } else if (identical(tdef$type, "Trialwise")) {
      spec <- trialwise(basis = hobj, lag = 0, id = tnm)
      if (!is.null(subset_expr)) spec$subset <- subset_expr
      res[[tnm]] <- spec
      next
    } else {
      stop(sprintf("Unknown term type '%s'", tdef$type), call. = FALSE)
    }

    expr <- rlang::expr(hrf(!!!rlang::syms(vars), basis = hobj, id = !!tnm, subset = !!subset_expr))
    spec <- rlang::eval_bare(expr, var_env)
    res[[tnm]] <- spec
  }

  res
}

#' Parse Baseline Basis Specification
#'
#' Helper for DSL2-501. Converts a baseline basis definition from the DSL
#' (either a shorthand string like "BSpline(3)" or a structured list)
#' into the arguments required by [baseline_model()].
#'
#' @param basis_def Baseline basis definition (character or list).
#' @return List with elements `basis` and `degree`.
#' @keywords internal
parse_baseline_basis <- function(basis_def) {
  if (is.null(basis_def)) {
    return(list(basis = "bs", degree = 3L))
  }
  if (is.list(basis_def)) {
    type <- basis_def$type %||% "BSpline"
    degree <- basis_def$parameters$degree %||% 3L
  } else if (is.character(basis_def)) {
    m <- regexec("^([A-Za-z]+)\\s*(?:\\((\\d+)\\))?$", basis_def)
    reg <- regmatches(basis_def, m)[[1]]
    if (length(reg) == 0) {
      stop(sprintf("Cannot parse basis specification '%s'", basis_def), call. = FALSE)
    }
    type <- reg[2]
    degree <- if (length(reg) >= 3 && nzchar(reg[3])) as.integer(reg[3]) else 3L
  } else {
    stop("Invalid basis specification", call. = FALSE)
  }

  type_map <- c(Polynomial = "poly", BSpline = "bs", NSpline = "ns", Constant = "constant")
  b <- type_map[[type]]
  if (is.null(b)) {
    stop(sprintf("Unknown baseline basis type '%s'", type), call. = FALSE)
  }
  if (identical(b, "constant")) degree <- 1L
  list(basis = b, degree = as.integer(degree))
}

#' Construct baseline_model from DSL model definition
#'
#' Implements DSL2-501. Uses the configuration, model specification,
#' and loaded subject data to create a [baseline_model()] object.
#'
#' @param fmri_config Validated `fmri_config` object.
#' @param model Single model definition from `fmri_config$spec$models`.
#' @param subject_data Result from [load_and_prepare_subject_data()].
#' @return A `baseline_model` object.
#' @keywords internal
build_baseline_model_from_dsl <- function(fmri_config, model, subject_data) {
  stopifnot(inherits(fmri_config, "fmri_config"))
  bl <- model$baseline %||% list()
  basis_info <- parse_baseline_basis(bl$basis %||% "BSpline(3)")
  intercept_map <- c(PerRun = "runwise", Global = "global", None = "none")
  intercept <- intercept_map[[bl$intercept %||% "PerRun"]]

  sframe <- sampling_frame(blocklens = subject_data$run_lengths,
                           TR = subject_data$TR)

  groups <- bl$include_confound_groups %||% list()
  group_cols <- unique(unlist(fmri_config$confounds_info$groups[groups], use.names = FALSE))

  nuisance_list <- NULL
  if (!is.null(subject_data$confounds_df) && length(group_cols) > 0) {
    cf_df <- subject_data$confounds_df
    cf_df <- dplyr::select(cf_df, dplyr::all_of(group_cols), run)
    nuisance_list <- split(as.matrix(cf_df[, group_cols, drop = FALSE]), cf_df$run)
  }

  baseline_model(basis = basis_info$basis,
                 degree = basis_info$degree,
                 sframe = sframe,
                 intercept = intercept,
                 nuisance_list = nuisance_list)
}

#' Construct event_model from DSL term specifications
#'
#' Implements DSL2-502. Given a list of `hrfspec` and `covariatespec`
#' objects produced by [convert_terms_to_specs()], this helper calls
#' [event_model()] with the appropriate data extracted from the
#' configuration and subject events.
#'
#' @param term_specs Named list of term specification objects.
#' @param fmri_config Validated `fmri_config` object.
#' @param subject_data List returned by [load_and_prepare_subject_data()].
#' @return An `event_model` object.
#' @keywords internal
build_event_model_from_dsl <- function(term_specs,
                                       fmri_config,
                                       subject_data) {
  stopifnot(is.list(term_specs))
  stopifnot(inherits(fmri_config, "fmri_config"))

  ev_defs <- fmri_config$spec$events
  on_col   <- ev_defs$onset_column
  dur_col  <- ev_defs$duration_column
  block_col <- ev_defs$block_column

  ev_df <- subject_data$events_df

  data_df <- data.frame(
    onset = ev_df[[on_col]],
    duration = ev_df[[dur_col]],
    block = ev_df[[block_col]],
    stringsAsFactors = FALSE
  )

  sframe <- sampling_frame(blocklens = subject_data$run_lengths,
                           TR = subject_data$TR)

  event_model(term_specs,
              data = data_df,
              block = data_df$block,
              sampling_frame = sframe,
              durations = data_df$duration)
}
