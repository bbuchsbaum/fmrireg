#' Parse a YAML Configuration File
#'
#' Reads a YAML file from disk and converts it into an R list. The
#' function performs basic error handling for missing files or malformed
#' YAML content.
#'
#' @param yaml_file_path Path to the YAML configuration file.
#'
#' @return A named list representing the contents of the YAML file.
#' @examples
#' \dontrun{
#' config_list <- parse_yaml_to_list("analysis.yaml")
#' }
#' @export
parse_yaml_to_list <- function(yaml_file_path) {
  if (!fs::file_exists(yaml_file_path)) {
    stop("YAML configuration file not found: ", yaml_file_path, call. = FALSE)
  }

  parsed <- NULL
  tryCatch({
    parsed <- yaml::read_yaml(yaml_file_path)
  }, error = function(e) {
    stop("Failed to parse YAML file '", yaml_file_path, "'. Error: ", e$message,
         call. = FALSE)
  })

  if (is.null(parsed) || !is.list(parsed)) {
    stop("YAML file '", yaml_file_path, "' could not be parsed into a list.",
         call. = FALSE)
  }

  parsed
}
#' ValidationErrors R6 class
#'
#' Collects validation errors encountered while parsing a configuration.
#' Each error is stored with the path to the offending field and a message.
#'
#' @section Methods:
#' \describe{
#'   \item{add_error(path, message)}{Record an error at a specific configuration path.}
#'   \item{is_valid()}{Return TRUE if no errors have been recorded.}
#'   \item{format_errors()}{Produce a human-readable summary of all errors.}
#'   \item{stop_if_invalid(context)}{Stop execution with formatted errors if any exist.}
#' }
#'
#' @importFrom R6 R6Class
#' @export
ValidationErrors <- R6::R6Class(
  "ValidationErrors",
  public = list(
    #' @field errors list of recorded errors
    errors = list(),

    #' @description
    #' Add an error entry
    #' @param path configuration path for the error
    #' @param message description of the error
    add_error = function(path, message) {
      self$errors <- c(self$errors, list(list(path = path, message = message)))
      invisible(self)
    },

    #' @description
    #' Determine if any errors have been recorded
    #' @return logical
    is_valid = function() {
      length(self$errors) == 0
    },

    #' @description
    #' Format errors for printing
    #' @return character string summarizing all errors
    format_errors = function() {
      if (self$is_valid()) {
        return("No validation errors.")
      }
      error_messages <- vapply(
        self$errors,
        function(err) {
          prefix <- if (is.null(err$path) || err$path == "") "[Top Level]" else paste0("At '", err$path, "'")
          paste0("  - ", prefix, ": ", err$message)
        },
        character(1)
      )
      paste(
        "Configuration validation failed with",
        length(error_messages), "error(s):\n",
        paste(error_messages, collapse = "\n")
      )
    },

    #' @description
    #' Stop with formatted errors if any exist
    #' @param context message prefix for the error
    stop_if_invalid = function(context = "Configuration is invalid") {
      if (!self$is_valid()) {
        stop(context, ":\n", self$format_errors(), call. = FALSE)
      }
      invisible(self)
    }
  )
)

#' Recursively Apply Default Values to a Configuration List
#'
#' Merges a list of defaults into a configuration list. For any element
#' present in `defaults` but missing from `data`, the default value is
#' inserted. If both `data` and `defaults` contain a given element and
#' both values are lists, the merging is performed recursively.
#'
#' @param data A list parsed from the user YAML configuration.
#' @param defaults A list of default values following the DSL specification.
#'
#' @return A list with defaults applied where values were missing.
#' @keywords internal
#' @noRd
apply_defaults <- function(data, defaults) {
  if (!is.list(data) || !is.list(defaults)) return(data)

  data_names     <- names(data) %||% character()
  defaults_names <- names(defaults) %||% character()

  # Add missing keys from defaults
  missing_keys <- setdiff(defaults_names, data_names)
  for (key in missing_keys) {
    data[[key]] <- defaults[[key]]
  }

  # Recurse into common list elements
  common_keys <- intersect(data_names, defaults_names)
  for (key in common_keys) {
    if (is.list(defaults[[key]]) && !is.null(data[[key]]) && is.list(data[[key]])) {
      if (length(defaults[[key]]) > 0 || length(data[[key]]) > 0) {
        data[[key]] <- apply_defaults(data[[key]], defaults[[key]])
      }
    }
  }

  data
}

#' Parse, Apply Defaults, and Validate a DSL Configuration
#'
#' Reads a YAML configuration file, merges in default values as
#' specified by the DSL, and performs basic schema validation for the
#' required top-level sections and the \code{dataset} block (DSL-105).
#'
#' Validation checks ensure required fields exist and have the correct
#' basic types. The \code{dataset$scan_params} sub-fields are also
#' flattened so that override objects become named vectors.
#'
#' @param yaml_file Path to the YAML configuration file.
#'
#' @return A list representing the validated configuration with
#'   defaults applied.  An error is thrown if validation fails.
#' @keywords internal
parse_and_validate_config <- function(yaml_file) {
  config_list <- parse_yaml_to_list(yaml_file)

  defaults <- list(
    dataset = list(
      relpath = "func",
      subjects = list(include = NULL, exclude = NULL),
      tasks = NULL,
      runs = NULL,
      scan_params = list(
        TR = NULL,
        TR_overrides = list(),
        run_lengths = list(),
        run_length_overrides = list()
      )
    ),
    events = list(
      onset_column = "onset",
      duration_column = "duration",
      block_column = "run"
    ),
    hrfs = list(
      canonical = list(type = "SPMCanonical")
    )
  )

  config_list <- apply_defaults(config_list, defaults)

  errors <- ValidationErrors$new()

  required_sections <- c("dataset", "events", "variables", "terms", "models")
  for (sec in required_sections) {
    check_required(config_list, sec, "", errors)
  }

  if (check_type(config_list, "dataset", "object", "", errors)) {
    ds <- config_list$dataset
    ds_path <- "dataset"

    if (check_required(ds, "path", ds_path, errors)) {
      check_type(ds, "path", "string", ds_path, errors)
    }

    if (exists("relpath", ds)) {
      rp <- ds$relpath
      if (is.character(rp)) {
        if (length(rp) > 1) {
          check_type(ds, "relpath", "array[string]", ds_path, errors)
        }
      } else {
        check_type(ds, "relpath", "array[string]", ds_path, errors)
        if (!is.null(ds$relpath)) ds$relpath <- unlist(ds$relpath)
      }
    }

    if (check_type(ds, "subjects", "object", ds_path, errors, allow_null = TRUE)) {
      if (!is.null(ds$subjects)) {
        sub_path <- paste0(ds_path, "$subjects")
        check_type(ds$subjects, "include", "array[string]", sub_path, errors, allow_null = TRUE)
        check_pattern(ds$subjects, "include", "^sub-[0-9A-Za-z]+$", sub_path, errors, allow_null = TRUE)
        check_type(ds$subjects, "exclude", "array[string]", sub_path, errors, allow_null = TRUE)
        check_pattern(ds$subjects, "exclude", "^sub-[0-9A-Za-z]+$", sub_path, errors, allow_null = TRUE)
      }
    }

    if (check_type(ds, "tasks", "array[string]", ds_path, errors, allow_null = TRUE)) {
      check_pattern(ds, "tasks", "^task-[a-zA-Z0-9]+$", ds_path, errors, allow_null = TRUE)
    }

    if (check_type(ds, "runs", "array[string]", ds_path, errors, allow_null = TRUE)) {
      check_pattern(ds, "runs", "^run-[0-9]+$", ds_path, errors, allow_null = TRUE)
    }

    if (check_type(ds, "scan_params", "object", ds_path, errors, allow_null = TRUE)) {
      sp <- ds$scan_params
      sp_path <- paste0(ds_path, "$scan_params")
      if (!is.null(sp)) {
        check_type(sp, "TR", "number", sp_path, errors, allow_null = TRUE)

        if (check_type(sp, "TR_overrides", "object", sp_path, errors, allow_null = TRUE)) {
          if (!is.null(sp$TR_overrides)) {
            for (nm in names(sp$TR_overrides)) {
              check_type(sp$TR_overrides, nm, "number", paste0(sp_path, "$TR_overrides$", nm), errors)
            }
            sp$TR_overrides <- unlist(sp$TR_overrides)
          }
        }

        if (check_type(sp, "run_lengths", "object", sp_path, errors, allow_null = TRUE)) {
          if (!is.null(sp$run_lengths)) {
            for (nm in names(sp$run_lengths)) {
              check_type(sp$run_lengths, nm, "integer", paste0(sp_path, "$run_lengths$", nm), errors)
            }
            sp$run_lengths <- unlist(sp$run_lengths)
          }
        }

        if (check_type(sp, "run_length_overrides", "object", sp_path, errors, allow_null = TRUE)) {
          if (!is.null(sp$run_length_overrides)) {
            for (nm in names(sp$run_length_overrides)) {
              check_type(sp$run_length_overrides, nm, "integer", paste0(sp_path, "$run_length_overrides$", nm), errors)
            }
            sp$run_length_overrides <- unlist(sp$run_length_overrides)
          }
        }
        ds$scan_params <- sp
      }
    }

    config_list$dataset <- ds
  }

  if (check_type(config_list, "events", "object", "", errors)) {
    ev <- config_list$events
    ev_path <- "events"

    if (check_required(ev, "onset_column", ev_path, errors)) {
      check_type(ev, "onset_column", "string", ev_path, errors)
    }

    if (check_required(ev, "duration_column", ev_path, errors)) {
      check_type(ev, "duration_column", "string", ev_path, errors)
    }

    if (check_required(ev, "block_column", ev_path, errors)) {
      check_type(ev, "block_column", "string", ev_path, errors)
    }

    config_list$events <- ev
  }

  if (exists("hrfs", config_list)) {
    if (check_type(config_list, "hrfs", "object", "", errors, allow_null = TRUE)) {
      hrfs <- config_list$hrfs
      hrfs_path <- "hrfs"
      if (!is.null(hrfs)) {
        for (nm in names(hrfs)) {
          if (check_type(hrfs, nm, "object", hrfs_path, errors)) {
            h <- hrfs[[nm]]
            h_path <- paste0(hrfs_path, "$", nm)

            if (check_required(h, "type", h_path, errors)) {
              check_enum(
                h,
                "type",
                c(
                  "SPMCanonical",
                  "SPMCanonicalDerivs",
                  "GammaFunction",
                  "Gaussian",
                  "BSplineBasisHRF",
                  "TentBasisHRF",
                  "FourierBasisHRF",
                  "DaguerreBasisHRF",
                  "CustomR"
                ),
                h_path,
                errors
              )
            }

            if (check_type(h, "derivatives", "array[string]", h_path, errors, allow_null = TRUE)) {
              check_enum(h, "derivatives", c("Temporal", "Dispersion"), h_path, errors, allow_null = TRUE)
            }

            check_type(h, "parameters", "object", h_path, errors, allow_null = TRUE)
            check_type(h, "definition", "string", h_path, errors, allow_null = TRUE)
          }
        }
      }
    }
  }

  if (check_type(config_list, "variables", "object", "", errors)) {
    vars <- config_list$variables
    vars_path <- "variables"
    for (nm in names(vars)) {
      if (check_type(vars, nm, "object", vars_path, errors)) {
        v <- vars[[nm]]
        v_path <- paste0(vars_path, "$", nm)
        if (check_required(v, "bids_column", v_path, errors)) {
          check_type(v, "bids_column", "string", v_path, errors)
        }
        if (check_required(v, "role", v_path, errors)) {
          check_enum(
            v,
            "role",
            c("Factor", "Numeric", "NuisanceSource", "TrialIndex", "GroupID"),
            v_path,
            errors
          )
        }
      }
    }
    config_list$variables <- vars
  }

  if (exists("transformations", config_list)) {
    if (check_type(config_list, "transformations", "object", "", errors, allow_null = TRUE)) {
      trans <- config_list$transformations
      trans_path <- "transformations"
      if (!is.null(trans)) {
        for (nm in names(trans)) {
          if (check_type(trans, nm, "object", trans_path, errors)) {
            tentry <- trans[[nm]]
            t_path <- paste0(trans_path, "$", nm)

            if (check_required(tentry, "source_variable", t_path, errors)) {
              check_type(tentry, "source_variable", "string", t_path, errors)
            }

            if (check_required(tentry, "ops", t_path, errors)) {
              ops <- tentry$ops
              op_field_path <- paste0(t_path, "$ops")
              if (is.null(ops) || !(is.list(ops) || (is.vector(ops) && !is.matrix(ops)))) {
                errors$add_error(op_field_path, "Field must be an array of operations.")
              } else {
                for (i in seq_along(ops)) {
                  op <- ops[[i]]
                  op_path <- paste0(op_field_path, "[", i, "]")
                  if (is.character(op) && length(op) == 1) {
                    check_enum(list(x = op), "x",
                              c("center", "scale-sd", "z-score", "log", "exp",
                                "factorize", "demean-by-group"),
                              op_path, errors)
                  } else if (is.list(op)) {
                    if (check_required(op, "type", op_path, errors)) {
                      check_enum(op, "type",
                                c("scale-within-group", "clip", "recode-levels"),
                                op_path, errors)
                    }
                    check_type(op, "group_by_variable", "string", op_path, errors, allow_null = TRUE)
                    check_type(op, "min", "number", op_path, errors, allow_null = TRUE)
                    check_type(op, "max", "number", op_path, errors, allow_null = TRUE)
                    check_type(op, "level_map", "object", op_path, errors, allow_null = TRUE)
                  } else {
                    errors$add_error(op_path, "Each op must be a string or object.")
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  if (exists("confound_groups", config_list)) {
    if (check_type(config_list, "confound_groups", "object", "", errors, allow_null = TRUE)) {
      cg <- config_list$confound_groups
      cg_path <- "confound_groups"
      if (!is.null(cg)) {
        for (nm in names(cg)) {
          if (check_type(cg, nm, "object", cg_path, errors)) {
            grp <- cg[[nm]]
            g_path <- paste0(cg_path, "$", nm)
            have_pat  <- exists("select_by_pattern", grp)
            have_col  <- exists("select_by_bids_column", grp)
            check_type(grp, "select_by_pattern", "array[string]", g_path, errors, allow_null = TRUE)
            check_type(grp, "select_by_bids_column", "array[string]", g_path, errors, allow_null = TRUE)
            if (!have_pat && !have_col) {
              errors$add_error(g_path, "Must specify 'select_by_pattern' or 'select_by_bids_column'.")
            }
          }
        }
      }
    }
  }

  if (check_type(config_list, "terms", "object", "", errors)) {
    terms <- config_list$terms
    terms_path <- "terms"
    for (nm in names(terms)) {
      if (check_type(terms, nm, "object", terms_path, errors)) {
        t <- terms[[nm]]
        t_path <- paste0(terms_path, "$", nm)

        if (check_required(t, "type", t_path, errors)) {
          check_enum(
            t,
            "type",
            c("EventRelated", "ParametricModulation", "Trialwise", "NuisanceRegressors"),
            t_path,
            errors
          )
        }

        # conditional required fields
        if (!is.null(t$type)) {
          if (t$type %in% c("EventRelated", "Trialwise")) {
            if (check_required(t, "event_variables", t_path, errors)) {
              check_type(t, "event_variables", "array[string]", t_path, errors)
            }
          } else if (t$type == "ParametricModulation") {
            if (check_required(t, "selector_vars", t_path, errors)) {
              check_type(t, "selector_vars", "array[string]", t_path, errors)
            }
            if (check_required(t, "mod_var", t_path, errors)) {
              check_type(t, "mod_var", "string", t_path, errors)
            }
          } else if (t$type == "NuisanceRegressors") {
            if (check_required(t, "nuisance_source_variables", t_path, errors)) {
              check_type(t, "nuisance_source_variables", "array[string]", t_path, errors)
            }
          }
        }

        if (exists("event_variables", t)) {
          check_type(t, "event_variables", "array[string]", t_path, errors, allow_null = TRUE)
        }
        if (exists("selector_vars", t)) {
          check_type(t, "selector_vars", "array[string]", t_path, errors, allow_null = TRUE)
        }
        if (exists("mod_var", t)) {
          check_type(t, "mod_var", "string", t_path, errors, allow_null = TRUE)
        }
        if (exists("nuisance_source_variables", t)) {
          check_type(t, "nuisance_source_variables", "array[string]", t_path, errors, allow_null = TRUE)
        }

        if (exists("hrf", t)) {
          check_type(t, "hrf", "string", t_path, errors, allow_null = TRUE)
        } else {
          t$hrf <- "canonical"
        }

        if (exists("subset", t)) {
          check_type(t, "subset", "string", t_path, errors, allow_null = TRUE)
        }

        if (exists("lag", t)) {
          check_type(t, "lag", "number", t_path, errors, allow_null = TRUE)
        } else {
          t$lag <- 0
        }

        if (exists("modulator_basis", t)) {
          if (check_type(t, "modulator_basis", "object", t_path, errors)) {
            mb <- t$modulator_basis
            mb_path <- paste0(t_path, "$modulator_basis")
            if (check_required(mb, "type", mb_path, errors)) {
              check_enum(
                mb,
                "type",
                c("Polynomial", "BSpline", "Standardized", "Identity", "NSpline"),
                mb_path,
                errors
              )
            }
            if (check_type(mb, "parameters", "object", mb_path, errors, allow_null = TRUE)) {
              if (!is.null(mb$type) && mb$type == "Polynomial") {
                params_path <- paste0(mb_path, "$parameters")
                if (check_required(mb$parameters, "degree", params_path, errors)) {
                  check_type(mb$parameters, "degree", "integer", params_path, errors)
                }
              }
            }
          }
        }

        terms[[nm]] <- t
      }
    }

    config_list$terms <- terms
  }

  if (exists("contrasts", config_list)) {
    if (check_type(config_list, "contrasts", "object", "", errors, allow_null = TRUE)) {
      cons <- config_list$contrasts
      cons_path <- "contrasts"
      if (!is.null(cons)) {
        for (nm in names(cons)) {
          if (check_type(cons, nm, "object", cons_path, errors)) {
            cdef <- cons[[nm]]
            c_path <- paste0(cons_path, "$", nm)

            if (check_required(cdef, "type", c_path, errors)) {
              check_enum(
                cdef,
                "type",
                c(
                  "Formula", "Pair", "OneAgainstAll",
                  "Unit", "Oneway", "Interaction",
                  "Polynomial", "ColumnRegex"
                ),
                c_path,
                errors
              )
            }

            check_type(cdef, "expression", "string", c_path, errors, allow_null = TRUE)
            check_type(cdef, "factors", "array[string]", c_path, errors, allow_null = TRUE)
            check_type(cdef, "where", "string", c_path, errors, allow_null = TRUE)
          }
        }
      }
    }
  }

  if (check_type(config_list, "models", "array[object]", "", errors)) {
    mods <- config_list$models
    mods_path <- "models"
    if (length(mods) == 0) {
      errors$add_error(mods_path, "At least one model must be specified.")
    }
    for (i in seq_along(mods)) {
      m <- mods[[i]]
      m_path <- paste0(mods_path, "[", i, "]")
      if (is.list(m)) {
        if (check_required(m, "name", m_path, errors)) {
          check_type(m, "name", "string", m_path, errors)
          check_pattern(list(name = m$name), "name", "^(?!(?:_meta|_generated)).*$", m_path, errors)
        }

        if (exists("baseline", m)) {
          if (check_type(m, "baseline", "object", m_path, errors, allow_null = TRUE)) {
            bl <- m$baseline
            bl_path <- paste0(m_path, "$baseline")
            if (exists("basis", bl)) {
              check_type(bl, "basis", "string", bl_path, errors, allow_null = TRUE)
            } else {
              bl$basis <- "BSpline(3)"
            }
            if (exists("intercept", bl)) {
              check_enum(bl, "intercept", c("PerRun", "Global", "None"), bl_path, errors, allow_null = TRUE)
            } else {
              bl$intercept <- "PerRun"
            }
            if (exists("include_confound_groups", bl)) {
              check_type(bl, "include_confound_groups", "array[string]", bl_path, errors, allow_null = TRUE)
            } else {
              bl$include_confound_groups <- list()
            }
            m$baseline <- bl
          }
        } else {
          m$baseline <- list(basis = "BSpline(3)", intercept = "PerRun", include_confound_groups = list())
        }

        if (check_required(m, "terms", m_path, errors)) {
          check_type(m, "terms", "array[string]", m_path, errors)
        }
        if (exists("contrasts", m)) {
          check_type(m, "contrasts", "array[string]", m_path, errors, allow_null = TRUE)
        } else {
          m$contrasts <- list()
        }
        mods[[i]] <- m
      } else {
        errors$add_error(m_path, "Each model must be an object.")
      }
    }
    config_list$models <- mods
  }

  if (exists("default_model", config_list)) {
    check_type(config_list, "default_model", "string", "", errors, allow_null = TRUE)
  }

  if (exists("validation_settings", config_list)) {
    if (check_type(config_list, "validation_settings", "object", "", errors, allow_null = TRUE)) {
      vs <- config_list$validation_settings
      vs_path <- "validation_settings"
      if (!is.null(vs)) {
        if (exists("cross_references", vs)) {
          check_enum(vs, "cross_references", c("Error", "Warn", "Off"), vs_path, errors, allow_null = TRUE)
        } else {
          vs$cross_references <- "Error"
        }
        if (exists("bids_content_checks", vs)) {
          check_enum(vs, "bids_content_checks", c("Error", "Warn", "Off"), vs_path, errors, allow_null = TRUE)
        } else {
          vs$bids_content_checks <- "Warn"
        }
        if (exists("allow_unknown_yaml_keys", vs)) {
          check_type(vs, "allow_unknown_yaml_keys", "boolean", vs_path, errors, allow_null = TRUE)
        } else {
          vs$allow_unknown_yaml_keys <- FALSE
        }
      }
      config_list$validation_settings <- vs
    }
  }

  errors$stop_if_invalid("Configuration validation failed")

  config_list
}

