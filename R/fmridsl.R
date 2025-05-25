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

  errors$stop_if_invalid("Configuration validation failed")

  config_list
}

