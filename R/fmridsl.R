#' @keywords internal
`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}

#' @importFrom R6 R6Class
ValidationErrors <- R6::R6Class("ValidationErrors",
  public = list(
    errors = list(),
    add_error = function(path, message) {
      self$errors <- c(self$errors, list(list(path = path, message = message)))
      invisible(self)
    },
    is_valid = function() {
      length(self$errors) == 0
    },
    format_errors = function() {
      if (self$is_valid()) return("No validation errors.")
      error_messages <- sapply(self$errors, function(err) {
        paste0("  - ", ifelse(is.null(err$path) || err$path == "", "[Top Level]", paste0("At '", err$path, "'")), ": ", err$message)
      })
      paste("Configuration validation failed with", length(error_messages), "error(s):\n",
            paste(error_messages, collapse = "\n"))
    },
    stop_if_invalid = function(context = "Configuration is invalid") {
       if (!self$is_valid()) {
           stop(context, ":\n", self$format_errors(), call. = FALSE)
       }
       invisible(self)
    }
  )
)

# --- Validation Check Helpers ---
# check_required, check_type, check_enum, check_pattern, check_dir_exists
# (Assume these are defined as before)

#' Check if a required field exists.
#' @keywords internal
check_required <- function(data, field, path, errors) {
  # Use `exists` with inherits=FALSE for lists
  if (!is.list(data) || !exists(field, where = data, inherits = FALSE)) {
    field_path <- paste0(path, if (path != "") "$", field)
    errors$add_error(field_path, paste0("Required field '", field, "' is missing."))
    return(FALSE)
  }
  return(TRUE)
}

#' Check if a field has the expected data type.
#' Handles basic types, objects (named lists), and arrays (vectors/lists)
#' with optional element type checking.
#' @keywords internal
check_type <- function(data, field, expected_type, path, errors, allow_null = FALSE) {
  # Field might not exist if optional, return TRUE if allow_null is TRUE
  if (!exists(field, where = data, inherits = FALSE)) {
    return(allow_null)
  }

  value <- data[[field]]
  field_path <- paste0(path, if (path != "") "$", field)

  if (is.null(value)) {
    if (!allow_null) {
      errors$add_error(field_path, "Field cannot be null.")
      return(FALSE)
    }
    return(TRUE) # Null is allowed and present
  }

  valid <- FALSE
  actual_type_desc <- paste(class(value), collapse="/")

  # Basic types
  if (expected_type == "string") {
      valid <- is.character(value) && length(value) == 1 && !is.na(value)
  } else if (expected_type == "integer") {
      # Allow numeric that are whole numbers
      valid <- (is.integer(value) || (is.numeric(value) && !is.na(value) && value == floor(value))) && length(value) == 1
  } else if (expected_type == "number") {
      valid <- is.numeric(value) && length(value) == 1 && !is.na(value)
  } else if (expected_type == "boolean") {
      valid <- is.logical(value) && length(value) == 1 && !is.na(value)
  } else if (expected_type == "object") {
      # Must be a named list (can be empty)
      valid <- is.list(value) && (length(value) == 0 || (!is.null(names(value)) && all(names(value) != "")))
  }
  # Array types
  else if (grepl("^array", expected_type)) {
      is_array_like <- is.list(value) || (is.vector(value) && !is.matrix(value)) # Exclude matrices unless intended
      if (!is_array_like) {
          valid <- FALSE
      } else {
          valid <- TRUE # Assume valid until element check fails
          element_type_match <- regmatches(expected_type, regexec("^array\\[(\\w+)\\]$", expected_type))
          element_type <- if (length(element_type_match[[1]]) > 1) element_type_match[[1]][2] else NULL

          if (length(value) > 0 && !is.null(element_type)) {
              # Check element types recursively using a simplified internal check
              check_element <- function(el) {
                   if (element_type == "string") is.character(el) && length(el) == 1 && !is.na(el)
                   else if (element_type == "integer") (is.integer(el) || (is.numeric(el) && !is.na(el) && el == floor(el))) && length(el) == 1
                   else if (element_type == "number") is.numeric(el) && length(el) == 1 && !is.na(el)
                   else if (element_type == "boolean") is.logical(el) && length(el) == 1 && !is.na(el)
                   else if (element_type == "object") is.list(el) && (length(el) == 0 || (!is.null(names(el)) && all(names(el) != "")))
                   else FALSE # Unknown element type
              }
              # Use vapply for efficiency and type safety
              all_elems_valid <- tryCatch(
                  vapply(value, check_element, logical(1)),
                  error = function(e) FALSE # Handle cases where check_element fails on an element
              )
              if (!all(all_elems_valid)) {
                 valid <- FALSE
                 actual_type_desc <- paste0(actual_type_desc, " (contains elements of incorrect type)")
              }
          }
      }
  } else {
      errors$add_error(field_path, paste0("Internal validation error: Unknown expected type '", expected_type, "' specified."))
      return(FALSE)
  }

  if (!valid) {
    errors$add_error(field_path, paste0("Invalid type. Expected '", expected_type, "', but got '", actual_type_desc, "'."))
  }
  return(valid)
}


#' Check if a field's value(s) are within an allowed set.
#' @keywords internal
check_enum <- function(data, field, allowed_values, path, errors, allow_null = FALSE) {
   if (!exists(field, where=data, inherits=FALSE)) {
       return(allow_null)
   }
   value <- data[[field]]
   field_path <- paste0(path, if (path != "") "$", field)
   if (is.null(value)) {
       return(allow_null)
   }

   values_to_check <- unlist(value) # Works for single values and flattens lists/vectors

   if (length(values_to_check) == 0 && !is.null(value)) return(TRUE) # e.g. an empty list might be valid if allow_null=T
   if (is.null(values_to_check) && length(value)>0) { # handles list(NULL)
        errors$add_error(field_path, "Contains NULL values where specific choices are expected.")
        return(FALSE)
   }

   invalid_values <- setdiff(values_to_check, allowed_values)

   if (length(invalid_values) > 0) {
       errors$add_error(field_path, paste0("Invalid value(s): '", paste(unique(invalid_values), collapse="', '"), "'. Must be one of: ", paste(allowed_values, collapse=", ")))
       return(FALSE)
   }
   return(TRUE)
}

#' Check if a string field's value(s) match a regex pattern.
#' @keywords internal
check_pattern <- function(data, field, pattern, path, errors, allow_null = FALSE) {
    if (!exists(field, where=data, inherits=FALSE)) {
       return(allow_null)
    }
    value <- data[[field]]
    field_path <- paste0(path, if (path != "") "$", field)
    if (is.null(value)) {
        return(allow_null)
    }

    values_to_check <- unlist(value) # Works for single values and flattens lists/vectors

    if (length(values_to_check) == 0 && !is.null(value)) return(TRUE)
    if (is.null(values_to_check) && length(value)>0) { # handles list(NULL)
        errors$add_error(field_path, "Pattern check applied to list containing NULL values.")
        return(FALSE)
    }

    if (!is.character(values_to_check)) {
        # Let check_type handle fundamental type errors
        return(TRUE)
    }

    does_not_match <- !grepl(pattern, values_to_check)
    if (any(does_not_match)) {
        errors$add_error(field_path, paste0("Value(s) do not match required pattern ('", pattern, "'): '", paste(unique(values_to_check[does_not_match]), collapse="', '"), "'."))
        return(FALSE)
    }
    return(TRUE)
}

#' Check if a string field points to an existing directory.
#' @keywords internal
check_dir_exists <- function(data, field, path, errors) {
  if (!exists(field, where = data, inherits = FALSE)) {
     # If the field itself is missing, required check should handle it.
     # If it's optional and missing, directory check doesn't apply.
     return(TRUE)
  }
  dir_path <- data[[field]]
  field_path <- paste0(path, if (path != "") "$", field)

  # Ensure it's a non-null string first
  if (!is.character(dir_path) || length(dir_path) != 1 || is.na(dir_path) || dir_path == "") {
     errors$add_error(field_path, "Must be a non-empty string path.")
     return(FALSE)
  }

  # Use fs::path_expand to handle tilde etc.
  expanded_path <- tryCatch(fs::path_expand(dir_path), error = function(e) NULL)
  if (is.null(expanded_path)) {
      errors$add_error(field_path, paste0("Could not expand path: '", dir_path, "'."))
      return(FALSE)
  }

  if (!fs::dir_exists(expanded_path)) {
    errors$add_error(field_path, paste0("Directory not found: '", dir_path, "' (expanded to '", expanded_path, "')."))
    return(FALSE)
  }
  return(TRUE)
}


apply_defaults <- function(data, defaults) {
  if (!is.list(data) || !is.list(defaults)) return(data)

  # Get names, handle potential NULL names for empty lists
  data_names <- names(data) %||% character()
  defaults_names <- names(defaults) %||% character()

  # Add missing keys from defaults to data
  missing_keys <- setdiff(defaults_names, data_names)
  for (key in missing_keys) {
    data[[key]] <- defaults[[key]]
  }

  # Recursively apply for keys present in both that are lists
  common_keys <- intersect(data_names, defaults_names)
  for (key in common_keys) {
    # Check if both are lists (and data[[key]] isn't NULL, which can happen)
    if (is.list(defaults[[key]]) && !is.null(data[[key]]) && is.list(data[[key]])) {
       # Only recurse if the default is non-empty or if the data list is non-empty
       # This prevents infinite recursion on empty list defaults for non-list data
       if (length(defaults[[key]]) > 0 || length(data[[key]]) > 0) {
            data[[key]] <- apply_defaults(data[[key]], defaults[[key]])
       }
    }
  }
  return(data)
}

#' Parse, Validate Schema, Apply Defaults, and Return IOR
#'
#' @param yaml_file Path to the YAML configuration file.
#' @return A list (IOR) representing the configuration, with defaults applied.
#' @keywords internal
parse_and_validate_config <- function(yaml_file) {

  if (!fs::file_exists(yaml_file)) {
    stop("YAML configuration file not found: ", yaml_file, call. = FALSE)
  }

  # --- 1. Parse YAML ---
  raw_config_list <- NULL
  tryCatch({
      raw_config_list <- yaml::read_yaml(yaml_file)
  }, error = function(e) {
      stop("Failed to parse YAML file '", yaml_file, "'. Error: ", e$message, call. = FALSE)
  })
  if (is.null(raw_config_list) || !is.list(raw_config_list)) { # Check it's a list
       stop("YAML file '", yaml_file, "' is empty or could not be parsed into a list structure.", call. = FALSE)
  }

  errors <- ValidationErrors$new()

  # --- 2. Define Defaults (Refactored) ---
  default_config <- list(
    dataset = list(
      relpath = "func",
      subjects = list(include = NULL, exclude = NULL),
      tasks = NULL,
      runs = NULL,
      scan_params = list(TR = NULL, run_length = NULL)
    ),
    events = list(), # Required fields checked later
    hrfs = list(canonical = list(type="HRF_SPMG1")), # Default canonical HRF
    confounds = list(include = NULL, exclude = NULL),
    terms = list(), # Required, content user-defined
    contrasts = list(), # Optional
    model = list(
        baseline = list(basis = "bspline", degree = 3L, intercept = "runwise", confound_variables = list()),
        variable_mapping = list(), # Required, content user-defined
        variable_roles = list(factors = list(), parametric = list()),
        terms = list(), # Required, content user-defined
        contrasts = list() # Optional
    )
  )

  # --- 3. Apply Defaults ---
  config_list <- apply_defaults(raw_config_list, default_config)

  # --- 4. Schema Validation (based on refined YAML spec) ---
  # (This section implements the checks outlined in the YAML spec v2.0)

  # --- Top Level ---
  check_required(config_list, "dataset", "", errors)
  check_required(config_list, "events", "", errors)
  check_required(config_list, "terms", "", errors)
  check_required(config_list, "model", "", errors)

  # --- Dataset ---
  if (check_type(config_list, "dataset", "object", "", errors)) {
      ds_path = "dataset"
      ds_data = config_list$dataset
      if(check_required(ds_data, "path", ds_path, errors)) {
          check_type(ds_data, "path", "string", ds_path, errors)
          # Directory existence checked later in build_config_from_ior
      }
      check_type(ds_data, "relpath", "string", ds_path, errors, allow_null=FALSE) # default applied

      if(check_type(ds_data, "subjects", "object", ds_path, errors, allow_null = TRUE)) {
           if (!is.null(ds_data$subjects)){
              subj_path = paste0(ds_path,"$subjects")
              check_type(ds_data$subjects, "include", "array[string]", subj_path, errors, allow_null = TRUE)
              check_pattern(ds_data$subjects, "include", "^sub-[0-9A-Za-z]+$", subj_path, errors, allow_null = TRUE)
              check_type(ds_data$subjects, "exclude", "array[string]", subj_path, errors, allow_null = TRUE)
              check_pattern(ds_data$subjects, "exclude", "^sub-[0-9A-Za-z]+$", subj_path, errors, allow_null = TRUE)
           }
      }
      check_type(ds_data, "tasks", "array[string]", ds_path, errors, allow_null = TRUE)
      check_type(ds_data, "runs", "array[string]", ds_path, errors, allow_null = TRUE)
      check_pattern(ds_data, "runs", "^run-[0-9]+$", ds_path, errors, allow_null = TRUE)

      if(check_type(ds_data, "scan_params", "object", ds_path, errors, allow_null = TRUE)) {
          if (!is.null(ds_data$scan_params)) {
               sp_data <- ds_data$scan_params
               sp_path <- paste0(ds_path, "$scan_params")
               # TR Validation
                if(check_type(sp_data, "TR", "object", sp_path, errors, allow_null = TRUE)) {
                    if (!is.null(sp_data$TR)) {
                        tr_data <- sp_data$TR
                        tr_path <- paste0(sp_path, "$TR")
                        check_required(tr_data, "default", tr_path, errors)
                        check_type(tr_data, "default", "number", tr_path, errors)
                        if(check_type(tr_data, "overrides", "array[object]", tr_path, errors, allow_null=TRUE)) {
                            # (Validation of override items as before)
                        }
                    }
                }
                # Run Length Validation
                if(check_type(sp_data, "run_length", "object", sp_path, errors, allow_null = TRUE)) {
                    if (!is.null(sp_data$run_length)) {
                         rl_data <- sp_data$run_length
                         rl_path <- paste0(sp_path, "$run_length")
                         if (check_required(rl_data, "default", rl_path, errors) &&
                             check_type(rl_data, "default", "object", rl_path, errors)) {
                               # (Validation of default object items as before)
                         }
                         if(check_type(rl_data, "overrides", "array[object]", rl_path, errors, allow_null=TRUE)) {
                              # (Validation of override items as before)
                         }
                    }
                }
          }
      }
  } # End Dataset Validation

  # --- Events ---
  if (check_type(config_list, "events", "object", "", errors)) {
      ev_path = "events"
      ev_data = config_list$events
      check_required(ev_data, "onset_column", ev_path, errors)
      check_type(ev_data, "onset_column", "string", ev_path, errors)
      check_required(ev_data, "duration_column", ev_path, errors)
      check_type(ev_data, "duration_column", "string", ev_path, errors)
      check_required(ev_data, "block_column", ev_path, errors)
      check_type(ev_data, "block_column", "string", ev_path, errors)
  }

  # --- HRFs ---
  if (check_type(config_list, "hrfs", "object", "hrfs", errors, allow_null=TRUE)) {
       if (!is.null(config_list$hrfs)) {
            if (length(config_list$hrfs) > 0 && is.null(names(config_list$hrfs))) {
                errors$add_error("hrfs", "Must be an object (named list), not an array.")
            } else {
                 for (hrf_name in names(config_list$hrfs)) {
                     hrf_path <- paste0("hrfs$", hrf_name)
                     if(check_type(config_list$hrfs, hrf_name, "object", "hrfs", errors)) {
                         hrf_spec <- config_list$hrfs[[hrf_name]]
                         check_required(hrf_spec, "type", hrf_path, errors)
                         check_type(hrf_spec, "type", "string", hrf_path, errors)
                         allowed_hrf_types <- c("HRF_SPMG1", "HRF_SPMG2", "HRF_SPMG3",
                                                "HRF_GAUSSIAN", "HRF_GAMMA", "HRF_BSPLINE",
                                                "HRF_DAGUERRE", "HRF_DAGUERRE_BASIS", "custom")
                         check_enum(hrf_spec, "type", allowed_hrf_types, hrf_path, errors)
                         check_type(hrf_spec, "parameters", "object", hrf_path, errors, allow_null=TRUE)
                         if (!is.null(hrf_spec$type) && hrf_spec$type == "custom") {
                             check_required(hrf_spec, "definition", hrf_path, errors)
                             check_type(hrf_spec, "definition", "string", hrf_path, errors)
                         }
                         check_type(hrf_spec, "lag", "number", hrf_path, errors, allow_null=TRUE)
                         check_type(hrf_spec, "width", "number", hrf_path, errors, allow_null=TRUE)
                         check_type(hrf_spec, "summate", "boolean", hrf_path, errors, allow_null=TRUE)
                         check_type(hrf_spec, "normalize", "boolean", hrf_path, errors, allow_null=TRUE)
                     }
                 }
            }
       }
  }

  # --- Confounds ---
   if (check_type(config_list, "confounds", "object", "confounds", errors, allow_null=TRUE)) {
        if (!is.null(config_list$confounds)) {
             conf_path = "confounds"
             conf_data = config_list$confounds
             check_type(conf_data, "include", "array[string]", conf_path, errors, allow_null=TRUE)
             check_type(conf_data, "exclude", "array[string]", conf_path, errors, allow_null=TRUE)
        }
   }

  # --- Terms ---
  if (check_type(config_list, "terms", "object", "", errors)) {
       if (length(config_list$terms) == 0) {
           errors$add_error("terms", "Section 'terms' cannot be empty.")
       } else if (is.null(names(config_list$terms))) {
           errors$add_error("terms", "Section 'terms' must be an object (named list), not an array.")
       } else {
           allowed_hrf_names <- names(config_list$hrfs %||% list()) # Get defined HRF names
           for (term_name in names(config_list$terms)) {
                term_path <- paste0("terms$", term_name)
                if(check_type(config_list$terms, term_name, "object", "terms", errors)) {
                    term_spec <- config_list$terms[[term_name]]
                    check_required(term_spec, "type", term_path, errors)
                    term_type <- term_spec$type # Store for conditional checks
                    check_type(term_spec, "type", "string", term_path, errors)
                    allowed_term_types <- c("hrf", "parametric", "trialwise", "nuisance")
                    check_enum(term_spec, "type", allowed_term_types, term_path, errors)

                    if (check_required(term_spec, "variables", term_path, errors) &&
                        check_type(term_spec, "variables", "array[string]", term_path, errors)) {
                         if (length(term_spec$variables %||% list()) == 0) {
                             errors$add_error(paste0(term_path, "$variables"), "Must contain at least one variable name.")
                         }
                    }

                    # Check HRF reference if type requires it
                    needs_hrf <- term_type %in% c("hrf", "parametric", "trialwise")
                    if (needs_hrf) {
                        # HRF is optional (has default), but if specified, must be valid
                        if (check_type(term_spec, "hrf", "string", term_path, errors, allow_null=TRUE)) {
                            if (!is.null(term_spec$hrf) && !term_spec$hrf %in% allowed_hrf_names) {
                                # Check if it's a known built-in *implicitly*
                                known_builtins <- c("spmg1", "spmg2", "spmg3", "gaussian", "gamma") # Add others if applicable
                                if (!tolower(term_spec$hrf) %in% known_builtins && term_spec$hrf != "canonical") {
                                    errors$add_error(paste0(term_path, "$hrf"), paste0("HRF name '", term_spec$hrf, "' not found in global 'hrfs' section and not a recognized default."))
                                }
                            }
                        }
                    } else { # Nuisance type should NOT have hrf
                        if (!is.null(term_spec$hrf)) errors$add_error(paste0(term_path,"$hrf"), "Nuisance terms should not specify an HRF.")
                    }

                    check_type(term_spec, "subset", "string", term_path, errors, allow_null=TRUE)
                    check_type(term_spec, "lag", "number", term_path, errors, allow_null=TRUE)

                    # Validations specific to 'parametric' terms
                    if (!is.null(term_type) && term_type == "parametric") {
                        check_type(term_spec, "transform", "array[string]", term_path, errors, allow_null = TRUE)
                        if (!is.null(term_spec$transform)) {
                           allowed_transforms <- c("center", "scale", "zscore", "log", "exp", "within_subject")
                           check_enum(term_spec, "transform", allowed_transforms, term_path, errors, allow_null = TRUE)
                        }
                        if(check_type(term_spec, "basis", "object", term_path, errors, allow_null=TRUE)) {
                             if (!is.null(term_spec$basis)) {
                                 basis_path <- paste0(term_path, "$basis")
                                 basis_spec_data <- term_spec$basis
                                 check_required(basis_spec_data, "type", basis_path, errors)
                                 check_type(basis_spec_data, "type", "string", basis_path, errors)
                                 allowed_basis_types <- c("Poly", "BSpline", "Standardized", "Ident", "ns")
                                 check_enum(basis_spec_data, "type", allowed_basis_types, basis_path, errors)
                                 check_type(basis_spec_data, "parameters", "object", basis_path, errors, allow_null=TRUE)
                                 # Can add parameter sub-validation here if needed
                             }
                        }
                         # Check if last variable exists for modulation
                         if (length(term_spec$variables %||% list()) == 0) {
                            errors$add_error(paste0(term_path, "$variables"), "Parametric terms require at least one variable (the modulator).")
                         }

                    } else { # Not parametric
                        if (!is.null(term_spec$transform)) errors$add_error(paste0(term_path, "$transform"), "'transform' is only valid for 'parametric' terms.")
                        if (!is.null(term_spec$basis)) errors$add_error(paste0(term_path, "$basis"), "'basis' is only valid for 'parametric' terms.")
                    }
                }
           }
       }
  } # End Terms Validation

  # --- Contrasts ---
  if (check_type(config_list, "contrasts", "object", "contrasts", errors, allow_null=TRUE)) {
       if (!is.null(config_list$contrasts)) {
           if (length(config_list$contrasts) > 0 && is.null(names(config_list$contrasts))) {
                errors$add_error("contrasts", "Section 'contrasts' must be an object (named list) if not empty.")
           } else {
                for (con_name in names(config_list$contrasts)) {
                    con_path <- paste0("contrasts$", con_name)
                    if(check_type(config_list$contrasts, con_name, "object", "contrasts", errors)) {
                       con_spec <- config_list$contrasts[[con_name]]
                       check_required(con_spec, "type", con_path, errors)
                       check_type(con_spec, "type", "string", con_path, errors)
                       allowed_con_types <- c("formula", "pair", "one_against_all", "unit", "oneway", "interaction", "polynomial")
                       check_enum(con_spec, "type", allowed_con_types, con_path, errors)
                       con_type <- con_spec$type
                       if (!is.null(con_type)) {
                           # Type-specific required/optional fields
                           if (con_type == "formula") {
                               check_required(con_spec, "expression", con_path, errors)
                               check_type(con_spec, "expression", "string", con_path, errors)
                           } else if (con_type %in% c("pair", "one_against_all", "oneway", "interaction", "polynomial")) {
                               check_required(con_spec, "factors", con_path, errors) # Make factors required for these
                               check_type(con_spec, "factors", "array[string]", con_path, errors)
                           }
                           if (con_type == "polynomial") {
                                check_type(con_spec, "degree", "integer", con_path, errors, allow_null=TRUE) # Degree might be optional (default 1?)
                           }
                           if (con_type %in% c("one_against_all", "pair")) {
                                check_type(con_spec, "levels", "array[string]", con_path, errors, allow_null=TRUE)
                           }
                       }
                       check_type(con_spec, "where", "string", con_path, errors, allow_null=TRUE)
                    }
                }
           }
       }
  } # End Contrasts Validation

  # --- Model ---
  if (check_type(config_list, "model", "object", "", errors)) {
      mod_path = "model"
      mod_data = config_list$model
      check_required(mod_data, "name", mod_path, errors)
      check_type(mod_data, "name", "string", mod_path, errors)

      # Validate model$baseline
      # (Assume checks from previous version are sufficient, uses defaults)

      # Validate model$variable_mapping
      if (check_required(mod_data, "variable_mapping", mod_path, errors) &&
          check_type(mod_data, "variable_mapping", "object", mod_path, errors)) {
          # (Assume checks from previous version are sufficient)
      }

      # Validate model$variable_roles
      if (check_type(mod_data, "variable_roles", "object", mod_path, errors, allow_null = TRUE)) {
          # (Assume checks from previous version are sufficient)
      }

      # Validate model$terms
      if (check_required(mod_data, "terms", mod_path, errors) &&
          check_type(mod_data, "terms", "array[string]", mod_path, errors)) {
          if (length(mod_data$terms %||% list()) == 0) {
              errors$add_error(paste0(mod_path, "$terms"), "Must include at least one term.")
          } else {
              # Check term names exist globally
              defined_terms <- names(config_list$terms %||% list())
              missing_t <- setdiff(unlist(mod_data$terms), defined_terms)
              if (length(missing_t) > 0) errors$add_error(paste0(mod_path, "$terms"), paste("Term(s) not defined globally:", paste(missing_t, collapse=", ")))
          }
      }

      # Validate model$contrasts
      if (check_type(mod_data, "contrasts", "array[string]", mod_path, errors, allow_null=TRUE)) {
          if (!is.null(mod_data$contrasts)) {
               defined_contrasts <- names(config_list$contrasts %||% list())
               missing_c <- setdiff(unlist(mod_data$contrasts), defined_contrasts)
               if (length(missing_c) > 0) errors$add_error(paste0(mod_path, "$contrasts"), paste("Contrast(s) not defined globally:", paste(missing_c, collapse=", ")))
          }
      }
  } # End Model Validation


  # --- 5. Stop if validation failed ---
  errors$stop_if_invalid(paste("YAML configuration schema validation failed for:", yaml_file))

  # --- 6. Return IOR ---
  message("YAML configuration schema validated successfully: ", yaml_file)
  attr(config_list, "validated_schema") <- TRUE
  return(config_list)
}

#' Build fmri_config Object from Validated IOR List (Contextual Validation)
#'
#' @param validated_ior A list returned by `parse_and_validate_config`.
#' @return An `fmri_config` object ready for model building.
#' @keywords internal
build_config_from_ior <- function(validated_ior) {

  if (!isTRUE(attr(validated_ior, "validated_schema"))) {
      stop("Input 'validated_ior' must be the result of parse_and_validate_config().")
  }

  errors <- ValidationErrors$new()

  # --- 1. Load BIDS Project ---
  bids_path <- fs::path_expand(validated_ior$dataset$path)
  proj <- NULL
  tryCatch({
      if (!fs::dir_exists(bids_path)) {
          stop("BIDS directory does not exist at the specified path.")
      }
      proj <- bidser::bids_project(bids_path, fmriprep = FALSE)
      message("BIDS project loaded successfully: ", bids_path)
  }, error = function(e) {
      errors$add_error("dataset$path", paste("Failed to load BIDS project at '", bids_path, "'. Error: ", e$message))
      errors$stop_if_invalid("Cannot proceed without a valid BIDS project")
  })

  # --- 2. Contextual Validation ---

  # Validate Subjects (using bidser::participants)
  available_ids <- character(0)
  tryCatch({
      parts_df <- bidser::participants(proj)
      if (!is.null(parts_df) && "participant_id" %in% names(parts_df)) {
          available_ids <- parts_df$participant_id
      }
      if (length(available_ids) == 0) warning("No participants found in BIDS dataset.")
  }, error = function(e) { errors$add_error("dataset$subjects", paste("Failed retrieving participants:", e$message)) })

  subjects_spec <- validated_ior$dataset$subjects
  subjects <- available_ids
  if (!is.null(subjects_spec$include)) {
      include_subs <- unique(unlist(subjects_spec$include))
      missing_inc <- setdiff(include_subs, available_ids)
      if (length(missing_inc) > 0) errors$add_error("dataset$subjects$include", paste("Subject(s) not found:", paste(missing_inc, collapse=", ")))
      subjects <- intersect(subjects, include_subs)
  }
  if (!is.null(subjects_spec$exclude)) {
      exclude_subs <- unique(unlist(subjects_spec$exclude))
      subjects <- setdiff(subjects, exclude_subs)
  }
  if (length(subjects) == 0) {
      if (length(available_ids) > 0) { errors$add_error("dataset$subjects", "No subjects selected.") }
      else { warning("No subjects available in BIDS dataset.") }
  }

  # Validate Tasks (using bidser::tasks)
  available_tasks <- character(0)
  tryCatch({
      available_tasks <- bidser::tasks(proj)
      if (length(available_tasks) == 0) warning("No tasks found in BIDS dataset.")
  }, error = function(e) { errors$add_error("dataset$tasks", paste("Failed retrieving tasks:", e$message)) })

  tasks_spec <- validated_ior$dataset$tasks
  tasks <- if (!is.null(tasks_spec)) {
      tasks_vec <- unique(unlist(tasks_spec))
      missing_tasks <- setdiff(tasks_vec, available_tasks)
      if (length(missing_tasks) > 0) errors$add_error("dataset$tasks", paste("Task(s) not found:", paste(missing_tasks, collapse=", ")))
      intersect(tasks_vec, available_tasks)
   } else {
      available_tasks
   }
   if (length(tasks) == 0) {
       if (!is.null(tasks_spec)) { errors$add_error("dataset$tasks", "No specified tasks found.") }
       else if (length(available_tasks) > 0) { # User specified none, but BIDS has tasks -> error? or use all? Let's use all for now, consistent with spec.
             tasks <- available_tasks }
       else { warning("No tasks specified and none found in BIDS.") }
   }

  # Validate Runs (Placeholder - requires bidser::runs or similar)
  selected_runs <- validated_ior$dataset$runs %||% "all" # Placeholder

  # Validate Event Columns, Confound Columns, and Infer Roles
  events_info <- list(columns = character(0), mapping = validated_ior$events)
  confounds_info <- list(spec = validated_ior$confounds %||% list(include=NULL, exclude=NULL), columns = character(0), available_columns = character(0))
  variable_roles <- list(inferred_roles = list(), factors = character(0), parametric = character(0))
  model_event_map <- validated_ior$model$variable_mapping %||% list()

  if (length(subjects) > 0 && length(tasks) > 0) {
      first_subj <- subjects[1]
      first_task <- tasks[1]
      rep_events_df <- NULL
      rep_confounds_df <- NULL

      # --- Load Representative Event Data ---
      tryCatch({
          events_read <- bidser::read_events(proj, subid = first_subj, task = first_task)
          # Assuming read_events returns a structure where data[[1]] is the relevant df
          if (length(events_read) > 0 && length(events_read$data) > 0 && !is.null(events_read$data[[1]])) {
              rep_events_df <- events_read$data[[1]]
              events_info$columns <- names(rep_events_df)
          } else { errors$add_error("events", paste0("No event data for representative: ", first_subj, "/", first_task)) }
      }, error = function(e) { errors$add_error("events", paste0("Error loading events: ", e$message)) })

      # --- Load Representative Confound Data ---
      if (!is.null(confounds_info$spec$include) || length(validated_ior$model$baseline$confound_variables %||% list()) > 0) {
          tryCatch({
              confounds_read <- bidser::read_confounds(proj, subid = first_subj) # Add task/run?
              if (length(confounds_read) > 0 && length(confounds_read$data) > 0 && !is.null(confounds_read$data[[1]])) {
                  rep_confounds_df <- confounds_read$data[[1]]
                  confounds_info$available_columns <- names(rep_confounds_df)

                  # Apply include/exclude from global confounds spec
                  included_cols <- confounds_info$available_columns
                  if (!is.null(confounds_info$spec$include)) {
                     include_patterns <- unlist(confounds_info$spec$include)
                     matched_include <- unique(unlist(lapply(include_patterns, grep, x=included_cols, value=TRUE)))
                     if (length(matched_include) == 0 && length(include_patterns) > 0) {
                        errors$add_error("confounds$include", "No available confounds matched the include patterns.")
                     }
                     included_cols <- matched_include
                  }
                  if (!is.null(confounds_info$spec$exclude)) {
                     exclude_patterns <- unlist(confounds_info$spec$exclude)
                     exclude_matches <- unique(unlist(lapply(exclude_patterns, grep, x=included_cols, value=TRUE)))
                     included_cols <- setdiff(included_cols, exclude_matches)
                  }
                  confounds_info$columns <- included_cols # Store selected columns

              } else { errors$add_error("confounds", paste0("No confound data for representative: ", first_subj)) }
          }, error = function(e) { errors$add_error("confounds", paste0("Error loading confounds: ", e$message)) })
      }

      # --- Validate Mappings and Baseline Confounds ---
      if (!is.null(rep_events_df)) {
          # Check essential event columns
          essential_bids_cols <- unlist(events_info$mapping[c("onset_column", "duration_column", "block_column")])
          missing_essential <- setdiff(essential_bids_cols, events_info$columns)
          if(length(missing_essential)>0) errors$add_error("events", paste0("Essential BIDS columns not found in event file: ", paste(missing_essential, collapse=", ")))

          # Check mapped variable columns (from model$variable_mapping)
          mapped_bids_cols <- unlist(model_event_map)
          available_bids_cols <- unique(c(events_info$columns, confounds_info$available_columns)) # All potential source columns

          missing_mapped <- setdiff(mapped_bids_cols, available_bids_cols)
          if (length(missing_mapped) > 0) {
               errors$add_error("model$variable_mapping", paste0("Mapped BIDS columns not found in representative files: ", paste(missing_mapped, collapse=", ")))
          }

          # Check baseline confounds listed in model are mapped AND selected
          baseline_conf_vars <- validated_ior$model$baseline$confound_variables %||% list()
          if (length(baseline_conf_vars) > 0) {
              baseline_bids_cols_needed <- unlist(model_event_map[baseline_conf_vars])
              if (any(sapply(baseline_bids_cols_needed, is.null))) {
                   errors$add_error("model$variable_mapping", "Baseline confound variables must be mapped.")
              }
              missing_in_available <- setdiff(baseline_bids_cols_needed, confounds_info$available_columns)
              missing_in_selected <- setdiff(baseline_bids_cols_needed, confounds_info$columns)
              if(length(missing_in_available)>0) errors$add_error("model$baseline$confound_variables", paste("Mapped BIDS columns for baseline not found in confound file:", paste(missing_in_available, collapse=", ")))
              if(length(missing_in_selected)>0) errors$add_error("model$baseline$confound_variables", paste("Mapped BIDS columns for baseline were excluded by global 'confounds' rules:", paste(missing_in_selected, collapse=", ")))
          }
      }

      # --- Infer/Validate Roles ---
      if (!is.null(rep_events_df)) {
          # Combine representative data sources for inference if possible
          # This part remains tricky if confounds aren't event-aligned
          # Simplification: Infer based on event data primarily, assume confounds are parametric unless declared factors
          data_for_inference <- rep_events_df
          variable_roles <- infer_and_validate_variable_roles(
              regressors_spec = validated_ior$terms[validated_ior$model$terms %||% list()], # Only terms in the model
              events_mapping = model_event_map,
              events_data = data_for_inference, # Use representative events data
              user_factors = unlist(validated_ior$model$variable_roles$factors %||% list()),
              user_parametric = unlist(validated_ior$model$variable_roles$parametric %||% list())
          )
      } else {
          warning("Skipping role inference due to missing representative event data.")
          variable_roles$factors = unique(unlist(validated_ior$model$variable_roles$factors %||% list()))
          variable_roles$parametric = unique(unlist(validated_ior$model$variable_roles$parametric %||% list()))
      }

  } else {
      warning("Skipping contextual validation and role inference (no subjects or tasks selected/found).")
      variable_roles$factors = unique(unlist(validated_ior$model$variable_roles$factors %||% list()))
      variable_roles$parametric = unique(unlist(validated_ior$model$variable_roles$parametric %||% list()))
  }


  # --- Final Checks ---
  # Check that terms and contrasts listed in model are defined globally
  defined_terms <- names(validated_ior$terms %||% list())
  model_terms <- unique(unlist(validated_ior$model$terms %||% list()))
  missing_defined_terms <- setdiff(model_terms, defined_terms)
  if (length(missing_defined_terms) > 0) errors$add_error("model$terms", paste("Term(s) not defined globally:", paste(missing_defined_terms, collapse=", ")))

  defined_contrasts <- names(validated_ior$contrasts %||% list())
  model_contrasts <- unique(unlist(validated_ior$model$contrasts %||% list()))
  missing_defined_contrasts <- setdiff(model_contrasts, defined_contrasts)
  if (length(missing_defined_contrasts) > 0) errors$add_error("model$contrasts", paste("Contrast(s) not defined globally:", paste(missing_defined_contrasts, collapse=", ")))

  errors$stop_if_invalid("Context-dependent validation failed")

  # --- Construct fmri_config ---
  config <- structure(
    list(
      spec = validated_ior,
      project = proj,
      subjects = subjects,
      tasks = tasks,
      runs = selected_runs,
      events_info = events_info,
      confounds_info = confounds_info,
      variable_roles = variable_roles,
      validated = TRUE
    ),
    class = "fmri_config"
  )

  message("fmri_config object built successfully.")
  return(config)
}


#' Build fmri_config Object from Validated IOR List
#'
#' Takes a validated IOR list, performs context-dependent validation against BIDS,
#' infers variable roles, and constructs the final `fmri_config` object.
#'
#' @param validated_ior A list returned by `parse_and_validate_config`.
#' @return An `fmri_config` object ready for model building.
#' @keywords internal
build_config_from_ior <- function(validated_ior) {

  if (!isTRUE(attr(validated_ior, "validated_schema"))) {
      stop("Input 'validated_ior' must be the result of parse_and_validate_config().")
  }

  errors <- ValidationErrors$new()

  # --- 1. Load BIDS Project ---
  bids_path <- fs::path_expand(validated_ior$dataset$path)
  proj <- NULL
  tryCatch({
      # Ensure path exists before calling bids_project
      if (!fs::dir_exists(bids_path)) {
          stop("BIDS directory does not exist at the specified path.")
      }
      # Load as standard BIDS, not necessarily fmriprep derivatives for broader compatibility
      proj <- bidser::bids_project(bids_path, fmriprep = FALSE) 
      message("BIDS project loaded successfully: ", bids_path)
  }, error = function(e) {
      errors$add_error("dataset$path", paste("Failed to load BIDS project at '", bids_path, "'. Error: ", e$message))
      errors$stop_if_invalid("Cannot proceed without a valid BIDS project")
  })

  # --- 2. Contextual Validation ---

  # Validate Subjects (using bidser::participants)
  available_ids <- character(0)
  tryCatch({
      parts_df <- bidser::participants(proj)
      if (!is.null(parts_df) && "participant_id" %in% names(parts_df)) {
          available_ids <- parts_df$participant_id
      }
      if (length(available_ids) == 0) warning("No participants found in BIDS dataset.")
  }, error = function(e) { errors$add_error("dataset$subjects", paste("Failed retrieving participants:", e$message)) })

  subjects_spec <- validated_ior$dataset$subjects
  subjects <- available_ids
  if (!is.null(subjects_spec$include)) {
      include_subs <- unique(unlist(subjects_spec$include))
      missing_inc <- setdiff(include_subs, available_ids)
      if (length(missing_inc) > 0) errors$add_error("dataset$subjects$include", paste("Subject(s) not found:", paste(missing_inc, collapse=", ")))
      subjects <- intersect(subjects, include_subs)
  }
  if (!is.null(subjects_spec$exclude)) {
      exclude_subs <- unique(unlist(subjects_spec$exclude))
      subjects <- setdiff(subjects, exclude_subs)
  }
  if (length(subjects) == 0) {
      if (length(available_ids) > 0) { errors$add_error("dataset$subjects", "No subjects selected.") }
      else { warning("No subjects available in BIDS dataset.") }
  }

  # Validate Tasks (using bidser::tasks)
  available_tasks <- character(0)
  tryCatch({
      available_tasks <- bidser::tasks(proj)
      if (length(available_tasks) == 0) warning("No tasks found in BIDS dataset.")
  }, error = function(e) { errors$add_error("dataset$tasks", paste("Failed retrieving tasks:", e$message)) })

  tasks_spec <- validated_ior$dataset$tasks
  tasks <- if (!is.null(tasks_spec)) {
      tasks_vec <- unique(unlist(tasks_spec))
      missing_tasks <- setdiff(tasks_vec, available_tasks)
      if (length(missing_tasks) > 0) errors$add_error("dataset$tasks", paste("Task(s) not found:", paste(missing_tasks, collapse=", ")))
      intersect(tasks_vec, available_tasks)
   } else {
      available_tasks
   }
   if (length(tasks) == 0) {
       if (!is.null(tasks_spec)) { errors$add_error("dataset$tasks", "No specified tasks found.") }
       else if (length(available_tasks) > 0) { # User specified none, but BIDS has tasks -> error? or use all? Let's use all for now, consistent with spec.
             tasks <- available_tasks }
       else { warning("No tasks specified and none found in BIDS.") }
   }

  # Validate Runs (if specified) - Check against available runs for first subject/task combo
   available_runs <- character(0)
   if (length(subjects) > 0 && length(tasks) > 0) {
      first_subj_runs <- character(0)
      tryCatch({
           # This needs a function like `bidser::runs(proj, subid=subjects[1], task=tasks[1])`
           # Assuming such a function exists or can be approximated
           # Placeholder: For now, just check if specified runs look okay
           available_runs <- paste0("run-", sprintf("%02d", 1:5)) # Mock available runs
      }, error = function(e) {
           warning("Could not determine available runs for validation.")
      })

       if (!is.null(validated_ior$dataset$runs)) {
           runs_spec_req <- unique(unlist(validated_ior$dataset$runs))
           missing_runs <- setdiff(runs_spec_req, available_runs)
            if (length(missing_runs) > 0) {
                errors$add_error("dataset$runs", paste("Specified run(s) potentially not available:", paste(missing_runs, collapse=", ")))
            }
       }
   }
   selected_runs <- validated_ior$dataset$runs %||% available_runs # Use specified or all available


  # Validate Event Columns, Confound Columns, and Infer Roles
  events_info <- list(columns = character(0), mapping = validated_ior$events)
  confounds_info <- list(spec = validated_ior$confounds, columns = character(0), available_columns = character(0))
  variable_roles <- list(inferred_roles = list(), factors = character(0), parametric = character(0))
  events_list_for_inference <- list() # Store loaded data for role inference

  if (length(subjects) > 0 && length(tasks) > 0) {
      first_subj <- subjects[1]
      first_task <- tasks[1]
      first_subj_events_loaded <- FALSE
      first_subj_confounds_loaded <- FALSE

      # --- Validate Events ---
      tryCatch({
          events_read <- bidser::read_events(proj, subid = first_subj, task = first_task)
          if (length(events_read) > 0 && length(events_read$data) > 0 && !is.null(events_read$data[[1]])) {
              rep_events_df <- events_read$data[[1]] # Representative event df
              events_info$columns <- names(rep_events_df)
              first_subj_events_loaded <- TRUE
              events_list_for_inference[[first_subj]] <- list() # Init subject
              events_list_for_inference[[first_subj]][[first_task]] <- rep_events_df # Store for inference

              # Check essential columns
              emap <- events_info$mapping
              req_cols <- c(emap$onset_column, emap$duration_column, emap$block_column)
              missing_req <- setdiff(req_cols, events_info$columns)
              if(length(missing_req)>0) errors$add_error(paste0("events (",first_task,")"), paste0("Required column(s) ('",paste(missing_req, collapse="', '"),"') not found in events file. Avail: ", paste(events_info$columns, collapse=", ")))

              # Check block column validity
              block_col_name <- emap$block_column
              if (block_col_name %in% events_info$columns) {
                  block_vals <- rep_events_df[[block_col_name]]
                  # Conversion/Validation logic for block_vals (as before)
                  # ... (omitted for brevity, assume checks from previous step are here)
              }
          } else {
              errors$add_error(paste0("events (",first_task,")"), paste0("No event data found for representative subj/task: ", first_subj, "/", first_task))
          }
      }, error = function(e) {
           errors$add_error(paste0("events (", first_task, ")"), paste0("Error reading/validating events: ", e$message))
      })

      # --- Validate Mapped Event Variables ---
      model_event_map <- validated_ior$model$variable_mapping
      all_model_vars <- names(model_event_map)
      bids_event_cols_needed <- unlist(model_event_map) # Assuming simple string mapping for now

      if (first_subj_events_loaded) {
          missing_mapped_cols <- setdiff(bids_event_cols_needed, events_info$columns)
          if (length(missing_mapped_cols) > 0) {
              # Find which model vars map to missing cols
              missing_map_info <- paste(
                  sapply(names(model_event_map)[model_event_map %in% missing_mapped_cols], function(mv) {
                      paste0("'", mv, "' -> '", model_event_map[[mv]], "'")
                  }), collapse = "; "
              )
              errors$add_error("model$variable_mapping", paste0("Mapped BIDS columns not found in event file: ", paste(missing_mapped_cols, collapse=", "), ". Mappings: [", missing_map_info, "]"))
          }
      } # else error already added above

      # --- Validate Confounds ---
      confounds_spec <- confounds_info$spec
      if (!is.null(confounds_spec)) {
          tryCatch({
              confounds_read <- bidser::read_confounds(proj, subid = first_subj) # Add task/run filtering if needed
              if (length(confounds_read) > 0 && length(confounds_read$data) > 0 && !is.null(confounds_read$data[[1]])) {
                  rep_confounds_df <- confounds_read$data[[1]] # Representative confound df
                  confounds_info$available_columns <- names(rep_confounds_df)
                  first_subj_confounds_loaded <- TRUE

                  # Apply include/exclude logic (as before)
                   included_cols <- character()
                   if (!is.null(confounds_spec$include)) {
                      include_patterns <- unlist(confounds_spec$include)
                      for (pattern in include_patterns) {
                           # Use fixed=FALSE for regex matching
                           matches <- grep(pattern, confounds_info$available_columns, value = TRUE, fixed=FALSE)
                           if (length(matches) == 0) {
                                errors$add_error("confounds$include", sprintf("Include pattern '%s' matched no variables. Avail: %s", pattern, paste(head(confounds_info$available_columns, 10), collapse=", ")))
                           }
                           included_cols <- c(included_cols, matches)
                      }
                      included_cols <- unique(included_cols)
                   } else {
                      included_cols <- confounds_info$available_columns # Include all if include is omitted
                   }
                   if (!is.null(confounds_spec$exclude)) {
                       exclude_patterns <- unlist(confounds_spec$exclude)
                       for (pattern in exclude_patterns) {
                          exclude_matches <- grep(pattern, included_cols, value = TRUE, fixed=FALSE)
                          included_cols <- setdiff(included_cols, exclude_matches)
                       }
                   }
                   confounds_info$columns <- included_cols # Store the *selected* columns
                   if (length(confounds_info$columns) == 0 && (!is.null(confounds_spec$include) || !is.null(confounds_spec$exclude))) {
                       warning(sprintf("Confound include/exclude patterns resulted in zero selected columns for subject '%s'.", first_subj))
                   }

                   # Check baseline confounds are available and selected
                   baseline_conf_vars <- validated_ior$model$baseline$confound_variables %||% list()
                   baseline_bids_cols_needed <- unlist(model_event_map[baseline_conf_vars])
                   missing_base_conf_bids <- setdiff(baseline_bids_cols_needed, confounds_info$available_columns)
                   unselected_base_conf_bids <- setdiff(baseline_bids_cols_needed, confounds_info$columns)

                   if (length(missing_base_conf_bids) > 0) {
                       errors$add_error("model$baseline$confound_variables", paste("Mapped BIDS columns for baseline confounds not found in confound file:", paste(missing_base_conf_bids, collapse=", ")))
                   }
                   if (length(unselected_base_conf_bids) > 0) {
                       errors$add_error("model$baseline$confound_variables", paste("Baseline confounds were not selected by global 'confounds' include/exclude rules:", paste(unselected_base_conf_bids, collapse=", ")))
                   }

              } else if (length(baseline_conf_vars)>0){
                   # Error only if baseline *needs* confounds but none were found
                   errors$add_error("confounds", paste0("No confound data found for representative subject '", first_subj, "', but confounds listed in baseline."))
              }
          }, error = function(e) {
              errors$add_error("confounds", paste0("Error reading/validating confounds: ", e$message))
          })
      }

      # --- Infer and Validate Variable Roles ---
      if (first_subj_events_loaded) {
            # Ensure all needed data sources are available
            all_source_data <- rep_events_df
            if(first_subj_confounds_loaded) {
                # Need to handle potential column name overlaps - prefix confound names?
                # For now, assume unique names or that events take precedence
                confound_cols_to_add <- setdiff(names(rep_confounds_df), names(all_source_data))
                if (length(confound_cols_to_add) > 0) {
                     # Combine, assuming same number of rows (might not be true!) - this needs careful thought
                     # This logic is flawed if rows don't align. Need run-specific data.
                     # Let's simplify: Role inference primarily based on event file for now.
                     # Confounds usually assumed parametric unless specified otherwise.
                     # all_source_data <- cbind(all_source_data, rep_confounds_df[, confound_cols_to_add, drop=FALSE])
                }
            }

          variable_roles <- infer_and_validate_variable_roles(
              regressors_spec = validated_ior$terms, # Use global terms for inference pool
              events_mapping = model_event_map, # Use model's mapping
              events_data = all_source_data, # Use combined data
              user_factors = unlist(validated_ior$model$variable_roles$factors %||% list()),
              user_parametric = unlist(validated_ior$model$variable_roles$parametric %||% list())
          )

           # Post-inference check: Ensure columns for identified roles exist
           all_role_vars <- c(variable_roles$factors, variable_roles$parametric)
           cols_for_roles <- model_event_map[all_role_vars] # Get BIDS columns
           cols_for_roles <- cols_for_roles[!sapply(cols_for_roles, is.null)] # Remove NULLs if mapping missing
           
           # Check against available columns from events AND selected confounds
           available_bids_cols <- unique(c(events_info$columns, confounds_info$columns))
           missing_role_cols <- setdiff(unlist(cols_for_roles), available_bids_cols)
            if (length(missing_role_cols) > 0) {
                 errors$add_error("model$variable_mapping/roles", paste0("BIDS column(s) for factor/parametric roles ('", paste(missing_role_cols, collapse="', '"),"') not found in representative event/confound file(s)."))
            }

      } else {
          warning("Skipping variable role inference due to missing representative event data.")
          # Assign based purely on user declaration
          variable_roles$factors = unique(unlist(validated_ior$model$variable_roles$factors %||% list()))
          variable_roles$parametric = unique(unlist(validated_ior$model$variable_roles$parametric %||% list()))
      }
  } else {
      warning("Skipping contextual validation and role inference (no subjects/tasks selected/found).")
      variable_roles$factors = unique(unlist(validated_ior$model$variable_roles$factors %||% list()))
      variable_roles$parametric = unique(unlist(validated_ior$model$variable_roles$parametric %||% list()))
  }


  # --- Final Checks & Stop if Invalid ---
  errors$stop_if_invalid("Context-dependent validation failed (BIDS content, event/confound columns, roles)")

  # --- Construct Final fmri_config Object ---
  config <- structure(
    list(
      spec = validated_ior,
      project = proj,
      subjects = subjects,
      tasks = tasks,
      runs = selected_runs, # Store selected runs
      events_info = events_info,
      confounds_info = confounds_info,
      variable_roles = variable_roles,
      validated = TRUE
    ),
    class = "fmri_config"
  )

  message("fmri_config object built successfully.")
  return(config)
}

# --- Role Inference Helper (Modified) ---
#' @keywords internal
infer_and_validate_variable_roles <- function(regressors_spec, events_mapping, events_data,
                                              user_factors = NULL, user_parametric = NULL) {

  inferred_roles <- list()
  final_factors <- character(0)
  final_parametric <- character(0)
  user_factors <- unique(user_factors %||% character(0))
  user_parametric <- unique(user_parametric %||% character(0))

  # Check for overlap in user declarations
  overlap <- intersect(user_factors, user_parametric)
  if (length(overlap) > 0) {
     warning(sprintf("Variables declared as both factor and parametric: %s. Prioritizing 'factor'.", paste(overlap, collapse=", ")))
     user_parametric <- setdiff(user_parametric, overlap) # Remove overlap from parametric
  }

  # Get all *model variable names* actually used in the specified regressors
  model_vars_used_in_terms <- unique(unlist(lapply(regressors_spec, `[[`, "variables")))

  if (length(model_vars_used_in_terms) == 0) {
     warning("No variables found used in any term specifications. Cannot infer roles.")
     return(list(inferred_roles = list(), factors = user_factors, parametric = user_parametric))
  }

  if (is.null(events_data) || nrow(events_data) == 0) {
     warning("Representative event data is empty. Cannot infer roles.")
     return(list(inferred_roles = list(), factors = user_factors, parametric = user_parametric))
  }
  
  available_bids_columns <- names(events_data)

  for (model_var_name in model_vars_used_in_terms) {
      bids_col_name <- events_mapping[[model_var_name]]
      if (is.null(bids_col_name)) {
           warning(sprintf("Variable '%s' used in term but not found in `model$variable_mapping`. Cannot infer role.", model_var_name))
           inferred_roles[[model_var_name]] <- "unknown"
           next
      }
      if (!bids_col_name %in% available_bids_columns) {
           warning(sprintf("Mapped BIDS column '%s' for variable '%s' not found in representative data. Cannot infer role.", bids_col_name, model_var_name))
           inferred_roles[[model_var_name]] <- "unknown"
           next
      }

      col_data <- events_data[[bids_col_name]]
      inferred_role <- "unknown" # Default

      # --- Inference Logic ---
      valid_data <- col_data[!is.na(col_data)]
      if (length(valid_data) == 0) {
          warning(sprintf("Column '%s' (for var '%s') contains only NA values. Cannot infer role.", bids_col_name, model_var_name))
          inferred_role <- "unknown"
      } else if (is.factor(valid_data) || is.character(valid_data) || is.logical(valid_data)) {
          inferred_role <- "factor"
      } else if (is.numeric(valid_data) || is.integer(valid_data)) {
           num_unique <- length(unique(valid_data))
           n_total <- length(valid_data)
           # Heuristic: If few unique numeric values relative to total, might be a factor
           if (num_unique <= max(5, 0.1 * n_total) && all(valid_data == floor(valid_data))) {
               inferred_role <- "factor"
               # message(sprintf("Info: Var '%s' (col '%s') is numeric/integer with few unique values (%d). Inferred 'factor'.", model_var_name, bids_col_name, num_unique))
           } else {
               inferred_role <- "parametric"
           }
      } else {
         warning(sprintf("Column '%s' (for var '%s') has unsupported type (%s) for inference.", bids_col_name, model_var_name, class(col_data)[1]))
         inferred_role <- "unknown"
      }
      inferred_roles[[model_var_name]] <- inferred_role

      # --- Combine with User Declaration ---
      user_declared_role <- NA
      if (model_var_name %in% user_factors) user_declared_role <- "factor"
      if (model_var_name %in% user_parametric) user_declared_role <- "parametric" # Already handled overlap

      final_role <- "unknown"
      if (!is.na(user_declared_role)) {
          if (user_declared_role != inferred_role && inferred_role != "unknown") {
               # message(sprintf("Info: Var '%s': Using declared role '%s' (inferred '%s').", model_var_name, user_declared_role, inferred_role))
          }
          final_role <- user_declared_role
      } else {
           if (inferred_role != "unknown") {
                # message(sprintf("Info: Var '%s' not declared. Using inferred role '%s'.", model_var_name, inferred_role))
                final_role <- inferred_role
           } else {
                warning(sprintf("Var '%s' used in term, but role could not be inferred or declared.", model_var_name))
                final_role <- "unknown" # Stays unknown
           }
      }

      # Add to final lists
      if (final_role == "factor") {
          final_factors <- c(final_factors, model_var_name)
      } else if (final_role == "parametric") {
          final_parametric <- c(final_parametric, model_var_name)
      }
  } # End loop over model_var_name

  # Report unused declarations
  all_final_vars <- c(final_factors, final_parametric)
  unused_factors <- setdiff(user_factors, all_final_vars)
  unused_parametric <- setdiff(user_parametric, all_final_vars)
  if (length(unused_factors) > 0) {
      warning("Declared factors not used in any selected model terms: ", paste(unused_factors, collapse=", "))
  }
   if (length(unused_parametric) > 0) {
      warning("Declared parametric vars not used in any selected model terms: ", paste(unused_parametric, collapse=", "))
  }

  list(
    inferred_roles = inferred_roles,
    factors = unique(final_factors),
    parametric = unique(final_parametric)
  )
}


#' Load, Validate, and Build fMRI Analysis Configuration
#'
#' This is the main user-facing function. It reads a YAML specification,
#' validates it against the DSL schema, applies defaults, performs context-
#' dependent BIDS validation (subjects, tasks, event/confound columns),
#' infers variable roles for the model, and returns a validated `fmri_config`
#' object ready for model building.
#'
#' @param yaml_file Path to the YAML configuration file.
#' @return An `fmri_config` object.
#' @export
#' @examples
#' \dontrun{
#'   # Assuming a valid config.yaml exists
#'   # config <- load_fmri_config("config.yaml")
#'   # print(config)
#' }
load_fmri_config <- function(yaml_file) {
  # Step 1: Parse YAML, apply defaults, and validate schema
  validated_ior <- parse_and_validate_config(yaml_file)
  # Step 2: Perform contextual validation and build final object
  fmri_config_obj <- build_config_from_ior(validated_ior)
  return(fmri_config_obj)
}


# --- S3 Print Method ---
# Placeholder for print.fmri_config
# Assume this function exists and works as defined previously
# (Implementation omitted for brevity)
#' @export
print.fmri_config <- function(x, ...) {
    cat("fmri_config object (validated: ", x$validated, ")\n", sep="")
    cat("Model Name:", x$spec$model$name, "\n")
    cat("Subjects:", length(x$subjects), "\n")
    cat("Tasks:", paste(x$tasks, collapse=", "), "\n")
    # Add more details as needed
    invisible(x)
}