###############################################################################
#                   Error Handling and Validation Helpers
###############################################################################

# The `%||%` operator
#' @keywords internal
`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}

# Simple error collector class (optional but helpful)
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
# (These will be used by parse_and_validate_config)

#' @keywords internal
check_required <- function(data, field, path, errors) {
  if (is.null(data[[field]])) {
    # Construct the full path for the error message
    field_path <- paste0(path, if (path != "") "$", field)
    errors$add_error(field_path, paste0("Required field '", field, "' is missing."))
    return(FALSE)
  }
  return(TRUE)
}

#' @keywords internal
check_type <- function(data, field, expected_type, path, errors, allow_null = FALSE) {
  value <- data[[field]]
  field_path <- paste0(path, if (path != "") "$", field)

  if (is.null(value)) {
    # If null is allowed, it's valid in terms of type check here.
    # If it was required, check_required would have caught it.
    return(allow_null)
  }

  valid <- FALSE
  actual_type_desc <- paste(class(value), collapse="/") # More informative for lists etc.

  if (expected_type == "string") {
      valid <- is.character(value) && length(value) == 1 && !is.na(value)
  } else if (expected_type == "integer") {
      valid <- (is.integer(value) || (is.numeric(value) && !is.na(value) && value == floor(value))) && length(value) == 1
  } else if (expected_type == "number") {
      valid <- is.numeric(value) && length(value) == 1 && !is.na(value)
  } else if (expected_type == "boolean") {
      valid <- is.logical(value) && length(value) == 1 && !is.na(value)
  } else if (expected_type == "object") {
      # Must be a named list (can be empty)
      valid <- is.list(value) && (length(value) == 0 || (!is.null(names(value)) && all(names(value) != "")))
  } else if (grepl("^array", expected_type)) {
      # Allows zero-length arrays (empty lists or vectors)
      element_type <- sub("^array\\[(\\w+)\\]$", "\\1", expected_type)
      if (element_type == expected_type) element_type <- NULL # Simple "array" type

      # Check if list or non-atomic vector (like list of lists)
      if (is.list(value) || (is.vector(value) && !is.atomic(value))) { 
          valid <- TRUE # Assume valid until proven otherwise
          if (length(value) > 0 && !is.null(element_type)) {
              # Check element types
               all_elems_valid <- all(sapply(value, function(el) {
                   # Check type of each element based on element_type
                   if (element_type == "string") is.character(el) && length(el) == 1 && !is.na(el)
                   else if (element_type == "integer") (is.integer(el) || (is.numeric(el) && !is.na(el) && el == floor(el))) && length(el) == 1
                   else if (element_type == "number") is.numeric(el) && length(el) == 1 && !is.na(el)
                   else if (element_type == "boolean") is.logical(el) && length(el) == 1 && !is.na(el)
                   else if (element_type == "object") is.list(el) && (length(el) == 0 || (!is.null(names(el)) && all(names(el) != "")))
                   else FALSE # Unknown element type
               }))
               if (!all_elems_valid) {
                   valid <- FALSE
                   actual_type_desc <- paste0(actual_type_desc, " (with invalid element types)")
               }
          }
      # Check for atomic vectors (simpler types like numeric(5), character(3)) 
      } else if (is.atomic(value)) { 
            valid <- TRUE
            if (length(value) > 0 && !is.null(element_type)) {
                 el_class_check <- switch(element_type,
                                           "string" = "character",
                                           "integer" = "integer", # or numeric check
                                           "number" = "numeric",
                                           "boolean" = "logical",
                                           NULL) # Unknown
                 if (!is.null(el_class_check)) {
                     # Basic check for atomic vectors
                     if (el_class_check == "integer") {
                          # Allow numeric vectors that are all integers
                          valid <- (is.integer(value) || (is.numeric(value) && all(!is.na(value)) && all(value == floor(value))))
                     } else {
                          valid <- inherits(value, el_class_check) && all(!is.na(value))
                     }

                     if (!valid) actual_type_desc <- paste0(actual_type_desc, " (invalid element types)")
                 } else {
                    valid <- FALSE # Unknown element type for atomic vector check
                 }
            }
      } else {
          valid <- FALSE # Not a list or known vector type
      }
  } else {
      # Unknown expected type
      errors$add_error(field_path, paste0("Internal validation error: Unknown expected type '", expected_type, "' specified."))
      return(FALSE)
  }

  if (!valid) {
    errors$add_error(field_path, paste0("Invalid type. Expected '", expected_type, "', but got '", actual_type_desc, "'."))
  }
  return(valid)
}

#' @keywords internal
check_enum <- function(data, field, allowed_values, path, errors, allow_null = FALSE) {
   value <- data[[field]]
   field_path <- paste0(path, if (path != "") "$", field)
   if (is.null(value)) {
       # Handled by allow_null. If it's required, check_required handles it.
       return(allow_null)
   }

   # Handle arrays/lists - check each element
   values_to_check <- unlist(value) # Flatten list/vector

   if (length(values_to_check) == 0 && allow_null) return(TRUE) # Empty list/vector might be ok if null allowed
   # Handle case where unlist produces NULL for list(NULL)?
   if (is.null(values_to_check) && length(value) > 0) {
        errors$add_error(field_path, "Contains NULL values where specific choices are expected.")
        return(FALSE)
   }

   invalid_values <- setdiff(values_to_check, allowed_values)

   if (length(invalid_values) > 0) {
       errors$add_error(field_path, paste0("Invalid value(s): '", paste(invalid_values, collapse="', '"), "'. Must be one of: ", paste(allowed_values, collapse=", ")))
       return(FALSE)
   }
   return(TRUE)
}

#' @keywords internal
check_pattern <- function(data, field, pattern, path, errors, allow_null = FALSE) {
    value <- data[[field]]
    field_path <- paste0(path, if (path != "") "$", field)
    if (is.null(value)) {
        return(allow_null)
    }

    # Handle arrays/lists - check each element
    values_to_check <- unlist(value) # Flatten list/vector

    if (length(values_to_check) == 0 && allow_null) return(TRUE)
    if (is.null(values_to_check) && length(value) > 0) {
        # If original value was a list containing NULLs
         errors$add_error(field_path, "Pattern check applied to list containing NULL values.")
         return(FALSE) 
    }

    if (!is.character(values_to_check)) {
        # Type check should catch this, but we can be robust
        # Let type check handle the fundamental type.
         return(TRUE) 
    }

    does_not_match <- !grepl(pattern, values_to_check)
    if (any(does_not_match)) {
        errors$add_error(field_path, paste0("Value(s) do not match required pattern ('", pattern, "'): '", paste(values_to_check[does_not_match], collapse="', '"), "'."))
        return(FALSE)
    }
    return(TRUE)
}

#' @keywords internal
check_dir_exists <- function(data, field, path, errors) {
  # Basic check during initial validation
  dir_path <- data[[field]]
  field_path <- paste0(path, if (path != "") "$", field)
  # Check only if the field is actually present
  if (!is.null(dir_path)) {
      # Ensure it's a string first (using check_type)
      if (!check_type(data, field, "string", path, errors)) return(FALSE)

      # Use fs::path_expand to handle tilde etc.
      expanded_path <- fs::path_expand(dir_path)
      if (!fs::dir_exists(expanded_path)) {
        errors$add_error(field_path, paste0("Directory not found: '", dir_path, "' (expanded to '", expanded_path, "')."))
        return(FALSE)
      }
  } else {
     # If the field itself is null, it might be okay if not required.
  }
  return(TRUE)
}

# --- Default Application Helper ---
#' @keywords internal
apply_defaults <- function(data, defaults) {
  if (!is.list(data) || !is.list(defaults)) return(data)

  # Ensure data is a named list if it's going to be treated like an object
  # Only apply defaults if 'defaults' has names
  if (!is.null(names(defaults))) {
      # Add default values for keys missing in data
      missing_keys <- setdiff(names(defaults), names(data))
      for (key in missing_keys) {
        data[[key]] <- defaults[[key]]
      }

      # Recursively apply defaults for nested lists/objects that are named lists in both
      for (key in intersect(names(data), names(defaults))) {
          # Check if both data[[key]] and defaults[[key]] are potential objects (named lists)
          is_data_obj <- is.list(data[[key]]) && (!is.null(names(data[[key]])) || length(data[[key]])==0) # Allow empty list
          is_defaults_obj <- is.list(defaults[[key]]) && !is.null(names(defaults[[key]]))

          if (is_data_obj && is_defaults_obj) {
             data[[key]] <- apply_defaults(data[[key]], defaults[[key]])
          }
      }
  }
  return(data)
}

###############################################################################
#             Step 1: Parse, Validate Schema, Apply Defaults (Create IOR)
###############################################################################

#' Parse, Validate, and Normalize fMRI Configuration from YAML
#'
#' Reads a YAML file, validates its structure and types against the DSL schema,
#' applies default values, and returns a validated intermediate list representation (IOR).
#' Stops with an informative error if validation fails.
#'
#' @param yaml_file Path to the YAML configuration file.
#' @return A nested list representing the validated and normalized configuration (IOR).
#' @keywords internal
#' @importFrom yaml read_yaml
#' @importFrom fs file_exists dir_exists path_expand
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
  # Handle completely empty YAML file case
  if (is.null(raw_config_list)) {
       stop("YAML file '", yaml_file, "' is empty or could not be parsed correctly.", call. = FALSE)
  }


  errors <- ValidationErrors$new()

  # --- 2. Define Defaults ---
  default_config <- list(
      dataset = list(
          relpath = "func",
          subjects = list(include = NULL, exclude = NULL),
          tasks = NULL,
          runs = NULL,
          scan_params = list(TR = NULL, run_length = NULL)
      ),
      confounds = list(include = NULL, exclude = NULL),
      baseline = list(
          basis = "bspline",
          degree = 3L, 
          intercept = "runwise",
          confounds = NULL 
      ),
      hrfs = list(), 
      parametric = list(), 
      events = list(), # Essential fields onset, duration, block are required, no defaults
      regressors = list(), 
      contrasts = list(), 
      model = list( 
          factors = list(), # Default to empty lists instead of NULL
          parametric = list(),
          events = NULL, # No default for model-specific event mapping
          regressors = NULL, # No default for model-specific regressors
          contrasts = NULL # No default for model-specific contrasts
      )
  )

  # --- 3. Apply Defaults --- 
  config_list <- apply_defaults(raw_config_list, default_config)


  # --- 4. Perform Schema Validation --- 

  # Top Level Required Fields
  check_required(config_list, "dataset", "", errors)
  check_required(config_list, "events", "", errors)
  check_required(config_list, "regressors", "", errors)
  check_required(config_list, "model", "", errors)

  # Dataset Validation
  if (check_type(config_list, "dataset", "object", "", errors)) {
      ds_path = "dataset"
      ds_data = config_list$dataset
      if(check_required(ds_data, "path", ds_path, errors)) {
          check_type(ds_data, "path", "string", ds_path, errors)
          check_dir_exists(ds_data, "path", ds_path, errors)
      }
      check_required(ds_data, "relpath", ds_path, errors)
      check_type(ds_data, "relpath", "string", ds_path, errors, allow_null=FALSE)

      # Subjects: Optional section, but if present, must be object. 
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

      # Scan Params: Optional section
       if(check_type(ds_data, "scan_params", "object", ds_path, errors, allow_null = TRUE)) {
            if (!is.null(ds_data$scan_params)) { # Only validate substructure if scan_params is present
               sp_data <- ds_data$scan_params
               sp_path <- paste0(ds_path, "$scan_params")
               # TR Validation: Optional section within scan_params
                if(check_type(sp_data, "TR", "object", sp_path, errors, allow_null = TRUE)) {
                    if (!is.null(sp_data$TR)) { # Validate only if TR section exists
                        tr_data <- sp_data$TR
                        tr_path <- paste0(sp_path, "$TR")
                        check_required(tr_data, "default", tr_path, errors)
                        check_type(tr_data, "default", "number", tr_path, errors)
                        if(check_type(tr_data, "overrides", "array[object]", tr_path, errors, allow_null=TRUE)) {
                            if (!is.null(tr_data$overrides)) {
                               for(i in seq_along(tr_data$overrides)){
                                  ovr_path <- paste0(tr_path, "$overrides[", i, "]")
                                  ovr_item <- tr_data$overrides[[i]]
                                  if(is.list(ovr_item) && !is.null(names(ovr_item))) {
                                      check_required(ovr_item, "pattern", ovr_path, errors)
                                      check_type(ovr_item, "pattern", "string", ovr_path, errors)
                                      check_required(ovr_item, "value", ovr_path, errors)
                                      check_type(ovr_item, "value", "number", ovr_path, errors)
                                  } else { errors$add_error(ovr_path, "Override item must be an object with 'pattern' and 'value'.")}
                               }
                            }
                        }
                    }
                }
                # Run Length Validation: Optional section within scan_params
                if(check_type(sp_data, "run_length", "object", sp_path, errors, allow_null = TRUE)) {
                    if (!is.null(sp_data$run_length)) { # Validate only if run_length section exists
                         rl_data <- sp_data$run_length
                         rl_path <- paste0(sp_path, "$run_length")
                         check_required(rl_data, "default", rl_path, errors)
                         if (check_type(rl_data, "default", "object", rl_path, errors)) {
                             def_obj = rl_data$default
                             all_num = all(sapply(def_obj, is.numeric)) && all(sapply(def_obj, length) == 1)
                             if (!all_num) errors$add_error(paste0(rl_path,"$default"), "All values within 'default' must be single numbers.")
                             # Check that keys are strings (task names)
                             if (!is.null(names(def_obj)) && !all(sapply(names(def_obj), is.character))) {
                                  errors$add_error(paste0(rl_path,"$default"), "Keys within 'default' must be strings (task names).")
                             }
                         }
                         if(check_type(rl_data, "overrides", "array[object]", rl_path, errors, allow_null=TRUE)) {
                              if (!is.null(rl_data$overrides)) {
                                for(i in seq_along(rl_data$overrides)){
                                    ovr_path <- paste0(rl_path, "$overrides[", i, "]")
                                    ovr_item <- rl_data$overrides[[i]]
                                    if(is.list(ovr_item) && !is.null(names(ovr_item))) {
                                        check_required(ovr_item, "pattern", ovr_path, errors)
                                        check_type(ovr_item, "pattern", "string", ovr_path, errors)
                                        check_required(ovr_item, "value", ovr_path, errors)
                                        check_type(ovr_item, "value", "number", ovr_path, errors) 
                                    } else { errors$add_error(ovr_path, "Override item must be an object with 'pattern' and 'value'.")}
                                }
                            }
                         }
                    }
                }
            }
       }
  }

  # Confounds (Optional section)
   # Check type only if the key exists after defaults (it might be NULL)
   if (check_type(config_list, "confounds", "object", "confounds", errors, allow_null=TRUE)) {
        if (!is.null(config_list$confounds)) {
             conf_path = "confounds"
             conf_data = config_list$confounds
             check_type(conf_data, "include", "array[string]", conf_path, errors, allow_null=TRUE)
             check_type(conf_data, "exclude", "array[string]", conf_path, errors, allow_null=TRUE)
        }
   } 

  # Baseline (Optional section)
   if (check_type(config_list, "baseline", "object", "baseline", errors, allow_null=TRUE)) {
       if (!is.null(config_list$baseline)) {
          bl_path = "baseline"
          bl_data = config_list$baseline 
          check_required(bl_data, "basis", bl_path, errors)
          check_type(bl_data, "basis", "string", bl_path, errors)
          check_enum(bl_data, "basis", c("constant", "poly", "bspline", "ns"), bl_path, errors)
          check_required(bl_data, "degree", bl_path, errors)
          check_type(bl_data, "degree", "integer", bl_path, errors)
          check_required(bl_data, "intercept", bl_path, errors)
          check_type(bl_data, "intercept", "string", bl_path, errors)
          check_enum(bl_data, "intercept", c("runwise", "global", "none"), bl_path, errors)
          # Validate baseline$confounds: Optional section
          if(check_type(bl_data, "confounds", "object", bl_path, errors, allow_null = TRUE)) {
              if(!is.null(bl_data$confounds)){
                 bl_conf_path <- paste0(bl_path, "$confounds")
                 bl_conf_data <- bl_data$confounds
                 check_type(bl_conf_data, "include", "array[string]", bl_conf_path, errors, allow_null=TRUE)
                 check_type(bl_conf_data, "exclude", "array[string]", bl_conf_path, errors, allow_null=TRUE)
              }
          }
      }
   } 

  # HRFs (Optional section, defaults to empty list)
   if (check_type(config_list, "hrfs", "object", "hrfs", errors)) {
        hrfs_path = "hrfs"
        hrfs_data <- config_list$hrfs
        if (!is.null(names(hrfs_data))) { # Check if it's actually a named list (object)
           if (length(hrfs_data) > 0) {
              for (hrf_name in names(hrfs_data)) {
                  hrf_path <- paste0(hrfs_path, "$", hrf_name)
                  if(check_type(hrfs_data, hrf_name, "object", hrfs_path, errors)) {
                      hrf_spec <- hrfs_data[[hrf_name]]
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
        } else if (length(hrfs_data) > 0 && is.null(names(hrfs_data))) {
            errors$add_error(hrfs_path, "'hrfs' must be an object (named list), not an array.")
        }
        # Empty hrfs: {} is valid and results in length=0, names=NULL, passes check_type("object")
   } 

   # Parametric Transformations (Optional section, defaults to empty list)
   if (check_type(config_list, "parametric", "object", "parametric", errors)) {
        par_path <- "parametric"
        par_data <- config_list$parametric
        if (!is.null(names(par_data))) { # Check if named list (object)
           if (length(par_data) > 0) {
               for (var_name in names(par_data)) {
                   var_path <- paste0(par_path, "$", var_name)
                   if(check_type(par_data, var_name, "object", par_path, errors)) {
                       var_spec <- par_data[[var_name]]
                       check_type(var_spec, "transform", "array[string]", var_path, errors, allow_null = TRUE)
                       if (!is.null(var_spec$transform)) {
                           allowed_transforms <- c("center", "scale", "within_subject", "log", "exp", "zscore")
                           check_enum(var_spec, "transform", allowed_transforms, var_path, errors, allow_null = TRUE)
                       }
                       if(check_type(var_spec, "basis", "object", var_path, errors, allow_null=TRUE)) {
                            if (!is.null(var_spec$basis)) {
                                basis_path <- paste0(var_path, "$basis")
                                basis_spec_data <- var_spec$basis
                                check_required(basis_spec_data, "type", basis_path, errors)
                                check_type(basis_spec_data, "type", "string", basis_path, errors)
                                allowed_basis_types <- c("polynomial", "bspline", "ns", "Poly", "BSpline", "Standardized", "Ident") 
                                check_enum(basis_spec_data, "type", allowed_basis_types, basis_path, errors)
                                check_type(basis_spec_data, "parameters", "object", basis_path, errors, allow_null=TRUE)
                            }
                       }
                   }
               }
            }
        } else if (length(par_data) > 0 && is.null(names(par_data))) {
             errors$add_error(par_path, "'parametric' must be an object (named list), not an array.")
        }
   } 

   # Events (Required Section, must be object)
   if (check_type(config_list, "events", "object", "", errors)) {
        ev_path = "events"
        ev_data = config_list$events
        check_required(ev_data, "onset", ev_path, errors)
        check_type(ev_data, "onset", "string", ev_path, errors)
        check_required(ev_data, "duration", ev_path, errors)
        check_type(ev_data, "duration", "string", ev_path, errors)
        check_required(ev_data, "block", ev_path, errors)
        check_type(ev_data, "block", "string", ev_path, errors)
        var_map_names <- setdiff(names(ev_data), c("onset", "duration", "block"))
         if (length(var_map_names) > 0) {
             for (var_name in var_map_names) {
                 var_path <- paste0(ev_path, "$", var_name)
                 if (check_type(ev_data, var_name, "object", ev_path, errors)) {
                    var_spec <- ev_data[[var_name]]
                    check_required(var_spec, "column", var_path, errors)
                    check_type(var_spec, "column", "string", var_path, errors)
                    check_type(var_spec, "onset", "string", var_path, errors, allow_null=TRUE)
                    check_type(var_spec, "duration", "string", var_path, errors, allow_null=TRUE)
                 }
             }
         }
   }

   # Regressors (Required Section, must be object)
   if (check_type(config_list, "regressors", "object", "", errors)) {
        regs_path <- "regressors"
        regs_data <- config_list$regressors
        if (!is.null(names(regs_data))) { # Check if named list (object)
           if (length(regs_data) > 0) {
                defined_hrf_names <- names(config_list$hrfs %||% list())
                for (reg_name in names(regs_data)) {
                    reg_path <- paste0(regs_path, "$", reg_name)
                    if(check_type(regs_data, reg_name, "object", regs_path, errors)) {
                        reg_spec <- regs_data[[reg_name]]
                        check_required(reg_spec, "type", reg_path, errors)
                        check_type(reg_spec, "type", "string", reg_path, errors)
                        allowed_reg_types <- c("hrf", "hrf_parametric", "trialwise") 
                        check_enum(reg_spec, "type", allowed_reg_types, reg_path, errors)
                        check_required(reg_spec, "variables", reg_path, errors)
                        valid_vars = check_type(reg_spec, "variables", "array[string]", reg_path, errors)
                        if (valid_vars && length(reg_spec$variables) == 0){
                             errors$add_error(paste0(reg_path, "$variables"), "Field 'variables' cannot be empty.")
                        }
                        check_type(reg_spec, "hrf", "string", reg_path, errors, allow_null=TRUE)
                        check_type(reg_spec, "lag", "number", reg_path, errors, allow_null=TRUE)
                        check_type(reg_spec, "subset", "string", reg_path, errors, allow_null=TRUE)
                         check_type(reg_spec, "transform", "array[string]", reg_path, errors, allow_null=TRUE)
                         if (!is.null(reg_spec$transform)) {
                           allowed_reg_transforms <- c("center", "scale", "zscore", "log", "exp") 
                           check_enum(reg_spec, "transform", allowed_reg_transforms, reg_path, errors, allow_null=TRUE)
                         }
                        if (!is.null(reg_spec$type) && reg_spec$type == "hrf_parametric") {
                            if (check_required(reg_spec, "basis", reg_path, errors) &&
                                check_type(reg_spec, "basis", "object", reg_path, errors)) {
                                   basis_path <- paste0(reg_path, "$basis")
                                   basis_spec_data <- reg_spec$basis
                                   check_required(basis_spec_data, "type", basis_path, errors)
                                   check_type(basis_spec_data, "type", "string", basis_path, errors)
                                   allowed_basis_types <- c("Poly", "BSpline", "Standardized", "Ident", "polynomial", "bspline", "ns") 
                                   check_enum(basis_spec_data, "type", allowed_basis_types, basis_path, errors)
                                   check_type(basis_spec_data, "parameters", "object", basis_path, errors, allow_null=TRUE)
                            }
                            if(valid_vars && length(reg_spec$variables) != 1) {
                                errors$add_error(paste0(reg_path, "$variables"), "'hrf_parametric' regressors require exactly one variable.")
                            }
                        } else {
                            if (!is.null(reg_spec$basis)) {
                                errors$add_error(paste0(reg_path, "$basis"), "'basis' section is only applicable for 'hrf_parametric' regressors.")
                            }
                        }
                    }
                }
           } else {
                errors$add_error(regs_path, "Required section 'regressors' cannot be empty.")
           }
        } else if (length(regs_data) > 0 && is.null(names(regs_data))) {
            errors$add_error(regs_path, "'regressors' must be an object (named list), not an array.")
        }
   }

   # Contrasts (Optional Section)
    if (check_type(config_list, "contrasts", "object", "contrasts", errors, allow_null=TRUE)) {
         if (!is.null(config_list$contrasts)) {
             cons_path <- "contrasts"
             cons_data <- config_list$contrasts
             if (!is.null(names(cons_data))) { # Check if named list
                if (length(cons_data) > 0) {
                    for (con_name in names(cons_data)) {
                        con_path <- paste0(cons_path, "$", con_name)
                        if(check_type(cons_data, con_name, "object", cons_path, errors)) {
                           con_spec <- cons_data[[con_name]]
                           check_required(con_spec, "type", con_path, errors)
                           check_type(con_spec, "type", "string", con_path, errors)
                           allowed_con_types <- c("formula", "pair", "one_against_all", "unit", "oneway", "interaction", "polynomial")
                           check_enum(con_spec, "type", allowed_con_types, con_path, errors)
                           con_type <- con_spec$type 
                           if (!is.null(con_type)) {
                               if (con_type == "formula") {
                                   check_required(con_spec, "expression", con_path, errors)
                                   check_type(con_spec, "expression", "string", con_path, errors)
                               } else if (con_type %in% c("pair", "one_against_all", "oneway", "interaction", "polynomial")) {
                                   check_type(con_spec, "factors", "array[string]", con_path, errors, allow_null=TRUE)
                               }
                               # Check degree for polynomial?
                               if (con_type == "polynomial") {
                                    check_type(con_spec, "degree", "integer", con_path, errors, allow_null=TRUE)
                               }
                               # Check levels for one_against_all?
                               if (con_type == "one_against_all") {
                                    check_type(con_spec, "levels", "array[string]", con_path, errors, allow_null=TRUE)
                               }
                           }
                           check_type(con_spec, "where", "string", con_path, errors, allow_null=TRUE)
                           check_type(con_spec, "include", "string", con_path, errors, allow_null=TRUE)
                           check_type(con_spec, "exclude", "string", con_path, errors, allow_null=TRUE)
                        }
                    }
                } # else: Empty contrasts object {} is valid
            } else if (length(cons_data) > 0 && is.null(names(cons_data))) {
                 errors$add_error(cons_path, "'contrasts' must be an object (named list), not an array.")
            }
         }
    } 

  # Model (Required Section)
  if (check_type(config_list, "model", "object", "", errors)) {
      mod_path = "model"
      mod_data = config_list$model
      check_required(mod_data, "name", mod_path, errors)
      check_type(mod_data, "name", "string", mod_path, errors)
      check_type(mod_data, "factors", "array[string]", mod_path, errors, allow_null=TRUE)
      check_type(mod_data, "parametric", "array[string]", mod_path, errors, allow_null=TRUE)
      
      # Model-specific event mapping: required, object
      if(check_required(mod_data, "events", mod_path, errors)) {
          if (check_type(mod_data, "events", "object", mod_path, errors)) {
              mod_ev_path <- paste0(mod_path, "$events")
              mod_ev_data <- mod_data$events
              check_required(mod_ev_data, "onset", mod_ev_path, errors)
              check_type(mod_ev_data, "onset", "string", mod_ev_path, errors)
              check_required(mod_ev_data, "duration", mod_ev_path, errors)
              check_type(mod_ev_data, "duration", "string", mod_ev_path, errors)
              check_required(mod_ev_data, "block", mod_ev_path, errors)
              check_type(mod_ev_data, "block", "string", mod_ev_path, errors)
              mod_var_names <- setdiff(names(mod_ev_data), c("onset", "duration", "block"))
              for(mod_var in mod_var_names) {
                  check_type(mod_ev_data, mod_var, "string", mod_ev_path, errors)
              }
          }
      }

      # Model-specific regressors: required, object
      if (check_required(mod_data, "regressors", mod_path, errors)) {
         if(check_type(mod_data, "regressors", "object", mod_path, errors)) {
             model_reg_names <- names(mod_data$regressors %||% list())
             top_level_reg_names <- names(config_list$regressors %||% list())
             missing_regs <- setdiff(model_reg_names, top_level_reg_names)
             if (length(missing_regs) > 0) {
                errors$add_error(paste0(mod_path, "$regressors"), paste("Regressor(s) used in model but not defined in top-level 'regressors':", paste(missing_regs, collapse=", ")))
             }
             if (any(!sapply(mod_data$regressors, is.null))) {
                 # Allow non-null values if schema permits overrides here, otherwise warn/error
                 # errors$add_error(paste0(mod_path, "$regressors"), "Model 'regressors' values should be null if just referencing top-level.")
             }
         }
      }

      # Model-specific contrasts: optional, object if present
      if (check_type(mod_data, "contrasts", "object", mod_path, errors, allow_null=TRUE)) {
           if (!is.null(mod_data$contrasts)) {
              model_con_names <- names(mod_data$contrasts %||% list())
              top_level_con_names <- names(config_list$contrasts %||% list())
              missing_cons <- setdiff(model_con_names, top_level_con_names)
              if (length(missing_cons) > 0) {
                 errors$add_error(paste0(mod_path, "$contrasts"), paste("Contrast(s) used in model but not defined in top-level 'contrasts':", paste(missing_cons, collapse=", ")))
              }
              if (any(!sapply(mod_data$contrasts, is.null))) {
                   # Allow non-null values if schema permits overrides here
                   # errors$add_error(paste0(mod_path, "$contrasts"), "Model 'contrasts' values should be null if just referencing top-level.")
               }
           }
      }
  }


  # --- 5. Stop if validation failed ---
  errors$stop_if_invalid(paste("YAML configuration schema validation failed for:", yaml_file))

  # --- 6. Return the validated and normalized list (IOR) ---
  message("YAML configuration schema validated successfully: ", yaml_file)
  attr(config_list, "validated_schema") <- TRUE 
  return(config_list) 
}


###############################################################################
#             Step 2: Build fmri_config Object from IOR
###############################################################################

#' Build fmri_config Object from Validated IOR List
#'
#' Takes a validated and normalized intermediate representation (IOR) list,
#' performs context-dependent validation (BIDS checks, event/confound column checks),
#' infers variable roles, and constructs the final `fmri_config` object.
#'
#' @param validated_ior A list returned by `parse_and_validate_config`.
#' @return An `fmri_config` object ready for model building.
#' @keywords internal
#' @importFrom bidser bids_project participants tasks read_events read_confounds
#' @importFrom stats na.omit
build_config_from_ior <- function(validated_ior) {

  if (!isTRUE(attr(validated_ior, "validated_schema"))) {
      stop("Input 'validated_ior' must be the result of parse_and_validate_config().")
  }

  errors <- ValidationErrors$new() 

  # --- 1. Load BIDS Project ---
  bids_path <- fs::path_expand(validated_ior$dataset$path) # Expand path
  proj <- NULL 
  tryCatch({
      proj <- bidser::bids_project(bids_path, fmriprep = TRUE)
      message("BIDS project loaded successfully: ", bids_path)
  }, error = function(e) {
      errors$add_error("dataset$path", paste("Failed to load BIDS project at '", bids_path, "'. Error: ", e$message))
      errors$stop_if_invalid("Cannot proceed without a valid BIDS project")
  })

  # --- 2. Contextual Validation (using BIDS project) ---

  # Validate Subjects 
  available_ids <- character(0)
  tryCatch({
      available_ids <- bidser::participants(proj)$participant_id
      if (length(available_ids) == 0) warning("No participants found in BIDS dataset at: ", bids_path)
  }, error = function(e) {
       errors$add_error("dataset$subjects", paste("Failed to retrieve participants from BIDS project:", e$message))
       # Allow continuing if subjects list might be empty, error out later if needed
  })

  subjects_spec <- validated_ior$dataset$subjects
  subjects <- available_ids 

  if (!is.null(subjects_spec) && !is.null(subjects_spec$include)) {
      include_subs <- unlist(subjects_spec$include)
      missing_inc <- setdiff(include_subs, available_ids)
      if (length(missing_inc) > 0) errors$add_error("dataset$subjects$include", paste("Subject(s) not found in BIDS dataset:", paste(missing_inc, collapse=", ")))
      subjects <- intersect(subjects, include_subs)
  }
   if (!is.null(subjects_spec) && !is.null(subjects_spec$exclude)) {
      exclude_subs <- unlist(subjects_spec$exclude)
      missing_exc <- setdiff(exclude_subs, available_ids)
       if (length(missing_exc) > 0) warning("Subject(s) listed in 'exclude' but not found in BIDS: ", paste(missing_exc, collapse=", "))
      subjects <- setdiff(subjects, exclude_subs)
   }
   if (length(subjects) == 0 && length(available_ids) > 0) {
       # Error if selection resulted in zero subjects, but BIDS had subjects
       errors$add_error("dataset$subjects", "No subjects selected after applying include/exclude criteria.")
   } else if (length(subjects) == 0 && length(available_ids) == 0) {
       # If BIDS had no subjects to begin with
       # This might be okay or an error depending on use case, warn for now
       warning("No subjects selected because no participants were found in the BIDS dataset.")
   }

  # Validate Tasks 
   available_tasks <- character(0)
   tryCatch({
      available_tasks <- bidser::tasks(proj)
      if (length(available_tasks) == 0) warning("No tasks found in BIDS dataset.")
   }, error = function(e) {
        errors$add_error("dataset$tasks", paste("Failed to retrieve tasks from BIDS project:", e$message))
   })

   # Default to all available tasks if dataset$tasks is NULL
   tasks_spec_req <- validated_ior$dataset$tasks 
   tasks_spec_vec <- if (!is.null(tasks_spec_req)) unlist(tasks_spec_req) else available_tasks
   
   missing_tasks <- setdiff(tasks_spec_vec, available_tasks)
   if (length(missing_tasks) > 0) {
       errors$add_error("dataset$tasks", paste("Specified task(s) not found in BIDS dataset:", paste(missing_tasks, collapse=", ")))
   }
   tasks <- intersect(tasks_spec_vec, available_tasks)

   if (length(tasks) == 0) {
        if (!is.null(tasks_spec_req)) {
           # Error if tasks were explicitly specified but none were found/matched
           errors$add_error("dataset$tasks", "None of the specified tasks were found in the BIDS dataset.")
        } else {
           # No tasks specified, none found in BIDS -> might be valid (e.g. resting state).
           # Warning already issued by tasks() call if available_tasks was empty.
           tasks <- character(0) # Ensure tasks is empty character vector
        }
   }

  # Validate Event Columns (Requires loading data for one subject/task)
   events_mapping <- validated_ior$events 
   events_list_for_inference <- list() 
   final_event_columns <- NULL 

   # Only proceed if we have subjects AND tasks to check against
   if (length(subjects) > 0 && length(tasks) > 0) {
        first_subj <- subjects[1]
        events_list_for_inference[[first_subj]] <- list() 
        some_event_file_found <- FALSE

        for (task in tasks) {
            task_events_data <- NULL 
            tryCatch({
               events_read <- bidser::read_events(proj, subid = first_subj, task = task)
               if (length(events_read) > 0 && length(events_read$data) > 0 && !is.null(events_read$data[[1]])) {
                  task_events_data <- events_read$data[[1]] 
                  some_event_file_found <- TRUE
                  events_list_for_inference[[first_subj]][[task]] <- task_events_data
                  current_cols <- names(task_events_data)
                  if (is.null(final_event_columns)) final_event_columns <- current_cols
                  
                  # Check required columns exist
                  req_col_names <- c(events_mapping$onset, events_mapping$duration, events_mapping$block)
                  missing_req <- setdiff(req_col_names, current_cols)
                  if(length(missing_req)>0) errors$add_error(paste0("events (",task,")"), paste0("Required column mapping(s) ('",paste(missing_req, collapse="', '"),"') point to columns not found in event file [subj ", first_subj, ", task ", task,"]. Available: ", paste(current_cols, collapse=", ")))

                  # Check mapped variable columns exist
                  var_map_names <- setdiff(names(events_mapping), c("onset", "duration", "block"))
                  for (var_name in var_map_names) {
                      col_to_check <- events_mapping[[var_name]]$column
                      if (!col_to_check %in% current_cols) {
                         errors$add_error(paste0("events$", var_name), paste0("Mapped column '", col_to_check, "' for variable '", var_name, "' not found in event file [subj ", first_subj, ", task ", task,"]. Available: ", paste(current_cols, collapse=", ")))
                      }
                      onset_override_col <- events_mapping[[var_name]]$onset
                      duration_override_col <- events_mapping[[var_name]]$duration
                      if (!is.null(onset_override_col) && !onset_override_col %in% current_cols) {
                           errors$add_error(paste0("events$", var_name), paste0("Onset override column '", onset_override_col, "' not found [subj ", first_subj, ", task ", task,"]."))
                      }
                      if (!is.null(duration_override_col) && !duration_override_col %in% current_cols) {
                           errors$add_error(paste0("events$", var_name), paste0("Duration override column '", duration_override_col, "' not found [subj ", first_subj, ", task ", task,"]."))
                      }
                  }
                  
                  # Validate block column type and order
                  block_col_name <- events_mapping$block
                  if (!is.null(block_col_name) && block_col_name %in% current_cols) {
                      block_vals <- task_events_data[[block_col_name]]
                      if (!is.numeric(block_vals) && !is.integer(block_vals)) {
                         if (is.factor(block_vals) || is.character(block_vals)) {
                             block_vals_num <- suppressWarnings(as.numeric(as.character(block_vals)))
                             if(any(is.na(block_vals_num))) {
                                  errors$add_error(paste0("events$block ('", block_col_name, "')"), paste0("Cannot convert block column to numeric [subj ", first_subj, ", task ", task,"]."))
                             } else {
                                 block_vals <- block_vals_num 
                             }
                         } else {
                             errors$add_error(paste0("events$block ('", block_col_name, "')"), paste0("Block column must be numeric/integer (or convertible) [subj ", first_subj, ", task ", task,"], got ", class(block_vals)[1], "."))
                         }
                      }
                      if ((is.numeric(block_vals) || is.integer(block_vals))) { 
                          if(any(is.na(block_vals))) {
                                errors$add_error(paste0("events$block ('", block_col_name, "')"), paste0("Block column contains NA values [subj ", first_subj, ", task ", task,"]."))
                          } else if (is.unsorted(block_vals, strictly = FALSE)) {
                              errors$add_error(paste0("events$block ('", block_col_name, "')"), paste0("Block values must be non-decreasing [subj ", first_subj, ", task ", task,"]."))
                          }
                      }
                  }
               } else {
                  # File not found for this task for the representative subject
                  errors$add_error(paste0("events (",task,")"), paste("No event file found or data is empty for representative subject ", first_subj, ", task ", task))
               }
            }, error = function(e) {
                errors$add_error(paste0("events (", task, ")"), paste("Failed reading/validating events [subj ", first_subj, ", task ", task, "]. Error: ", e$message))
            })
        } # End task loop
        
        if (!some_event_file_found && length(tasks) > 0) {
             # If we looped through all specified tasks and found no event files for the first subject
             errors$add_error("events", paste0("No event files could be found or read for any specified task for representative subject ", first_subj, ". Cannot proceed with event-related validation."))
        }
        
   } else {
      # Handle cases where no subjects or no tasks were selected/found
      if (length(subjects) == 0) warning("No subjects selected, skipping event file validation.")
      if (length(tasks) == 0 && length(available_tasks) > 0) { 
            # This means tasks were likely specified but none matched BIDS 
            # Error should have been added earlier. We just skip validation here.
      } else if (length(tasks) == 0 && length(available_tasks) == 0) {
            warning("No tasks found in BIDS, skipping event file validation.")
      }
   }

   events_info <- list(
       columns = final_event_columns %||% character(0), 
       mapping = events_mapping
   )

  # Validate Confound Columns (Requires loading data for one subject)
   confounds_info <- NULL
   selected_confound_cols <- character(0) 
   confounds_spec <- validated_ior$confounds %||% validated_ior$baseline$confounds
    if (is.null(confounds_spec) && !is.null(validated_ior$confounds)) {
       confounds_spec <- validated_ior$confounds
    }

   if (!is.null(confounds_spec) && length(subjects) > 0) {
       first_subj <- subjects[1]
       available_conf_cols <- NULL
       tryCatch({
           # Reading confounds might depend on runs, this is a simplified check
           confounds_read <- bidser::read_confounds(proj, subid = first_subj) 

           if (length(confounds_read) > 0 && length(confounds_read$data) > 0 && !is.null(confounds_read$data[[1]])) {
                confound_data <- confounds_read$data[[1]] # Assume first is representative
                available_conf_cols <- names(confound_data)
                included_cols <- character()
                if (!is.null(confounds_spec$include)) {
                   include_patterns <- unlist(confounds_spec$include)
                   for (pattern in include_patterns) {
                       matches <- grep(pattern, available_conf_cols, value = TRUE)
                       if (length(matches) == 0) {
                           errors$add_error("confounds$include", sprintf("Include pattern '%s' matched no variables [subj '%s']. Avail: %s", pattern, first_subj, paste(head(available_conf_cols), collapse=", ")))
                       }
                       included_cols <- c(included_cols, matches)
                   }
                   included_cols <- unique(included_cols)
                } else {
                   included_cols <- available_conf_cols
                }
                if (!is.null(confounds_spec$exclude)) {
                    exclude_patterns <- unlist(confounds_spec$exclude)
                    for (pattern in exclude_patterns) {
                       exclude_matches <- grep(pattern, included_cols, value = TRUE)
                       included_cols <- setdiff(included_cols, exclude_matches)
                    }
                }
                selected_confound_cols <- included_cols
                if (length(selected_confound_cols) == 0) {
                    warning(sprintf("Confound include/exclude patterns resulted in zero selected columns for subject '%s'.", first_subj))
                }
           } else {
                # Only error if confounds were specified but none found
                 errors$add_error("confounds", paste0("No confound data found for representative subject '", first_subj, "', but confounds were specified in YAML."))
           }
       }, error = function(e) {
           errors$add_error("confounds", paste0("Failed reading/validating confounds [subj ", first_subj, "]. Error: ", e$message))
       })

       confounds_info <- list(
           spec = confounds_spec,          
           columns = selected_confound_cols,
           available_columns = available_conf_cols # Store available columns for reference
       )
   }

  # --- 3. Infer and Validate Variable Roles ---
  variable_roles <- list(inferred_roles = list(), factors = character(0), parametric = character(0)) # Default
  # Check if we have the necessary components for inference
  can_infer <- length(subjects) > 0 && 
               length(tasks) > 0 && 
               !is.null(events_list_for_inference[[subjects[1]]]) && 
               length(events_list_for_inference[[subjects[1]]]) > 0 &&
               !is.null(validated_ior$regressors) &&
               length(validated_ior$regressors) > 0 &&
               !is.null(validated_ior$model$events)
               
  if (can_infer) {
       variable_roles <- infer_and_validate_variable_roles(
           regressors_spec = validated_ior$regressors,
           events_mapping = validated_ior$model$events, 
           events_list = events_list_for_inference, 
           user_factors = unlist(validated_ior$model$factors %||% list()), # Ensure vectors
           user_parametric = unlist(validated_ior$model$parametric %||% list()) 
       )
       
       # Post-inference check: Ensure columns for identified factors/parametrics exist
       all_role_vars <- c(variable_roles$factors, variable_roles$parametric)
       model_event_map_vars <- setdiff(names(validated_ior$model$events), c("onset", "duration", "block"))
       cols_for_roles <- character(0)
        for (model_var in intersect(all_role_vars, model_event_map_vars)) {
            cols_for_roles <- c(cols_for_roles, validated_ior$model$events[[model_var]])
        }
        cols_for_roles <- unique(cols_for_roles)
        missing_role_cols <- setdiff(cols_for_roles, events_info$columns)
        if (length(missing_role_cols) > 0) {
             errors$add_error("model$events", paste0("Column(s) used for factor/parametric roles ('", paste(missing_role_cols, collapse="', '"),"') not found in representative event file(s)."))
        }
        
  } else {
     warning("Skipping variable role inference due to missing subjects, tasks, event data, regressors, or model event mappings.")
     # Assign roles based only on user declaration if available
     variable_roles$factors = unique(unlist(validated_ior$model$factors %||% list()))
     variable_roles$parametric = unique(unlist(validated_ior$model$parametric %||% list()))
  }

  # --- 4. Final Checks & Stop if Invalid ---
  errors$stop_if_invalid("Context-dependent validation failed (BIDS content, event/confound columns, roles)")

  # --- 5. Construct Final fmri_config Object ---
  config <- structure(
    list(
      spec = validated_ior,      
      project = proj,            
      subjects = subjects,         
      tasks = tasks,             
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


###############################################################################
#             Step 3: Orchestrator Function (User-Facing API)
###############################################################################

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
#'   # Create a dummy config.yaml and BIDS structure for example
#'   bids_dir <- tempfile("bids_")
#'   dir.create(file.path(bids_dir, "sub-01", "func"), recursive = TRUE)
#'   writeLines(c("participant_id\tage", "sub-01\t25"), file.path(bids_dir, "participants.tsv"))
#'   cat("onset\tduration\tstim_type\tblock\n0\t1\tA\t1\n2\t1\tB\t1", 
#'       file = file.path(bids_dir, "sub-01", "func", "sub-01_task-test_events.tsv"))
#'   
#'   yaml_content <- "
#'   dataset:
#'     path: """ # Placeholder for path
#'   events:
#'     onset: onset
#'     duration: duration
#'     block: block
#'     stim:
#'       column: stim_type
#'   regressors:
#'     stim_A:
#'       type: hrf
#'       variables: [stim]
#'       subset: \"stim_type == 'A'\" # Ensure quotes are handled if needed
#'     stim_B:
#'       type: hrf
#'       variables: [stim]
#'       subset: \"stim_type == 'B'\" # Ensure quotes are handled if needed
#'   model:
#'     name: basic_model
#'     factors: [stim] # Declare stim (variable name) as factor
#'     events: # Model-specific mapping
#'       onset: onset 
#'       duration: duration
#'       block: block
#'       stim: stim_type # Map model variable 'stim' to column 'stim_type'
#'     regressors:
#'       stim_A: # Use regressor defined above
#'       stim_B:
#'   "
#'   # Inject the temp BIDS path into the YAML content
#'   # Use shQuote for safety, especially on Windows
#'   # Ensure backslashes are doubled on Windows if not using shQuote or normalizePath
#'   safe_bids_dir <- shQuote(normalizePath(bids_dir, winslash = "/"), type = "cmd")
#'   yaml_content <- sub('path: ""', paste0("path: ", safe_bids_dir), yaml_content)
#'   
#'   yaml_file <- tempfile(fileext = ".yaml")
#'   writeLines(yaml_content, yaml_file)
#'   message("Example YAML file created at: ", yaml_file)
#'   message("Example BIDS directory at: ", bids_dir)
#'   message("YAML content:\n", yaml_content)
#'   
#'   # Load the configuration
#'   # Wrap in tryCatch for cleanup even if loading fails
#'   config <- NULL
#'   tryCatch({
#'     config <- load_fmri_config(yaml_file)
#'     print(config)
#'   }, error = function(e) {
#'     message("Error loading config: ", e$message)
#'   }, finally = {
#'     # Clean up temp files
#'     unlink(bids_dir, recursive=TRUE, force=TRUE)
#'     unlink(yaml_file, force=TRUE)
#'     message("Cleaned up temporary files.")
#'   })
#' }
load_fmri_config <- function(yaml_file) {
  validated_ior <- parse_and_validate_config(yaml_file)
  fmri_config_obj <- build_config_from_ior(validated_ior)
  return(fmri_config_obj)
}


###############################################################################
#        Helper Function to Infer and Validate Variable Roles (Used by Step 2)
###############################################################################

#' @keywords internal
infer_and_validate_variable_roles <- function(regressors_spec, events_mapping, events_list,
                                              user_factors = NULL, user_parametric = NULL) {

  inferred_roles <- list()
  final_factors <- character(0)
  final_parametric <- character(0)
  user_factors <- user_factors %||% character(0)
  user_parametric <- user_parametric %||% character(0)

  all_regressor_model_vars <- unique(unlist(lapply(regressors_spec, `[[`, "variables")))

  if (length(all_regressor_model_vars) == 0) {
      # warning("No variables found referenced in any regressor specifications. Cannot infer roles.")
      return(list(inferred_roles = list(), factors = unique(user_factors), parametric = unique(user_parametric)))
  }

  representative_events_df <- NULL
  if (length(events_list) > 0 && length(events_list[[1]]) > 0) {
      # Find the first non-null data frame
      first_subj_tasks <- events_list[[1]]
      non_null_tasks <- Filter(Negate(is.null), first_subj_tasks)
       if (length(non_null_tasks) > 0) {
            representative_events_df <- non_null_tasks[[1]]
       } # else remains NULL
  }

  if (is.null(representative_events_df)) {
       warning("Could not find representative event data for role inference. Check loading errors.")
       return(list(inferred_roles = list(), factors = unique(user_factors), parametric = unique(user_parametric)))
  }

   model_var_to_column_map <- list()
   event_colnames <- names(representative_events_df)
   model_event_var_names <- setdiff(names(events_mapping), c("onset", "duration", "block"))
   for (model_var in model_event_var_names) {
       column_name <- events_mapping[[model_var]] 
       if (is.null(column_name)) {
           warning(sprintf("Model event mapping for '%s' is NULL. Skipping role inference.", model_var))
           next
       }
       if (!column_name %in% event_colnames) {
            warning(sprintf("Column '%s' (mapped from model var '%s') not found in representative event data. Skipping role inference.", column_name, model_var))
            next
        }
       model_var_to_column_map[[model_var]] <- column_name
   }

  for (model_var_name in all_regressor_model_vars) {
    col_name <- model_var_to_column_map[[model_var_name]]
    if (is.null(col_name)) {
         warning(sprintf("Variable '%s' used in regressor but not found in `model$events` mapping. Cannot infer role.", model_var_name))
         next
    }
    col_data <- representative_events_df[[col_name]]
    inferred_role <- "unknown" 
    if (is.factor(col_data)) {
      inferred_role <- "factor"
    } else if (is.character(col_data)) {
       unique_vals <- unique(stats::na.omit(col_data))
       num_unique <- length(unique_vals)
       vals_num <- suppressWarnings(as.numeric(unique_vals))
       is_potentially_numeric <- !any(is.na(vals_num))
       if (is_potentially_numeric) {
            inferred_role <- "parametric"
            # message(sprintf("Info: Var '%s' (col '%s') is character but looks numeric. Inferred 'parametric'.", model_var_name, col_name))
       } else if (num_unique <= max(5, 0.1 * length(stats::na.omit(col_data)))) { # Adjusted heuristic
             inferred_role <- "factor"
             # message(sprintf("Info: Var '%s' (col '%s') is character with few unique values (%d). Inferred 'factor'.", model_var_name, col_name, num_unique))
       } else {
            inferred_role <- "factor" 
            # warning(sprintf("Var '%s' (col '%s') is character with many unique values (%d). Defaulting to 'factor'.", model_var_name, col_name, num_unique))
       }
    } else if (is.numeric(col_data) || is.integer(col_data)) {
       unique_vals <- unique(stats::na.omit(col_data))
       num_unique <- length(unique_vals)
       is_integer_like <- all(unique_vals == floor(unique_vals), na.rm = TRUE) # Handle potential NAs from unique
       if (is_integer_like && num_unique < max(5, 0.05 * length(stats::na.omit(col_data)))) { 
           inferred_role <- "factor"
           # message(sprintf("Info: Var '%s' (col '%s') is numeric/integer with few unique values (%d). Inferred 'factor'.", model_var_name, col_name, num_unique))
       } else {
           inferred_role <- "parametric" 
       }
    } else if (is.logical(col_data)) {
        inferred_role <- "factor" 
    } else {
       warning(sprintf("Var '%s' (col '%s') has unsupported type (%s). Declare role explicitly.",
                       model_var_name, col_name, class(col_data)[1]))
       inferred_role <- "unknown"
    }
    inferred_roles[[model_var_name]] <- inferred_role
    user_declared_role <- NA
    if (model_var_name %in% user_factors) user_declared_role <- "factor"
    if (model_var_name %in% user_parametric) {
        if (!is.na(user_declared_role)) {
             warning(sprintf("Var '%s' declared as both factor and parametric. Prioritizing 'factor'.", model_var_name))
             user_declared_role <- "factor"
        } else {
            user_declared_role <- "parametric"
        }
    }
    final_role <- "unknown" 
    if (!is.na(user_declared_role)) {
        if (user_declared_role != inferred_role && inferred_role != "unknown") {
             warning(sprintf("Var '%s': Inferred '%s' conflicts with declared '%s'. Using declared role.",
                            model_var_name, inferred_role, user_declared_role))
        } 
        final_role <- user_declared_role 
    } else {
         if (inferred_role != "unknown") {
              # message(sprintf("Info: Var '%s' not declared. Using inferred role '%s'.", model_var_name, inferred_role))
              final_role <- inferred_role
         } else {
              warning(sprintf("Var '%s' used in regressors, role unknown/undeclared.", model_var_name))
              final_role <- "unknown"
         }
    }
    if (final_role == "factor") {
        final_factors <- c(final_factors, model_var_name)
    } else if (final_role == "parametric") {
        final_parametric <- c(final_parametric, model_var_name)
    }
  } 
  all_final_vars <- c(final_factors, final_parametric)
  unused_factors <- setdiff(user_factors, all_final_vars)
  unused_parametric <- setdiff(user_parametric, all_final_vars)
  if (length(unused_factors) > 0) {
      warning("Declared factors not used in model regressors: ", paste(unused_factors, collapse=", "))
  }
   if (length(unused_parametric) > 0) {
      warning("Declared parametric vars not used in model regressors: ", paste(unused_parametric, collapse=", "))
  }
   unknown_role_vars <- names(inferred_roles)[sapply(inferred_roles, `==`, "unknown")]
   unknown_role_vars <- intersect(unknown_role_vars, all_regressor_model_vars) 
   unknown_role_vars <- setdiff(unknown_role_vars, c(final_factors, final_parametric))
   if (length(unknown_role_vars) > 0) {
        warning("Could not determine role for regressor variables: ", paste(unknown_role_vars, collapse=", "), ". Errors likely during model construction.")
   }
  list(
    inferred_roles = inferred_roles,
    factors = unique(final_factors),
    parametric = unique(final_parametric)
  )
}


###############################################################################
#                   S3 Print Method for fmri_config (Updated)
###############################################################################

#' @export
#' @importFrom crayon bold green blue cyan red italic yellow
print.fmri_config <- function(x, ...) {
    format_list <- function(lst, indent = "  ", none_msg = crayon::yellow("none")) {
        if (is.null(lst) || length(lst) == 0) return(none_msg)
        if (!is.null(names(lst))) {
             paste(paste0(indent, "- ", crayon::italic(names(lst)), ": ", lst), collapse = "\n")
        } else { 
             paste(paste0(indent, "- ", lst), collapse = "\n")
        }
    }
    section_header <- function(title) paste0("\n", crayon::bold(crayon::blue("=== ", title, " ===")))
    subsection_header <- function(title) paste0(crayon::bold(crayon::cyan("--- ", title, " ---")))
    cat(crayon::bold(crayon::green("fMRI Analysis Configuration\n")))
    cat(section_header("Dataset"))
    cat("\nBIDS Project Path:", x$project$path)
    cat("\nRelative Path (func):", x$spec$dataset$relpath)
    cat("\nSelected Subjects:", length(x$subjects), "total")
    cat("\n", format_list(x$subjects))
    cat("\nSelected Tasks:")
    cat("\n", format_list(x$tasks))
     if (!is.null(x$spec$dataset$runs)) {
         cat("\nSelected Runs:")
         cat("\n", format_list(unlist(x$spec$dataset$runs)))
     }
    if (!is.null(x$spec$dataset$scan_params)) {
       cat("\n", subsection_header("Scan Parameters"))
       sp <- x$spec$dataset$scan_params
       if(!is.null(sp$TR)) cat("\n  TR Default:", sp$TR$default %||% crayon::yellow("N/A"), " (Overrides:", length(sp$TR$overrides %||% list()), ")")
       if(!is.null(sp$run_length)) cat("\n  Run Length Defaults:", length(sp$run_length$default %||% list()), "tasks specified", " (Overrides:", length(sp$run_length$overrides %||% list()), ")")
    }
    cat(section_header("Events"))
    cat("\n", subsection_header("Top-Level Event Column Mappings"))
    ev_map <- x$events_info$mapping
    cat(sprintf("\n  Onset Column: '%s'", ev_map$onset %||% crayon::red("Not Specified")))
    cat(sprintf("\n  Duration Column: '%s'", ev_map$duration %||% crayon::red("Not Specified")))
    cat(sprintf("\n  Block Column: '%s'", ev_map$block %||% crayon::red("Not Specified")))
    var_maps <- setdiff(names(ev_map), c("onset", "duration", "block"))
    if (length(var_maps) > 0) {
        cat("\n", subsection_header("Variable Mappings (Top Level)"))
        for (vm in var_maps) {
             map_info <- ev_map[[vm]]
             cat(sprintf("\n  %s -> column: '%s'%s%s",
                        crayon::italic(vm),
                        map_info$column %||% crayon::red("N/A"),
                        if(!is.null(map_info$onset)) sprintf(" (onset: '%s')", map_info$onset) else "",
                        if(!is.null(map_info$duration)) sprintf(" (duration: '%s')", map_info$duration) else "" ))
        }
    }
     cat("\nValidated Event Columns (Rep. File):", if(length(x$events_info$columns)>0) paste(x$events_info$columns, collapse=", ") else crayon::yellow("Not validated/found"))
    if (!is.null(x$confounds_info)) {
        cat(section_header("Confounds"))
        conf_spec <- x$confounds_info$spec
        if (!is.null(conf_spec$include)) {
          cat("\n", subsection_header("Include Patterns"))
          cat("\n", format_list(unlist(conf_spec$include)))
        }
        if (!is.null(conf_spec$exclude)) {
          cat("\n", subsection_header("Exclude Patterns"))
          cat("\n", format_list(unlist(conf_spec$exclude)))
        }
        cat("\nSelected Columns (Rep. File):")
        cat("\n", format_list(x$confounds_info$columns, none_msg = crayon::yellow("None selected or validated")))
        # Optionally show available confound columns
        # cat("\nAvailable Columns (Rep. File):")
        # cat("\n", format_list(x$confounds_info$available_columns, none_msg = crayon::yellow("Not available")))
    }
    cat(section_header(paste("Model:", x$spec$model$name %||% crayon::red("Unnamed"))))
    cat("\n", subsection_header("Model Variable Roles"))
    cat("\nFactors:", format_list(x$variable_roles$factors))
    cat("\nParametric:", format_list(x$variable_roles$parametric))
     cat("\n", subsection_header("Model-Specific Event Variable -> Column Maps"))
     mod_ev_map <- x$spec$model$events
     cat(sprintf("\n  Onset Column (model): '%s'", mod_ev_map$onset %||% crayon::yellow("Using top-level")))
     cat(sprintf("\n  Duration Column (model): '%s'", mod_ev_map$duration %||% crayon::yellow("Using top-level")))
     cat(sprintf("\n  Block Column (model): '%s'", mod_ev_map$block %||% crayon::yellow("Using top-level")))
     mod_var_maps <- setdiff(names(mod_ev_map), c("onset", "duration", "block"))
     if (length(mod_var_maps) > 0) {
         for (vm in mod_var_maps) {
              cat(sprintf("\n  %s -> column: '%s'", crayon::italic(vm), mod_ev_map[[vm]] %||% crayon::red("N/A")))
         }
     } else { cat("\n  (No model-specific variable mappings)") }
    model_reg_names <- names(x$spec$model$regressors %||% list())
    hrfs_used_in_model <- unique(unlist(lapply(model_reg_names, function(name) x$spec$regressors[[name]]$hrf)))
    if (length(hrfs_used_in_model) > 0) {
        cat("\n", subsection_header("HRF Specifications (Used in Model)"))
        for (hrf_name in hrfs_used_in_model) {
           hrf_spec <- x$spec$hrfs[[hrf_name]]
            if (!is.null(hrf_spec)) {
                 cat(sprintf("\n  %s: type=%s", crayon::italic(hrf_name), hrf_spec$type %||% "N/A"))
                 if (!is.null(hrf_spec$parameters)) {
                    params <- paste(names(hrf_spec$parameters), hrf_spec$parameters, sep = "=", collapse = ", ")
                    cat(sprintf(" (params: %s)", params))
                 }
                  if (!is.null(hrf_spec$type) && hrf_spec$type == "custom") cat(sprintf(" (definition: %s)", hrf_spec$definition %||% "N/A"))
            } else {
                cat(sprintf("\n  %s: %s", crayon::italic(hrf_name), crayon::yellow("Definition not found in 'hrfs' (using built-in?)")))
            }
        }
    }
     cat("\n", subsection_header("Regressors (Used in Model)"))
     if (length(model_reg_names) > 0) {
        for (reg_name in model_reg_names) {
            reg_spec <- x$spec$regressors[[reg_name]]
            cat(sprintf("\n  %s: type=%s, vars=[%s]", crayon::italic(reg_name), reg_spec$type %||% "N/A", paste(unlist(reg_spec$variables), collapse=", ")))
            if (!is.null(reg_spec$hrf)) cat(sprintf(", hrf=%s", reg_spec$hrf))
            if (!is.null(reg_spec$basis)) cat(sprintf(", basis=%s", reg_spec$basis$type))
            if (!is.null(reg_spec$subset)) cat(crayon::yellow(" (subset)"))
        }
     } else { cat("\n", format_list(NULL)) }
    model_con_names <- names(x$spec$model$contrasts %||% list())
    if (length(model_con_names) > 0) {
        cat("\n", subsection_header("Contrasts (Used in Model)"))
        for (con_name in model_con_names) {
            con_spec <- x$spec$contrasts[[con_name]]
             cat(sprintf("\n  %s: type=%s", crayon::italic(con_name), con_spec$type %||% "N/A"))
              if (!is.null(con_spec$expression)) cat(sprintf(", expr='%s'", con_spec$expression))
              if (!is.null(con_spec$factors)) cat(sprintf(", factors=[%s]", paste(unlist(con_spec$factors), collapse=", ")))
               if (!is.null(con_spec$where)) cat(crayon::yellow(" (where)"))
        }
    }
    cat(section_header("Build Status"))
    cat("\nStatus:", if (x$validated) crayon::green("✓ Successfully Built") else crayon::red("✗ Build Failed or Incomplete"))
    cat("\n")
    invisible(x)
}

# Note: Older standalone validation/utility functions like validate_subjects,
# validate_tasks, validate_events, validate_confounds, get_scan_tr,
# get_scan_length are removed as their logic is integrated or superseded.
# Keeping %||%, infer_and_validate_variable_roles, print.fmri_config,
# load_fmri_config, parse_and_validate_config, build_config_from_ior,
# and the validation helper functions. 