#' Schema Validation Helper Functions
#'
#' These internal helpers validate fields in a configuration list
#' against simplified JSON-Schema like constraints. They are used
#' by the DSL parser to perform low level schema checks.
#'
#' @keywords internal
#' @noRd
check_required <- function(data, field, path, errors) {
  if (!is.list(data) || !exists(field, where = data, inherits = FALSE)) {
    field_path <- paste0(path, if (nzchar(path)) "$", field)
    errors$add_error(field_path,
                     sprintf("Required field '%s' is missing.", field))
    return(FALSE)
  }
  TRUE
}

#' @keywords internal
#' @noRd
check_type <- function(data, field, expected_type, path, errors,
                       allow_null = FALSE) {
  if (!exists(field, where = data, inherits = FALSE)) {
    return(allow_null)
  }
  value <- data[[field]]
  field_path <- paste0(path, if (nzchar(path)) "$", field)

  if (is.null(value)) {
    if (!allow_null) {
      errors$add_error(field_path, "Field cannot be null.")
      return(FALSE)
    }
    return(TRUE)
  }

  valid <- FALSE
  actual_type <- paste(class(value), collapse = "/")

  if (expected_type == "string") {
    valid <- is.character(value) && length(value) == 1 && !is.na(value)
  } else if (expected_type == "integer") {
    valid <- (is.integer(value) ||
               (is.numeric(value) && !is.na(value) && value == floor(value))) &&
             length(value) == 1
  } else if (expected_type == "number") {
    valid <- is.numeric(value) && length(value) == 1 && !is.na(value)
  } else if (expected_type == "boolean") {
    valid <- is.logical(value) && length(value) == 1 && !is.na(value)
  } else if (expected_type == "object") {
    valid <- is.list(value) &&
      (length(value) == 0 || (!is.null(names(value)) && all(names(value) != "")))
  } else if (grepl("^array", expected_type)) {
    is_array_like <- is.list(value) || (is.vector(value) && !is.matrix(value))
    if (is_array_like) {
      valid <- TRUE
      m <- regexec("^array\\[(\\w+)\\]$", expected_type)
      mt <- regmatches(expected_type, m)
      element_type <- if (length(mt[[1]]) > 1) mt[[1]][2] else NULL
      if (length(value) > 0 && !is.null(element_type)) {
        check_el <- function(el) {
          if (element_type == "string")
            is.character(el) && length(el) == 1 && !is.na(el)
          else if (element_type == "integer")
            (is.integer(el) || (is.numeric(el) && !is.na(el) && el == floor(el))) &&
              length(el) == 1
          else if (element_type == "number")
            is.numeric(el) && length(el) == 1 && !is.na(el)
          else if (element_type == "boolean")
            is.logical(el) && length(el) == 1 && !is.na(el)
          else if (element_type == "object")
            is.list(el) &&
              (length(el) == 0 || (!is.null(names(el)) && all(names(el) != "")))
          else FALSE
        }
        ok <- tryCatch(vapply(value, check_el, logical(1)),
                        error = function(e) FALSE)
        if (!all(ok)) {
          valid <- FALSE
          actual_type <- paste0(actual_type,
                                " (contains elements of incorrect type)")
        }
      }
    } else {
      valid <- FALSE
    }
  } else {
    errors$add_error(field_path,
                     sprintf("Internal validation error: Unknown expected type '%s' specified.",
                             expected_type))
    return(FALSE)
  }

  if (!valid) {
    errors$add_error(field_path,
                     sprintf("Invalid type. Expected '%s', but got '%s'.",
                             expected_type, actual_type))
  }
  valid
}

#' @keywords internal
#' @noRd
check_enum <- function(data, field, allowed_values, path, errors,
                       allow_null = FALSE) {
  if (!exists(field, where = data, inherits = FALSE)) {
    return(allow_null)
  }
  value <- data[[field]]
  field_path <- paste0(path, if (nzchar(path)) "$", field)
  if (is.null(value)) {
    return(allow_null)
  }

  values <- unlist(value)
  if (length(values) == 0 && !is.null(value)) return(TRUE)
  if (is.null(values) && length(value) > 0) {
    errors$add_error(field_path,
                     "Contains NULL values where specific choices are expected.")
    return(FALSE)
  }
  bad <- setdiff(values, allowed_values)
  if (length(bad) > 0) {
    errors$add_error(field_path,
                     sprintf("Invalid value(s): '%s'. Must be one of: %s",
                             paste(unique(bad), collapse = "', '"),
                             paste(allowed_values, collapse = ", ")))
    return(FALSE)
  }
  TRUE
}

#' @keywords internal
#' @noRd
check_pattern <- function(data, field, pattern, path, errors,
                          allow_null = FALSE) {
  if (!exists(field, where = data, inherits = FALSE)) {
    return(allow_null)
  }
  value <- data[[field]]
  field_path <- paste0(path, if (nzchar(path)) "$", field)
  if (is.null(value)) {
    return(allow_null)
  }

  values <- unlist(value)
  if (length(values) == 0 && !is.null(value)) return(TRUE)
  if (is.null(values) && length(value) > 0) {
    errors$add_error(field_path,
                     "Pattern check applied to list containing NULL values.")
    return(FALSE)
  }
  if (!is.character(values)) {
    return(TRUE)
  }
  bad <- !grepl(pattern, values)
  if (any(bad)) {
    errors$add_error(field_path,
                     sprintf("Value(s) do not match required pattern ('%s'): '%s'",
                             pattern,
                             paste(unique(values[bad]), collapse = "', '")))
    return(FALSE)
  }
  TRUE
}
