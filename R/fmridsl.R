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
