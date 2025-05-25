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
