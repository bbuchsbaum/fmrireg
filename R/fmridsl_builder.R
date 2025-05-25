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

  config <- list(spec = validated_ior, project = project, validated = TRUE)
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
