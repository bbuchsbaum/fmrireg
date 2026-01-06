#' Error handling utilities for fmrireg
#'
#' This module provides structured error reporting with informative messages
#' and catchable error classes. All errors inherit from `fmrireg_error`.
#'
#' @name fmrireg_errors
#' @keywords internal
NULL

#' Abort with a fmrireg error
#'
#' @param message Character vector of error messages. First element is the
#'   main message, subsequent elements are bullet points.
#' @param class Additional error classes (prepended to "fmrireg_error")
#' @param ... Additional arguments passed to [cli::cli_abort()]
#' @param call The execution environment of the error
#' @return Never returns; always throws an error
#' @keywords internal
fmrireg_abort <- function(message, class = NULL, ..., call = rlang::caller_env()) {

cli::cli_abort(
    message,
    class = c(class, "fmrireg_error"),
    ...,
    call = call
  )
}

#' Signal a fmrireg warning
#'
#' @inheritParams fmrireg_abort
#' @return NULL invisibly
#' @keywords internal
fmrireg_warn <- function(message, class = NULL, ..., call = rlang::caller_env()) {
  cli::cli_warn(
    message,
    class = c(class, "fmrireg_warning"),
    ...,
    call = call
  )
}

#' Abort due to dimension mismatch
#'
#' @param expected Expected dimension value
#' @param actual Actual dimension value
#' @param context Optional context string describing where mismatch occurred
#' @param expected_name Name for expected dimension (default "expected")
#' @param actual_name Name for actual dimension (default "got")
#' @param call The execution environment
#' @return Never returns; always throws an error
#' @keywords internal
fmrireg_abort_dimension <- function(expected, actual, context = NULL,
                                     expected_name = "expected",
                                     actual_name = "got",
                                     call = rlang::caller_env()) {
  context_msg <- if (!is.null(context)) paste0(" in ", context) else ""
  fmrireg_abort(
    c(
      paste0("Dimension mismatch", context_msg),
      "x" = paste0(expected_name, ": {.val {expected}}"),
      "x" = paste0(actual_name, ": {.val {actual}}")
    ),
    class = "fmrireg_error_dimension",
    call = call
  )
}

#' Abort due to invalid input
#'
#' @param arg Argument name that has invalid input
#' @param must Description of what the argument must be
#' @param not Optional description of what was actually provided
#' @param call The execution environment
#' @return Never returns; always throws an error
#' @keywords internal
fmrireg_abort_input <- function(arg, must, not = NULL, call = rlang::caller_env()) {
  msg <- c("{.arg {arg}} must be {must}")
  if (!is.null(not)) {
    msg <- c(msg, "x" = "Got {.cls {not}}")
  }

  fmrireg_abort(msg, class = "fmrireg_error_input", call = call)
}
#' Abort due to configuration error
#'
#' @param param Configuration parameter name
#' @param message Description of the configuration problem
#' @param suggestion Optional suggestion for fixing the problem
#' @param call The execution environment
#' @return Never returns; always throws an error
#' @keywords internal
fmrireg_abort_config <- function(param, message, suggestion = NULL,
                                  call = rlang::caller_env()) {
  msg <- c(
    "Configuration error for {.arg {param}}",
    "x" = message
)
  if (!is.null(suggestion)) {
    msg <- c(msg, ">" = suggestion)
  }
  fmrireg_abort(msg, class = "fmrireg_error_config", call = call)
}

#' Abort due to missing file
#'
#' @param path File path that was not found
#' @param arg Optional argument name that contained the path
#' @param call The execution environment
#' @return Never returns; always throws an error
#' @keywords internal
fmrireg_abort_file <- function(path, arg = NULL, call = rlang::caller_env()) {
  msg <- "File not found"
  if (!is.null(arg)) {
    msg <- paste0("File not found for {.arg {arg}}")
  }
  fmrireg_abort(
    c(msg, "x" = "Path: {.file {path}}"),
    class = "fmrireg_error_file",
    call = call
  )
}

#' Abort because feature is not implemented
#'
#' @param feature Description of the unimplemented feature
#' @param call The execution environment
#' @return Never returns; always throws an error
#' @keywords internal
fmrireg_abort_not_implemented <- function(feature, call = rlang::caller_env()) {
  fmrireg_abort(
    c(
      "Feature not implemented",
      "i" = feature
    ),
    class = "fmrireg_error_not_implemented",
    call = call
  )
}

#' Abort due to computation failure
#'
#' @param operation Description of the operation that failed
#' @param reason Reason for the failure
#' @param suggestion Optional suggestion for resolving the issue
#' @param call The execution environment
#' @return Never returns; always throws an error
#' @keywords internal
fmrireg_abort_computation <- function(operation, reason, suggestion = NULL,
                                       call = rlang::caller_env()) {
  msg <- c(
    paste0("Failed to ", operation),
    "x" = reason
  )
  if (!is.null(suggestion)) {
    msg <- c(msg, ">" = suggestion)
  }
  fmrireg_abort(msg, class = "fmrireg_error_computation", call = call)
}
