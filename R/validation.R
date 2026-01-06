#' Input validation utilities for fmrireg
#'
#' Centralized validation functions that provide consistent, informative
#' error messages across the package. These functions check conditions and
#' abort with structured errors if validation fails.
#'
#' @name validation
#' @keywords internal
NULL

#' Check that an object inherits from specified class(es)
#'
#' @param x Object to check
#' @param class Character vector of class names (any match is OK)
#' @param arg Argument name for error messages
#' @param call The execution environment
#' @return `x` invisibly if valid
#' @keywords internal
check_inherits <- function(x, class, arg = NULL, call = rlang::caller_env()) {
  if (is.null(arg)) {
    arg <- deparse(substitute(x))
  }
  if (!inherits(x, class)) {
    class_str <- paste(class, collapse = " or ")
    fmrireg_abort_input(
      arg = arg,
      must = paste0("a {.cls ", class_str, "} object"),
      not = class(x)[1],
      call = call
    )
  }
  invisible(x)
}

#' Check that a value is numeric
#'
#' @inheritParams check_inherits
#' @return `x` invisibly if valid
#' @keywords internal
check_numeric <- function(x, arg = NULL, call = rlang::caller_env()) {
  if (is.null(arg)) {
    arg <- deparse(substitute(x))
  }
  if (!is.numeric(x)) {
    fmrireg_abort_input(
      arg = arg,
      must = "numeric",
      not = class(x)[1],
      call = call
    )
  }
  invisible(x)
}

#' Check that a value is a positive number
#'
#' @inheritParams check_inherits
#' @param allow_zero Logical; if TRUE, zero is allowed
#' @return `x` invisibly if valid
#' @keywords internal
check_positive <- function(x, arg = NULL, allow_zero = FALSE,
                           call = rlang::caller_env()) {
  if (is.null(arg)) {
    arg <- deparse(substitute(x))
  }
  check_numeric(x, arg, call)

  condition <- if (allow_zero) any(x < 0) else any(x <= 0)
  if (condition) {
    qualifier <- if (allow_zero) "non-negative" else "positive"
    fmrireg_abort_input(
      arg = arg,
      must = qualifier,
      call = call
    )
  }
  invisible(x)
}

#' Check that a value is a logical scalar
#'
#' @inheritParams check_inherits
#' @return `x` invisibly if valid
#' @keywords internal
check_logical_scalar <- function(x, arg = NULL, call = rlang::caller_env()) {
  if (is.null(arg)) {
    arg <- deparse(substitute(x))
  }
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    fmrireg_abort_input(
      arg = arg,
      must = "a single TRUE or FALSE value",
      not = if (!is.logical(x)) class(x)[1] else NULL,
      call = call
    )
  }
  invisible(x)
}

#' Check that a value is one of allowed choices
#'
#' @inheritParams check_inherits
#' @param choices Character vector of allowed values
#' @return `x` invisibly if valid
#' @keywords internal
check_one_of <- function(x, choices, arg = NULL, call = rlang::caller_env()) {
  if (is.null(arg)) {
    arg <- deparse(substitute(x))
  }
  if (!x %in% choices) {
    choices_str <- paste0("{.val ", choices, "}", collapse = ", ")
    fmrireg_abort(
      c(
        paste0("{.arg ", arg, "} must be one of: ", choices_str),
        "x" = "Got {.val {x}}"
      ),
      class = "fmrireg_error_input",
      call = call
    )
  }
  invisible(x)
}

#' Check that dimensions match
#'
#' @param x First matrix/vector
#' @param y Second matrix/vector
#' @param x_dim Which dimension of x to check ("rows" or "cols")
#' @param y_dim Which dimension of y to check ("rows" or "cols")
#' @param context Optional context for error message
#' @param call The execution environment
#' @return TRUE invisibly if valid
#' @keywords internal
check_dimensions_match <- function(x, y, x_dim = "rows", y_dim = "rows",
                                   context = NULL, call = rlang::caller_env()) {
  nx <- if (x_dim == "rows") NROW(x) else NCOL(x)
  ny <- if (y_dim == "rows") NROW(y) else NCOL(y)

  if (nx != ny) {
    x_name <- paste(deparse(substitute(x)), x_dim)
    y_name <- paste(deparse(substitute(y)), y_dim)
    fmrireg_abort_dimension(
      expected = nx,
      actual = ny,
      context = context,
      expected_name = x_name,
      actual_name = y_name,
      call = call
    )
  }
  invisible(TRUE)
}

#' Check that a file exists
#'
#' @param path File path to check
#' @param arg Optional argument name for error messages
#' @param call The execution environment
#' @return `path` invisibly if valid
#' @keywords internal
check_file_exists <- function(path, arg = NULL, call = rlang::caller_env()) {
  if (!file.exists(path)) {
    fmrireg_abort_file(path, arg, call)
  }
  invisible(path)
}

#' Check that a value is a formula
#'
#' @inheritParams check_inherits
#' @return `x` invisibly if valid
#' @keywords internal
check_formula <- function(x, arg = NULL, call = rlang::caller_env()) {
  if (is.null(arg)) {
    arg <- deparse(substitute(x))
  }
  if (!inherits(x, "formula")) {
    fmrireg_abort_input(
      arg = arg,
      must = "a formula",
      not = class(x)[1],
      call = call
    )
  }
  invisible(x)
}

#' Check that a matrix is not rank deficient
#'
#' @param x Matrix to check
#' @param tolerance Tolerance for rank determination
#' @param context Optional context for error/warning
#' @param warn_only If TRUE, warn instead of error
#' @param call The execution environment
#' @return Rank of the matrix invisibly
#' @keywords internal
check_matrix_rank <- function(x, tolerance = .Machine$double.eps * 100,
                              context = NULL, warn_only = FALSE,
                              call = rlang::caller_env()) {
  rank <- Matrix::rankMatrix(x, tol = tolerance)
  expected <- min(nrow(x), ncol(x))

  if (rank < expected) {
    context_msg <- if (!is.null(context)) paste0(" ", context) else ""
    msg <- c(
      paste0("Design matrix is rank deficient", context_msg),
      "i" = "Rank: {.val {rank}} / {.val {expected}}",
      ">" = "Consider removing collinear regressors"
    )

    if (warn_only) {
      fmrireg_warn(msg, call = call)
    } else {
      fmrireg_abort(msg, class = "fmrireg_error_rank", call = call)
    }
  }
  invisible(as.integer(rank))
}

#' Check for NA values in numeric data
#'
#' @param x Numeric vector or matrix
#' @param arg Argument name for error messages
#' @param context Optional context for error message
#' @param call The execution environment
#' @return `x` invisibly if valid
#' @keywords internal
check_no_na <- function(x, arg = NULL, context = NULL,
                        call = rlang::caller_env()) {
  if (is.null(arg)) {
    arg <- deparse(substitute(x))
  }
  if (anyNA(x)) {
    context_msg <- if (!is.null(context)) paste0(" in ", context) else ""
    fmrireg_abort(
      c(
        paste0("NA values detected", context_msg),
        "x" = "{.arg {arg}} contains {.val {sum(is.na(x))}} NA value(
s)"
      ),
      class = "fmrireg_error_na",
      call = call
    )
  }
  invisible(x)
}

#' Check that a value has expected length
#'
#' @inheritParams check_inherits
#' @param expected Expected length
#' @return `x` invisibly if valid
#' @keywords internal
check_length <- function(x, expected, arg = NULL, call = rlang::caller_env()) {
  if (is.null(arg)) {
    arg <- deparse(substitute(x))
  }
  if (length(x) != expected) {
    fmrireg_abort(
      c(
        "{.arg {arg}} must have length {.val {expected}}",
        "x" = "Got length {.val {length(x)}}"
      ),
      class = "fmrireg_error_input",
      call = call
    )
  }
  invisible(x)
}
