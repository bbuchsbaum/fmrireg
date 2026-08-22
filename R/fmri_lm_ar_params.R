#' Extract Estimated AR Parameters from fmri_lm Fit
#'
#' Retrieves the estimated autoregressive parameters from a fitted
#' fMRI linear model that used AR error modeling.
#'
#' @param object An object of class \code{fmri_lm}
#' @param scope Character; \code{"average"} (default) returns the pooled
#'   average AR coefficients, \code{"per_run"} returns a list of the run-level
#'   estimates, and \code{"raw"} returns the stored structure without
#'   post-processing.
#' @param ... Additional arguments (currently unused)
#'
#' @return Depending on \code{scope}, either a numeric vector of averaged AR
#'   coefficients, a list of per-run coefficient vectors, or the raw stored
#'   structure. Returns \code{NULL} when no AR modeling was performed.
#'
#' @export
#' @examples
#' \dontrun{
#' # Fit model with AR(1) errors
#' fit <- fmri_lm(onset ~ hrf(cond), dataset = dset, cor_struct = "ar1")
#' ar_parameters(fit)  # Extract estimated AR(1) coefficient
#' }
ar_parameters <- function(object, ...) {
  UseMethod("ar_parameters")
}

#' @rdname ar_parameters
#' @export
ar_parameters.fmri_lm <- function(object, scope = c("average", "per_run", "raw"), ...) {
  scope <- match.arg(scope)

  # Pull stored coefficients from the fit
  ar_data <- object$ar_coef %||% object$result$ar_coef %||% attr(object, "ar_coef")

  if (is.null(ar_data)) {
    cfg <- attr(object, "config")
    if (is.null(cfg) || cfg$ar$struct %in% c("iid", "none")) {
      return(NULL)
    }
    return(NULL)
  }

  if (scope == "raw") {
    return(ar_data)
  }

  phi_vecs <- .collect_ar_vectors(ar_data)
  phi_vecs <- Filter(function(x) length(x) > 0, phi_vecs)

  if (!length(phi_vecs)) {
    return(NULL)
  }

  if (scope == "per_run") {
    return(phi_vecs)
  }

  max_len <- max(lengths(phi_vecs))
  phi_mat <- vapply(phi_vecs, function(phi) {
    c(phi, rep(NA_real_, max_len - length(phi)))
  }, numeric(max_len))

  if (is.null(dim(phi_mat))) {
    phi_mat <- matrix(phi_mat, nrow = max_len, ncol = length(phi_vecs))
  }

  rowMeans(phi_mat, na.rm = TRUE)
}

#' @keywords internal
#' @noRd
.collect_ar_vectors <- function(x) {
  out <- list()
  idx <- 0L

  append_vec <- function(v) {
    idx <<- idx + 1L
    out[[idx]] <<- as.numeric(v)
  }

  recurse <- function(item) {
    if (is.null(item)) {
      return()
    }
    if (is.numeric(item)) {
      append_vec(item)
    } else if (is.list(item)) {
      lapply(item, recurse)
    }
  }

  recurse(x)
  out
}

#' Extract Estimated MA Parameters from fmri_lm Fit
#'
#' Retrieves the estimated moving-average parameters from a fitted fMRI linear
#' model that used ARMA or MA temporal-noise modeling.
#'
#' @param object An object of class \code{fmri_lm}.
#' @param scope Character; \code{"average"} (default) returns coefficients
#'   averaged across runs, \code{"per_run"} returns the run-level estimates,
#'   and \code{"raw"} returns the stored structure without post-processing.
#' @param ... Additional arguments (currently unused).
#'
#' @return Depending on \code{scope}, either a numeric vector of averaged MA
#'   coefficients, a list of per-run coefficient vectors, or the raw stored
#'   structure. Returns \code{NULL} when no MA modeling was performed.
#' @export
ma_parameters <- function(object, ...) {
  UseMethod("ma_parameters")
}

#' @rdname ma_parameters
#' @export
ma_parameters.fmri_lm <- function(object,
                                  scope = c("average", "per_run", "raw"), ...) {
  scope <- match.arg(scope)
  ma_data <- object$ma_coef %||% object$result$ma_coef %||% attr(object, "ma_coef")

  if (is.null(ma_data)) {
    return(NULL)
  }
  if (scope == "raw") {
    return(ma_data)
  }

  theta_vecs <- .collect_ar_vectors(ma_data)
  theta_vecs <- Filter(function(x) length(x) > 0L, theta_vecs)
  if (!length(theta_vecs)) {
    return(NULL)
  }
  if (scope == "per_run") {
    return(theta_vecs)
  }

  max_len <- max(lengths(theta_vecs))
  theta_mat <- vapply(theta_vecs, function(theta) {
    c(theta, rep(NA_real_, max_len - length(theta)))
  }, numeric(max_len))
  if (is.null(dim(theta_mat))) {
    theta_mat <- matrix(theta_mat, nrow = max_len, ncol = length(theta_vecs))
  }
  rowMeans(theta_mat, na.rm = TRUE)
}
