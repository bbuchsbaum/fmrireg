#' Configuration for fmri_lm fitting
#'
#' `fmri_lm_control()` creates an `fmri_lm_config` object collecting all
#' options for robust and autoregressive modelling. It validates inputs and
#' applies defaults so downstream functions receive a single structured list.
#'
#' For common use cases, `fmri_lm()` provides convenience parameters that
#' are easier to use than these detailed option lists:
#' \itemize{
#'   \item \code{volume_weights = TRUE} enables volume weighting with defaults
#'   \item \code{volume_weights = "tukey"} enables with Tukey method
#'   \item \code{nuisance_projection = N} enables soft projection with matrix N
#'   \item \code{nuisance_projection = "mask.nii"} enables with mask file
#' }
#'
#' Use the `*_options` lists below only when you need fine-grained control.
#'
#' @param robust_options list of robust fitting options. See Details.
#' @param ar_options list of autoregressive modelling options. See Details.
#' @param volume_weights_options list of volume weighting options. See Details.
#'   For simple cases, use the \code{volume_weights} parameter in \code{fmri_lm()} instead.
#' @param soft_subspace_options list of soft subspace projection options. See Details.
#'   For simple cases, use the \code{nuisance_projection} parameter in \code{fmri_lm()} instead.
#' @return An object of class `fmri_lm_config`.
#' @details
#' `robust_options` may contain:
#'   * `type` (`FALSE`, "huber", "bisquare")
#'   * `k_huber`
#'   * `c_tukey`
#'   * `max_iter`
#'   * `scale_scope` ("run", "global", "voxel")
#'   * `reestimate_phi` (logical)
#'
#' `ar_options` may contain:
#'   * `struct` ("iid", "ar1", "ar2", "arp")
#'   * `p` (order for "arp")
#'   * `iter_gls` (integer number of GLS iterations)
#'   * `global` (logical, use global phi)
#'   * `voxelwise` (logical)
#'   * `exact_first` (logical)
#'   * `censor` (integer vector of timepoints to exclude from AR estimation,
#'     logical vector where TRUE = censored, or "auto" to extract from dataset)
#'
#' `volume_weights_options` may contain:
#'   * `enabled` (logical, whether to compute and apply volume weights)
#'   * `method` ("inverse_squared", "soft_threshold", "tukey")
#'   * `threshold` (numeric, DVARS threshold for weighting)
#'   * `weights` (optional pre-computed weight vector)
#'
#' `soft_subspace_options` may contain:
#'   * `enabled` (logical, whether to apply soft subspace projection)
#'   * `nuisance_mask` (path to NIfTI mask or logical vector)
#'   * `nuisance_matrix` (pre-computed nuisance timeseries matrix)
#'   * `lambda` (numeric, "auto", or "gcv")
#'   * `warn_redundant` (logical, warn if baseline has nuisance terms)
#'
#' @export
fmri_lm_control <- function(robust_options = list(),
                            ar_options = list(),
                            volume_weights_options = list(),
                            soft_subspace_options = list()) {
  # Handle NULL inputs
  if (is.null(robust_options)) robust_options <- list()
  if (is.null(ar_options)) ar_options <- list()
  if (is.null(volume_weights_options)) volume_weights_options <- list()
  if (is.null(soft_subspace_options)) soft_subspace_options <- list()
  
  # defaults for robust fitting
  default_robust <- list(
    type = FALSE,
    k_huber = 1.345,
    c_tukey = 4.685,
    max_iter = 2L,
    scale_scope = "run",
    reestimate_phi = FALSE
  )
  robust <- utils::modifyList(default_robust, robust_options)

  tryCatch({
    robust$type <- match.arg(as.character(robust$type),
                             choices = c("FALSE", "huber", "bisquare"))
  }, error = function(e) {
    stop("Invalid robust_psi/type. Must be one of: FALSE, 'huber', 'bisquare'")
  })
  if (identical(robust$type, "FALSE")) robust$type <- FALSE
  robust$scale_scope <- match.arg(robust$scale_scope, c("run", "global", "voxel", "local"))
  if (identical(robust$scale_scope, "local")) robust$scale_scope <- "voxel"
  stopifnot(is.numeric(robust$k_huber), is.numeric(robust$c_tukey))
  stopifnot(is.numeric(robust$max_iter))
  # Check max_iter even if robust is FALSE - parameter validation should always occur
  if (robust$max_iter < 1) {
    stop("robust_max_iter must be at least 1")
  }
  stopifnot(is.logical(robust$reestimate_phi), length(robust$reestimate_phi) == 1)

  # defaults for autoregressive modelling
  default_ar <- list(
    struct = "iid",
    p = NULL,
    iter_gls = 1L,
    global = FALSE,
    voxelwise = FALSE,
    exact_first = FALSE,
    censor = NULL
  )
  ar <- utils::modifyList(default_ar, ar_options)

  ar$struct <- match.arg(ar$struct, c("iid", "ar1", "ar2", "arp"))
  if (!is.null(ar$p)) stopifnot(is.numeric(ar$p), ar$p >= 1)
  # Validate that p is provided when struct is "arp"
  if (ar$struct == "arp" && is.null(ar$p)) {
    stop("p must be specified in ar_options when struct is 'arp'")
  }
  stopifnot(is.numeric(ar$iter_gls), ar$iter_gls >= 0)
  stopifnot(is.logical(ar$global), length(ar$global) == 1)
  stopifnot(is.logical(ar$voxelwise), length(ar$voxelwise) == 1)
  stopifnot(is.logical(ar$exact_first), length(ar$exact_first) == 1)
  # Validate censor: NULL, "auto", integer vector, or logical vector

  if (!is.null(ar$censor)) {
    if (is.character(ar$censor)) {
      ar$censor <- match.arg(ar$censor, "auto")
    } else if (!is.numeric(ar$censor) && !is.logical(ar$censor)) {
      stop("ar_options$censor must be NULL, 'auto', an integer vector, or a logical vector")
    }
  }

  # defaults for volume weighting
  default_volume_weights <- list(
    enabled = FALSE,
    method = "inverse_squared",
    threshold = 1.5,
    weights = NULL
  )
  volume_weights <- utils::modifyList(default_volume_weights, volume_weights_options)

  stopifnot(is.logical(volume_weights$enabled), length(volume_weights$enabled) == 1)
  volume_weights$method <- match.arg(volume_weights$method,
                                      c("inverse_squared", "soft_threshold", "tukey"))
  stopifnot(is.numeric(volume_weights$threshold), volume_weights$threshold > 0)
  if (!is.null(volume_weights$weights)) {
    stopifnot(is.numeric(volume_weights$weights))
  }

  # defaults for soft subspace projection
  default_soft_subspace <- list(
    enabled = FALSE,
    nuisance_mask = NULL,
    nuisance_matrix = NULL,
    lambda = "auto",
    warn_redundant = TRUE
  )
  soft_subspace <- utils::modifyList(default_soft_subspace, soft_subspace_options)

  stopifnot(is.logical(soft_subspace$enabled), length(soft_subspace$enabled) == 1)
  stopifnot(is.logical(soft_subspace$warn_redundant), length(soft_subspace$warn_redundant) == 1)
  if (soft_subspace$enabled) {
    if (is.null(soft_subspace$nuisance_mask) && is.null(soft_subspace$nuisance_matrix)) {
      stop("soft_subspace requires either nuisance_mask or nuisance_matrix when enabled")
    }
  }
  # Validate lambda
  if (!is.null(soft_subspace$lambda)) {
    if (is.character(soft_subspace$lambda)) {
      soft_subspace$lambda <- match.arg(soft_subspace$lambda, c("auto", "gcv"))
    } else {
      stopifnot(is.numeric(soft_subspace$lambda), soft_subspace$lambda >= 0)
    }
  }

  cfg <- list(robust = robust, ar = ar,
              volume_weights = volume_weights, soft_subspace = soft_subspace)
  class(cfg) <- "fmri_lm_config"
  cfg
}

#' @export
print.fmri_lm_config <- function(x, ...) {
  cat("<fmri_lm_config>\n")
  str(list(robust = x$robust, ar = x$ar,
           volume_weights = x$volume_weights, soft_subspace = x$soft_subspace),
      give.attr = FALSE)
  invisible(x)
}
