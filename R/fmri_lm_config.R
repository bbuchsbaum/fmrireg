#' Configuration for fmri_lm fitting
#'
#' `fmri_lm_control()` creates an `fmri_lm_config` object collecting all
#' options for robust and autoregressive modelling. It validates inputs and
#' applies defaults so downstream functions receive a single structured list.
#'
#' @param robust_options list of robust fitting options. See Details.
#' @param ar_options list of autoregressive modelling options. See Details.
#' @return An object of class `fmri_lm_config`.
#' @details
#' `robust_options` may contain:
#'   * `type` (`FALSE`, "huber", "bisquare")
#'   * `k_huber`
#'   * `c_tukey`
#'   * `max_iter`
#'   * `scale_scope` ("run", "global")
#'   * `reestimate_phi` (logical)
#'
#' `ar_options` may contain:
#'   * `struct` ("iid", "ar1", "ar2", "arp")
#'   * `p` (order for "arp")
#'   * `iter_gls` (integer number of GLS iterations)
#'   * `global` (logical, use global phi)
#'   * `voxelwise` (logical)
#'   * `exact_first` (logical)
#'
#' @export
fmri_lm_control <- function(robust_options = list(),
                            ar_options = list()) {
  # Handle NULL inputs
  if (is.null(robust_options)) robust_options <- list()
  if (is.null(ar_options)) ar_options <- list()

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

  # Validate robust$type
  type_str <- as.character(robust$type)
  valid_types <- c("FALSE", "huber", "bisquare")
  if (!type_str %in% valid_types) {
    fmrireg_abort_config(
      param = "robust_options$type",
      message = "Must be one of: FALSE, 'huber', or 'bisquare'",
      suggestion = paste0("Got '", type_str, "'")
    )
  }
  robust$type <- if (type_str == "FALSE") FALSE else type_str

  # Validate scale_scope
  if (!robust$scale_scope %in% c("run", "global")) {
    fmrireg_abort_config(
      param = "robust_options$scale_scope",
      message = "Must be 'run' or 'global'",
      suggestion = paste0("Got '", robust$scale_scope, "'")
    )
  }

  # Validate numeric parameters
  check_numeric(robust$k_huber, arg = "robust_options$k_huber")
  check_numeric(robust$c_tukey, arg = "robust_options$c_tukey")
  check_numeric(robust$max_iter, arg = "robust_options$max_iter")

  # Check max_iter even if robust is FALSE - parameter validation should always occur
  if (robust$max_iter < 1) {
    fmrireg_abort_config(
      param = "robust_options$max_iter",
      message = "Must be at least 1",
      suggestion = paste0("Got ", robust$max_iter)
    )
  }

  check_logical_scalar(robust$reestimate_phi, arg = "robust_options$reestimate_phi")

  # defaults for autoregressive modelling
  default_ar <- list(
    struct = "iid",
    p = NULL,
    iter_gls = 1L,
    global = FALSE,
    voxelwise = FALSE,
    exact_first = FALSE
  )
  ar <- utils::modifyList(default_ar, ar_options)

  # Validate ar$struct
  valid_structs <- c("iid", "ar1", "ar2", "arp")
  if (!ar$struct %in% valid_structs) {
    fmrireg_abort_config(
      param = "ar_options$struct",
      message = paste0("Must be one of: ", paste(valid_structs, collapse = ", ")),
      suggestion = paste0("Got '", ar$struct, "'")
    )
  }

  # Validate p if provided
  if (!is.null(ar$p)) {
    check_numeric(ar$p, arg = "ar_options$p")
    if (ar$p < 1) {
      fmrireg_abort_config(
        param = "ar_options$p",
        message = "AR order must be at least 1",
        suggestion = paste0("Got ", ar$p)
      )
    }
  }

  # Validate that p is provided when struct is "arp"
  if (ar$struct == "arp" && is.null(ar$p)) {
    fmrireg_abort_config(
      param = "ar_options$p",
      message = "Required when struct is 'arp'",
      suggestion = "Specify the AR order, e.g., ar_options = list(struct = 'arp', p = 3)"
    )
  }

  # Validate iter_gls
  check_numeric(ar$iter_gls, arg = "ar_options$iter_gls")
  if (ar$iter_gls < 0) {
    fmrireg_abort_config(
      param = "ar_options$iter_gls",
      message = "Must be non-negative",
      suggestion = paste0("Got ", ar$iter_gls)
    )
  }

  # Validate logical parameters
  check_logical_scalar(ar$global, arg = "ar_options$global")
  check_logical_scalar(ar$voxelwise, arg = "ar_options$voxelwise")
  check_logical_scalar(ar$exact_first, arg = "ar_options$exact_first")

  cfg <- list(robust = robust, ar = ar)
  class(cfg) <- "fmri_lm_config"
  cfg
}

#' @export
print.fmri_lm_config <- function(x, ...) {
  cat("<fmri_lm_config>\n")
  str(list(robust = x$robust, ar = x$ar), give.attr = FALSE)
  invisible(x)
}
