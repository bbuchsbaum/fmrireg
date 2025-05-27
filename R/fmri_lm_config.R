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

  robust$type <- match.arg(as.character(robust$type),
                           choices = c("FALSE", "huber", "bisquare"))
  if (identical(robust$type, "FALSE")) robust$type <- FALSE
  robust$scale_scope <- match.arg(robust$scale_scope, c("run", "global"))
  stopifnot(is.numeric(robust$k_huber), is.numeric(robust$c_tukey))
  stopifnot(is.numeric(robust$max_iter), robust$max_iter >= 0)
  stopifnot(is.logical(robust$reestimate_phi), length(robust$reestimate_phi) == 1)

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

  ar$struct <- match.arg(ar$struct, c("iid", "ar1", "ar2", "arp"))
  if (!is.null(ar$p)) stopifnot(is.numeric(ar$p), ar$p >= 1)
  stopifnot(is.numeric(ar$iter_gls), ar$iter_gls >= 0)
  stopifnot(is.logical(ar$global), length(ar$global) == 1)
  stopifnot(is.logical(ar$voxelwise), length(ar$voxelwise) == 1)
  stopifnot(is.logical(ar$exact_first), length(ar$exact_first) == 1)

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
