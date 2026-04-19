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

  tryCatch({
    robust$type <- match.arg(as.character(robust$type),
                             choices = c("FALSE", "huber", "bisquare"))
  }, error = function(e) {
    stop("Invalid robust_psi/type. Must be one of: FALSE, 'huber', 'bisquare'")
  })
  if (identical(robust$type, "FALSE")) robust$type <- FALSE
  robust$scale_scope <- match.arg(robust$scale_scope, c("run", "global"))
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
    exact_first = FALSE
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

#' Get AR order from config
#'
#' Extracts the autoregressive order from an fmri_lm_config object.
#' This centralizes the mapping from struct names to numeric order.
#'
#' @param cfg An `fmri_lm_config` object or a list with `$ar$struct` and optionally `$ar$p`
#' @return Integer AR order (0L for iid/none, 1L for ar1, 2L for ar2, etc.)
#' @keywords internal
#' @noRd
get_ar_order <- function(cfg) {

  struct <- cfg$ar$struct %||% "iid"
  switch(as.character(struct),
         iid = 0L,
         none = 0L,
         ar1 = 1L,
         ar2 = 2L,
         ar3 = 3L,
         ar4 = 4L,
         arp = as.integer(cfg$ar$p %||% 0L),
         0L)
}

#' Build fmri_lm config from mixed parameter sources
#'
#' Consolidates config building logic used by fmri_lm methods.
#' Handles backward compatibility parameters and option merging.
#'
#' @param robust Logical or character ("huber", "bisquare")
#' @param robust_options List of robust fitting options
#' @param ar_options List of AR modeling options
#' @param engine_robust_options Override robust options from engine
#' @param engine_ar_options Override AR options from engine
#' @param engine_cfg Pre-built config from engine
#' @param cor_struct Shorthand for ar_options$struct
#' @param cor_iter Shorthand for ar_options$iter_gls
#' @param cor_global Shorthand for ar_options$global
#' @param ar1_exact_first Shorthand for ar_options$exact_first
#' @param ar_p Shorthand for ar_options$p
#' @param ar_voxelwise Shorthand for ar_options$voxelwise
#' @param robust_psi Shorthand for robust_options$type
#' @param robust_max_iter Shorthand for robust_options$max_iter
#' @param robust_scale_scope Shorthand for robust_options$scale_scope
#' @return An fmri_lm_config object
#' @keywords internal
#' @noRd
build_fmri_lm_cfg <- function(robust = FALSE, robust_options = NULL, ar_options = NULL,
                               engine_robust_options = NULL, engine_ar_options = NULL,
                               engine_cfg = NULL,
                               cor_struct = NULL, cor_iter = NULL, cor_global = NULL,
                               ar1_exact_first = NULL, ar_p = NULL, ar_voxelwise = FALSE,
                               robust_psi = NULL, robust_max_iter = NULL,
                               robust_scale_scope = NULL) {

  # Convert robust parameter to type
  robust_type <- if (is.logical(robust)) {
    if (robust) "huber" else FALSE
  } else {
    robust
  }

  # Initialize options lists
  if (is.null(robust_options)) robust_options <- list()
  if (is.null(ar_options)) ar_options <- list()

  # Merge engine-supplied options
  if (!is.null(engine_robust_options)) {
    robust_options <- utils::modifyList(robust_options, engine_robust_options)
  }
  if (!is.null(engine_ar_options)) {
    ar_options <- utils::modifyList(ar_options, engine_ar_options)
  }

  # Apply robust type if not already set

  if (!is.null(robust_type) && !("type" %in% names(robust_options))) {
    robust_options$type <- robust_type
  }

  # Merge individual robust parameters (backward compatibility)
  if (!is.null(robust_psi) && !("type" %in% names(robust_options))) {
    robust_options$type <- robust_psi
  }
  if (!is.null(robust_max_iter) && !("max_iter" %in% names(robust_options))) {
    robust_options$max_iter <- robust_max_iter
  }
  if (!is.null(robust_scale_scope) && !("scale_scope" %in% names(robust_options))) {
    robust_options$scale_scope <- robust_scale_scope
  }

  # Handle low-rank shorthand order -> struct/p mapping
  if (!is.null(ar_options$order) && is.null(ar_options$struct)) {
    ar_order_tmp <- as.integer(ar_options$order[1])
    if (!is.finite(ar_order_tmp) || ar_order_tmp <= 0L) {
      ar_options$struct <- "iid"
    } else if (ar_order_tmp == 1L) {
      ar_options$struct <- "ar1"
    } else if (ar_order_tmp == 2L) {
      ar_options$struct <- "ar2"
    } else {
      ar_options$struct <- "arp"
      ar_options$p <- ar_options$p %||% ar_order_tmp
    }
    ar_options$order <- NULL
  }

  # Merge individual AR parameters (backward compatibility)
  if (!is.null(cor_struct) && !("struct" %in% names(ar_options))) {
    ar_options$struct <- cor_struct
  }
  if (!is.null(cor_iter) && !("iter_gls" %in% names(ar_options))) {
    ar_options$iter_gls <- cor_iter
  }
  if (!is.null(cor_global) && !("global" %in% names(ar_options))) {
    ar_options$global <- cor_global
  }
  if (!is.null(ar1_exact_first) && !("exact_first" %in% names(ar_options))) {
    ar_options$exact_first <- ar1_exact_first
  }
  if (!is.null(ar_p) && !("p" %in% names(ar_options))) {
    ar_options$p <- ar_p
  }
  if (!("voxelwise" %in% names(ar_options))) {
    ar_options$voxelwise <- ar_voxelwise
  }

  # Create config object
  cfg <- if (!is.null(engine_cfg) && inherits(engine_cfg, "fmri_lm_config")) {
    engine_cfg
  } else {
    fmri_lm_control(robust_options = robust_options, ar_options = ar_options)
  }

  # Merge with engine_cfg if both provided
  if (!is.null(engine_cfg) && inherits(engine_cfg, "fmri_lm_config")) {
    cfg$robust <- engine_cfg$robust
    cfg$ar <- utils::modifyList(cfg$ar, engine_cfg$ar)
  }

  cfg
}
