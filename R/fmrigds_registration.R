# fmrireg's spatial FDR (Efron central-matching empirical null, see
# spatial_fdr()) exposed as an fmrigds post-hoc step. It is registered under a
# namespaced name so fmrigds's own "fdr:spatial" (Simes + weighted BH) and its
# reducers are never replaced; fmrireg does not register any reducer.

.fmrigds_spatial_posthoc_compute <- function(arrays, group, alpha) {
  # Get p or z from either top-level arrays or arrays$assays
  p_arr <- arrays[["p"]]
  if (is.null(p_arr) && !is.null(arrays$assays)) p_arr <- arrays$assays[["p"]]
  z_arr <- arrays[["z"]]
  if (is.null(z_arr) && !is.null(arrays$assays)) z_arr <- arrays$assays[["z"]]

  # Normalize to [P x S x K]
  if (is.null(p_arr) && is.null(z_arr)) {
    stop("fmrireg:fdr:spatial requires p (or z to derive p)", call. = FALSE)
  }
  if (!is.null(z_arr)) {
    z_dims <- dim(z_arr)
    if (is.null(z_dims)) z_arr <- array(z_arr, dim = c(length(z_arr), 1L, 1L))
    if (length(dim(z_arr)) == 2L) z_arr <- array(z_arr, dim = c(dim(z_arr)[1L], 1L, dim(z_arr)[2L]))
    # keep z for empirical-null handling inside spatial_fdr(z=...)
  } else {
    p_dims <- dim(p_arr)
    if (is.null(p_dims)) p_arr <- array(p_arr, dim = c(length(p_arr), 1L, 1L))
    if (length(dim(p_arr)) == 2L) p_arr <- array(p_arr, dim = c(dim(p_arr)[1L], 1L, dim(p_arr)[2L]))
  }

  dims <- if (!is.null(z_arr)) dim(z_arr) else dim(p_arr)
  if (is.null(dims) || length(dims) != 3L) {
    stop("fmrireg:fdr:spatial expected assay dimensions [features x subjects x contrasts]", call. = FALSE)
  }
  Q <- array(NA_real_, dim = dims)
  for (j in seq_len(dims[2L])) {
    for (k in seq_len(dims[3L])) {
      if (!is.null(z_arr)) {
        zv <- z_arr[, j, k]
        sres <- fmrireg::spatial_fdr(z = zv, group = group, alpha = alpha)
      } else {
        pv <- p_arr[, j, k]
        sres <- fmrireg::spatial_fdr(p = pv, group = group, alpha = alpha)
      }
      Q[, j, k] <- sres$q
    }
  }
  Q
}

.fmrigds_spatial_posthoc <- function(arrays, opts) {
  group <- opts$group
  if (is.null(group)) stop("fmrireg:fdr:spatial requires options$group (length n_features)", call. = FALSE)
  alpha <- opts$alpha %||% 0.05
  q <- .fmrigds_spatial_posthoc_compute(arrays, group, alpha)
  list(q = q)
}

.fmrireg_posthoc_name <- "fmrireg:fdr:spatial"

.register_fmrireg_posthoc <- function() {
  if (!is.null(fmrigds::get_posthoc(.fmrireg_posthoc_name))) {
    return(invisible(FALSE))
  }
  fmrigds::register_posthoc(
    name = .fmrireg_posthoc_name,
    fun = .fmrigds_spatial_posthoc,
    requires = c("p"),
    provides = c("q"),
    overwrite = FALSE
  )
  invisible(TRUE)
}
