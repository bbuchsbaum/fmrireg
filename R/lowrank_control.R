#' Low-rank / sketch controls for fast GLM
#'
#' Control object to enable the optional sketched GLM engine.
#'
#' @param parcels Optional parceling, e.g., a neuroim2::ClusteredNeuroVol or
#'   integer vector (length = number of voxels in mask).
#' @param landmarks Optional integer; number of landmark voxels for optional
#'   Nyström extension (NULL = off).
#' @param k_neighbors Integer; k for k-NN in Nyström extension.
#' @param time_sketch List controlling the temporal sketch, with elements
#'   \describe{
#'     \item{`method`}{`"gaussian"` (default), `"countsketch"`, `"srht"`
#'       (subsampled randomized Hadamard transform) or `"ihs"` (iterative
#'       Hessian sketch).}
#'     \item{`m`}{Number of sketch rows; `NULL` uses `min(8 * p, T)` for a
#'       design with `p` columns and `T` scans. The sketch-and-solve methods
#'       require `m > p`.}
#'     \item{`iters`}{`"ihs"` only: maximum number of iterations (default
#'       100). Must be at least 1.}
#'     \item{`tol`}{`"ihs"` only: stop once every coefficient moves by less
#'       than `tol` OLS standard errors (default `1e-3`); `0` runs exactly
#'       `iters` iterations. A fit that stops at `iters` before reaching
#'       `tol` warns.}
#'   }
#'   All sketches are normalised (\eqn{E\|Sr\|^2 = \|r\|^2}), so residual
#'   variances, `rss` and covariances are reported on the data scale.
#'   `"gaussian"`, `"countsketch"` and `"srht"` are sketch-and-solve
#'   estimators: they fit the model to `m` sketched rows, and their standard
#'   errors are the conditional-on-sketch standard errors of that estimator
#'   (roughly `sqrt(T / m)` times the full-data OLS ones), with Satterthwaite
#'   residual degrees of freedom (about `m - p`). `"ihs"` iterates to the
#'   full-data least-squares solution and reports exact OLS covariance,
#'   residual variance and `T - p` degrees of freedom.
#'
#'   Estimability is judged on the full (whitened) design with the same
#'   pivoted-QR rule that exact fits use. Aliased coefficients are reported
#'   as `NA`, and contrasts that load on them are `NA` with a warning. If the
#'   sketched design restricted to the estimable columns is numerically rank
#'   deficient, `m` is too small for the design and the fit stops with an
#'   error rather than regularising.
#' @param ncomp Optional integer; number of latent components within parcels (PCA).
#' @param noise_pcs Integer; optional GLMdenoise-style PCs from low-R2 parcels.
#' @return A list with class "lowrank_control".
#' @export
lowrank_control <- function(parcels = NULL,
                            landmarks = NULL,
                            k_neighbors = 16L,
                            time_sketch = list(method = "gaussian", m = NULL),
                            ncomp = NULL,
                            noise_pcs = 0L) {
  if (!is.null(landmarks)) {
    stopifnot(is.numeric(landmarks), all(landmarks > 0))
    stopifnot(is.numeric(k_neighbors), all(k_neighbors > 0))
  }
  if (!is.null(time_sketch)) {
    stopifnot(is.list(time_sketch))
    if (is.null(time_sketch[["method"]])) time_sketch$method <- "gaussian"
    method <- time_sketch[["method"]]
    if (!is.character(method) || length(method) != 1L ||
        !method %in% c("gaussian", "countsketch", "srht", "ihs")) {
      stop("`time_sketch$method` must be one of \"gaussian\", \"countsketch\", ",
           "\"srht\" or \"ihs\"", call. = FALSE)
    }
    iters <- time_sketch[["iters"]]
    if (identical(method, "ihs") && !is.null(iters) &&
        (length(iters) != 1L || is.na(iters) || iters < 1L)) {
      stop("`time_sketch$iters` must be >= 1 for method = \"ihs\"", call. = FALSE)
    }
  }
  structure(list(parcels = parcels,
                 landmarks = landmarks,
                 k_neighbors = as.integer(k_neighbors),
                 time_sketch = time_sketch,
                 ncomp = ncomp,
                 noise_pcs = as.integer(noise_pcs)),
            class = "lowrank_control")
}

