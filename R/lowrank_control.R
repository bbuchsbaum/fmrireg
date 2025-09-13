#' Low-rank / sketch controls for fast GLM
#'
#' Control object to enable the optional sketched GLM engine.
#'
#' @param parcels Optional parceling, e.g., a neuroim2::ClusteredNeuroVol or
#'   integer vector (length = number of voxels in mask).
#' @param landmarks Optional integer; number of landmark voxels for optional
#'   Nyström extension (NULL = off).
#' @param k_neighbors Integer; k for k-NN in Nyström extension.
#' @param time_sketch List(method = "gaussian" | "countsketch", m = NULL, iters = 0L).
#' @param ncomp Optional integer; number of latent components within parcels (PCA).
#' @param noise_pcs Integer; optional GLMdenoise-style PCs from low-R2 parcels.
#' @return A list with class "lowrank_control".
#' @export
lowrank_control <- function(parcels = NULL,
                            landmarks = NULL,
                            k_neighbors = 16L,
                            time_sketch = list(method = "gaussian", m = NULL, iters = 0L),
                            ncomp = NULL,
                            noise_pcs = 0L) {
  if (!is.null(landmarks)) {
    stopifnot(is.numeric(landmarks), all(landmarks > 0))
    stopifnot(is.numeric(k_neighbors), all(k_neighbors > 0))
  }
  if (!is.null(time_sketch)) {
    stopifnot(is.list(time_sketch))
    if (is.null(time_sketch$method)) time_sketch$method <- "gaussian"
  }
  structure(list(parcels = parcels,
                 landmarks = landmarks,
                 k_neighbors = as.integer(k_neighbors),
                 time_sketch = time_sketch,
                 ncomp = ncomp,
                 noise_pcs = as.integer(noise_pcs)),
            class = "lowrank_control")
}

