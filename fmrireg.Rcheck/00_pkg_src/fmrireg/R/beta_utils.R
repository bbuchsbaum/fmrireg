#' @keywords internal
#' @noRd
build_design_data <- function(bdes) {
  X <- if (is.null(bdes$dmat_fixed)) bdes$dmat_ran else cbind(bdes$dmat_ran, bdes$dmat_fixed)
  Base <- as.matrix(bdes$dmat_base)
  X[is.na(X)] <- 0
  list(Base = Base, X = X)
}

#' @keywords internal
#' @noRd
masked_vectors <- function(dset) {
  neuroim2::vectors(get_data(dset), subset = which(get_mask(dset) > 0))
}

#' Apply a function to voxel vectors with optional progress bar
#'
#' @keywords internal
#' @noRd
map_voxels <- function(vecs, FUN, ..., .progress = TRUE) {
  worker <- function(v) FUN(v, ...)
  if (.progress) {
    with_package("progressr")
    progressr::with_progress({
      p <- progressr::progressor(along = vecs)
      res <- furrr::future_map(vecs, function(v) { p(); worker(v) })
    })
  } else {
    res <- furrr::future_map(vecs, worker)
  }
  do.call(cbind, res)
}
