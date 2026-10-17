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
  data_matrix <- .dset_data_matrix(dset)
  neuroim2::vectors(data_matrix, subset = seq_len(ncol(data_matrix)))
}

#' Apply a function to voxel vectors through the active future plan
#'
#' @keywords internal
#' @noRd
map_voxels <- function(vecs, FUN, ..., .progress = TRUE) {
  res <- .fmrireg_lapply(vecs, function(v) FUN(v, ...), parallel = TRUE,
                         progress = .progress, label = "Fitting voxels")
  do.call(cbind, res)
}
