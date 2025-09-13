#' Build sparse landmark weight matrix W (V x L) using k-NN heat-kernel
#' @param coords V x d matrix of voxel coordinates (d = 2 or 3)
#' @param lcoords L x d matrix of landmark coordinates (subset of coords or arbitrary within space)
#' @param k integer; number of nearest landmarks used per voxel (default 16)
#' @param h optional bandwidth; if NULL, uses median of k-NN distances
#' @return dgCMatrix W with rows summing to 1 (each voxel is a convex combination of k landmarks)
#' @keywords internal
build_landmark_weights <- function(coords, lcoords, k = 16L, h = NULL) {
  stopifnot(is.matrix(coords), is.matrix(lcoords))
  stopifnot(nrow(coords) > 0, nrow(lcoords) > 0)
  k <- as.integer(k)
  k <- max(1L, min(k, nrow(lcoords)))

  nn <- RANN::nn2(data = lcoords, query = coords, k = k)
  idx <- as.integer(nn$nn.idx)
  dst <- as.numeric(nn$nn.dists)

  if (is.null(h)) {
    # median distance to the kth neighbor across all voxels
    # guard against zeros / degeneracies
    kth <- nn$nn.dists[, k]
    h <- stats::median(kth)
    if (!is.finite(h) || h <= 0) h <- 1.0
  }

  # Heat-kernel weights per voxel over its k nearest landmarks
  w <- exp(-(dst * dst) / (2 * h * h))

  V <- nrow(coords); L <- nrow(lcoords)
  i <- rep.int(seq_len(V), k)
  j <- idx
  W <- Matrix::sparseMatrix(i = i, j = j, x = w, dims = c(V, L))

  # Normalize rows to sum to 1
  rs <- Matrix::rowSums(W)
  rs[rs == 0] <- 1
  W <- Matrix::Diagonal(x = 1 / rs) %*% W
  Matrix::drop0(W)
}

#' Apply Nyström / barycentric extension: B (p x V) from B_L (p x L)
#' @param BL p x L dense matrix of betas on landmarks
#' @param W  V x L sparse weight matrix (rows sum to 1)
#' @return p x V dense matrix of extended betas
#' @keywords internal
extend_betas_landmarks <- function(BL, W) {
  stopifnot(is.matrix(BL), inherits(W, "sparseMatrix") || inherits(W, "dgCMatrix"))
  as.matrix(BL %*% Matrix::t(W))
}

