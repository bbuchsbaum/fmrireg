# tests/testthat/test-nystrom-identity.R

test_that("Nyström extension reduces to identity when landmarks == voxels (k = 1)", {
  skip_on_cran()
  # Matrix and RANN are in Imports; no skip needed in this repo

  # Small 2D grid of voxel coordinates
  nx <- 6L; ny <- 6L; nz <- 1L
  coords <- as.matrix(expand.grid(x = seq_len(nx), y = seq_len(ny), z = seq_len(nz)))
  V <- nrow(coords)

  # Landmarks are all voxels, in the same order
  lcoords <- coords

  # Build W with k = 1 → each voxel maps to itself with weight 1
  W <- fmrireg:::build_landmark_weights(coords, lcoords, k = 1L)

  # Check row sums are exactly 1
  rs <- Matrix::rowSums(W)
  expect_true(max(abs(rs - 1)) < 1e-12)

  # Check W is (numerically) identity
  I <- Matrix::Diagonal(V)
  expect_true(max(abs(W - I)) < 1e-12)

  # Extension preserves betas exactly
  set.seed(7)
  p <- 3L
  BL <- matrix(rnorm(p * V), nrow = p, ncol = V)
  B  <- fmrireg:::extend_betas_landmarks(BL, W)
  expect_equal(B, BL, tolerance = 1e-12)
})
