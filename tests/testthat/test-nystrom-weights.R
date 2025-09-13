test_that("Landmark weight matrix has k-sparsity and row sums = 1", {
  skip_on_cran()

  set.seed(5)
  nx <- 7L; ny <- 5L; nz <- 1L
  coords <- as.matrix(expand.grid(x = seq_len(nx), y = seq_len(ny), z = seq_len(nz)))
  V <- nrow(coords)
  # Choose L random landmarks from coords
  L <- 12L
  lcoords <- coords[sample.int(V, L), , drop = FALSE]
  k <- 4L

  W <- fmrireg:::build_landmark_weights(coords, lcoords, k = k)

  # Row sums are ~1
  rs <- Matrix::rowSums(W)
  expect_true(max(abs(rs - 1)) < 1e-10)

  # Each row has at most k nonzeros
  nnz_per_row <- Matrix::rowSums(W != 0)
  expect_true(all(nnz_per_row <= k))
})

