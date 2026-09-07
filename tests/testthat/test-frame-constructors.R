# The frame constructors generate the axis IDs, so any dimnames the caller's
# matrix happens to carry are a second, unchecked naming of the same axes.
# fmridataset rejects that as ambiguous rather than silently dropping it, so
# the constructors must resolve it: adopt the labels as IDs where that is
# meaningful, and drop them otherwise. Matrices carrying names are ordinary
# input (prcomp writes PC1..PCr, neuroim2::series() labels voxel columns), so
# without this every such caller hit an alignment error.

test_that("matrix_frame accepts a matrix carrying dimnames", {
  Y <- matrix(rnorm(20), 10, 2,
    dimnames = list(paste0("t", 1:10), c("roi_a", "roi_b"))
  )
  frame <- matrix_frame(Y, TR = 2, run_length = 10)

  # Column names are meaningful feature labels and are kept as the IDs.
  expect_identical(fmridataset::feature_ids(frame), c("roi_a", "roi_b"))
  # Row names describe nothing the frame models; volume IDs are generated.
  expect_identical(
    fmridataset::observation_ids(frame)[1:2],
    c("run-1-vol-00001", "run-1-vol-00002")
  )
  expect_equal(unname(frame_data(frame)), unname(Y))
})

test_that("matrix_frame keeps generated feature IDs when column names are unusable", {
  Y <- matrix(rnorm(20), 10, 2, dimnames = list(NULL, c("dup", "dup")))
  frame <- matrix_frame(Y, TR = 2, run_length = 10)

  expect_identical(fmridataset::feature_ids(frame), c("feature-1", "feature-2"))
  expect_equal(unname(frame_data(frame)), unname(Y))
})

test_that("latent_frame adopts component names from a decomposition", {
  set.seed(1)
  Y <- matrix(rnorm(40 * 8), 40, 8)
  pr <- prcomp(Y, center = TRUE, rank. = 3)
  mask <- neuroim2::LogicalNeuroVol(
    array(TRUE, c(2L, 2L, 2L)), neuroim2::NeuroSpace(c(2L, 2L, 2L))
  )

  frame <- latent_frame(
    list(basis = pr$x, loadings = pr$rotation, mask = mask),
    TR = 2, run_length = 40
  )

  # prcomp names its columns PC1..PC3; those are the component identities.
  expect_identical(fmridataset::feature_ids(frame), colnames(pr$x))
  expect_equal(unname(frame_data(frame)), unname(pr$x))
  # The loadings remain reachable as the basis synthesis operator.
  expect_equal(
    dim(fmridataset::basis_synthesis(fmridataset::space(frame))),
    c(8L, 3L)
  )
})

test_that("latent_frame generates component IDs when the basis is unnamed", {
  set.seed(2)
  scores <- matrix(rnorm(40 * 3), 40, 3)
  loadings <- matrix(rnorm(8 * 3), 8, 3)
  mask <- neuroim2::LogicalNeuroVol(
    array(TRUE, c(2L, 2L, 2L)), neuroim2::NeuroSpace(c(2L, 2L, 2L))
  )

  frame <- latent_frame(
    list(basis = scores, loadings = loadings, mask = mask),
    TR = 2, run_length = 40
  )

  expect_identical(
    fmridataset::feature_ids(frame),
    c("component-1", "component-2", "component-3")
  )
})

test_that("neurovec_frame builds a frame from a named series matrix", {
  dim3 <- c(2L, 2L, 2L)
  arr <- array(rnorm(prod(dim3) * 6), c(dim3, 6L))
  vec <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(c(dim3, 6L)))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dim3), neuroim2::NeuroSpace(dim3))

  frame <- neurovec_frame(vec, mask = mask, TR = 2)

  expect_equal(dim(frame), c(6L, 8L))
  expect_s3_class(fmridataset::space(frame), "volume_space")
})
