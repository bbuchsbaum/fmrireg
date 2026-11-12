# matrix_frame(), neurovec_frame(), nifti_frame() and latent_frame() are
# deprecated forwarders to fmridataset::as_fmri_frame() / read_nifti_bold().
# They must warn, and return exactly the frame the replacement returns.
#
# The constructors generate the axis IDs, so any dimnames the caller's matrix
# carries are resolved: labels are adopted as IDs where that is meaningful and
# dropped otherwise (prcomp writes PC1..PCr, neuroim2::series() labels voxel
# columns). The list form of latent_frame() is built in fmrireg itself.

frame_digest <- function(frame) {
  fmridataset::identity_descriptor(frame)$digest
}

expect_deprecated_frame <- function(expr) {
  expect_warning(expr, class = "deprecatedWarning")
}

test_that("matrix_frame accepts a matrix carrying dimnames", {
  Y <- matrix(rnorm(20), 10, 2,
    dimnames = list(paste0("t", 1:10), c("roi_a", "roi_b"))
  )
  frame <- as_fmri_frame(Y, TR = 2, run_length = 10)

  # Column names are meaningful feature labels and are kept as the IDs.
  expect_identical(fmridataset::feature_ids(frame), c("roi_a", "roi_b"))
  # Row names describe nothing the frame models; volume IDs are generated.
  expect_identical(
    fmridataset::observation_ids(frame)[1:2],
    c("run-1_vol-0001", "run-1_vol-0002")
  )
  expect_equal(unname(frame_data(frame)), unname(Y))
})

test_that("matrix_frame keeps generated feature IDs when column names are unusable", {
  Y <- matrix(rnorm(20), 10, 2, dimnames = list(NULL, c("dup", "dup")))
  frame <- as_fmri_frame(Y, TR = 2, run_length = 10)

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

  expect_deprecated_frame(frame <- latent_frame(
    list(basis = pr$x, loadings = pr$rotation, mask = mask),
    TR = 2, run_length = 40
  ))

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

  expect_deprecated_frame(frame <- latent_frame(
    list(basis = scores, loadings = loadings, mask = mask),
    TR = 2, run_length = 40
  ))

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

  frame <- as_fmri_frame(vec, mask = mask, TR = 2)

  expect_equal(dim(frame), c(6L, 8L))
  expect_s3_class(fmridataset::space(frame), "volume_space")
})

test_that("each deprecated constructor warns and returns the replacement's frame", {
  set.seed(3)
  Y <- matrix(rnorm(24), 12, 2)
  events <- data.frame(onset = c(2, 10), run = c(1, 2))
  censor <- c(rep(FALSE, 11), TRUE)
  expect_deprecated_frame(
    old <- matrix_frame(Y, TR = 2, run_length = c(6, 6), event_table = events,
                        censor = censor)
  )
  new <- as_fmri_frame(Y, TR = 2, run_length = c(6, 6), events = events,
                       censor = censor)
  expect_identical(frame_digest(old), frame_digest(new))

  # a vector is one feature, as before
  expect_deprecated_frame(v <- matrix_frame(Y[, 1], TR = 2, run_length = 12))
  expect_equal(dim(v), c(12L, 1L))

  dim3 <- c(2L, 2L, 2L)
  sp4 <- neuroim2::NeuroSpace(c(dim3, 6L))
  vecs <- lapply(1:2, function(i) neuroim2::NeuroVec(array(rnorm(8 * 6), c(dim3, 6L)), sp4))
  mask <- neuroim2::LogicalNeuroVol(array(c(TRUE, FALSE), dim3), neuroim2::NeuroSpace(dim3))
  expect_deprecated_frame(old <- neurovec_frame(vecs, mask, TR = 2))
  new <- as_fmri_frame(vecs, mask = mask, TR = 2)
  expect_identical(frame_digest(old), frame_digest(new))
  expect_equal(unname(fmridataset::collect_assay(old)), unname(fmridataset::collect_assay(new)))

  dir <- withr::local_tempdir()
  files <- vapply(1:2, function(i) {
    f <- file.path(dir, sprintf("run-%d.nii", i))
    neuroim2::write_vec(vecs[[i]], f)
    f
  }, character(1))
  neuroim2::write_vol(mask, file.path(dir, "mask.nii"))
  expect_deprecated_frame(
    old <- nifti_frame(basename(files), "mask.nii", TR = 2, base_path = dir)
  )
  new <- read_nifti_bold(files, file.path(dir, "mask.nii"), TR = 2)
  expect_equal(dim(old), dim(new))
  expect_identical(fmridataset::observation_ids(old), fmridataset::observation_ids(new))
  expect_equal(unname(fmridataset::collect_assay(old)),
               unname(fmridataset::collect_assay(new)), tolerance = 1e-6)
})

test_that("saved jobs naming the old constructors realize without a warning", {
  Y <- matrix(rnorm(20), 10, 2)
  spec <- dataset_spec("matrix_frame",
                       args = list(datamat = Y, TR = 2, run_length = 10),
                       source = "inline")
  job <- structure(list(id = "s1", dataset_spec = spec), class = "fmri_job")
  expect_no_warning(frame <- realize_dataset(job))
  expect_identical(frame_digest(frame),
                   frame_digest(as_fmri_frame(Y, TR = 2, run_length = 10)))
})
