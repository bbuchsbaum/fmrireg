# group_data_nifti.R: NIfTI IO helpers, dimension validation, deprecation.

make_tiny_nifti <- function(vals = rnorm(8), dims = c(2, 2, 2)) {
  stopifnot(length(vals) == prod(dims))
  sp <- neuroim2::NeuroSpace(dims)
  vol <- neuroim2::NeuroVol(array(vals, dim = dims), sp)
  path <- tempfile(fileext = ".nii")
  neuroim2::write_vol(vol, path)
  path
}

test_that("group_data_from_nifti emits deprecation when not suppressed", {
  old <- getOption("fmrireg.suppress_deprecation")
  on.exit(options(fmrireg.suppress_deprecation = old), add = TRUE)
  options(fmrireg.suppress_deprecation = FALSE)

  expect_warning(
    tryCatch(
      group_data_from_nifti(beta_paths = "a.nii", se_paths = "b.nii", validate = FALSE),
      error = function(e) e
    ),
    "deprecated|group_data"
  )
})

test_that("read_nifti_header / validate_nifti_dimensions round-trip tiny volumes", {
  skip_if_not_installed("neuroim2")
  set.seed(3)
  path <- make_tiny_nifti(1:8, c(2, 2, 2))
  on.exit(unlink(path), add = TRUE)

  hdr <- fmrireg:::read_nifti_header(path)
  expect_equal(as.integer(hdr$dim), c(2L, 2L, 2L))
  expect_true(length(hdr$voxel_size) >= 3L)

  expect_true(fmrireg:::validate_nifti_dimensions(path, c(2, 2, 2)))
  expect_error(
    fmrireg:::validate_nifti_dimensions(path, c(3, 3, 3)),
    "Dimension mismatch"
  )
  # Non-existent paths are skipped inside the loop
  expect_true(
    fmrireg:::validate_nifti_dimensions("/no/such/file.nii", c(2, 2, 2))
  )
})

test_that("read_nifti_data returns full array or masked vector", {
  skip_if_not_installed("neuroim2")
  set.seed(4)
  vals <- rnorm(8)
  path <- make_tiny_nifti(vals, c(2, 2, 2))
  on.exit(unlink(path), add = TRUE)

  arr <- fmrireg:::read_nifti_data(path)
  expect_true(is.atomic(arr) || inherits(arr, "array") || inherits(arr, "NeuroVol"))
  flat <- as.vector(as.array(arr))
  expect_equal(length(flat), 8L)

  mask <- array(c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE), dim = c(2, 2, 2))
  masked <- fmrireg:::read_nifti_data(path, mask = mask)
  expect_true(is.atomic(masked))
  expect_equal(length(masked), sum(mask))

  expect_error(
    suppressWarnings(fmrireg:::read_nifti_data("/definitely/missing/file.nii")),
    "Failed to read NIfTI"
  )
})

test_that("read_nifti_chunk / read_nifti_full validate object class and read data", {
  skip_if_not_installed("neuroim2")
  set.seed(5)
  beta1 <- make_tiny_nifti(rnorm(8), c(2, 2, 2))
  beta2 <- make_tiny_nifti(rnorm(8), c(2, 2, 2))
  se1 <- make_tiny_nifti(abs(rnorm(8)) + 0.1, c(2, 2, 2))
  se2 <- make_tiny_nifti(abs(rnorm(8)) + 0.1, c(2, 2, 2))
  on.exit(unlink(c(beta1, beta2, se1, se2)), add = TRUE)

  gd <- fmrireg:::.group_data_from_nifti_impl(
    beta_paths = c(beta1, beta2),
    se_paths = c(se1, se2),
    subjects = c("s1", "s2"),
    validate = TRUE
  )
  expect_s3_class(gd, "group_data_nifti")
  expect_equal(gd$n_subjects, 2L)

  chunk <- fmrireg:::read_nifti_chunk(gd, voxel_indices = c(1L, 3L, 5L))
  expect_true(is.list(chunk) || is.matrix(chunk$beta) || !is.null(chunk$beta))
  if (is.list(chunk) && !is.null(chunk$beta)) {
    expect_equal(nrow(chunk$beta), 2L)
    expect_equal(ncol(chunk$beta), 3L)
  }

  full <- fmrireg:::read_nifti_full(gd, use_mask = FALSE)
  expect_true(is.list(full))
  expect_true(!is.null(full$beta) || !is.null(full$betas))
})

test_that("t_to_beta_se covers n=NULL and n-provided SE scaling", {
  out <- t_to_beta_se(t = c(2, -1), df = 20)
  expect_equal(out$beta, c(2, -1))
  expect_equal(as.numeric(out$se), 1)

  out_n <- t_to_beta_se(t = c(2, -1), df = 20, n = 25)
  expect_equal(as.numeric(out_n$se), sqrt(1 / 25))
  expect_equal(out_n$beta, c(2, -1) * sqrt(1 / 25))
})
