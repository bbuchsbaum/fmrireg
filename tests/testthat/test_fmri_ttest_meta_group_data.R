test_that("fmri_ttest(meta) works with group_data_nifti and weights", {
  skip_on_cran()
  skip_if_not_installed("neuroim2")

  library(neuroim2)

  # Tiny synthetic 3D space, 1 slice
  space <- NeuroSpace(c(4, 4, 1), spacing = c(2, 2, 2))
  n_vox <- prod(dim(space))

  # Subjects and groups
  ids <- sprintf("sub-%02d", 1:6)
  grp <- factor(rep(c("A", "B"), each = 3))

  # Create per-subject beta and SE images
  tmpdir <- tempdir()
  beta_paths <- character(length(ids))
  se_paths   <- character(length(ids))

  set.seed(1)
  for (i in seq_along(ids)) {
    b <- array(0, dim = c(4, 4, 1))
    b[] <- if (grp[i] == "A") 0.5 else 1.5
    b <- b + array(rnorm(length(b), sd = 0.01), dim = c(4, 4, 1))
    s <- array(0.25, dim = c(4, 4, 1))
    write_vol(NeuroVol(b, space), file.path(tmpdir, sprintf("%s_beta.nii.gz", ids[i])))
    write_vol(NeuroVol(s, space), file.path(tmpdir, sprintf("%s_se.nii.gz", ids[i])))
    beta_paths[i] <- file.path(tmpdir, sprintf("%s_beta.nii.gz", ids[i]))
    se_paths[i]   <- file.path(tmpdir, sprintf("%s_se.nii.gz", ids[i]))
  }

  # Mask (all voxels)
  mask_path <- file.path(tmpdir, "mask.nii.gz")
  write_vol(NeuroVol(array(1, dim = c(4, 4, 1)), space), mask_path)

  gd <- group_data_from_nifti(
    beta_paths = beta_paths,
    se_paths   = se_paths,
    subjects   = ids,
    covariates = data.frame(group = grp),
    mask       = mask_path
  )

  # Meta engine with equal weights
  fit_eq <- fmri_ttest(gd, formula = ~ 1 + group, engine = "meta", weights = "equal")
  expect_s3_class(fit_eq, "fmri_ttest_fit")
  expect_true(all(dim(fit_eq$beta) == c(2, n_vox)))
  expect_true(all(dim(fit_eq$se) == c(2, n_vox)))
  expect_true(all(is.finite(fit_eq$z)))

  # Meta engine with custom weights
  w_subj <- rep(1, length(ids))
  fit_cu <- fmri_ttest(gd, formula = ~ 1 + group, engine = "meta",
                       weights = "custom", weights_custom = w_subj)
  expect_true(all(dim(fit_cu$beta) == c(2, n_vox)))
  expect_true(all(is.finite(fit_cu$z)))
})

