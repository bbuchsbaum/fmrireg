# Pass 3, stream R2-4: write_results.fmri_meta keeps the user's `desc` in the
# per-statistic desc- label (desc-<desc><stat>, sanitised to alphanumerics),
# as pass 1 (E6) did for fmri_lm outputs. The default desc ("Meta") keeps the
# historical names (desc-<stat>).

make_meta_r24 <- function() {
  dims <- c(2L, 2L, 1L)
  n_vox <- prod(dims)
  coef_mat <- cbind("(Intercept)" = seq_len(n_vox) / 10, groupB = seq_len(n_vox) / 20)
  se_mat <- matrix(0.1, nrow = n_vox, ncol = 2L, dimnames = list(NULL, colnames(coef_mat)))
  structure(
    list(
      coefficients = coef_mat, se = se_mat,
      tau2 = rep(0.01, n_vox), I2 = rep(0.2, n_vox), Q = rep(3, n_vox), Q_df = rep(2, n_vox),
      model = list(X = cbind("(Intercept)" = rep(1, 4))),
      method = "fe", robust = "none", weights = "ivw",
      data = structure(
        list(dim = dims, voxel_size = c(2, 2, 2), format = "nifti",
             n_subjects = 4L, n_voxels = n_vox),
        class = c("group_data_nifti", "group_data")
      ),
      formula = ~ 1, n_voxels = n_vox, n_subjects = 4L, voxel_indices = seq_len(n_vox)
    ),
    class = "fmri_meta"
  )
}

test_that("R2-4: by_stat filenames keep a custom desc (sanitised) and default names stay", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("neuroim2")
  fit <- make_meta_r24()
  out_dir <- tempfile()
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  res_a <- write_results(fit, path = out_dir, task = "meta", desc = "run-A",
                         format = "nifti", strategy = "by_stat",
                         coefficient_stats = c("beta", "z"), heterogeneity = TRUE)
  expect_match(basename(res_a$beta$nifti), "^task-meta_desc-runAbeta_bold\\.nii\\.gz$")
  expect_match(basename(res_a$z$nifti), "^task-meta_desc-runAz_bold\\.nii\\.gz$")
  expect_match(basename(res_a$heterogeneity$nifti),
               "^task-meta_desc-runAheterogeneity_bold\\.nii\\.gz$")

  # A second analysis label into the same directory does not collide
  res_b <- write_results(fit, path = out_dir, task = "meta", desc = "run_B",
                         format = "nifti", strategy = "by_stat",
                         coefficient_stats = c("beta", "z"), heterogeneity = TRUE)
  expect_match(basename(res_b$beta$nifti), "^task-meta_desc-runBbeta_bold\\.nii\\.gz$")
  expect_true(file.exists(res_a$beta$nifti))

  # The default desc keeps the historical names
  res_d <- write_results(fit, path = out_dir, task = "meta",
                         format = "nifti", strategy = "by_stat",
                         coefficient_stats = "beta", heterogeneity = TRUE)
  expect_match(basename(res_d$beta$nifti), "^task-meta_desc-beta_bold\\.nii\\.gz$")
  expect_match(basename(res_d$heterogeneity$nifti),
               "^task-meta_desc-heterogeneity_bold\\.nii\\.gz$")

  # Predicted files match written files: rewriting without overwrite refuses
  expect_error(
    write_results(fit, path = out_dir, task = "meta", desc = "run-A",
                  format = "nifti", strategy = "by_stat",
                  coefficient_stats = c("beta", "z"), heterogeneity = TRUE),
    "already exist"
  )
})

test_that("R2-4: by_coefficient with one statistic keeps the custom desc", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("neuroim2")
  fit <- make_meta_r24()
  out_dir <- tempfile()
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  res <- write_results(fit, path = out_dir, task = "meta", desc = "Study2",
                       format = "nifti", strategy = "by_coefficient",
                       coefficients = "groupB", coefficient_stats = "z",
                       heterogeneity = FALSE)
  expect_match(basename(res$groupB$nifti),
               "^task-meta_contrast-groupB_desc-Study2z_bold\\.nii\\.gz$")

  res_d <- write_results(fit, path = out_dir, task = "meta",
                         format = "nifti", strategy = "by_coefficient",
                         coefficients = "groupB", coefficient_stats = "z",
                         heterogeneity = FALSE)
  expect_match(basename(res_d$groupB$nifti),
               "^task-meta_contrast-groupB_desc-z_bold\\.nii\\.gz$")
})

test_that("R2-4: meta desc helper", {
  expect_identical(fmrireg:::.meta_stat_desc("Meta", "beta"), "beta")
  expect_identical(fmrireg:::.meta_stat_desc(NULL, "pval"), "pval")
  expect_identical(fmrireg:::.meta_stat_desc("my_study-1", "se"), "mystudy1se")
})
