library(testthat)
library(fmrireg)

make_test_fmri_meta <- function() {
  dims <- c(2L, 2L, 1L)
  n_vox <- prod(dims)
  coef_mat <- cbind(
    "(Intercept)" = seq_len(n_vox) / 10,
    groupB = seq_len(n_vox) / 20
  )
  se_mat <- matrix(0.1, nrow = n_vox, ncol = 2L)
  colnames(se_mat) <- colnames(coef_mat)

  structure(
    list(
      coefficients = coef_mat,
      se = se_mat,
      tau2 = rep(0.01, n_vox),
      I2 = rep(0.2, n_vox),
      Q = rep(3, n_vox),
      Q_df = rep(2, n_vox),
      model = list(X = cbind("(Intercept)" = rep(1, 4))),
      method = "fe",
      robust = "none",
      weights = "ivw",
      data = structure(
        list(
          dim = dims,
          voxel_size = c(2, 2, 2),
          format = "nifti",
          n_subjects = 4L,
          n_voxels = n_vox
        ),
        class = c("group_data_nifti", "group_data")
      ),
      formula = ~ 1,
      n_voxels = n_vox,
      n_subjects = 4L,
      voxel_indices = seq_len(n_vox)
    ),
    class = "fmri_meta"
  )
}

test_that("write_results.fmri_meta writes regular by-coefficient h5 and nifti outputs", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("neuroim2")

  fit <- make_test_fmri_meta()
  out_dir <- tempfile()
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  result <- write_results(
    fit,
    path = out_dir,
    task = "meta",
    space = "MNI152",
    desc = "Meta",
    format = c("h5", "nifti"),
    strategy = "by_coefficient",
    coefficients = "(Intercept)",
    coefficient_stats = c("beta", "se", "z", "pval")
  )

  expect_named(result, c("(Intercept)", "heterogeneity"))
  expect_true(file.exists(result[["(Intercept)"]]$h5))
  expect_true(file.exists(result[["(Intercept)"]]$nifti))
  expect_true(file.exists(result[["(Intercept)"]]$json))
  expect_true(file.exists(result$heterogeneity$h5))
  expect_true(file.exists(result$heterogeneity$nifti))
  expect_true(file.exists(result$heterogeneity$json))

  expect_match(
    basename(result[["(Intercept)"]]$h5),
    "task-meta_space-MNI152_contrast-Intercept_desc-Meta_bold\\.h5$"
  )
  expect_match(
    basename(result[["(Intercept)"]]$nifti),
    "task-meta_space-MNI152_contrast-Intercept_desc-Meta_bold\\.nii\\.gz$"
  )

  h5_meta <- read_h5_metadata(result[["(Intercept)"]]$h5)
  expect_equal(h5_meta$labels, c("beta", "se", "z", "pval"))

  json_meta <- jsonlite::read_json(result[["(Intercept)"]]$json)
  expect_equal(unlist(json_meta$StatisticOrder), c("beta", "se", "z", "pval"))
  expect_equal(unlist(json_meta$CoefficientOrder), "(Intercept)")

  expect_error(
    write_results(
      fit,
      path = out_dir,
      task = "meta",
      space = "MNI152",
      desc = "Meta",
      format = c("h5", "nifti"),
      strategy = "by_coefficient",
      coefficients = "(Intercept)",
      coefficient_stats = c("beta", "se", "z", "pval")
    ),
    "Output files already exist"
  )
})

test_that("write_results.fmri_meta writes by-stat outputs with coefficient labels", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("neuroim2")

  fit <- make_test_fmri_meta()
  out_dir <- tempfile()
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  result <- write_results(
    fit,
    path = out_dir,
    task = "meta",
    desc = "Meta",
    format = "nifti",
    strategy = "by_stat",
    coefficient_stats = c("beta", "z"),
    heterogeneity = FALSE
  )

  expect_named(result, c("beta", "z"))
  expect_true(file.exists(result$beta$nifti))
  expect_true(file.exists(result$beta$json))
  expect_false("h5" %in% names(result$beta))
  expect_match(basename(result$beta$nifti), "task-meta_desc-beta_bold\\.nii\\.gz$")

  json_meta <- jsonlite::read_json(result$beta$json)
  expect_equal(unlist(json_meta$CoefficientOrder), c("(Intercept)", "groupB"))
  expect_equal(unlist(json_meta$StatisticOrder), "beta")
})

test_that("coef_image.fmri_meta reconstructs voxel-space group_data_gds results", {
  skip_if_not_installed("fmrigds")
  skip_if_not_installed("neuroim2")

  fit <- make_test_fmri_meta()
  dims <- c(2L, 2L, 1L)
  gds <- fmrigds::new_gds(
    assays = list(
      beta = array(seq_len(prod(dims)), dim = c(prod(dims), 1L, 1L)),
      se = array(1, dim = c(prod(dims), 1L, 1L))
    ),
    space = fmrigds::space_voxel(dims, diag(4), mask_bitmap = array(TRUE, dim = dims)),
    subjects = "s1",
    contrasts = "c1"
  )
  class(gds) <- c("group_data_gds", "group_data", class(gds))
  fit$data <- gds

  img <- coef_image(fit, coef = "(Intercept)", statistic = "estimate")
  expect_s4_class(img, "NeuroVol")
  expect_equal(dim(img), dims)
})
