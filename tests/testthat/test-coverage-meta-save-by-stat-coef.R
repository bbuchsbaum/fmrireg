# Fifteenth wave: fmri_meta by_stat / by_coefficient NIfTI writers via group_data_nifti.

make_nifti_meta <- function(seed = 421L) {
  skip_if_not_installed("neuroim2")
  set.seed(seed)
  dims <- c(2L, 2L, 1L)
  nvox <- prod(dims)
  coefs <- matrix(
    rnorm(nvox * 2), nvox, 2,
    dimnames = list(NULL, c("Intercept", "age"))
  )
  se <- matrix(
    runif(nvox * 2, 0.1, 0.4), nvox, 2,
    dimnames = dimnames(coefs)
  )
  structure(
    list(
      coefficients = coefs,
      se = se,
      tau2 = runif(nvox),
      I2 = runif(nvox, 0, 40),
      Q = rchisq(nvox, df = 2),
      Q_df = rep(2, nvox),
      method = "DL",
      robust = "none",
      formula = ~ 1 + age,
      n_subjects = 6L,
      n_voxels = nvox,
      data = structure(
        list(dim = dims, voxel_size = c(1, 1, 1)),
        class = "group_data_nifti"
      )
    ),
    class = c("fmri_meta", "list")
  )
}

test_that(".fmri_meta_stat_image reconstructs NeuroVols for nifti group data", {
  meta <- make_nifti_meta()
  beta <- fmrireg:::.fmri_meta_stat_image(meta, 1L, "beta")
  expect_true(inherits(beta, "NeuroVol"))
  se_img <- fmrireg:::.fmri_meta_stat_image(meta, "age", "se")
  expect_true(inherits(se_img, "NeuroVol"))
  z_img <- fmrireg:::.fmri_meta_stat_image(meta, 2L, "z")
  expect_true(inherits(z_img, "NeuroVol"))

  expect_error(
    fmrireg:::.fmri_meta_stat_image(meta, 1L, "not_a_stat"),
    "Unsupported"
  )
})

test_that(".save_fmri_meta_by_stat writes nifti+json per statistic", {
  meta <- make_nifti_meta(seed = 422L)
  td <- tempfile("meta-by-stat-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ent <- list(subject = "01", task = "meta", space = "MNI")

  created <- fmrireg:::.save_fmri_meta_by_stat(
    x = meta,
    path = td,
    entities = ent,
    desc = "META",
    coefficient_indices = 1:2,
    coefficient_names = c("Intercept", "age"),
    coefficient_stats = c("beta", "se"),
    format = "nifti",
    output_formats = "nifti"
  )
  expect_true(is.list(created))
  expect_true(all(c("beta", "se") %in% names(created)))
  expect_true(file.exists(created$beta$nifti))
  expect_true(file.exists(created$beta$json))
  expect_true(file.exists(created$se$nifti))
  expect_true(file.exists(created$se$json))
  expect_match(basename(created$beta$nifti), "desc-beta")
  expect_match(basename(created$se$nifti), "desc-se")
})

test_that(".save_fmri_meta_by_coefficient writes nifti+json per coefficient", {
  meta <- make_nifti_meta(seed = 423L)
  td <- tempfile("meta-by-coef-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ent <- list(subject = "02", task = "meta", space = "T1w")

  created <- fmrireg:::.save_fmri_meta_by_coefficient(
    x = meta,
    path = td,
    entities = ent,
    desc = "META",
    coefficient_indices = 1:2,
    coefficient_names = c("Intercept", "age"),
    coefficient_stats = c("beta", "se"),
    format = "nifti",
    output_formats = "nifti"
  )
  expect_true(is.list(created))
  expect_true(all(c("Intercept", "age") %in% names(created)))
  expect_true(file.exists(created$Intercept$nifti))
  expect_true(file.exists(created$Intercept$json))
  expect_true(file.exists(created$age$nifti))
  expect_match(basename(created$age$nifti), "contrast-age|age")
})
