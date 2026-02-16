test_that("fmri_meta.gds equals legacy NIfTI path (FE/PM)", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("fmrigds")

  # Create tiny synthetic NIfTI volumes for 3 subjects
  set.seed(123)
  dims <- c(3,3,3)
  P <- prod(dims)
  S <- 3
  subjects <- paste0("sub-", 1:S)

  tmpdir <- tempdir()
  beta_paths <- character(S)
  se_paths   <- character(S)

  for (i in seq_len(S)) {
    beta_vol <- array(rnorm(P, mean = 0.1 * i, sd = 0.01), dim = dims)
    se_vol   <- array(0.1 + 0.01 * i, dim = dims)
    bp <- file.path(tmpdir, sprintf("beta_s%02d.nii", i))
    sp <- file.path(tmpdir, sprintf("se_s%02d.nii", i))
    RNifti::writeNifti(beta_vol, bp)
    RNifti::writeNifti(se_vol, sp)
    beta_paths[i] <- bp; se_paths[i] <- sp
  }

  # Mask (all ones)
  mask <- array(1L, dim = dims)
  mask_path <- file.path(tmpdir, "mask.nii")
  RNifti::writeNifti(mask, mask_path)

  # Legacy path
  gd_legacy <- group_data_from_nifti(beta_paths = beta_paths,
                                     se_paths   = se_paths,
                                     subjects   = subjects,
                                     mask       = mask_path,
                                     validate   = TRUE)

  # gds path via group_data (handles list structure)
  gd_gds <- group_data(data = list(beta = beta_paths, se = se_paths),
                       format = "nifti",
                       subjects = subjects,
                       mask = mask_path)

  # FE method
  fit_old_fe <- fmri_meta(gd_legacy, formula = ~ 1, method = "fe", robust = "none", verbose = FALSE)
  fit_new_fe <- fmri_meta(gd_gds,    formula = ~ 1, method = "fe", robust = "none", verbose = FALSE)

  expect_equal(dim(fit_new_fe$coefficients), dim(fit_old_fe$coefficients))
  expect_equal(dim(fit_new_fe$se),           dim(fit_old_fe$se))
  expect_equal(fit_new_fe$coefficients, fit_old_fe$coefficients, tolerance = 1e-8)
  expect_equal(fit_new_fe$se,           fit_old_fe$se,           tolerance = 1e-6)

  # PM method
  fit_old_pm <- fmri_meta(gd_legacy, formula = ~ 1, method = "pm", robust = "none", verbose = FALSE)
  fit_new_pm <- fmri_meta(gd_gds,    formula = ~ 1, method = "pm", robust = "none", verbose = FALSE)

  expect_equal(dim(fit_new_pm$coefficients), dim(fit_old_pm$coefficients))
  expect_equal(dim(fit_new_pm$se),           dim(fit_old_pm$se))
  expect_equal(fit_new_pm$coefficients, fit_old_pm$coefficients, tolerance = 1e-8)
  expect_equal(fit_new_pm$se,           fit_old_pm$se,           tolerance = 1e-6)
  # tau2 parity for PM
  expect_equal(fit_new_pm$tau2, fit_old_pm$tau2, tolerance = 1e-6)
})
