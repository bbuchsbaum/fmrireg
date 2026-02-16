test_that("fmri_meta.gds equals legacy HDF5 path (FE/PM)", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmrigds")

  set.seed(123)
  dims <- c(3,3,3); P <- prod(dims); S <- 3
  subjects <- paste0("sub-", 1:S)

  tmpdir <- tempdir()
  h5_paths <- character(S)

  # Create a simple logical mask and space
  mask <- array(TRUE, dim = dims)
  space <- neuroim2::NeuroSpace(dims)
  mask_vol <- neuroim2::LogicalNeuroVol(mask, space)

  for (i in seq_len(S)) {
    # 4D array: last dim = labels (beta, se)
    beta_vol <- array(rnorm(P, mean = 0.05 * i, sd = 0.005), dim = dims)
    se_vol   <- array(0.1 + 0.01 * i, dim = dims)
    arr <- array(0, dim = c(dims, 2))
    arr[,,,1] <- beta_vol
    arr[,,,2] <- se_vol
    labels <- c("beta","se")
    h5p <- file.path(tmpdir, sprintf("sub%02d_statmaps.h5", i))
    ospace <- neuroim2::add_dim(space, length(labels))
    vec <- neuroim2::NeuroVec(arr, ospace)
    h <- fmristore::write_labeled_vec(vec, mask_vol, labels, file = h5p)
    h$close_all()
    h5_paths[i] <- h5p
  }

  # Skip if legacy metadata reader cannot handle current fmristore handle
  meta_try <- try(fmrireg:::read_h5_metadata(h5_paths[1]), silent = TRUE)
  if (inherits(meta_try, "try-error")) {
    skip("fmristore handle incompatible with legacy readers; skipping HDF5 parity test")
  }

  # Legacy path
  gd_legacy <- group_data_from_h5(paths = h5_paths,
                                  subjects = subjects,
                                  stat = c("beta","se"),
                                  validate = TRUE)

  # gds path via fmrigds
  gd_gds <- fmrigds::gds(data = h5_paths, format = "h5", subjects = subjects)

  # FE
  fit_old_fe <- fmri_meta(gd_legacy, formula = ~ 1, method = "fe", robust = "none", verbose = FALSE)
  fit_new_fe <- fmri_meta(gd_gds,    formula = ~ 1, method = "fe", robust = "none", verbose = FALSE)
  expect_equal(dim(fit_new_fe$coefficients), dim(fit_old_fe$coefficients))
  expect_equal(dim(fit_new_fe$se),           dim(fit_old_fe$se))
  expect_equal(fit_new_fe$coefficients, fit_old_fe$coefficients, tolerance = 1e-8)
  expect_equal(fit_new_fe$se,           fit_old_fe$se,           tolerance = 1e-6)

  # PM
  fit_old_pm <- fmri_meta(gd_legacy, formula = ~ 1, method = "pm", robust = "none", verbose = FALSE)
  fit_new_pm <- fmri_meta(gd_gds,    formula = ~ 1, method = "pm", robust = "none", verbose = FALSE)
  expect_equal(fit_new_pm$coefficients, fit_old_pm$coefficients, tolerance = 1e-8)
  expect_equal(fit_new_pm$se,           fit_old_pm$se,           tolerance = 1e-6)
  expect_equal(fit_new_pm$tau2,         fit_old_pm$tau2,         tolerance = 1e-6)
})
