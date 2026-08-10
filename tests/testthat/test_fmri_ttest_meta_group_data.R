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
  t_paths    <- character(length(ids))

  set.seed(1)
  for (i in seq_along(ids)) {
    b <- array(0, dim = c(4, 4, 1))
    b[] <- if (grp[i] == "A") 0.5 else 1.5
    b <- b + array(rnorm(length(b), sd = 0.01), dim = c(4, 4, 1))
    s <- array(0.25, dim = c(4, 4, 1))
    write_vol(NeuroVol(b, space), file.path(tmpdir, sprintf("%s_beta.nii.gz", ids[i])))
    write_vol(NeuroVol(s, space), file.path(tmpdir, sprintf("%s_se.nii.gz", ids[i])))
    write_vol(NeuroVol(b / s, space), file.path(tmpdir, sprintf("%s_t.nii.gz", ids[i])))
    beta_paths[i] <- file.path(tmpdir, sprintf("%s_beta.nii.gz", ids[i]))
    se_paths[i]   <- file.path(tmpdir, sprintf("%s_se.nii.gz", ids[i]))
    t_paths[i]    <- file.path(tmpdir, sprintf("%s_t.nii.gz", ids[i]))
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
  expect_identical(fit_eq$method, "pm")
  expect_identical(fit_eq$weights, "equal")
  expect_null(fit_eq$combine)
  expect_equal(fit_eq$n_features, n_vox)

  # Meta engine with custom weights
  w_subj <- rep(1, length(ids))
  fit_cu <- fmri_ttest(gd, formula = ~ 1 + group, engine = "meta",
                       weights = "custom", weights_custom = w_subj)
  expect_true(all(dim(fit_cu$beta) == c(2, n_vox)))
  expect_true(all(is.finite(fit_cu$z)))

  fit_cov <- fmri_meta(
    gd, formula = ~ 1 + group, method = "pm", return_cov = "tri",
    verbose = FALSE
  )
  group_coef <- setdiff(colnames(fit_cov$coefficients), "(Intercept)")
  expect_identical(group_coef, "groupB")
  group_weights <- stats::setNames(
    as.numeric(colnames(fit_cov$coefficients) == group_coef),
    colnames(fit_cov$coefficients)
  )
  group_contrast <- contrast(fit_cov, group_weights)
  expect_equal(mean(group_contrast$estimate), 1, tolerance = 0.03)
  expect_true(all(is.finite(group_contrast$se)))
  expect_error(
    contrast(fit_cov, c(group = 1)),
    "Coefficient\\(s\\) not found.*groupB"
  )

  gd_t <- group_data_from_nifti(
    t_paths = t_paths,
    df = 60,
    subjects = ids,
    mask = mask_path
  )
  fit_st <- fmri_meta(
    gd_t, formula = ~ 1, combine = "stouffer", weights = "equal",
    verbose = FALSE
  )
  expect_identical(fit_st$method, "combine:stouffer")
  expect_identical(fit_st$weights, "equal")

  had_global_out <- exists("out", envir = .GlobalEnv, inherits = FALSE)
  if (had_global_out) old_global_out <- get("out", envir = .GlobalEnv)
  on.exit({
    if (had_global_out) {
      assign("out", old_global_out, envir = .GlobalEnv)
    } else if (exists("out", envir = .GlobalEnv, inherits = FALSE)) {
      rm("out", envir = .GlobalEnv)
    }
  }, add = TRUE)
  assign("out", "unrelated render output path", envir = .GlobalEnv)

  fit_tt_la <- fmri_ttest(
    gd_t, formula = ~ 1, engine = "meta", combine = "lancaster",
    weights = "custom", weights_custom = w_subj
  )
  expect_identical(fit_tt_la$method, "combine:lancaster")
  expect_identical(fit_tt_la$combine, "lancaster")
  expect_identical(fit_tt_la$weights, "custom")
  expect_equal(fit_tt_la$n_features, n_vox)
})

test_that("public ROI t-only fits return the requested inferential result", {
  ids <- sprintf("s%02d", 1:8)
  roi_t <- data.frame(
    subject = ids,
    roi = "ExampleROI",
    t = seq(1.5, 2.2, length.out = length(ids)),
    df = 40
  )
  expect_no_warning(gd <- group_data(
    roi_t, format = "csv",
    effect_cols = c(t = "t", df = "df"),
    subject_col = "subject",
    roi_col = "roi"
  ))

  fit_st <- fmri_meta(
    gd, formula = ~ 1, combine = "stouffer", weights = "equal",
    verbose = FALSE
  )
  fit_la <- fmri_meta(
    gd, formula = ~ 1, combine = "lancaster", weights = "custom",
    weights_custom = rep(1, length(ids)), verbose = FALSE
  )
  fit_tt <- fmri_ttest(
    gd, formula = ~ 1, engine = "meta", combine = "stouffer",
    weights = "equal"
  )
  stouffer_oracle <- combine_t_statistics(
    matrix(roi_t$t, ncol = 1L), df = roi_t$df,
    method = "stouffer", weights = "equal"
  )[[1L]]

  expect_s3_class(gd, "group_data_gds")
  expect_identical(fit_st$method, "combine:stouffer")
  expect_identical(fit_la$method, "combine:lancaster")
  expect_equal(as.numeric(fit_st$coefficients), stouffer_oracle, tolerance = 1e-12)
  expect_equal(as.numeric(fit_st$se), 1)
  expect_equal(as.numeric(fit_tt$z), stouffer_oracle, tolerance = 1e-12)
  expect_true(all(is.finite(fit_tt$p)))
})

test_that("fmri_ttest classic and welch materialize group_data_nifti betas", {
  skip_on_cran()
  skip_if_not_installed("neuroim2")

  library(neuroim2)

  space <- NeuroSpace(c(3, 3, 1), spacing = c(2, 2, 2))
  n_vox <- prod(dim(space))
  ids <- sprintf("sub-%02d", 1:6)
  grp <- factor(rep(c("A", "B"), each = 3))

  tmpdir <- tempdir()
  beta_paths <- character(length(ids))
  se_paths <- character(length(ids))

  set.seed(2)
  for (i in seq_along(ids)) {
    b <- array(rnorm(n_vox, sd = 0.02), dim = c(3, 3, 1))
    b <- b + if (grp[i] == "A") 1 else 2
    s <- array(0.2, dim = c(3, 3, 1))
    beta_paths[i] <- file.path(tmpdir, sprintf("%s_classic_beta.nii.gz", ids[i]))
    se_paths[i] <- file.path(tmpdir, sprintf("%s_classic_se.nii.gz", ids[i]))
    write_vol(NeuroVol(b, space), beta_paths[i])
    write_vol(NeuroVol(s, space), se_paths[i])
  }

  mask_path <- file.path(tmpdir, "classic_mask.nii.gz")
  write_vol(NeuroVol(array(1, dim = c(3, 3, 1)), space), mask_path)

  gd <- group_data_from_nifti(
    beta_paths = beta_paths,
    se_paths = se_paths,
    subjects = ids,
    covariates = data.frame(group = grp),
    mask = mask_path
  )

  fit_classic <- fmri_ttest(gd, formula = ~ 1 + group, engine = "classic")
  fit_welch <- fmri_ttest(gd, formula = ~ 1 + group, engine = "welch")

  expect_equal(dim(fit_classic$beta), c(2, n_vox))
  expect_equal(dim(fit_welch$beta), c(2, n_vox))
  expect_equal(rownames(fit_classic$beta), c("(Intercept)", "group"))
  expect_equal(rownames(fit_welch$beta), c("(Intercept)", "group"))
  expect_true(mean(fit_classic$beta["group", ]) < 0)
  expect_true(mean(fit_welch$beta["group", ]) < 0)
})
