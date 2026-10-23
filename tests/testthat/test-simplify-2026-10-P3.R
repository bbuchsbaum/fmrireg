# Simplification 2026-10, step P3 (routing): group_data(format = "nifti")
# reads variance maps and t-only maps through fmrigds, like beta/se maps, so
# the legacy group_data_nifti reader is reachable only from the deprecated
# group_data_from_nifti(). Results equal the legacy reader's.

skip_if_not_installed("fmrigds")
skip_if_not_installed("neuroim2")

p3_maps <- function(dir, S = 8L, dims = c(3L, 3L, 2L), seed = 41) {
  set.seed(seed)
  sp <- neuroim2::NeuroSpace(dims, spacing = c(2, 2, 2))
  grp <- factor(rep(c("A", "B"), each = S / 2))
  write <- function(kind, f) vapply(seq_len(S), function(i) {
    p <- file.path(dir, sprintf("sub-%02d_%s.nii.gz", i, kind))
    neuroim2::write_vol(neuroim2::NeuroVol(array(f(i), dims), sp), p)
    p
  }, "")
  list(
    beta = write("beta", function(i) stats::rnorm(prod(dims), 0.3 + 0.5 * (grp[i] == "B"))),
    var = write("var", function(i) stats::runif(prod(dims), 0.05, 0.5)),
    t = write("stat-t", function(i) stats::rnorm(prod(dims), 1)),
    subjects = sprintf("s%02d", seq_len(S)),
    covariates = data.frame(group = grp)
  )
}

test_that("P3: beta + var and t-only NIfTI maps are read by fmrigds", {
  m <- p3_maps(withr::local_tempdir())
  gd <- group_data(list(beta = m$beta, var = m$var), format = "nifti",
                   subjects = m$subjects, covariates = m$covariates)
  expect_s3_class(gd, "group_data_gds")
  df <- 15 + seq_along(m$t)
  gd_t <- group_data(list(t = m$t, df = df), format = "nifti", subjects = m$subjects)
  expect_s3_class(gd_t, "group_data_gds")

  # The default weights = "ivw" cannot use t statistics; like the legacy
  # reader, the combination warns and uses equal weights.
  expect_warning(
    fit_ivw <- fmri_meta(gd_t, ~ 1, combine = "stouffer", verbose = FALSE),
    "no SE available; using equal weights"
  )
  fit_eq <- fmri_meta(gd_t, ~ 1, combine = "stouffer", weights = "equal", verbose = FALSE)
  expect_equal(fit_ivw$coefficients, fit_eq$coefficients)
})

# LEGACY-GROUP-DELETE (next minor): parity test against the legacy reader (group_data_from_nifti() + fmri_meta.group_data_nifti()).
test_that("P3: beta + var NIfTI maps are read by fmrigds and fit as before", {
  m <- p3_maps(withr::local_tempdir())
  gd <- group_data(list(beta = m$beta, var = m$var), format = "nifti",
                   subjects = m$subjects, covariates = m$covariates)
  expect_s3_class(gd, "group_data_gds")
  legacy <- suppressWarnings(group_data_from_nifti(
    beta_paths = m$beta, var_paths = m$var, subjects = m$subjects,
    covariates = m$covariates
  ))
  for (method in c("fe", "pm")) {
    fit <- fmri_meta(gd, ~ 1 + group, method = method, verbose = FALSE)
    ref <- suppressWarnings(fmri_meta(legacy, ~ 1 + group, method = method, verbose = FALSE))
    expect_equal(unname(fit$coefficients), unname(ref$coefficients), tolerance = 1e-10)
    expect_equal(unname(fit$se), unname(ref$se), tolerance = 1e-10)
  }
})

# LEGACY-GROUP-DELETE (next minor): parity test against the legacy reader (group_data_from_nifti() + fmri_meta.group_data_nifti()); its gds-only assertions are duplicated above.
test_that("P3: t-only NIfTI maps are read by fmrigds and combine as before", {
  m <- p3_maps(withr::local_tempdir())
  df <- 15 + seq_along(m$t)
  gd <- group_data(list(t = m$t, df = df), format = "nifti", subjects = m$subjects)
  expect_s3_class(gd, "group_data_gds")
  legacy <- suppressWarnings(group_data_from_nifti(t_paths = m$t, df = df, subjects = m$subjects))
  for (cm in c("stouffer", "lancaster")) {
    fit <- fmri_meta(gd, ~ 1, combine = cm, weights = "equal", verbose = FALSE)
    ref <- suppressWarnings(fmri_meta(legacy, ~ 1, combine = cm, weights = "equal", verbose = FALSE))
    expect_equal(unname(fit$coefficients), unname(ref$coefficients), tolerance = 1e-10)
  }

  # The default weights = "ivw" cannot use t statistics; like the legacy
  # reader, the combination warns and uses equal weights.
  expect_warning(
    fit_ivw <- fmri_meta(gd, ~ 1, combine = "stouffer", verbose = FALSE),
    "no SE available; using equal weights"
  )
  fit_eq <- fmri_meta(gd, ~ 1, combine = "stouffer", weights = "equal", verbose = FALSE)
  expect_equal(fit_ivw$coefficients, fit_eq$coefficients)
})

# LEGACY-GROUP-DELETE (next minor): tests the legacy group_data_nifti -> gds conversion in .group_data_gds() and parity with fmri_meta.group_data_nifti().
test_that("P3: a legacy group_data_nifti object converts to gds with the same fit", {
  m <- p3_maps(withr::local_tempdir())
  legacy_var <- suppressWarnings(group_data_from_nifti(
    beta_paths = m$beta, var_paths = m$var, subjects = m$subjects
  ))
  gd_var <- group_data(legacy_var, format = "nifti", subjects = m$subjects)
  expect_s3_class(gd_var, "group_data_gds")
  expect_equal(
    unname(fmri_meta(gd_var, ~ 1, verbose = FALSE)$coefficients),
    unname(suppressWarnings(fmri_meta(legacy_var, ~ 1, verbose = FALSE))$coefficients),
    tolerance = 1e-10
  )

  legacy_t <- suppressWarnings(group_data_from_nifti(t_paths = m$t, df = 20, subjects = m$subjects))
  gd_t <- group_data(legacy_t, format = "nifti", subjects = m$subjects)
  expect_s3_class(gd_t, "group_data_gds")
  expect_equal(
    unname(fmri_meta(gd_t, ~ 1, combine = "stouffer", weights = "equal", verbose = FALSE)$coefficients),
    unname(suppressWarnings(fmri_meta(legacy_t, ~ 1, combine = "stouffer", weights = "equal",
                                      verbose = FALSE))$coefficients),
    tolerance = 1e-10
  )
})
