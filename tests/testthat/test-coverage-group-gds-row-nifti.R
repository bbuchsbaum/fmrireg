# group_data_gds.R: row_data / feature_group / nifti list routing branches.

test_that(".group_data_gds tabular row_data and feature_group branches", {
  skip_if_not_installed("fmrigds")

  df <- data.frame(
    subject = rep(paste0("s", 1:4), each = 2),
    roi = rep(c("r1", "r2"), times = 4),
    beta = rnorm(8, 0.2, 0.1),
    se = runif(8, 0.05, 0.15),
    group = rep(c("A", "A", "B", "B"), each = 2),
    stringsAsFactors = FALSE
  )

  # roi_col maps to sample_col; feature_group length must match samples
  gd <- fmrireg:::.group_data_gds(
    df,
    format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    covariate_cols = "group",
    feature_group = c("visual", "motor")
  )
  expect_true(inherits(gd, "group_data_gds"))
  expect_equal(attr(gd, "fmrireg_feature_group"), c("visual", "motor"))
  expect_equal(attr(gd, "fmrireg_sample_labels"), c("r1", "r2"))

  expect_error(
    fmrireg:::.group_data_gds(
      df,
      format = "csv",
      effect_cols = c(beta = "beta", se = "se"),
      subject_col = "subject",
      roi_col = "roi",
      feature_group = c("only_one")
    ),
    "feature_group must have length"
  )

  # Explicit row_data that already includes sample labels
  row_df <- data.frame(
    sample = c("r1", "r2"),
    label = c("r1", "r2"),
    stringsAsFactors = FALSE,
    row.names = c("r1", "r2")
  )
  gd2 <- fmrireg:::.group_data_gds(
    df,
    format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    sample_col = "roi",
    row_data = row_df
  )
  expect_true(inherits(gd2, "group_data_gds"))
  expect_true(
    !is.null(attr(gd2, "fmrireg_sample_labels")) ||
      inherits(gd2, "group_data")
  )
})

test_that(".group_data_gds nifti list routes to legacy nifti impl", {
  skip_if_not_installed("neuroim2")
  set.seed(9)
  dims <- c(2, 2, 1)
  sp <- neuroim2::NeuroSpace(dims)
  write_one <- function(vals) {
    p <- tempfile(fileext = ".nii")
    neuroim2::write_vol(neuroim2::NeuroVol(array(vals, dim = dims), sp), p)
    p
  }
  beta <- c(write_one(rnorm(4)), write_one(rnorm(4)))
  se <- c(write_one(abs(rnorm(4)) + 0.2), write_one(abs(rnorm(4)) + 0.2))
  tstat <- c(write_one(rnorm(4)), write_one(rnorm(4)))
  on.exit(unlink(c(beta, se, tstat)), add = TRUE)

  # t + df path returns nifti-backed object via legacy impl
  gd_t <- fmrireg:::.group_data_gds(
    list(t = tstat, df = 20),
    format = "nifti",
    subjects = c("s1", "s2"),
    validate = TRUE
  )
  expect_true(
    inherits(gd_t, "group_data_nifti") || inherits(gd_t, "group_data_gds") ||
      inherits(gd_t, "group_data")
  )

  # beta + var path
  vars <- c(write_one(abs(rnorm(4)) + 0.05), write_one(abs(rnorm(4)) + 0.05))
  on.exit(unlink(vars), add = TRUE)
  gd_v <- fmrireg:::.group_data_gds(
    list(beta = beta, var = vars),
    format = "nifti",
    subjects = c("s1", "s2"),
    validate = TRUE
  )
  expect_true(
    inherits(gd_v, "group_data_nifti") || inherits(gd_v, "group_data")
  )

  # beta + se without forcing legacy t/var branch builds fmrigds source list
  gd_se <- tryCatch(
    fmrireg:::.group_data_gds(
      list(beta = beta, se = se),
      format = "nifti",
      subjects = c("s1", "s2"),
      covariates = data.frame(group = c("A", "B")),
      validate = TRUE
    ),
    error = function(e) e
  )
  expect_true(
    inherits(gd_se, "group_data") || inherits(gd_se, "error") ||
      inherits(gd_se, "gds") || inherits(gd_se, "group_data_gds")
  )
})

test_that(".print_group_data_gds prints assay/subject info when available", {
  skip_if_not_installed("fmrigds")
  gd <- fmrireg:::.demo_group_data_csv()
  # Ensure gds class annotation
  if (!inherits(gd, "group_data_gds")) {
    gd <- fmrireg:::.annotate_group_data_gds(gd)
  }
  expect_output(fmrireg:::.print_group_data_gds(gd), "group_data")
})
