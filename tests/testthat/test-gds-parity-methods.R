test_that("fmri_meta.gds equals legacy CSV path across methods", {
  skip_if_not_installed("fmrigds")

  methods <- c("fe","dl","pm","reml")

  gd_csv <- fmrireg:::.demo_group_data_csv()
  df <- gd_csv$data
  gd_gds <- fmrireg::group_data(
    df, format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    contrast_col = "contrast",
    covariate_cols = "age"
  )

  for (m in methods) {
    fit_old <- fmri_meta(gd_csv, formula = ~ 1 + age, method = m, robust = "none", verbose = FALSE)
    fit_new <- fmri_meta(gd_gds, formula = ~ 1 + age, method = m, robust = "none", verbose = FALSE)

    expect_equal(dim(fit_new$coefficients), dim(fit_old$coefficients), info = m)
    expect_equal(dim(fit_new$se),           dim(fit_old$se),           info = m)
    expect_equal(colnames(fit_new$coefficients), colnames(fit_old$coefficients), info = m)

    expect_equal(fit_new$coefficients, fit_old$coefficients, tolerance = 1e-8, info = m)
    expect_equal(fit_new$se,           fit_old$se,           tolerance = 1e-6, info = m)

    # Compare heterogeneity/diagnostics where applicable (vectors per feature)
    if (m %in% c("dl","pm","reml")) {
      expect_equal(fit_new$tau2, fit_old$tau2, tolerance = 1e-6, info = paste(m, "tau2"))
    }
    # I2/Q exist for FE baseline as well; compare when available
    if (!is.null(fit_old$I2)) {
      expect_equal(fit_new$I2, fit_old$I2, tolerance = 1e-6, info = paste(m, "I2"))
    }
    if (!is.null(fit_old$Q)) {
      expect_equal(fit_new$Q, fit_old$Q, tolerance = 1e-6, info = paste(m, "Q"))
    }
  }
})

test_that("fmri_meta.gds returns cov_tri with correct size for multi-coef", {
  skip_if_not_installed("fmrigds")
  gd_csv <- fmrireg:::.demo_group_data_csv()
  df <- gd_csv$data
  gd_gds <- fmrireg::group_data(
    df, format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    contrast_col = "contrast",
    covariate_cols = "age"
  )
  fit <- fmri_meta(gd_gds, formula = ~ 1 + age, method = "pm", return_cov = "tri", verbose = FALSE)
  expect_true(!is.null(fit$cov))
  K <- ncol(fit$coefficients)
  tsize <- K*(K+1)/2
  expect_equal(nrow(fit$cov$tri), tsize)
})

test_that("fmri_meta.gds custom unit weights match equal weights", {
  skip_if_not_installed("fmrigds")
  gd_csv <- fmrireg:::.demo_group_data_csv()
  df <- gd_csv$data
  gd_gds <- fmrireg::group_data(
    df, format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi"
  )
  n_subjects <- length(unique(df$subject))
  fit_equal <- fmri_meta(
    gd_gds, weights = "equal", method = "fe", verbose = FALSE
  )
  fit_custom <- fmri_meta(
    gd_gds, weights = "custom", weights_custom = rep(1, n_subjects),
    method = "fe", verbose = FALSE
  )
  expect_identical(fit_custom$weights, "custom")
  expect_equal(fit_custom$coefficients, fit_equal$coefficients, tolerance = 1e-10)
  expect_equal(fit_custom$se, fit_equal$se, tolerance = 1e-10)
  expect_error(
    fmri_meta(gd_gds, weights = "custom", weights_custom = rep(1, n_subjects + 1L)),
    "weights_custom must be length S or SxP"
  )
})

test_that("fmri_ttest.gds custom weights match the declared meta estimator", {
  skip_if_not_installed("fmrigds")
  gd_csv <- fmrireg:::.demo_group_data_csv()
  df <- gd_csv$data
  gd_gds <- fmrireg::group_data(
    df, format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi"
  )
  n_subjects <- length(unique(df$subject))
  fit_meta <- fmri_meta(
    gd_gds, method = getOption("fmrireg.meta.method", "pm"),
    weights = "custom", weights_custom = rep(1, n_subjects),
    verbose = FALSE
  )
  fit_custom <- fmri_ttest(
    gd_gds, engine = "meta", weights = "custom",
    weights_custom = rep(1, n_subjects)
  )
  expect_identical(fit_custom$weights, "custom")
  expect_identical(
    fit_custom$method, getOption("fmrireg.meta.method", "pm")
  )
  expect_equal(fit_custom$beta, t(fit_meta$coefficients), tolerance = 1e-10)
  expect_equal(fit_custom$se, t(fit_meta$se), tolerance = 1e-10)
  expect_equal(
    fit_custom$z,
    t(fit_meta$coefficients / fit_meta$se),
    tolerance = 1e-10
  )
})
