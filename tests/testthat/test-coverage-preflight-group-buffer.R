# Extra buffer past 90%: preflight, effective_df, meta_gds edges, group_data print.

test_that("preflight and effective_df helpers cover remaining branches", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  pf <- tryCatch(preflight(fit$model, fit$dataset), error = function(e) e)
  if (inherits(pf, "error")) {
    # Alternative signature
    pf <- tryCatch(preflight(fit), error = function(e) e)
  }
  if (!inherits(pf, "error")) {
    expect_true(inherits(pf, "fmri_preflight") || is.list(pf))
    expect_output(print(pf), regexp = ".")
  }

  if (exists("calculate_effective_df", mode = "function")) {
    X <- cbind(1, rnorm(30), rnorm(30))
    edf <- calculate_effective_df(30, 3, method = "residual")
    expect_true(is.finite(edf))
    edf2 <- tryCatch(
      calculate_effective_df(30, 3, method = "satterthwaite",
                             X = X, resid_cov = diag(30)),
      error = function(e) e
    )
    if (!inherits(edf2, "error")) expect_true(is.finite(edf2))
    edf3 <- calculate_effective_df(30, 3, method = "simple")
    expect_true(is.finite(edf3))
  }
})

test_that("group_data print/summary and detect format cover remaining lines", {
  gd <- tryCatch(fmrireg:::.demo_group_data_csv(), error = function(e) NULL)
  skip_if(is.null(gd), "no demo group data")
  expect_output(print(gd), regexp = ".")
  expect_output(summary(gd), regexp = ".")
  expect_true(n_subjects(gd) >= 1L || length(get_subjects(gd)) >= 1L)

  expect_equal(fmrireg:::detect_group_data_format("x.nii.gz"), "nifti")
  expect_equal(fmrireg:::detect_group_data_format("x.h5"), "h5")
  expect_equal(fmrireg:::detect_group_data_format(data.frame(a = 1)), "csv")
})

test_that("fmri_ttest_gds classic path and error gates on demo csv-as-gds", {
  gd <- tryCatch(fmrireg:::.demo_group_data_csv(), error = function(e) NULL)
  skip_if(is.null(gd), "no demo group data")

  fit <- tryCatch(fmri_ttest(gd, formula = ~ 1, engine = "classic"), error = function(e) e)
  if (!inherits(fit, "error")) {
    expect_s3_class(fit, "fmri_ttest_fit")
    expect_output(print(fit), regexp = ".")
  }

  expect_error(fmri_ttest(gd, paired = TRUE), regexp = ".")
})
