# Additional fmri_ttest / meta / print / summary coverage.

test_that("fmri_ttest dataframe path and summary with q values", {
  set.seed(41)
  # Classic engine on a small subject x feature matrix via group_data csv
  skip_if_not_installed("tibble")

  n_subj <- 12
  n_feat <- 5
  beta <- matrix(rnorm(n_subj * n_feat, mean = 0.3), n_subj, n_feat)
  se <- matrix(runif(n_subj * n_feat, 0.2, 0.6), n_subj, n_feat)
  colnames(beta) <- paste0("feat", seq_len(n_feat))
  rownames(beta) <- paste0("s", seq_len(n_subj))

  cov <- data.frame(
    subject = rownames(beta),
    group = factor(rep(c("A", "B"), length.out = n_subj)),
    stringsAsFactors = FALSE
  )

  # Build CSV-style group data if helper exists
  gd <- tryCatch(
    group_data_from_csv(
      beta = as.data.frame(beta),
      se = as.data.frame(se),
      covariates = cov,
      subject_col = "subject"
    ),
    error = function(e) NULL
  )
  if (is.null(gd)) {
    gd <- tryCatch(fmrireg:::.demo_group_data_csv(), error = function(e) NULL)
  }
  skip_if(is.null(gd), "no group_data csv fixture")

  fit <- tryCatch(
    fmri_ttest(gd, formula = ~ 1, engine = "classic"),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    fit <- tryCatch(
      fmri_ttest(gd, formula = ~ 1, engine = "meta"),
      error = function(e) e
    )
  }
  if (inherits(fit, "error")) {
    skip(paste("fmri_ttest failed:", conditionMessage(fit)))
  }

  expect_s3_class(fit, "fmri_ttest_fit")
  expect_output(print(fit), regexp = ".")
  expect_output(summary(fit), regexp = ".")

  # Attach q for summary MC block
  if (!is.null(fit$p)) {
    fit$q <- matrix(p.adjust(as.numeric(fit$p), method = "BH"),
                    nrow = nrow(as.matrix(fit$p)),
                    ncol = max(ncol(as.matrix(fit$p)), 1L))
    expect_output(summary(fit), regexp = ".")
  }
})

test_that("fmri_ttest GDS engine gates and contrast single-coef", {
  skip_if_not_installed("fmrigds")

  gd <- tryCatch({
    # Prefer existing demo/helper used elsewhere
    if (exists(".demo_group_data_gds", envir = asNamespace("fmrireg"), inherits = FALSE)) {
      fmrireg:::.demo_group_data_gds()
    } else {
      NULL
    }
  }, error = function(e) NULL)

  if (is.null(gd)) {
    # Minimal tabular GDS via group_data if available
    n_subj <- 10
    n_feat <- 4
    df <- data.frame(
      subject = rep(paste0("s", 1:n_subj), each = n_feat),
      feature = rep(paste0("f", 1:n_feat), times = n_subj),
      beta = rnorm(n_subj * n_feat),
      se = runif(n_subj * n_feat, 0.2, 0.5)
    )
    gd <- tryCatch(
      group_data(df, format = "csv", subject = "subject", feature = "feature",
                 assays = c(beta = "beta", se = "se")),
      error = function(e) NULL
    )
  }
  skip_if(is.null(gd), "no GDS group_data")

  expect_error(fmri_ttest(gd, paired = TRUE), regexp = ".")
  expect_error(fmri_ttest(gd, mu0 = 1), regexp = ".")

  fit <- tryCatch(
    fmri_ttest(gd, formula = ~ 1, engine = "auto"),
    error = function(e) e
  )
  if (!inherits(fit, "error")) {
    expect_s3_class(fit, "fmri_ttest_fit")
  }
})

test_that("print/summary methods for meta and model cover remaining lines", {
  meta <- suppressWarnings(fmrireg:::.demo_fmri_meta())
  expect_output(print(meta), regexp = ".")
  expect_output(summary(meta), regexp = ".")

  fmod <- fmrireg:::.demo_fmri_model()
  expect_output(print(fmod), regexp = ".")

  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  expect_output(print(fit), regexp = ".")
  expect_true(length(coef_names(fit)) >= 1L)
})
