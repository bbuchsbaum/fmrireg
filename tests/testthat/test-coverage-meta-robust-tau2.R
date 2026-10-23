# meta_gds robust/reml leftovers.

test_that(".fmri_meta_group_data_gds_impl robust and reml branches", {
  skip_if_not_installed("fmrigds")
  set.seed(41)
  n_subj <- 10L
  n_roi <- 3L
  df <- expand.grid(
    subject = paste0("s", seq_len(n_subj)),
    roi = paste0("ROI", seq_len(n_roi)),
    contrast = "A_vs_B",
    stringsAsFactors = FALSE
  )
  df$group <- factor(ifelse(as.integer(sub("s", "", df$subject)) <= 5L, "A", "B"))
  df$beta <- rnorm(nrow(df), 0.2, 0.25)
  df$se <- runif(nrow(df), 0.05, 0.2)
  df$beta[1:2] <- df$beta[1:2] + 5

  gd <- fmrireg:::.group_data_gds(
    df, format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject", roi_col = "roi", contrast_col = "contrast",
    covariate_cols = "group"
  )

  hub <- fmrireg:::.fmri_meta_group_data_gds_impl(
    gd, formula = ~ 1, method = "fe", robust = "huber",
    weights = "ivw", verbose = FALSE
  )
  expect_s3_class(hub, "fmri_meta")

  reml <- tryCatch(
    fmrireg:::.fmri_meta_group_data_gds_impl(
      gd, formula = ~ group, method = "reml", robust = "none",
      weights = "ivw", verbose = FALSE
    ),
    error = function(e) e
  )
  expect_true(inherits(reml, "error") || inherits(reml, "fmri_meta"))

  expect_error(
    fmrireg:::.fmri_meta_group_data_gds_impl(
      gd, formula = ~ 1, method = "fe",
      weights = "custom", weights_custom = NULL, verbose = FALSE
    ),
    "weights_custom"
  )
})

# LEGACY-GROUP-DELETE (next minor): tests fit_meta_single() and fmri_meta.group_data_csv() on .demo_group_data_csv().
test_that("fit_meta_single and fmri_meta.group_data_csv edges", {
  set.seed(43)
  y <- rnorm(9, 0.2, 0.3)
  se <- runif(9, 0.05, 0.2)
  X <- cbind(1, rep(0:1, length.out = 9))

  single <- fmrireg:::fit_meta_single(
    y, se, X, method = "fe", robust = "none",
    weights = "ivw", weights_custom = NULL
  )
  expect_true(is.list(single))
  expect_true(!is.null(single$coefficients))

  single_eq <- fmrireg:::fit_meta_single(
    y, se, X, method = "dl", robust = "huber",
    weights = "equal", weights_custom = NULL
  )
  expect_true(is.list(single_eq))

  single_cw <- fmrireg:::fit_meta_single(
    y, se, X, method = "fe", robust = "none",
    weights = "custom", weights_custom = rep(1, 9)
  )
  expect_true(is.list(single_cw))

  gd <- fmrireg:::.demo_group_data_csv()
  fit <- tryCatch(
    fmri_meta(gd, formula = ~ 1, method = "fe", verbose = FALSE),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    expect_match(conditionMessage(fit), ".")
  } else {
    expect_s3_class(fit, "fmri_meta")
  }

  expect_error(
    fmri_meta(gd, formula = ~ 1, method = "not-a-method", verbose = FALSE)
  )
})
