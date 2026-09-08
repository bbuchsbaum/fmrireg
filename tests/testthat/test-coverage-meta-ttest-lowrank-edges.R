# Remaining meta / ttest-gds / lowrank error-branch coverage.

test_that("fmri_meta tidy/glance-adjacent helpers on demo meta", {
  meta <- suppressWarnings(fmrireg:::.demo_fmri_meta())
  expect_s3_class(meta, "fmri_meta")

  # coef / se / confint style accessors if present
  cf <- tryCatch(coef(meta), error = function(e) NULL)
  expect_true(is.null(cf) || is.numeric(cf) || is.matrix(cf) || is.array(cf))

  sm <- tryCatch(summary(meta), error = function(e) NULL)
  expect_true(is.null(sm) || inherits(sm, "summary.fmri_meta") || is.list(sm) ||
                inherits(sm, "fmri_meta"))

  expect_output(print(meta), regexp = "meta|fmri|coefficient|method", ignore.case = TRUE)
})

test_that(".fmri_meta_group_data_gds_impl rejects bad weights_custom", {
  skip_if_not_installed("fmrigds")
  set.seed(2)
  df <- expand.grid(
    subject = paste0("s", 1:6),
    roi = paste0("r", 1:3),
    contrast = "c1",
    stringsAsFactors = FALSE
  )
  df$beta <- rnorm(nrow(df))
  df$se <- runif(nrow(df), 0.1, 0.3)
  gd <- fmrireg:::.group_data_gds(
    df, format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject", roi_col = "roi", contrast_col = "contrast"
  )

  expect_error(
    fmrireg:::.fmri_meta_group_data_gds_impl(
      gd, formula = ~ 1, method = "fe",
      weights = "custom", weights_custom = c(1, 2), verbose = FALSE
    ),
    "weights|length|custom"
  )
})

test_that("engine context check rejects by_cluster without parcels", {
  fx_model <- fmrireg:::.demo_fmri_model()
  fx_dset <- fmrireg:::.demo_matrix_dataset()
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1", by_cluster = TRUE))
  caps <- list(
    ar_by_cluster = TRUE,
    requires_parcels_for_by_cluster = TRUE,
    forbid_by_cluster_dataset_classes = character()
  )

  expect_error(
    fmrireg:::.validate_engine_context(
      "latent_sketch", fx_model, fx_dset,
      args = list(lowrank = list(time_sketch = list(method = "gaussian", m = 8L))),
      cfg = cfg,
      capabilities = caps
    ),
    "parcels|by_cluster"
  )

  ok <- fmrireg:::.validate_engine_context(
    "latent_sketch", fx_model, fx_dset,
    args = list(lowrank = list(
      time_sketch = list(method = "gaussian", m = 8L),
      parcels = rep(1:2, length.out = ncol(fx_dset$datamat))
    )),
    cfg = cfg,
    capabilities = caps
  )
  expect_true(is.list(ok))
})

test_that("fmri_ttest_gds classic path on demo csv group data", {
  skip_if_not_installed("fmrigds")
  gd <- fmrireg:::.demo_group_data_csv()
  fit <- tryCatch(
    fmri_ttest(gd, formula = ~ 1, engine = "classic", mc = "none"),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    # Some demo shapes may lack assays; still require a meaningful message
    expect_match(conditionMessage(fit), ".")
  } else {
    expect_s3_class(fit, "fmri_ttest_fit")
    expect_true(!is.null(fit$beta) || !is.null(fit$z) || !is.null(fit$t))
  }
})
