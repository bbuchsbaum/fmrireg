# fmri_meta_gds.R: exercise real group_data_gds meta paths and combine gates.

test_that(".fmri_meta_group_data_gds_impl FE/DL paths on tabular GDS", {
  skip_if_not_installed("fmrigds")

  set.seed(5)
  n_subj <- 8L
  n_roi <- 4L
  df <- expand.grid(
    subject = paste0("s", seq_len(n_subj)),
    roi = paste0("ROI", seq_len(n_roi)),
    contrast = "A_vs_B",
    stringsAsFactors = FALSE
  )
  df$group <- factor(ifelse(as.integer(sub("s", "", df$subject)) <= 4L, "A", "B"))
  df$beta <- rnorm(nrow(df), mean = 0.3, sd = 0.2)
  df$se <- runif(nrow(df), 0.08, 0.2)

  gd <- fmrireg:::.group_data_gds(
    df,
    format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    contrast_col = "contrast",
    covariate_cols = "group"
  )

  fe <- fmrireg:::.fmri_meta_group_data_gds_impl(
    gd, formula = ~ 1, method = "fe", robust = "none",
    weights = "ivw", verbose = FALSE
  )
  expect_s3_class(fe, "fmri_meta")
  expect_true(!is.null(fe$coefficients))
  expect_true(nrow(fe$coefficients) >= 1L)
  expect_equal(fe$method, "fe")

  dl <- fmrireg:::.fmri_meta_group_data_gds_impl(
    gd, formula = ~ group, method = "dl", robust = "none",
    weights = "equal", verbose = FALSE
  )
  expect_true(ncol(dl$coefficients) >= 1L)
  expect_true(!is.null(dl$tau2) || !is.null(dl$I2) || !is.null(dl$coefficients))

  # Custom weights vector
  w <- rep(1, n_subj)
  custom <- fmrireg:::.fmri_meta_group_data_gds_impl(
    gd, formula = ~ 1, method = "fe",
    weights = "custom", weights_custom = w, verbose = FALSE
  )
  expect_true(!is.null(custom$coefficients))

  # return_cov tri on pm when supported
  pm <- fmrireg:::.fmri_meta_group_data_gds_impl(
    gd, formula = ~ 1, method = "pm",
    weights = "ivw", return_cov = "tri", verbose = FALSE
  )
  expect_true(!is.null(pm$coefficients))
})

test_that(".fmri_meta_group_data_gds_impl combine validation branches", {
  skip_if_not_installed("fmrigds")

  set.seed(6)
  n_subj <- 6L
  n_roi <- 3L
  df <- expand.grid(
    subject = paste0("s", seq_len(n_subj)),
    roi = paste0("ROI", seq_len(n_roi)),
    contrast = "A_vs_B",
    stringsAsFactors = FALSE
  )
  df$t <- rnorm(nrow(df))
  df$df <- 20
  df$beta <- df$t * 0.1
  df$se <- rep(0.1, nrow(df))

  # Prefer a GDS that exposes t+df if the builder supports it
  gd <- tryCatch(
    fmrireg:::.group_data_gds(
      df,
      format = "csv",
      effect_cols = c(t = "t", df = "df", beta = "beta", se = "se"),
      subject_col = "subject",
      roi_col = "roi",
      contrast_col = "contrast"
    ),
    error = function(e) {
      fmrireg:::.group_data_gds(
        df,
        format = "csv",
        effect_cols = c(beta = "beta", se = "se"),
        subject_col = "subject",
        roi_col = "roi",
        contrast_col = "contrast"
      )
    }
  )

  # Early gates that don't require t assay to be present still exercise match.arg
  expect_error(
    fmrireg:::.fmri_meta_group_data_gds_impl(
      gd, weights = "custom", weights_custom = NULL, verbose = FALSE
    ),
    "weights_custom"
  )

  an_try <- try(fmrigds::compute(gd), silent = TRUE)
  anames <- fmrireg:::.gds_safe_assay_names(
    if (inherits(an_try, "try-error")) gd else an_try
  )
  if ("t" %in% anames) {
    expect_error(
      fmrireg:::.fmri_meta_group_data_gds_impl(
        gd, formula = ~ group, combine = "stouffer",
        weights = "equal", verbose = FALSE
      ),
      "intercept-only"
    )
    expect_error(
      fmrireg:::.fmri_meta_group_data_gds_impl(
        gd, formula = ~ 1, combine = "stouffer",
        weights = "ivw", verbose = FALSE
      ),
      "Inverse-variance|equal or custom"
    )
    expect_error(
      fmrireg:::.fmri_meta_group_data_gds_impl(
        gd, formula = ~ 1, combine = "fisher",
        weights = "custom", weights_custom = rep(1, n_subj),
        verbose = FALSE
      ),
      "Fisher|equal"
    )

    comb <- fmrireg:::.fmri_meta_group_data_gds_impl(
      gd, formula = ~ 1, combine = "stouffer",
      weights = "equal", verbose = FALSE
    )
    expect_true(!is.null(comb$coefficients))
    expect_match(comb$method, "combine")
  } else {
    # Still cover a successful intercept FE path as fallback
    fe <- fmrireg:::.fmri_meta_group_data_gds_impl(
      gd, formula = ~ 1, method = "fe", verbose = FALSE
    )
    expect_true(!is.null(fe$coefficients))
  }
})

test_that(".fmri_meta_group_data_gds_impl robust huber branch", {
  skip_if_not_installed("fmrigds")
  set.seed(8)
  df <- expand.grid(
    subject = paste0("s", 1:10),
    roi = paste0("R", 1:3),
    contrast = "c1",
    stringsAsFactors = FALSE
  )
  df$beta <- rnorm(nrow(df), 0.2, 0.25)
  df$beta[1] <- 5
  df$se <- runif(nrow(df), 0.1, 0.25)
  gd <- fmrireg:::.group_data_gds(
    df, format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject", roi_col = "roi", contrast_col = "contrast"
  )
  out <- fmrireg:::.fmri_meta_group_data_gds_impl(
    gd, formula = ~ 1, method = "fe", robust = "huber",
    weights = "ivw", verbose = FALSE
  )
  expect_equal(out$robust, "huber")
  expect_true(all(is.finite(out$coefficients) | is.na(out$coefficients)))
})
