test_that("fmri_meta ROI CSV with group covariate fits FE correctly", {
  skip_on_cran()

  set.seed(123)
  n_per_group <- 4
  subjects <- sprintf("s%02d", 1:(2 * n_per_group))
  group <- factor(rep(c("A", "B"), each = n_per_group))

  # Deterministic effects by group, constant SE
  beta <- ifelse(group == "A", 0.5, 1.5)
  se <- rep(0.2, length(beta))

  # One ROI, named "overall" (the name the legacy reader gave a table without
  # an ROI column).
  df <- data.frame(subject = subjects, roi = "overall", beta = beta, se = se,
                   group = group, stringsAsFactors = FALSE)

  gd <- group_data(
    df, format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    covariate_cols = c("group")
  )

  fit <- fmrireg:::fmri_meta(gd, formula = ~ 1 + group, method = "fe", verbose = FALSE)

  # One ROI named "overall"
  expect_s3_class(gd, "group_data_gds")
  expect_true(inherits(fit, "fmri_meta_roi"))
  expect_identical(fit$roi_names, "overall")
  expect_equal(fit$n_subjects, length(subjects))
  expect_equal(fit$n_rois, 1)

  coefs <- as.numeric(fit$coefficients[1, ])
  names(coefs) <- colnames(fit$coefficients)

  # Intercept equals group A mean; groupB coefficient equals difference B-A
  expect_equal(unname(coefs["(Intercept)"]), 0.5, tolerance = 1e-10)
  expect_equal(unname(coefs[grep("group", names(coefs))]), 1.0, tolerance = 1e-10)

  # SE under IVW with se=0.2 => w=25 each
  # Var(beta0) = 0.01 (SE=0.1); Var(groupB)=0.02 (SE~0.1414)
  se_hat <- as.numeric(fit$se[1, ])
  expect_equal(se_hat[1], 0.1, tolerance = 1e-6)
  expect_equal(se_hat[2], sqrt(0.02), tolerance = 1e-6)

  # Q df = n - p = 8 - 2 = 6
  expect_equal(fit$Q_df[1], 6)
  expect_equal(fit$Q[1], 0)
  expect_equal(fit$tau2[1], 0)
})

test_that("equal weights in ROI CSV give the OLS group fit (relative weights)", {
  skip_on_cran()

  # Equal weights are relative weights with an estimated residual scale, so
  # the fit is ordinary least squares on the subject effects. (This test used
  # to assert that the equal-weights SE exceeds the IVW SE, which only held
  # because equal weights were treated as unit sampling variances.)
  n_per_group <- 4
  subjects <- sprintf("s%02d", 1:(2 * n_per_group))
  group <- factor(rep(c("A", "B"), each = n_per_group))
  beta <- ifelse(group == "A", 0.5, 1.5) + c(0.1, -0.2, 0.05, 0.05, -0.1, 0.3, -0.15, -0.05)
  se <- rep(0.2, length(beta))
  # One ROI, named "overall" (the name the legacy reader gave a table without
  # an ROI column).
  df <- data.frame(subject = subjects, roi = "overall", beta = beta, se = se,
                   group = group, stringsAsFactors = FALSE)

  gd <- group_data(
    df, format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    covariate_cols = c("group")
  )

  fit_ivw <- fmrireg:::fmri_meta(gd, formula = ~ 1 + group, method = "fe", weights = "ivw", verbose = FALSE)
  fit_eq  <- fmrireg:::fmri_meta(gd, formula = ~ 1 + group, method = "fe", weights = "equal", verbose = FALSE)

  # Coefficients match (all SE equal so weighting doesn't change beta)
  expect_equal(as.numeric(fit_eq$coefficients), as.numeric(fit_ivw$coefficients), tolerance = 1e-10)

  # Equal-weights SE is the OLS standard error
  ols <- summary(stats::lm(beta ~ group, data = df))$coefficients
  expect_equal(as.numeric(fit_eq$se[1, ]), unname(ols[, "Std. Error"]), tolerance = 1e-10)
})

test_that("fmri_meta intercept-only matches grand mean in ROI CSV (PM)", {
  skip_on_cran()

  subjects <- sprintf("s%02d", 1:6)
  beta <- c(0.5, 0.6, 0.4, 1.5, 1.6, 1.4)
  se <- rep(0.3, length(beta))
  df <- data.frame(subject = subjects, roi = "overall", beta = beta, se = se,
                   stringsAsFactors = FALSE)

  gd <- group_data(
    df, format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi"
  )

  fit <- fmrireg:::fmri_meta(gd, formula = ~ 1, method = "pm", verbose = FALSE)
  grand_mean <- mean(beta)

  expect_equal(as.numeric(fit$coefficients[1,1]), grand_mean, tolerance = 1e-10)
  expect_equal(fit$Q_df[1], length(beta) - 1)
  expect_gte(fit$tau2[1], 0)
})
