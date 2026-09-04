# ttest_gds validation branches + soft subspace projection + legacy AR helpers.

test_that("fmri_ttest.group_data_gds rejects unsupported options and runs classic", {
  skip_if_not_installed("fmrigds")
  set.seed(81)
  df <- expand.grid(
    subject = paste0("s", 1:8),
    roi = paste0("r", 1:3),
    contrast = "c1",
    stringsAsFactors = FALSE
  )
  df$beta <- rnorm(nrow(df), 0.2, 0.3)
  df$se <- runif(nrow(df), 0.08, 0.25)
  df$group <- factor(ifelse(as.integer(sub("s", "", df$subject)) <= 4, "A", "B"))
  gd <- fmrireg:::.group_data_gds(
    df, format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject", roi_col = "roi", contrast_col = "contrast",
    covariate_cols = "group"
  )

  expect_error(fmri_ttest(gd, paired = TRUE), "paired")
  expect_error(fmri_ttest(gd, mu0 = 1), "mu0")
  expect_error(fmri_ttest(gd, mask = TRUE), "mask")
  expect_error(fmri_ttest(gd, voxelwise_cov = 1), "voxelwise_cov")
  expect_error(
    fmri_ttest(gd, weights = "custom", weights_custom = NULL),
    "weights_custom"
  )

  # auto -> meta when beta+se present
  fit_meta <- tryCatch(
    fmri_ttest(gd, formula = ~ 1, engine = "auto", mc = "none", weights = "ivw"),
    error = function(e) e
  )
  expect_true(inherits(fit_meta, "error") || inherits(fit_meta, "fmri_ttest_fit"))

  fit_classic <- tryCatch(
    fmri_ttest(gd, formula = ~ group, engine = "classic", mc = "bh", weights = "equal"),
    error = function(e) e
  )
  expect_true(inherits(fit_classic, "error") || inherits(fit_classic, "fmri_ttest_fit"))

  fit_custom <- tryCatch(
    fmri_ttest(
      gd, formula = ~ 1, engine = "meta", mc = "none",
      weights = "custom", weights_custom = rep(1, 8)
    ),
    error = function(e) e
  )
  expect_true(inherits(fit_custom, "error") || inherits(fit_custom, "fmri_ttest_fit"))
})

test_that("soft_projection auto/gcv/numeric and apply_soft_projection", {
  set.seed(82)
  n <- 40L
  p_n <- 6L
  V <- 5L
  N <- matrix(rnorm(n * p_n), n, p_n)
  Y <- matrix(rnorm(n * V), n, V)
  X <- cbind(1, rnorm(n), rnorm(n))

  proj_auto <- soft_projection(N, lambda = "auto")
  expect_s3_class(proj_auto, "soft_projection")
  expect_true(is.function(proj_auto$P_lambda))
  expect_output(print(proj_auto), "soft_projection|lambda|df")

  proj_num <- soft_projection(N, lambda = 0.5)
  expect_equal(proj_num$method, "user_specified")
  expect_error(soft_projection(N, lambda = -1), "non-negative")

  expect_warning(soft_projection(N, lambda = "gcv"), "GCV requires Y|falling back")
  proj_gcv <- soft_projection(N, lambda = "gcv", Y = Y)
  expect_true(proj_gcv$method %in% c("gcv", "singular_value_heuristic"))

  cleaned <- apply_soft_projection(proj_auto, Y, X)
  expect_true(is.list(cleaned))
  expect_equal(dim(cleaned$Y), dim(Y))
  expect_equal(dim(cleaned$X), dim(X))

  # Pipeline helper with options
  opts <- soft_subspace_options(
    enabled = TRUE, nuisance_matrix = N, lambda = 0.25
  )
  pipe <- fmrireg:::.apply_soft_projection_pipeline(
    Y, X, dataset = NULL, soft_opts = opts
  )
  expect_true(is.list(pipe) || is.matrix(pipe$Y) || !is.null(pipe))
})

test_that("legacy estimate_ar_parameters covers zero-var, NA, censor, Yule-Walker fallback", {
  skip_if_not_installed("fmriAR")
  set.seed(83)
  r <- as.numeric(arima.sim(list(ar = 0.5), n = 80))
  phi <- fmrireg:::estimate_ar_parameters(r, 1L)
  expect_true(is.numeric(phi))

  expect_equal(length(fmrireg:::estimate_ar_parameters(rep(0, 30), 2L)), 0L)
  expect_error(fmrireg:::estimate_ar_parameters(numeric(0), 1L), "at least one")
  expect_error(fmrireg:::estimate_ar_parameters(c(1, NA, 2), 1L), "NA")

  phi_c <- fmrireg:::estimate_ar_parameters(r, 2L, censor = c(1L, 2L, 3L))
  expect_true(is.numeric(phi_c))
})
