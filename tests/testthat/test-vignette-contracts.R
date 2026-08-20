test_that("factorial vignette contrasts resolve sanitized design term tags", {
  design <- expand.grid(
    category = factor(c("face", "scene")),
    attention = factor(c("attend", "ignore")),
    replication = 1:2
  )
  design$onset <- seq(1, 100, length.out = nrow(design))
  design$block <- 1L

  sframe <- sampling_frame(blocklens = 120, TR = 2)
  event <- event_model(
    onset ~ hrf(category, attention), block = ~ block,
    data = design, sampling_frame = sframe
  )
  X <- as.matrix(design_matrix(event))
  truth <- cbind(
    voxel1 = c(1.2, 0.4, 0.8, 0.2),
    voxel2 = c(1.0, 0.5, 0.7, 0.3),
    voxel3 = c(0.8, 0.3, 0.6, 0.2)
  )
  Y <- X %*% truth

  contrasts <- contrast_set(
    pair_contrast(
      ~ category == "face", ~ category == "scene",
      name = "face_vs_scene"
    ),
    pair_contrast(
      ~ attention == "attend", ~ attention == "ignore",
      name = "attend_vs_ignore"
    ),
    contrast(
      ~ (`face:attend` - `face:ignore`) -
        (`scene:attend` - `scene:ignore`),
      name = "category_X_attention"
    )
  )
  dataset <- matrix_dataset(
    Y, TR = 2, run_length = 120, event_table = design
  )

  expect_no_warning(
    fit <- fmri_lm(
      onset ~ hrf(category, attention, contrasts = contrasts),
      block = ~ block, dataset = dataset
    )
  )

  observed <- as.matrix(coef(fit, type = "contrasts"))
  expected <- cbind(
    face_vs_scene = c(0.70, 0.45, 0.45),
    attend_vs_ignore = c(0.30, 0.25, 0.15),
    category_X_attention = c(0.20, 0.10, 0.10)
  )
  expect_equal(unname(observed), unname(expected), tolerance = 1e-8)
})

test_that("vignette omnibus contrast executes as a two-df F test", {
  events <- data.frame(
    onset = seq(5, 80, length.out = 12),
    condition = factor(rep(c("faces", "objects", "scenes"), 4)),
    run = 1L
  )
  sframe <- sampling_frame(100, TR = 2)
  event <- event_model(
    onset ~ hrf(condition), block = ~ run,
    data = events, sampling_frame = sframe
  )
  X <- as.matrix(design_matrix(event))
  set.seed(404)
  Y <- X %*% matrix(c(1.0, 0.4, -0.2), ncol = 1) +
    matrix(rnorm(100, sd = 0.15), ncol = 1)
  dataset <- matrix_dataset(
    Y, TR = 2, run_length = 100, event_table = events
  )
  omnibus <- oneway_contrast(~ condition, name = "condition_F")
  weights <- contrast_weights(omnibus, terms(event)[[1]])$weights

  fit <- fmri_lm(
    onset ~ hrf(condition, contrasts = contrast_set(omnibus)),
    block = ~ run, dataset = dataset
  )
  f_stat <- as.numeric(stats(fit, type = "F")$condition_F)

  expect_equal(ncol(weights), 2L)
  expect_true(is.finite(f_stat) && f_stat > 0)
})

test_that("vignette tables preserve scientific labels and nonzero p-value text", {
  cell_order <- expand.grid(
    category = factor(c("face", "scene")),
    attention = factor(c("attend", "ignore"))
  )
  cell_labels <- paste(cell_order$category, cell_order$attention, sep = " / ")
  displayed_p <- format.pval(
    c(2.57e-32, 3.92e-13, 2.62e-5, 5.39e-53, 9.98e-47),
    digits = 3,
    eps = 1e-4
  )

  expect_equal(
    cell_labels,
    c("face / attend", "scene / attend", "face / ignore", "scene / ignore")
  )
  expect_false(any(grepl("category[._]|attention[._]", cell_labels)))
  expect_false(any(displayed_p == "0"))
  expect_true(all(grepl("^<1e-04$", displayed_p)))
})

test_that("multisubject teaching results cross the public GDS group boundary", {
  skip_if_not_installed("fmrigds")

  teaching_effects <- data.frame(
    job_id = c("sub-01", "sub-02"),
    estimate = c(-0.1, -0.02),
    se = c(0.2, 0.25),
    roi = "teaching_voxel_1",
    contrast = "condition_B"
  )
  gd <- fmrireg::group_data(
    teaching_effects,
    format = "csv",
    effect_cols = c(beta = "estimate", se = "se"),
    subject_col = "job_id",
    roi_col = "roi",
    contrast_col = "contrast"
  )
  fit <- fmri_meta(gd, ~ 1, method = "fe", verbose = FALSE)

  expect_s3_class(gd, "group_data_gds")
  expect_s3_class(fit, "fmri_meta")
  expect_true(all(is.finite(fit$coefficients)))
  expect_true(all(is.finite(fit$se)))
})

test_that("public group data executes named fit-time group contrasts", {
  skip_if_not_installed("fmrigds")

  group_df <- data.frame(
    subject = sprintf("s%02d", 1:10),
    roi = "teaching_roi",
    contrast = "task",
    beta = c(0.40, 0.45, 0.50, 0.55, 0.60,
             1.40, 1.45, 1.50, 1.55, 1.60),
    se = 0.25,
    group = factor(rep(c("A", "B"), each = 5))
  )
  expect_no_warning(
    gd <- fmrireg::group_data(
      group_df,
      format = "csv",
      effect_cols = c(beta = "beta", se = "se"),
      subject_col = "subject",
      roi_col = "roi",
      contrast_col = "contrast",
      covariate_cols = "group"
    )
  )
  weights <- c("(Intercept)" = 0, groupB = 1)
  fit_cov <- fmri_meta(
    gd, ~ 1 + group, method = "pm", return_cov = "tri", verbose = FALSE
  )
  post_hoc <- contrast(fit_cov, weights)
  fit_time <- fmri_meta(
    gd, ~ 1 + group, method = "pm", contrasts = weights, verbose = FALSE
  )

  expect_s3_class(gd, "group_data_gds")
  expect_false(is.null(fit_time$contrasts))
  expect_equal(
    as.numeric(fit_time$contrasts$estimate),
    as.numeric(post_hoc$estimate),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(fit_time$contrasts$se),
    as.numeric(post_hoc$se),
    tolerance = 1e-10
  )
  expect_error(
    fmri_meta(
      gd, ~ 1 + group, method = "pm",
      contrasts = c("(Intercept)" = 0, missing_term = 1),
      verbose = FALSE
    ),
    "Unknown contrast coefficient"
  )

  fit_welch <- fmri_ttest(
    gd, ~ 1 + group, engine = "welch", sign = "BminusA"
  )
  fit_classic <- fmri_ttest(
    gd, ~ 1 + group, engine = "classic", sign = "BminusA"
  )
  fit_custom <- fmri_ttest(
    gd, ~ 1 + group, engine = "meta",
    weights = "custom", weights_custom = rep(1, 10), sign = "BminusA"
  )
  fit_custom_meta <- fmri_meta(
    gd, ~ 1 + group, method = "pm",
    weights = "custom", weights_custom = rep(1, 10),
    verbose = FALSE
  )
  expect_gt(as.numeric(fit_welch$t["group", 1]), 0)
  expect_gt(as.numeric(fit_classic$t["group", 1]), 0)
  expect_identical(fit_custom$weights, "custom")
  expect_identical(fit_custom$method, "pm")
  expect_equal(
    unname(fit_custom$beta), unname(t(fit_custom_meta$coefficients)),
    tolerance = 1e-10
  )
  expect_equal(
    unname(fit_custom$se), unname(t(fit_custom_meta$se)), tolerance = 1e-10
  )
  expect_identical(rownames(fit_custom$beta), c("(Intercept)", "group"))
  expect_true(all(is.finite(fit_custom$z)))
})
