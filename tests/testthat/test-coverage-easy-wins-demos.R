# Easy coverage wins: zero-hit exports/aliases, demo helpers, rank utilities.

test_that("estimate() deprecation and fast_preproject alias behave", {
  expect_error(estimate(), "deprecated")
  expect_error(estimate(1, 2, a = 3), "estimate_betas")

  set.seed(1)
  X <- cbind(1, rnorm(20), rnorm(20))
  alias_proj <- fmrireg:::fast_preproject(X)
  direct_proj <- fmrireg:::.fast_preproject(X)
  expect_equal(alias_proj$dfres, direct_proj$dfres)
  expect_equal(alias_proj$XtXinv, direct_proj$XtXinv, tolerance = 1e-12)
  expect_true(isTRUE(alias_proj$is_full_rank) || isFALSE(alias_proj$is_full_rank))
})

test_that("design_matrix.fmri_lm delegates to the model design matrix", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  expect_s3_class(fit, "fmri_lm")
  dm_fit <- design_matrix(fit)
  dm_model <- design_matrix(fit$model)
  expect_equal(as.matrix(dm_fit), as.matrix(dm_model))
  expect_equal(nrow(dm_fit), sum(fit$model$event_model$sampling_frame$blocklens))
  expect_true(ncol(dm_fit) >= 1L)
})

test_that(".demo_* helpers construct coherent tiny objects", {
  ev <- fmrireg:::.demo_event_data()
  expect_equal(nrow(ev), 4L)
  expect_true(all(c("onsets", "condition", "run") %in% names(ev)))

  sframe <- fmrireg:::.demo_sampling_frame()
  expect_equal(sum(sframe$blocklens), 8L)

  emod <- fmrireg:::.demo_event_model()
  expect_s3_class(emod, "event_model")
  expect_true(ncol(design_matrix(emod)) >= 1L)

  dset <- fmrireg:::.demo_matrix_dataset()
  expect_s3_class(dset, "matrix_dataset")
  expect_equal(dim(dset$datamat), c(8L, 2L))

  fmod <- fmrireg:::.demo_fmri_model()
  expect_s3_class(fmod, "fmri_model")

  gd <- fmrireg:::.demo_group_data_csv()
  expect_true(inherits(gd, "group_data") || inherits(gd, "group_data_csv"))
  expect_true(nrow(gd$data) >= 3L)

  meta <- suppressWarnings(fmrireg:::.demo_fmri_meta())
  expect_s3_class(meta, "fmri_meta")
  expect_true(is.matrix(meta$coefficients) || is.matrix(meta$beta) ||
                !is.null(meta$result) || length(coef(meta)) > 0)

  tfit <- fmrireg:::.demo_fmri_ttest_fit()
  expect_s3_class(tfit, "fmri_ttest_fit")
  expect_equal(dim(tfit$beta), c(2L, 1L))
  expect_equal(rownames(tfit$beta), c("conditionA", "conditionB"))
})

test_that("fmri_lm_internal rank helpers cover full and deficient designs", {
  expect_true(fmrireg:::is.formula(~ x))
  expect_false(fmrireg:::is.formula("x"))

  X_full <- cbind(1, c(0, 1, 0, 1), c(1, 2, 3, 4))
  info <- fmrireg:::.design_rank_info(X_full)
  expect_true(info$is_full_rank)
  expect_equal(info$rank, 3L)
  expect_equal(length(info$aliased), 0L)

  X_def <- cbind(1, c(0, 1, 0, 1), c(0, 1, 0, 1))
  colnames(X_def) <- c("Intercept", "A", "A_dup")
  info_def <- fmrireg:::.design_rank_info(X_def)
  expect_false(info_def$is_full_rank)
  expect_equal(info_def$rank, 2L)
  expect_true(length(info_def$aliased) >= 1L)

  msg <- fmrireg:::.rank_deficiency_message(
    rank = 2L, p = 3L, aliased = 3L,
    varnames = colnames(X_def), context = "Design matrix"
  )
  expect_match(msg, "rank deficient")
  expect_match(msg, "A_dup")

  long_msg <- fmrireg:::.rank_deficiency_message(
    rank = 1L, p = 12L, aliased = 2:12, varnames = paste0("v", 1:12)
  )
  expect_match(long_msg, "\\.\\.\\.")

  attached <- fmrireg:::.attach_rank_attrs(list(ok = TRUE), info_def)
  expect_equal(attr(attached, "rank"), 2L)
  expect_false(attr(attached, "is_full_rank"))
  expect_equal(attr(attached, "aliased"), info_def$aliased)

  expect_error(
    fmrireg:::.design_rank_info(matrix(numeric(), 0, 0)),
    "at least one row"
  )

  status_ok <- fmrireg:::.fmri_lm_voxel_status(rnorm(10))
  expect_true(is.list(status_ok) || is.character(status_ok) || is.logical(status_ok))
  status_bad <- fmrireg:::.fmri_lm_voxel_status(c(1, NA, Inf))
  expect_false(identical(status_ok, status_bad))
})
