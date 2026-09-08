# Coverage for group_data validation error paths, ttest helpers, meta_gds edges.

test_that("group_data_h5 / nifti / csv validation error paths", {
  expect_error(
    fmrireg:::group_data_from_h5(1:3),
    "character vector"
  )
  expect_error(
    fmrireg:::group_data_from_h5(c("/no/such/a.h5", "/no/such/b.h5")),
    "do not exist"
  )
  expect_error(
    fmrireg:::validate_group_data_h5(list(paths = "a")),
    "Missing required fields"
  )
  expect_error(
    fmrireg:::validate_group_data_h5(list(
      paths = c("a", "b"), subjects = "s1", format = "h5",
      dim = 1:3, labels = "x"
    )),
    "paths and subjects"
  )
  expect_true(
    fmrireg:::validate_group_data_h5(list(
      paths = c("a", "b"), subjects = c("s1", "s2"), format = "h5",
      dim = 1:3, labels = "x"
    ))
  )
  expect_error(
    fmrireg:::read_h5_full(list()),
    "group_data_h5"
  )

  # NIfTI constructor validation (no files needed for early stops)
  expect_error(
    fmrireg:::.group_data_from_nifti_impl(),
    "beta_paths' or 't_paths'"
  )
  expect_error(
    fmrireg:::.group_data_from_nifti_impl(beta_paths = c("a.nii", "b.nii")),
    "se_paths' or 'var_paths'"
  )
  expect_error(
    fmrireg:::.group_data_from_nifti_impl(
      beta_paths = "a.nii", se_paths = "s.nii", var_paths = "v.nii"
    ),
    "not both"
  )
  expect_error(
    fmrireg:::.group_data_from_nifti_impl(t_paths = c("t1.nii", "t2.nii")),
    "must also provide 'df'"
  )
  expect_error(
    fmrireg:::.group_data_from_nifti_impl(
      beta_paths = c("a.nii", "b.nii"), se_paths = "s.nii"
    ),
    "Length of 'se_paths'"
  )
  expect_error(
    fmrireg:::.group_data_from_nifti_impl(
      beta_paths = c("a.nii", "b.nii"),
      se_paths = c("s1.nii", "s2.nii"),
      subjects = "only_one",
      validate = FALSE
    ),
    "subjects' must match"
  )
  expect_error(
    fmrireg:::.group_data_from_nifti_impl(
      beta_paths = c("a.nii", "b.nii"),
      se_paths = c("s1.nii", "s2.nii"),
      covariates = list(age = 1:2),
      validate = FALSE
    ),
    "data frame"
  )
  expect_error(
    fmrireg:::.group_data_from_nifti_impl(
      beta_paths = c("a.nii", "b.nii"),
      se_paths = c("s1.nii", "s2.nii"),
      covariates = data.frame(age = 1),
      validate = FALSE
    ),
    "rows in 'covariates'"
  )
  expect_error(
    fmrireg:::.group_data_from_nifti_impl(
      t_paths = c("t1.nii", "t2.nii"),
      df = c(10, 20, 30),
      validate = FALSE
    ),
    "scalar or have length"
  )
  expect_error(
    fmrireg:::.group_data_from_nifti_impl(
      beta_paths = c("/no/a.nii", "/no/b.nii"),
      se_paths = c("/no/s1.nii", "/no/s2.nii")
    ),
    "do not exist"
  )

  expect_error(
    fmrireg:::validate_group_data_nifti(list()),
    "beta_paths or t_paths"
  )
  expect_error(
    fmrireg:::validate_group_data_nifti(list(beta_paths = "a.nii")),
    "se_paths or var_paths"
  )
  expect_error(
    fmrireg:::validate_group_data_nifti(list(t_paths = "t.nii")),
    "must also have df"
  )
  expect_error(fmrireg:::read_nifti_full(list()), "group_data_nifti")

  # CSV validation paths
  expect_error(
    fmrireg:::.group_data_from_csv_impl("/no/such.csv", effect_cols = c(beta = "b")),
    "does not exist"
  )
  expect_error(
    fmrireg:::.group_data_from_csv_impl(1:3, effect_cols = c(beta = "b")),
    "file path or data frame"
  )
  expect_error(
    fmrireg:::.group_data_from_csv_impl(
      data.frame(id = 1, beta = 1, se = 0.1),
      effect_cols = c(beta = "beta", se = "se"),
      subject_col = "subject"
    ),
    "Subject column"
  )
  expect_error(
    fmrireg:::.group_data_from_csv_impl(
      data.frame(subject = "s1", beta = 1, se = 0.1),
      effect_cols = c(beta = "beta", se = "se"),
      roi_col = "roi"
    ),
    "ROI column"
  )
  expect_error(
    fmrireg:::.group_data_from_csv_impl(
      data.frame(subject = "s1", beta = 1, se = 0.1),
      effect_cols = c(beta = "beta", se = "se"),
      contrast_col = "contrast"
    ),
    "Contrast column"
  )
  expect_error(
    fmrireg:::.group_data_from_csv_impl(
      data.frame(subject = "s1", beta = 1, se = 0.1),
      effect_cols = c(beta = "beta", se = "se"),
      covariate_cols = "age"
    ),
    "Covariate columns not found"
  )
  expect_error(
    fmrireg:::validate_effect_cols(c("beta"), data.frame(beta = 1)),
    "named vector"
  )
  expect_error(
    fmrireg:::validate_effect_cols(c(beta = "missing_col"), data.frame(beta = 1)),
    "not found"
  )
  expect_error(
    fmrireg:::validate_effect_cols(c(beta = "beta"), data.frame(beta = 1)),
    "se' or 'var'"
  )
  expect_error(
    fmrireg:::wide_to_long_format(data.frame(x = 1), list(), "subject"),
    "not yet implemented"
  )

  expect_error(
    fmrireg:::validate_group_data_csv(list(data = 1)),
    "Missing required fields"
  )
  expect_error(
    fmrireg:::validate_group_data_csv(list(
      data = 1, subjects = "s1", subject_col = "subject",
      effect_cols = list(beta = "b"), format = "csv"
    )),
    "data frame"
  )
  expect_error(
    fmrireg:::validate_group_data_csv(list(
      data = data.frame(x = 1), subjects = "s1", subject_col = "subject",
      effect_cols = c(beta = "x"), format = "csv"
    )),
    "named list"
  )
  expect_error(fmrireg:::extract_csv_data(list()), "group_data_csv")
  expect_error(fmrireg:::get_rois(list()), "group_data_csv")
  expect_error(fmrireg:::get_contrasts(list()), "group_data_csv")

  # Happy-path CSV extract/filter errors on synthetic object
  gd_csv <- structure(
    list(
      data = data.frame(
        subject = c("s1", "s2"),
        beta = c(0.5, 0.7),
        se = c(0.1, 0.2),
        roi = c("r1", "r1"),
        contrast = c("A", "A"),
        stringsAsFactors = FALSE
      ),
      subjects = c("s1", "s2"),
      subject_col = "subject",
      effect_cols = list(beta = "beta", se = "se"),
      roi_col = "roi",
      contrast_col = "contrast",
      rois = "r1",
      contrasts = "A",
      format = "csv"
    ),
    class = c("group_data_csv", "group_data")
  )
  expect_true(fmrireg:::validate_group_data_csv(gd_csv))
  got <- fmrireg:::extract_csv_data(gd_csv, roi = "r1", contrast = "A")
  expect_equal(length(got$beta), 2L)
  expect_error(fmrireg:::extract_csv_data(gd_csv, roi = "missing"), "No data found for ROI")
  expect_error(
    fmrireg:::extract_csv_data(gd_csv, contrast = "missing"),
    "No data found for contrast"
  )
  expect_equal(fmrireg:::get_rois(gd_csv), "r1")
  expect_equal(fmrireg:::get_contrasts(gd_csv), "A")

  gd_no_roi <- gd_csv
  gd_no_roi$roi_col <- NULL
  expect_error(fmrireg:::extract_csv_data(gd_no_roi, roi = "r1"), "no ROI column")
  gd_no_con <- gd_csv
  gd_no_con$contrast_col <- NULL
  expect_error(
    fmrireg:::extract_csv_data(gd_no_con, contrast = "A"),
    "no contrast column"
  )
})

test_that("fmri_ttest helper materialize/weight/contrast branches", {
  expect_error(
    fmrireg:::.fmri_ttest_materialize_effects(list()),
    "cannot be materialized"
  )
  expect_error(
    fmrireg:::.fmri_ttest_materialize_effects(
      structure(list(n_voxels = 0L), class = "group_data_nifti")
    ),
    "missing beta|beta paths|group_data_nifti"
  )

  expect_error(
    fmrireg:::.fmri_ttest_validate_weights_custom("bad", 3, 4),
    "length S or SxP"
  )
  expect_error(
    fmrireg:::.fmri_ttest_validate_weights_custom(c(1, 2), 3, 4),
    "length S or SxP"
  )
  expect_error(
    fmrireg:::.fmri_ttest_validate_weights_custom(matrix(1, 2, 2), 3, 4),
    "length S or SxP"
  )
  expect_error(
    fmrireg:::.fmri_ttest_validate_weights_custom(matrix(c(1, -1, 1), 3, 1), 3, 1),
    "positive"
  )

  X <- cbind(`(Intercept)` = 1, groupB = c(0, 0, 1, 1))
  covars <- data.frame(group = factor(c("A", "A", "B", "B")))
  ginfo <- fmrireg:::.fmri_ttest_group_term(X, covars)
  expect_equal(ginfo$canonical_name, "group")
  expect_equal(ginfo$raw_name, "groupB")
  expect_null(fmrireg:::.fmri_ttest_group_term(X, data.frame(age = 1:4)))

  mat <- matrix(1:4, 2, 2, dimnames = list(c("groupB", "age"), NULL))
  renamed <- fmrireg:::.fmri_ttest_rename_rows(mat, "groupB", "group")
  expect_equal(rownames(renamed)[1], "group")
  expect_null(fmrireg:::.fmri_ttest_rename_rows(NULL, "a", "b"))

  res <- list(
    beta = mat,
    se = mat,
    t = mat,
    z = mat,
    p = mat,
    q = mat,
    df = mat
  )
  normed <- fmrireg:::.fmri_ttest_normalize_group_rows(res, ginfo)
  expect_equal(rownames(normed$beta)[1], "group")
  expect_identical(
    fmrireg:::.fmri_ttest_normalize_group_rows(res, NULL),
    res
  )

  expect_equal(
    fmrireg:::.fmri_ttest_canonical_coef_names(c("groupB", "age"), ginfo),
    c("group", "age")
  )

  signed <- fmrireg:::.fmri_ttest_apply_group_sign(
    list(
      beta = matrix(c(2, 1), 2, 1, dimnames = list(c("group", "age"), NULL)),
      t = matrix(c(3, 1), 2, 1, dimnames = list(c("group", "age"), NULL)),
      z = matrix(c(2.5, 1), 2, 1, dimnames = list(c("group", "age"), NULL)),
      z_contrast = 1.5,
      t_contrast = 2
    ),
    ginfo,
    target_sign = -1,
    source_sign = 1
  )
  expect_equal(unname(signed$beta["group", 1]), -2)
  expect_equal(signed$z_contrast, -1.5)

  w_unnamed <- fmrireg:::.fmri_ttest_resolve_contrast(c(0, 1), c("a", "b"))
  expect_equal(as.numeric(w_unnamed), c(0, 1))
  expect_error(
    fmrireg:::.fmri_ttest_resolve_contrast(c(1, 0, 0), c("a", "b")),
    "length equal"
  )
  # resolve_contrast expects canonical coef names (as fmri_ttest passes)
  w_named <- fmrireg:::.fmri_ttest_resolve_contrast(
    c(group = 1), c("group", "age"), ginfo
  )
  expect_equal(unname(w_named["group"]), 1)
  w_alias <- fmrireg:::.fmri_ttest_resolve_contrast(
    c(groupB = 1), c("group", "age"), ginfo
  )
  expect_equal(unname(w_alias["group"]), 1)
  expect_error(
    fmrireg:::.fmri_ttest_resolve_contrast(c(missing = 1), c("a", "b")),
    "Unknown contrast"
  )

  raw_w <- fmrireg:::.fmri_ttest_raw_contrast_weights(
    c(group = 1, age = 0),
    coef_names = c("groupB", "age"),
    group_info = ginfo,
    target_sign = -1,
    source_sign = 1
  )
  expect_equal(unname(raw_w["groupB"]), -1)

  single <- fmrireg:::.fmri_ttest_single_coef_contrast(
    list(
      beta = rbind(group = c(2, 3), age = c(1, 1)),
      se = rbind(group = c(0.5, 0.4), age = c(0.2, 0.2)),
      t = rbind(group = c(4, 5), age = c(1, 1)),
      z = rbind(group = c(3, 4), age = c(1, 1)),
      p = rbind(group = c(0.01, 0.02), age = c(0.1, 0.1)),
      df = rbind(group = c(10, 10), age = c(10, 10))
    ),
    weights = c(group = -1, age = 0)
  )
  expect_equal(unname(single$estimate), c(-2, -3))
  expect_equal(unname(single$se), c(0.5, 0.4))
  expect_equal(unname(single$t), c(-4, -5))
  expect_null(
    fmrireg:::.fmri_ttest_single_coef_contrast(
      list(beta = matrix(1, 2, 1)),
      weights = c(1, -1)
    )
  )

  Xols <- cbind(1, c(0, 0, 1, 1))
  Yols <- matrix(rnorm(4 * 3), 4, 3)
  ols <- fmrireg:::fmri_ols_fit(Yols, Xols)
  exact <- fmrireg:::.fmri_ttest_exact_ols_contrast(ols, Xols, c(0, 1))
  expect_equal(length(exact$estimate), 3L)
  expect_true(all(is.finite(exact$se)))

  stored <- fmrireg:::.fmri_ttest_store_contrast(list(beta = 1), exact)
  expect_equal(stored$beta_contrast, exact$estimate)
  expect_identical(fmrireg:::.fmri_ttest_store_contrast(list(a = 1), NULL), list(a = 1))
})

test_that("fmri_meta_gds early validation branches", {
  expect_error(
    fmrireg:::.fmri_meta_group_data_gds_impl(
      structure(list(), class = c("group_data_gds", "group_data")),
      weights = "custom",
      weights_custom = NULL,
      verbose = FALSE
    ),
    "weights_custom"
  )

  # combine + non-intercept formula fails before needing real assays when
  # t assay path is entered; construct a minimal stub that short-circuits
  # via missing compute by using combine with formula error on fake gd.
  # The formula check runs only when combine is set and "t" is in assay names.
  # Without fmrigds assays, we still cover the custom-weights gate above and
  # the method/robust match.arg paths via a dry call that errors later.
  skip_if_not_installed("fmrigds")

  # If fmrigds is present but data is empty, expect a downstream error rather
  # than a crash; custom-weights path already asserted.
  expect_error(
    fmrireg:::.fmri_meta_group_data_gds_impl(
      structure(list(), class = c("group_data_gds", "group_data")),
      formula = ~ 1,
      method = "fe",
      verbose = FALSE
    )
  )
})

test_that("detect_h5_file_type and subject id extraction", {
  expect_equal(
    fmrireg:::detect_h5_file_type(list(labels = c("beta", "se", "tstat"))),
    "by_stat"
  )
  expect_equal(
    fmrireg:::detect_h5_file_type(list(labels = c("beta_1", "beta_2"))),
    "betas"
  )
  expect_equal(
    fmrireg:::detect_h5_file_type(list(labels = c("FacesVsPlaces", "GoVsNoGo"))),
    "by_contrast"
  )
  expect_equal(
    fmrireg:::detect_h5_file_type(list()),
    "by_stat"
  )

  ids <- fmrireg:::extract_subject_ids_from_paths(c(
    "/data/sub-01_task-nback_stat-beta.h5",
    "/data/sub-02_task-nback_stat-beta.h5"
  ))
  expect_equal(ids, c("sub-01", "sub-02"))

  ids2 <- fmrireg:::extract_subject_ids_from_paths(c(
    "/tmp/subjA_stat.h5",
    "/tmp/subjB_stat.h5"
  ))
  expect_equal(ids2, c("subjA_stat", "subjB_stat"))
})
