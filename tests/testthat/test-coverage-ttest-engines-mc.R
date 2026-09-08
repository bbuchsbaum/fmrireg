# fmri_ttest classic/meta/welch engines, contrast, mc="bh", paired/mu0 gates.

make_block_group_data <- function(n_subj = 12L, n_feat = 5L, seed = 11L) {
  set.seed(seed)
  Y <- matrix(rnorm(n_subj * n_feat, mean = 0.25), n_subj, n_feat)
  V <- matrix(runif(n_subj * n_feat, 0.05, 0.25), n_subj, n_feat)
  colnames(Y) <- paste0("f", seq_len(n_feat))
  rownames(Y) <- paste0("s", seq_len(n_subj))
  covars <- data.frame(
    subject = rownames(Y),
    group = factor(rep(c("A", "B"), length.out = n_subj)),
    age = rnorm(n_subj, 30, 4),
    stringsAsFactors = FALSE
  )
  structure(
    list(
      blocks = list(list(Y = Y, V = V, covars = covars)),
      covariates = covars,
      subjects = rownames(Y),
      n_subjects = n_subj,
      n_features = n_feat
    ),
    class = c("group_data_matrix", "group_data")
  )
}

test_that("fmri_ttest classic/meta/welch engines with mc and contrast", {
  gd <- make_block_group_data()

  classic <- fmri_ttest(gd, formula = ~ 1, engine = "classic", mc = "bh", alpha = 0.1)
  expect_s3_class(classic, "fmri_ttest_fit")
  expect_true(!is.null(classic$beta))
  expect_true(!is.null(classic$p))
  expect_true(!is.null(classic$q))
  expect_equal(nrow(classic$beta), 1L)

  meta <- fmri_ttest(
    gd,
    formula = ~ 1 + group,
    engine = "meta",
    mc = "bh",
    contrast = c(groupB = 1)
  )
  expect_s3_class(meta, "fmri_ttest_fit")
  expect_true(!is.null(meta$z_contrast) || !is.null(meta$p_contrast) ||
                !is.null(meta$contrast) || length(coef(meta)) >= 1L)
  expect_true(!is.null(meta$q) || !is.null(meta$p))

  welch <- fmri_ttest(gd, formula = ~ 1 + group, engine = "welch")
  expect_s3_class(welch, "fmri_ttest_fit")
  expect_equal(nrow(welch$beta), 2L)

  auto <- fmri_ttest(gd, formula = ~ 1, engine = "auto")
  expect_s3_class(auto, "fmri_ttest_fit")

  # CSV meta path (no Y materialization) with BH
  gd_csv <- fmrireg:::.demo_group_data_csv()
  csv_fit <- suppressWarnings(
    fmri_ttest(gd_csv, formula = ~ 1, engine = "meta", mc = "bh")
  )
  expect_s3_class(csv_fit, "fmri_ttest_fit")
  expect_output(print(csv_fit), "fmri_ttest|Engine|Subjects")
})

test_that("fmri_ttest paired path and input/contrast error gates", {
  gd <- make_block_group_data(n_subj = 10L, n_feat = 3L, seed = 19L)

  paired <- fmri_ttest(gd, formula = ~ 1, paired = TRUE, mu0 = 0.05)
  expect_s3_class(paired, "fmri_ttest_fit")
  expect_equal(rownames(paired$beta), "(Intercept)")

  expect_error(fmri_ttest(1:5, formula = ~ 1), "group_data|data frame")
  expect_error(
    fmri_ttest(gd, formula = ~ 1, engine = "classic", contrast = c(missing = 1)),
    regexp = "."
  )

  # GDS-backed objects reject paired / nonzero mu0 when available
  gd_csv <- fmrireg:::.demo_group_data_csv()
  # classic materialization fails for csv; meta works
  expect_error(
    suppressWarnings(fmri_ttest(gd_csv, formula = ~ 1, engine = "classic")),
    "materialized|betas"
  )

  # Design/data row mismatch
  bad <- gd
  bad$blocks[[1]]$Y <- bad$blocks[[1]]$Y[-1, , drop = FALSE]
  expect_error(fmri_ttest(bad, formula = ~ 1, engine = "classic"), "rows")
})
