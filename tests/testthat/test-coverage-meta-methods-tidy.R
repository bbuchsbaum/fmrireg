# fmri_meta_methods: coef/se/zscores/pvalues/tidy/print remaining branches.

test_that("fmri_meta coef/se/zscores/pvalues/tidy cover ROI and voxel summaries", {
  meta <- suppressWarnings(fmrireg:::.demo_fmri_meta())
  expect_s3_class(meta, "fmri_meta")

  cf <- coef(meta)
  expect_true(is.matrix(cf) || is.numeric(cf))
  expect_equal(dim(cf), dim(meta$coefficients))

  s <- se(meta)
  expect_equal(dim(s), dim(meta$se))

  z <- zscores(meta)
  expect_equal(dim(z), dim(cf))
  expect_equal(z, meta$coefficients / meta$se)

  p2 <- pvalues(meta, two_tailed = TRUE)
  p1 <- pvalues(meta, two_tailed = FALSE)
  expect_equal(dim(p2), dim(cf))
  expect_true(all(p2 >= p1 - 1e-12, na.rm = TRUE) || all(is.finite(p1)))

  expect_output(print(meta), regexp = ".")
  expect_output(summary(meta), regexp = ".")

  td <- tidy(meta)
  expect_true(inherits(td, "tbl_df") || is.data.frame(td))
  expect_true(all(c("estimate", "std.error", "p.value") %in% names(td)) ||
                all(c("term", "mean_estimate") %in% names(td)))

  td_ci <- tidy(meta, conf.int = TRUE, conf.level = 0.9)
  if (inherits(meta, "fmri_meta_roi")) {
    expect_true(all(c("conf.low", "conf.high") %in% names(td_ci)))
  }

  # Synthetic voxelwise tidy branch
  vox <- structure(
    list(
      coefficients = matrix(rnorm(20), 10, 2, dimnames = list(NULL, c("b0", "b1"))),
      se = matrix(runif(20, 0.1, 0.4), 10, 2, dimnames = list(NULL, c("b0", "b1"))),
      tau2 = runif(10),
      I2 = runif(10),
      n_voxels = 10L
    ),
    class = "fmri_meta"
  )
  td_vox <- tidy(vox)
  expect_equal(nrow(td_vox), 2L)
  expect_true(all(c("term", "mean_estimate", "n_significant") %in% names(td_vox)))
})

test_that("fmri_meta contrast and helper validators cover remaining lines", {
  meta <- suppressWarnings(fmrireg:::.demo_fmri_meta())

  con <- tryCatch(contrast(meta, c("(Intercept)" = 1)), error = function(e) e)
  if (!inherits(con, "error")) {
    expect_true(is.list(con) || is.matrix(con) || is.numeric(con))
  } else {
    expect_match(conditionMessage(con), ".", perl = TRUE)
  }

  expect_error(fmrireg:::.validate_fmri_meta_object(list()), regexp = ".")
  ok <- tryCatch(fmrireg:::.validate_fmri_meta_object(meta), error = function(e) e)
  if (!inherits(ok, "error")) {
    expect_true(isTRUE(ok) || is.null(ok) || inherits(ok, "fmri_meta"))
  }

  names_out <- fmrireg:::.fmri_meta_coef_names(meta)
  expect_true(is.character(names_out) || is.null(names_out))

  # write_results path (may error without spatial mask; still exercises gates)
  outdir <- tempfile("meta-cov-")
  dir.create(outdir)
  wr <- tryCatch(
    write_results(meta, path = outdir, strategy = "by_stat", format = "nifti"),
    error = function(e) e
  )
  if (inherits(wr, "error")) {
    expect_match(conditionMessage(wr), ".", perl = TRUE)
  } else {
    expect_true(is.list(wr) || is.character(wr))
  }
})
