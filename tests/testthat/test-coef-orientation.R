# Orientation and naming contract for the fmri_lm accessor family
# (bbuchsbaum/fmrireg#227 and #217).
#
# Every accessor returns voxels x terms: one row per voxel, one column per
# coefficient/contrast, term names on the column margin. Values are checked
# against betas that are known exactly because the data are generated from the
# fitted design matrix itself.

.orientation_fit <- function(n_vox, noise_sd = 1e-6, seed = 227) {
  set.seed(seed)
  ev <- data.frame(
    onset = seq(10, 190, by = 20),
    cond = factor(rep(c("A", "B"), 5)),
    run = 1
  )
  form <- onset ~ hrf(cond, contrasts = pair_contrast(~ cond == "A", ~ cond == "B", name = "AvB"))

  # First pass only to obtain the design matrix that fmri_lm() will use.
  probe <- fmri_lm(form, block = ~run,
                   dataset = matrix_frame(matrix(rnorm(100 * n_vox), 100, n_vox),
                                          TR = 2, run_length = 100, event_table = ev))
  X <- as.matrix(design_matrix(probe$model))
  p <- ncol(X)
  ev_idx <- probe$result$event_indices

  # Known voxels x coefficients truth; event betas differ per voxel so that a
  # transposition cannot go unnoticed.
  B <- matrix(rnorm(n_vox * p, sd = 0.5), n_vox, p)
  B[, ev_idx] <- outer(seq_len(n_vox), seq_along(ev_idx), function(v, k) v + 10 * k)
  colnames(B) <- colnames(X)
  Y <- X %*% t(B) + matrix(rnorm(100 * n_vox, sd = noise_sd), 100, n_vox)

  fit <- fmri_lm(form, block = ~run,
                 dataset = matrix_frame(Y, TR = 2, run_length = 100, event_table = ev))
  list(fit = fit, B = B, X = X, ev_idx = ev_idx)
}

for (n_vox in c(3L, 2L, 1L)) {
  # n_vox = 2 equals the number of event coefficients: orientation cannot be
  # inferred from dim() there, so it must be right by construction.
  test_that(sprintf("coef.fmri_lm() is voxels x terms by default (%d voxel(s))", n_vox), {
    s <- .orientation_fit(n_vox)
    b <- coef(s$fit)
    expect_true(is.matrix(b))
    expect_identical(dim(b), c(n_vox, length(s$ev_idx)))
    expect_identical(colnames(b), colnames(s$X)[s$ev_idx])
    expect_equal(unname(b), unname(s$B[, s$ev_idx, drop = FALSE]), tolerance = 1e-4)
  })

  test_that(sprintf("coef.fmri_lm(include_baseline = TRUE) is voxels x terms (%d voxel(s))", n_vox), {
    s <- .orientation_fit(n_vox)
    bb <- coef(s$fit, include_baseline = TRUE)
    expect_true(is.matrix(bb))
    expect_identical(dim(bb), c(n_vox, ncol(s$X)))
    expect_identical(colnames(bb), colnames(s$X))
    expect_equal(unname(bb), unname(s$B), tolerance = 1e-4)
    # Event block of the full matrix is exactly the default result.
    expect_equal(bb[, s$ev_idx, drop = FALSE], coef(s$fit))
  })

  test_that(sprintf("coef.fmri_lm(type = 'contrasts') is voxels x contrasts (%d voxel(s))", n_vox), {
    s <- .orientation_fit(n_vox)
    cc <- as.matrix(coef(s$fit, type = "contrasts"))
    expect_identical(dim(cc), c(n_vox, 1L))
    expect_identical(colnames(cc), "AvB")
    truth <- s$B[, s$ev_idx[1]] - s$B[, s$ev_idx[2]]
    expect_equal(unname(cc[, 1]), unname(truth), tolerance = 1e-4)
  })

  test_that(sprintf("stats/standard_error/p_values share coef()'s orientation (%d voxel(s))", n_vox), {
    s <- .orientation_fit(n_vox, noise_sd = 0.1)
    b <- coef(s$fit)
    se <- as.matrix(standard_error(s$fit, type = "betas"))
    tt <- as.matrix(stats(s$fit, type = "betas"))
    pp <- as.matrix(p_values(s$fit, type = "betas"))
    for (m in list(se, tt, pp)) {
      expect_identical(dim(m), dim(b))
      expect_identical(colnames(m), colnames(b))
    }
    expect_equal(unname(tt), unname(b / se), tolerance = 1e-6)
    expect_identical(colnames(b), coef_names(s$fit))
  })
}

test_that("stats(type = 'betas') returns t-statistics, not estimates (#217)", {
  s <- .orientation_fit(3L, noise_sd = 0.1)
  b <- coef(s$fit)
  tt <- as.matrix(stats(s$fit, type = "betas"))
  se <- as.matrix(standard_error(s$fit, type = "betas"))
  expect_equal(unname(tt), unname(b / se), tolerance = 1e-6)
  expect_false(isTRUE(all.equal(unname(tt), unname(b))))
  # The default family is "betas".
  expect_equal(as.matrix(stats(s$fit)), tt)
  expect_equal(as.matrix(standard_error(s$fit)), se)
})

test_that("stats(type = 'estimates') is deprecated in favour of 'betas' (#217)", {
  s <- .orientation_fit(3L, noise_sd = 0.1)
  withr::local_options(fmrireg.suppress_deprecation = FALSE)
  expect_warning(old <- stats(s$fit, type = "estimates"),
                 class = "fmrireg_deprecated_estimates_type")
  expect_equal(old, stats(s$fit, type = "betas"))
  # standard_error()/p_values() keep "estimates" as a silent synonym: the
  # quantity they return is named by the function, so it cannot be misread.
  expect_no_warning(se_old <- standard_error(s$fit, type = "estimates"))
  expect_equal(se_old, standard_error(s$fit, type = "betas"))
  expect_no_warning(p_old <- p_values(s$fit, type = "estimates"))
  expect_equal(p_old, p_values(s$fit, type = "betas"))

  withr::local_options(fmrireg.suppress_deprecation = TRUE)
  expect_no_warning(stats(s$fit, type = "estimates"))
})

test_that("fit_contrasts() is correct when voxel count equals coefficient count", {
  # With 2 voxels and 2 event coefficients the old terms x voxels default of
  # coef() had the same dim() as voxels x terms, and fit_contrasts() silently
  # used the transpose.
  s <- .orientation_fit(2L)
  con <- pair_contrast(~ cond == "A", ~ cond == "B", name = "AvB2")
  res <- fit_contrasts(s$fit, list(con))
  truth <- s$B[, s$ev_idx[1]] - s$B[, s$ev_idx[2]]
  expect_equal(as.numeric(res[[1]]$estimate), unname(truth), tolerance = 1e-4)
})
