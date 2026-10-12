# Review 2026-10, pass 3, stream R3-4.
#
# * compute_lm_contrasts(drop_failed = TRUE) dropped a misspelled contrast
#   without a word; it now warns and names every dropped contrast.
# * estimate_betas() stored the user's `basemod` argument (usually NULL)
#   instead of the baseline model it actually regressed out.
# * glm_lss() reported "Cholesky factorization failed" whenever any trial beta
#   was non-finite (`all(...) || any(...)`), which is not what happened.

.r3_suffstats <- function() {
  set.seed(31)
  X <- cbind(1, rnorm(40), rnorm(40))
  colnames(X) <- c("Intercept", "A", "B")
  Y <- drop(X %*% c(1, 2, -1)) + matrix(rnorm(40 * 3), 40, 3)
  qrx <- qr(X)
  B <- qr.coef(qrx, Y)
  df <- nrow(X) - qrx$rank
  list(X = X, Y = Y, B = B, df = df,
       sigma2 = colSums(qr.resid(qrx, Y)^2) / df,
       XtXinv = chol2inv(qr.R(qrx)))
}

test_that("R3-4: drop_failed = TRUE warns and names each dropped contrast", {
  s <- .r3_suffstats()
  expect_warning(
    out <- compute_lm_contrasts(
      B = s$B, XtXinv = s$XtXinv, df = s$df, sigma2 = s$sigma2,
      contrasts = list(good = c(A = 1, B = -1), typo = c(A = 1, BB = -1)),
      f_contrasts = list(badF = matrix(1, 2, 2, dimnames = list(NULL, c("A", "Z"))))
    ),
    "typo.*badF|badF.*typo"
  )
  expect_identical(unique(out$name), "good")
  # the surviving contrast is correct (vs lm)
  m <- lm(s$Y[, 1] ~ 0 + s$X)
  cv <- c(0, 1, -1)
  expect_equal(out$data[[1]]$estimate[1], sum(cv * coef(m)), tolerance = 1e-10)

  # no warning when nothing fails
  expect_no_warning(compute_lm_contrasts(
    B = s$B, XtXinv = s$XtXinv, df = s$df, sigma2 = s$sigma2,
    contrasts = list(good = c(A = 1, B = -1))
  ))
  # drop_failed = FALSE still errors
  expect_error(compute_lm_contrasts(
    B = s$B, XtXinv = s$XtXinv, df = s$df, sigma2 = s$sigma2,
    contrasts = list(typo = c(BB = 1)), drop_failed = FALSE
  ), "typo")
})

.r3_lss_setup <- function() {
  set.seed(7)
  etab <- data.frame(onset = seq(5, 70, by = 8), cond = factor("a"), run = 1L)
  Y <- matrix(rnorm(80 * 3), 80, 3)
  dset <- matrix_frame(Y, TR = 1, run_length = 80, event_table = etab)
  sframe <- fmrihrf::sampling_frame(blocklens = 80, TR = 1)
  emod <- event_model(onset ~ hrf(cond), data = etab, block = ~ run,
                      sampling_frame = sframe)
  list(dset = dset, emod = emod, sframe = sframe, etab = etab)
}

test_that("R3-4: estimate_betas returns the baseline model it used", {
  s <- .r3_lss_setup()
  res <- estimate_betas(s$dset, ran = onset ~ trialwise(), block = ~ run,
                        method = "lss", progress = FALSE)
  expect_s3_class(res$basemod, "baseline_model")
  expect_equal(unname(as.matrix(design_matrix(res$basemod))),
               unname(as.matrix(res$design_base)))

  bmod <- baseline_model(basis = "poly", degree = 2, sframe = s$sframe)
  res2 <- estimate_betas(s$dset, ran = onset ~ trialwise(), block = ~ run,
                         method = "lss", basemod = bmod, progress = FALSE)
  expect_identical(res2$basemod, bmod)
  expect_equal(unname(as.matrix(design_matrix(res2$basemod))),
               unname(as.matrix(res2$design_base)))
})

test_that("R3-4: glm_lss reports non-finite trial betas accurately", {
  s <- .r3_lss_setup()
  bad <- matrix(c(1, NaN, 2, Inf, 3, 4), 2, 3)
  local_mocked_bindings(
    estimate_betas = function(...) structure(list(betas_ran = bad),
                                             class = "fmri_betas"),
    .package = "fmrireg"
  )
  err <- tryCatch(glm_lss(s$dset, s$emod, fmrihrf::HRF_SPMG1, progress = FALSE),
                  error = conditionMessage)
  expect_true(is.character(err))
  expect_no_match(err, "Cholesky", fixed = TRUE)
  expect_match(err, "non-finite", fixed = TRUE)
  expect_match(err, "2 of 6", fixed = TRUE)

  good <- matrix(1:6 + 0, 2, 3)
  local_mocked_bindings(
    estimate_betas = function(...) structure(list(betas_ran = good),
                                             class = "fmri_betas"),
    .package = "fmrireg"
  )
  expect_s3_class(glm_lss(s$dset, s$emod, fmrihrf::HRF_SPMG1, progress = FALSE),
                  "fmri_betas")
})
