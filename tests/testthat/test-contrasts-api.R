pivoted_qr_t_oracle <- function(X, Y, contrast) {
  fit <- stats::lm.fit(X, Y)
  p <- ncol(X)
  R <- qr.R(fit$qr)[seq_len(p), seq_len(p), drop = FALSE]
  cov_pivoted <- chol2inv(R)
  pivot <- fit$qr$pivot[seq_len(p)]
  cov_unscaled <- matrix(0, nrow = p, ncol = p)
  cov_unscaled[pivot, pivot] <- cov_pivoted

  estimate <- drop(crossprod(contrast, fit$coefficients))
  df <- nrow(X) - fit$rank
  sigma2 <- colSums(fit$residuals^2) / df
  variance_factor <- drop(crossprod(contrast, cov_unscaled %*% contrast))
  se <- sqrt(sigma2 * variance_factor)
  statistic <- estimate / se

  list(
    estimate = estimate,
    se = se,
    stat = statistic,
    prob = 2 * stats::pt(-abs(statistic), df = df),
    df = df
  )
}

expect_close_with_tolerance <- function(actual, expected,
                                        atol = 1e-12, rtol = 1e-10) {
  expect_true(all(is.finite(actual)))
  expect_true(all(is.finite(expected)))
  expect_true(all(abs(actual - expected) <= atol + rtol * abs(expected)))
}

subset_contrast_fixture <- function() {
  x <- c(-2.0, -1.5, -0.5, 0.25, 1.0, 1.75, 2.5)
  X <- cbind(intercept = 1, x = x)
  Y <- cbind(
    voxel_a = 1.25 + 2.0 * x + c(0.10, -0.08, 0.03, -0.04, 0.07, -0.06, -0.02),
    voxel_b = -0.75 + 0.4 * x + c(-0.05, 0.09, -0.02, 0.06, -0.08, 0.01, -0.01)
  )
  fit <- fmri_ols_fit(Y, X)
  sigma2 <- colSums((Y - X %*% fit$beta)^2) / fit$df

  list(X = X, Y = Y, fit = fit, sigma2 = sigma2)
}

expect_subset_contrasts_match <- function(result, oracle) {
  expect_setequal(result$name, c("named", "indexed", "full"))

  rows <- lapply(c("named", "indexed", "full"), function(nm) {
    result$data[[match(nm, result$name)]]
  })
  names(rows) <- c("named", "indexed", "full")

  for (field in c("estimate", "se", "stat", "prob")) {
    expect_close_with_tolerance(rows$named[[field]], rows$full[[field]])
    expect_close_with_tolerance(rows$indexed[[field]], rows$full[[field]])
    expect_close_with_tolerance(rows$named[[field]], oracle[[field]])
  }

  expect_equal(result$df.residual, rep(oracle$df, 3L))
  expect_equal(result$colind[[match("named", result$name)]], 2L)
  expect_equal(result$colind[[match("indexed", result$name)]], 2L)
}

test_that("subset t contrasts survive normalization and match a pivoted-QR oracle", {
  fixture <- subset_contrast_fixture()
  XtXinv <- solve(crossprod(fixture$X))
  oracle <- pivoted_qr_t_oracle(fixture$X, fixture$Y, c(0, 1))

  result <- compute_lm_contrasts(
    B = fixture$fit$beta,
    XtXinv = XtXinv,
    df = fixture$fit$df,
    sigma2 = fixture$sigma2,
    t_contrasts = list(
      named = c(x = 1),
      indexed = structure(1, colind = 2L),
      full = c(0, 1)
    ),
    columns = colnames(fixture$X)
  )

  expect_subset_contrasts_match(result, oracle)
})

test_that("sufficient-statistics API preserves subset t contrasts", {
  fixture <- subset_contrast_fixture()
  oracle <- pivoted_qr_t_oracle(fixture$X, fixture$Y, c(0, 1))

  result <- compute_lm_contrasts_from_suffstats(
    XtX = crossprod(fixture$X),
    XtS = crossprod(fixture$X, fixture$Y),
    StS = colSums(fixture$Y^2),
    df = fixture$fit$df,
    t_contrasts = list(
      named = c(x = 1),
      indexed = structure(1, colind = 2L),
      full = c(0, 1)
    ),
    columns = colnames(fixture$X)
  )

  expect_subset_contrasts_match(result, oracle)
})

test_that("invalid subset t-contrast metadata fails with precise errors", {
  fixture <- subset_contrast_fixture()
  common <- list(
    B = fixture$fit$beta,
    XtXinv = solve(crossprod(fixture$X)),
    df = fixture$fit$df,
    sigma2 = fixture$sigma2,
    columns = colnames(fixture$X),
    drop_failed = FALSE
  )

  expect_error(
    do.call(
      compute_lm_contrasts,
      c(common, list(t_contrasts = list(bad_name = c(unknown = 1))))
    ),
    "Names in t-contrast 'bad_name' not found in `columns`",
    fixed = TRUE
  )
  expect_error(
    do.call(
      compute_lm_contrasts,
      c(common, list(t_contrasts = list(bad_index = structure(1, colind = 3L))))
    ),
    "`colind` for t-contrast 'bad_index' must contain indices between 1 and p (2)",
    fixed = TRUE
  )
  expect_error(
    do.call(
      compute_lm_contrasts,
      c(common, list(t_contrasts = list(bad_length = structure(1, colind = 1:2))))
    ),
    "`colind` for t-contrast 'bad_length' must have one index per weight",
    fixed = TRUE
  )
})
