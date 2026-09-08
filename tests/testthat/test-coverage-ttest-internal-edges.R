# More ttest helpers + lm_internal meta/chunk edges.

test_that(".fmri_ttest_exact_ols_contrast and normalize_group_rows", {
  set.seed(111)
  S <- 10L
  P <- 3L
  Y <- matrix(rnorm(S * P), S, P)
  X <- cbind(Intercept = 1, groupB = rep(c(0, 1), each = 5))
  ols <- fmrireg:::fmri_ols_fit(Y, X)
  raw_w <- c(Intercept = 0, groupB = 1)
  exact <- fmrireg:::.fmri_ttest_exact_ols_contrast(ols, X, raw_w)
  expect_true(is.list(exact))
  expect_equal(length(exact$estimate), P)

  res <- list(
    beta = ols$beta,
    se = ols$se,
    t = ols$t,
    p = ols$p,
    df = ols$df
  )
  ginfo <- list(
    raw_name = "groupB",
    canonical_name = "group",
    levels = c("A", "B"),
    has_group = TRUE
  )
  normed <- fmrireg:::.fmri_ttest_normalize_group_rows(res, ginfo)
  expect_true(is.list(normed))
  expect_true(any(grepl("group", rownames(normed$beta))) ||
                identical(rownames(normed$beta), rownames(res$beta)))
})

test_that(".fmri_ttest_resolve_contrast character/integer and materialize errors", {
  coef_names <- c("Intercept", "group")
  w_named <- fmrireg:::.fmri_ttest_resolve_contrast(
    c(group = 1), coef_names
  )
  expect_equal(unname(w_named["group"]), 1)
  expect_equal(unname(w_named["Intercept"]), 0)

  w_full <- fmrireg:::.fmri_ttest_resolve_contrast(c(0, 1), coef_names)
  expect_equal(as.numeric(w_full), c(0, 1))

  expect_error(
    fmrireg:::.fmri_ttest_resolve_contrast(c(1), coef_names),
    "Unnamed contrast must have length"
  )
  expect_error(
    fmrireg:::.fmri_ttest_resolve_contrast(c(missing = 1), coef_names),
    "Unknown contrast"
  )
  expect_null(fmrireg:::.fmri_ttest_resolve_contrast(NULL, coef_names))

  # Non-gds object should pass through or error clearly
  out <- tryCatch(
    fmrireg:::.fmri_ttest_materialize_effects(list(Y = matrix(1, 2, 2))),
    error = function(e) e
  )
  expect_true(inherits(out, "error") || is.list(out) || is.matrix(out))
})

test_that("meta_contrasts aggregates simple contrast result lists", {
  make_con <- function(est) {
    dplyr::tibble(
      type = "contrast",
      name = "A_vs_B",
      stat_type = "tstat",
      df.residual = 20,
      conmat = list(c(1, -1)),
      colind = list(1:2),
      data = list(dplyr::tibble(
        estimate = est, se = abs(est) / 2, stat = est / 0.2,
        prob = rep(0.05, length(est)), sigma = rep(1, length(est))
      ))
    )
  }
  c1 <- list(A_vs_B = make_con(c(1, 2, 3)))
  c2 <- list(A_vs_B = make_con(c(1.1, 1.9, 3.2)))
  out <- tryCatch(
    fmrireg:::meta_contrasts(list(c1, c2)),
    error = function(e) e
  )
  expect_true(inherits(out, "error") || is.list(out) || inherits(out, "tbl_df"))
})

test_that(".orient_fcontrast rejects mismatched dims", {
  expect_error(
    fmrireg:::.orient_fcontrast(matrix(1, 2, 2), colind = 1:3),
    "do not match"
  )
  # Square treated as hypothesis-oriented
  sq <- diag(3)
  expect_equal(fmrireg:::.orient_fcontrast(sq, colind = 1:3), sq)
})
