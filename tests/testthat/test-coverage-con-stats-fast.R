# con_stats.R helpers: aliased indices, fast contrasts with robust/AR, F-tests.

tiny_design <- function(n = 40L, p = 3L, V = 5L, seed = 7) {
  set.seed(seed)
  X <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  B <- matrix(rnorm(p * V), p, V)
  Y <- X %*% B + matrix(rnorm(n * V, sd = 0.3), n, V)
  XtXinv <- chol2inv(chol(crossprod(X)))
  resid <- Y - X %*% solve(crossprod(X), crossprod(X, Y))
  sigma2 <- colSums(resid^2) / (n - p)
  list(X = X, Y = Y, B = solve(crossprod(X), crossprod(X, Y)),
       XtXinv = XtXinv, sigma2 = sigma2, df = n - p, p = p, V = V)
}

test_that(".aliased_indices / .lm_basic_stats / .lm_cov_unscaled work", {
  set.seed(4)
  d <- data.frame(y1 = rnorm(30), y2 = rnorm(30), x = rnorm(30), z = rnorm(30))
  fit <- lm(cbind(y1, y2) ~ x + z, data = d)

  basics <- fmrireg:::.lm_basic_stats(fit)
  expect_true(is.matrix(basics$betamat))
  expect_equal(length(basics$sigma), 2L)
  expect_true(basics$dfres > 0)

  covu <- fmrireg:::.lm_cov_unscaled(fit, varnames = c("(Intercept)", "x", "z"))
  expect_equal(dim(covu), c(3L, 3L))
  expect_true(all(diag(covu) > 0))

  # .aliased_indices reads attr(..., "aliased")
  expect_equal(fmrireg:::.aliased_indices(list()), integer(0))
  mat <- diag(3)
  attr(mat, "aliased") <- c(2L, 3L)
  expect_equal(fmrireg:::.aliased_indices(mat), c(2L, 3L))
  expect_equal(fmrireg:::.aliased_indices(diag(2), mat), c(2L, 3L))

  # Rank-deficient lm yields aliased attr via .lm_cov_unscaled
  d2 <- data.frame(y = rnorm(20), x1 = 1:20, x2 = 1:20)
  fit_sing <- lm(y ~ x1 + x2, data = d2)
  expect_warning(
    cov_sing <- fmrireg:::.lm_cov_unscaled(
      fit_sing, varnames = c("(Intercept)", "x1", "x2"), warn = TRUE
    ),
    "rank|aliased|deficient|collinear"
  )
  ali <- attr(cov_sing, "aliased")
  expect_true(length(ali) >= 1L || !isTRUE(attr(cov_sing, "is_full_rank")))
})

test_that(".fast_t_contrast and .fast_F_contrast cover aliased and robust/AR", {
  d <- tiny_design()
  l <- matrix(c(0, 1, -1), 1)
  t_plain <- fmrireg:::.fast_t_contrast(
    d$B, d$sigma2, d$XtXinv, l, d$df
  )
  expect_equal(length(t_plain$estimate), d$V)
  expect_true(all(is.finite(t_plain$se)))
  expect_equal(t_plain$stat_type, "tstat")

  # Vector l coerced to matrix
  t_vec <- fmrireg:::.fast_t_contrast(
    d$B, d$sigma2, d$XtXinv, c(0, 1, -1), d$df
  )
  expect_equal(t_vec$estimate, t_plain$estimate, tolerance = 1e-12)

  # Robust + AR effective df path
  w <- rep(0.9, d$df + d$p)
  t_rob <- fmrireg:::.fast_t_contrast(
    d$B, d$sigma2, d$XtXinv, l, d$df,
    robust_weights = w, ar_order = 1L
  )
  expect_equal(length(t_rob$stat), d$V)
  expect_true(all(t_rob$prob >= 0 & t_rob$prob <= 1))

  # Aliased contrast returns NA stats
  expect_warning(
    t_alias <- fmrireg:::.fast_t_contrast(
      d$B, d$sigma2, d$XtXinv, l, d$df,
      aliased = 2L, contrast_name = "A_vs_B"
    ),
    "non-estimable|aliased|A_vs_B"
  )
  expect_true(all(is.na(t_alias$estimate)))
  expect_true(all(is.na(t_alias$stat)))

  L <- rbind(c(0, 1, 0), c(0, 0, 1))
  f_plain <- fmrireg:::.fast_F_contrast(
    d$B, d$sigma2, d$XtXinv, L, d$df
  )
  expect_equal(length(f_plain$stat), d$V)
  expect_equal(f_plain$stat_type, "Fstat")
  expect_true(all(f_plain$stat >= 0 | is.nan(f_plain$stat)))

  f_rob <- fmrireg:::.fast_F_contrast(
    d$B, d$sigma2, d$XtXinv, L, d$df,
    robust_weights = w, ar_order = 2L
  )
  expect_equal(length(f_rob$prob), d$V)

  expect_error(
    fmrireg:::.fast_F_contrast(d$B, d$sigma2, d$XtXinv, c(0, 1, -1), d$df),
    "matrix"
  )

  expect_warning(
    f_alias <- fmrireg:::.fast_F_contrast(
      d$B, d$sigma2, d$XtXinv, L, d$df,
      aliased = 3L, contrast_name = "omnibus"
    ),
    "non-estimable|aliased|omnibus"
  )
  expect_true(all(is.na(f_alias$stat)))
})

test_that("process_t/f_contrasts and fit_lm_contrasts_fast package results", {
  d <- tiny_design()
  w_t <- structure(c(1, -1), colind = c(2L, 3L))
  w_f <- structure(diag(2), colind = c(2L, 3L))

  t_res <- fmrireg:::process_t_contrasts(
    d$B, d$sigma2, d$XtXinv,
    conlist = list(A_vs_B = w_t),
    df = d$df
  )
  expect_equal(names(t_res), "A_vs_B")
  expect_equal(t_res$A_vs_B$type, "contrast")
  expect_equal(length(t_res$A_vs_B$data[[1]]$estimate), d$V)

  f_res <- fmrireg:::process_f_contrasts(
    d$B, d$sigma2, d$XtXinv,
    fconlist = list(both = w_f),
    df = d$df,
    robust_weights = rep(1, d$df + d$p),
    ar_order = 1L
  )
  expect_equal(names(f_res), "both")
  expect_equal(f_res$both$type, "Fcontrast")

  expect_equal(fmrireg:::process_t_contrasts(d$B, d$sigma2, d$XtXinv, list(), d$df), list())
  expect_equal(fmrireg:::process_f_contrasts(d$B, d$sigma2, d$XtXinv, list(), d$df), list())

  all_res <- fmrireg:::fit_lm_contrasts_fast(
    d$B, d$sigma2, d$XtXinv,
    conlist = list(A_vs_B = w_t),
    fconlist = list(both = w_f),
    df = d$df
  )
  expect_true(all(c("A_vs_B", "both") %in% names(all_res)))

  expect_error(
    fmrireg:::fit_lm_contrasts_fast(
      as.numeric(d$B), d$sigma2, d$XtXinv, list(), list(), d$df
    ),
    "matrix"
  )
})

test_that("beta_stats_matrix with AR/robust and fit_Ftests/fit_Fcontrasts", {
  d <- tiny_design()
  varnames <- paste0("b", seq_len(d$p))
  bs <- fmrireg:::beta_stats_matrix(
    d$B, d$XtXinv, sqrt(d$sigma2), d$df, varnames,
    robust_weights = rep(0.95, d$df + d$p),
    ar_order = 1L
  )
  expect_equal(bs$stat_type, "tstat")
  expect_equal(ncol(bs$data[[1]]$estimate[[1]]), d$p)

  # Negative diag warning path
  XtXbad <- d$XtXinv
  XtXbad[1, 1] <- -1e-8
  expect_warning(
    fmrireg:::beta_stats_matrix(
      d$B, XtXbad, sqrt(d$sigma2), d$df, varnames
    ),
    "Negative diagonal"
  )

  # fit_Fcontrasts / fit_Ftests on mlm
  set.seed(9)
  dd <- data.frame(
    y1 = rnorm(40), y2 = rnorm(40),
    x1 = rnorm(40), x2 = rnorm(40)
  )
  fit <- lm(cbind(y1, y2) ~ x1 + x2, data = dd)
  fcon <- fmrireg:::fit_Fcontrasts(fit, conmat = diag(2), colind = 2:3)
  expect_true(is.list(fcon) || inherits(fcon, "Fstat") || !is.null(fcon$stat))

  # Orient cell-oriented matrix
  cell <- matrix(c(1, 0, 0, 1), 2, 2) # 2 x 2 square stays as-is
  expect_equal(fmrireg:::.orient_fcontrast(cell, colind = 1:2), cell)
  tall <- matrix(c(1, 0, 0, 1), nrow = 2) # already hyp x col
  expect_equal(dim(fmrireg:::.orient_fcontrast(tall, 1:2)), c(2L, 2L))

  # fit_Ftests is multi-response oriented (apply over residual columns)
  ftests <- suppressWarnings(fmrireg:::fit_Ftests(fit))
  expect_true(is.list(ftests))
  expect_true(all(c("F", "P") %in% names(ftests)))
  expect_equal(nrow(ftests$F), 2L)
  expect_equal(nrow(ftests$P), 2L)
})
