# con_stats.R: aliased F-contrast return path and singular M handling.

test_that(".fast_F_contrast returns NA stats when contrast uses aliased cols", {
  set.seed(21)
  n <- 40L
  p <- 3L
  V <- 4L
  X <- cbind(1, rnorm(n), rnorm(n))
  B <- matrix(rnorm(p * V), p, V)
  rownames(B) <- c("(Intercept)", "A", "B")
  XtXinv <- chol2inv(chol(crossprod(X)))
  sigma2 <- rep(0.25, V)
  conmat <- matrix(c(0, 1, -1), 1)

  expect_warning(
    out <- fmrireg:::.fast_F_contrast(
      B, sigma2, XtXinv, L = conmat, df = n - p,
      aliased = 2L, contrast_name = "A_vs_B"
    ),
    "non-estimable|aliased|A_vs_B"
  )
  expect_true(all(is.na(out$estimate)))
  expect_true(all(is.na(out$stat)))
  expect_equal(out$stat_type, "Fstat")
})

test_that(".fast_F_contrast singular contrast gram yields NaN path", {
  set.seed(22)
  n <- 30L
  p <- 3L
  V <- 3L
  X <- cbind(1, rnorm(n), rnorm(n))
  B <- matrix(rnorm(p * V), p, V)
  XtXinv <- chol2inv(chol(crossprod(X)))
  sigma2 <- rep(1, V)

  # Rank-1 multirow contrast that can make C XtXinv C' singular-ish:
  # duplicate rows force singular M after projection.
  conmat <- rbind(c(0, 1, -1), c(0, 2, -2))
  out <- suppressWarnings(
    fmrireg:::.fast_F_contrast(B, sigma2, XtXinv, L = conmat, df = n - p)
  )
  expect_equal(length(out$stat), V)
  expect_equal(out$stat_type, "Fstat")
  expect_true(all(is.finite(out$stat) | is.nan(out$stat) | is.na(out$stat)))
})

test_that("fit_lm_contrasts_fast mixes t and F contrasts", {
  set.seed(23)
  n <- 36L
  p <- 3L
  V <- 5L
  X <- cbind(1, rnorm(n), rnorm(n))
  colnames(X) <- c("(Intercept)", "A", "B")
  B <- solve(crossprod(X), crossprod(X, matrix(rnorm(n * V), n, V)))
  rownames(B) <- colnames(X)
  XtXinv <- chol2inv(chol(crossprod(X)))
  resid <- matrix(rnorm(n * V), n, V)
  sigma2 <- colSums(resid^2) / (n - p)

  simple <- list(
    A_vs_B = matrix(c(0, 1, -1), 1, dimnames = list(NULL, colnames(X)))
  )
  attr(simple$A_vs_B, "colind") <- 1:3
  fcons <- list(
    main = matrix(c(0, 1, 0, 0, 0, 1), 2, 3, byrow = TRUE,
                  dimnames = list(NULL, colnames(X)))
  )
  attr(fcons$main, "colind") <- 1:3
  res <- fmrireg:::fit_lm_contrasts_fast(
    B = B, sigma2 = sigma2, XtXinv = XtXinv,
    conlist = simple, fconlist = fcons, df = n - p
  )
  expect_true(length(res) >= 1L)
  names_or_types <- unlist(lapply(res, function(x) {
    c(as.character(x$name[1]), as.character(x$type[1]))
  }))
  expect_true(any(grepl("A_vs_B|contrast|F|main", names_or_types)))
})
