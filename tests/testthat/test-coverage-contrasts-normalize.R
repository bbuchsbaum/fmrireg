# contrasts_api .normalize_contrasts_for_engine name/colind branches.

test_that(".normalize_contrasts_for_engine covers named, colind, and F paths", {
  p <- 4L
  cols <- c("Intercept", "A", "B", "C")

  norm <- fmrireg:::.normalize_contrasts_for_engine(
    p = p,
    contrasts = list(
      named_t = c(A = 1, B = -1),
      full_t = c(0, 1, -1, 0),
      matrix_t = matrix(c(0, 1, 0, 0), nrow = 1),
      F_ab = matrix(
        c(0, 1, 0, 0, 0, 0, 1, 0),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(NULL, cols)
      )
    ),
    t_contrasts = list(
      indexed = structure(c(1, -1), colind = c(2L, 3L))
    ),
    f_contrasts = list(
      F_named = matrix(
        c(1, -1, 0),
        nrow = 1,
        dimnames = list(NULL, c("A", "B", "C"))
      ),
      F_full = diag(4)
    ),
    columns = cols,
    drop_failed = TRUE
  )

  expect_true(length(norm$t_list) >= 3L)
  expect_true(length(norm$f_list) >= 2L)
  expect_equal(attr(norm$t_list$named_t, "colind"), c(2L, 3L))
  expect_equal(attr(norm$t_list$indexed, "colind"), c(2L, 3L))
  expect_equal(attr(norm$t_list$full_t, "colind"), 1:4)
  expect_equal(length(attr(norm$f_list$F_named, "colind")), 3L)
})

test_that(".normalize_contrasts_for_engine errors and drop_failed branches", {
  p <- 3L
  cols <- c("a", "b", "c")

  expect_error(
    fmrireg:::.normalize_contrasts_for_engine(
      p, contrasts = c(1, 0, 0), t_contrasts = NULL, f_contrasts = NULL,
      columns = cols, drop_failed = FALSE
    ),
    "named list"
  )

  # Unnamed contrasts get synthetic names
  unnamed <- fmrireg:::.normalize_contrasts_for_engine(
    p,
    contrasts = list(c(a = 1), c(b = 1)),
    t_contrasts = NULL,
    f_contrasts = NULL,
    columns = cols,
    drop_failed = TRUE
  )
  # Force empty names path
  cons <- list(c(1, 0, 0), c(0, 1, 0))
  names(cons) <- NULL
  auto_named <- fmrireg:::.normalize_contrasts_for_engine(
    p, contrasts = cons, t_contrasts = NULL, f_contrasts = NULL,
    columns = cols, drop_failed = TRUE
  )
  expect_true(all(nzchar(names(auto_named$t_list))))

  expect_error(
    fmrireg:::.normalize_contrasts_for_engine(
      p, NULL,
      t_contrasts = list(bad = c(missing = 1)),
      f_contrasts = NULL,
      columns = cols,
      drop_failed = FALSE
    ),
    "not found|Names"
  )

  expect_error(
    fmrireg:::.normalize_contrasts_for_engine(
      p, NULL,
      t_contrasts = list(bad = structure(1, colind = 99L)),
      f_contrasts = NULL,
      columns = cols,
      drop_failed = FALSE
    ),
    "between 1 and p|colind"
  )

  expect_error(
    fmrireg:::.normalize_contrasts_for_engine(
      p, NULL,
      t_contrasts = list(bad = structure(c(1, -1), colind = c(1L, 1L))),
      f_contrasts = NULL,
      columns = cols,
      drop_failed = FALSE
    ),
    "duplicates"
  )

  expect_error(
    fmrireg:::.normalize_contrasts_for_engine(
      p, NULL,
      t_contrasts = list(bad = structure(1, colind = 1.5)),
      f_contrasts = NULL,
      columns = cols,
      drop_failed = FALSE
    ),
    "integer"
  )

  expect_error(
    fmrireg:::.normalize_contrasts_for_engine(
      p, NULL, NULL,
      f_contrasts = list(bad = matrix(1:2, nrow = 1)),
      columns = cols,
      drop_failed = FALSE
    ),
    "colnames|ncol"
  )

  # drop_failed swallows bad contrasts
  dropped <- fmrireg:::.normalize_contrasts_for_engine(
    p,
    contrasts = list(ok = c(0, 1, 0), bad = c(missing = 1)),
    t_contrasts = NULL,
    f_contrasts = NULL,
    columns = cols,
    drop_failed = TRUE
  )
  expect_true("ok" %in% names(dropped$t_list))
  expect_false("bad" %in% names(dropped$t_list))
})

test_that("compute_lm_contrasts exercises normalize through public API", {
  set.seed(3)
  X <- cbind(1, rnorm(30), rnorm(30))
  colnames(X) <- c("Intercept", "A", "B")
  Y <- as.vector(X %*% c(0.5, 1, -0.5)) + matrix(rnorm(30 * 2), 30, 2)
  fit <- fmri_ols_fit(Y, X)
  sigma2 <- colSums((Y - X %*% fit$beta)^2) / fit$df
  XtXinv <- solve(crossprod(X))

  out <- compute_lm_contrasts(
    B = fit$beta,
    XtXinv = XtXinv,
    df = fit$df,
    sigma2 = sigma2,
    contrasts = list(
      A_vs_0 = c(A = 1),
      F_AB = matrix(c(0, 1, 0, 0, 0, 1), nrow = 2, byrow = TRUE)
    ),
    columns = colnames(X)
  )
  expect_true(is.data.frame(out) || inherits(out, "tbl_df") || is.list(out))
})
