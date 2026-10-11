# Review 2026-10, stream C, C2: tidy(fit, type = "contrasts") returned only
# the first contrast.

.tidy_review_fit <- function() {
  set.seed(11)
  n <- 90L
  etab <- data.frame(onset = seq(4, 78, by = 6),
                     cond = factor(rep(c("A", "B", "C"), length.out = 13)),
                     run = 1L)
  Y <- matrix(rnorm(n * 4), n, 4)
  dset <- matrix_frame(Y, TR = 1, run_length = n, event_table = etab)
  cons <- contrast_set(
    pair_contrast(~ cond == "A", ~ cond == "B", name = "AvB"),
    pair_contrast(~ cond == "B", ~ cond == "C", name = "BvC"),
    oneway_contrast(~ cond, name = "main")
  )
  fit <- fmri_lm(onset ~ hrf(cond, contrasts = cons), block = ~ run, dataset = dset)
  list(fit = fit, Y = Y)
}

test_that("C2: tidy(type = 'contrasts') returns every t- and F-contrast", {
  s <- .tidy_review_fit()
  fit <- s$fit
  td <- tidy(fit, type = "contrasts")
  V <- ncol(s$Y)

  expect_setequal(unique(td$term), c("AvB", "BvC", "main"))
  expect_equal(nrow(td), 3L * V)

  for (nm in c("AvB", "BvC")) {
    rows <- td[td$term == nm, ]
    rows <- rows[order(rows$voxel), ]
    expect_equal(rows$voxel, seq_len(V))
    expect_equal(rows$estimate, coef(fit, type = "contrasts")[[nm]])
    expect_equal(rows$std_error, standard_error(fit, type = "contrasts")[[nm]])
    expect_equal(rows$statistic, stats(fit, type = "contrasts")[[nm]])
    expect_equal(rows$p_value, p_values(fit, type = "contrasts")[[nm]])
  }
  Frows <- td[td$term == "main", ]
  Frows <- Frows[order(Frows$voxel), ]
  expect_equal(Frows$statistic, stats(fit, type = "F")[["main"]])

  # Reference values from stats::lm on the same design.
  X <- as.matrix(design_matrix(fit$model))
  ev <- fit$result$event_indices
  for (v in seq_len(V)) {
    m <- lm(s$Y[, v] ~ 0 + X)
    l <- numeric(ncol(X)); l[ev[1]] <- 1; l[ev[2]] <- -1
    est <- sum(l * coef(m))
    se <- sqrt(drop(t(l) %*% vcov(m) %*% l))
    r <- td[td$term == "AvB" & td$voxel == v, ]
    expect_equal(r$estimate, est, tolerance = 1e-6)
    expect_equal(r$statistic, est / se, tolerance = 1e-6)
    m0 <- lm(s$Y[, v] ~ 0 + X[, -ev] + I(X[, ev[1]] + X[, ev[2]] + X[, ev[3]]))
    a <- anova(m0, m)
    r <- td[td$term == "main" & td$voxel == v, ]
    expect_equal(r$statistic, a$F[2], tolerance = 1e-6)
    expect_equal(r$p_value, a$`Pr(>F)`[2], tolerance = 1e-6)
  }
  expect_true("type" %in% names(td))
  expect_equal(unique(td$type[td$term == "main"]), "Fcontrast")
})
