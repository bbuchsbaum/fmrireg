# Review 2026-10, stream B: p-value / z combination (B1, B3, B4) and the
# fixed-effects run pooling helper.

# --- B1: run-pooled F-contrasts use a real Stouffer z ----------------------

test_that("B1: meta_stouffer with k identical p-values gives sqrt(k) * qnorm(1 - p)", {
  for (k in c(2L, 4L, 9L)) {
    p <- c(0.2, 0.01, 1e-6)
    pval <- matrix(rep(p, k), nrow = length(p))
    se <- matrix(1, nrow = length(p), ncol = k)
    res <- fmrireg:::meta_stouffer(pval, se)
    z_expected <- sqrt(k) * stats::qnorm(p, lower.tail = FALSE)
    expect_equal(res$stat, z_expected, tolerance = 1e-12, info = paste("k =", k))
    expect_equal(res$prob, stats::pnorm(z_expected, lower.tail = FALSE),
                 tolerance = 1e-12, info = paste("k =", k))
  }
})

test_that("B1: meta_stouffer keeps very small p-values finite", {
  pval <- matrix(c(1e-20, 1e-30), nrow = 1)
  res <- fmrireg:::meta_stouffer(pval, matrix(1, 1, 2))
  z_expected <- sum(stats::qnorm(c(1e-20, 1e-30), lower.tail = FALSE)) / sqrt(2)
  expect_true(is.finite(res$stat))
  expect_equal(res$stat, z_expected, tolerance = 1e-12)
})

test_that("B1: run-pooled F p-values are calibrated under the null", {
  set.seed(20261011)
  n_vox <- 4000L
  k <- 4L
  make_row <- function(p) {
    tibble::tibble(
      type = "Fcontrast", name = "F1", conmat = list(NULL), colind = list(1L),
      data = list(tibble::tibble(prob = p))
    )
  }
  ftab <- dplyr::bind_rows(lapply(seq_len(k), function(i) make_row(stats::runif(n_vox))))
  out <- fmrireg:::meta_Fcontrasts(ftab)
  prob <- out$data[[1]]$prob
  rate <- mean(prob < 0.05)
  # Binomial band: 0.05 +/- 4 * sqrt(0.05 * 0.95 / n)
  band <- 4 * sqrt(0.05 * 0.95 / n_vox)
  expect_gt(rate, 0.05 - band)
  expect_lt(rate, 0.05 + band)
  # pooled z is standard normal under the null
  expect_equal(stats::sd(out$data[[1]]$stat), 1, tolerance = 0.05)
})

test_that("do_fixef keeps upper-tail precision for large z", {
  res <- fmrireg:::do_fixef(se = matrix(1, 1, 1), beta = matrix(30, 1, 1), "inv_var")
  expect_gt(res$prob, 0)
  expect_equal(res$prob, 2 * stats::pnorm(-30), tolerance = 1e-10)
})

# --- B3: one t -> z helper that does not overflow --------------------------

test_that("B3: t to z conversion is finite for extreme t and matches log-scale reference", {
  t <- c(-40, -3, 0, 2.5, 40)
  df <- 20
  ref <- -stats::qnorm(stats::pt(-abs(t), df = df, log.p = TRUE), log.p = TRUE) * sign(t)
  z <- fmrireg:::.t_to_z(t, df)
  expect_true(all(is.finite(z)))
  expect_equal(z, ref, tolerance = 1e-12)
  # Moderate values agree with the naive formula
  expect_equal(z[2:4], stats::qnorm(stats::pt(t[2:4], df)), tolerance = 1e-10)
})

test_that("B3: classic and Welch fmri_ttest engines give finite z for huge effects", {
  set.seed(1)
  S <- 22L
  Y <- cbind(100 + stats::rnorm(S, sd = 0.01), stats::rnorm(S))
  gd <- structure(list(blocks = list(list(Y = Y, V = NULL,
                                           covars = data.frame(subject = seq_len(S))))),
                  class = "group_data")
  fit <- fmri_ttest(gd, formula = ~ 1, engine = "classic")
  ref <- stats::t.test(Y[, 1])
  expect_equal(unname(fit$t[1, 1]), unname(ref$statistic), tolerance = 1e-8)
  expect_true(is.finite(fit$z[1, 1]))
  expect_equal(unname(fit$z[1, 1]),
               -stats::qnorm(stats::pt(-abs(unname(ref$statistic)), ref$parameter, log.p = TRUE),
                             log.p = TRUE), tolerance = 1e-10)

  group <- factor(rep(c("A", "B"), each = S / 2))
  Yw <- cbind(c(100 + stats::rnorm(S / 2, sd = 0.01), stats::rnorm(S / 2, sd = 0.01)),
              stats::rnorm(S))
  gdw <- structure(list(blocks = list(list(Y = Yw, V = NULL,
                                            covars = data.frame(group = group)))),
                   class = "group_data")
  fitw <- fmri_ttest(gdw, formula = ~ 1 + group, engine = "welch")
  expect_true(all(is.finite(fitw$z["group", ])))
  refw <- stats::t.test(Yw[group == "A", 1], Yw[group == "B", 1])
  expect_equal(unname(fitw$t["group", 1]), unname(refw$statistic), tolerance = 1e-8)
})

# --- B4: t-combination ignores missing subjects ----------------------------

test_that("B4: Stouffer/Lancaster/Fisher with NA subjects equal the result with them removed", {
  set.seed(42)
  S <- 8L
  P <- 5L
  tmat <- matrix(stats::rnorm(S * P, mean = 1), S, P)
  df <- seq(10, 24, length.out = S)
  w <- seq(0.5, 2, length.out = S)
  miss <- cbind(c(2, 5, 7, 1), c(1, 1, 3, 4))
  tna <- tmat
  tna[miss] <- NA

  for (method in c("stouffer", "lancaster", "fisher")) {
    for (wt in c("equal", "custom")) {
      if (method == "fisher" && wt == "custom") next
      full <- fmrireg:::combine_t_statistics(
        tna, df, method = method, weights = wt,
        weights_custom = if (wt == "custom") w else NULL
      )
      ref <- vapply(seq_len(P), function(j) {
        keep <- !is.na(tna[, j])
        fmrireg:::combine_t_statistics(
          tmat[keep, j, drop = FALSE], df[keep], method = method, weights = wt,
          weights_custom = if (wt == "custom") w[keep] else NULL
        )
      }, numeric(1))
      expect_equal(full, ref, tolerance = 1e-12, info = paste(method, wt))
    }
  }
})

test_that("B4: equal-weight Stouffer equals sum(z) / sqrt(k)", {
  tmat <- matrix(c(2, 1, -0.5, 3), ncol = 1)
  df <- c(15, 20, 30, 12)
  z_i <- stats::qnorm(stats::pt(tmat[, 1], df))
  expect_equal(
    fmrireg:::combine_t_statistics(tmat, df, method = "stouffer", weights = "equal"),
    sum(z_i) / 2, tolerance = 1e-12
  )
})
