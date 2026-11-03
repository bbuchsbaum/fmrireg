# Simplification 2026-10, step P2: fmri_ttest() has one input normaliser and
# estimates with fmrigds's "ols:voxelwise" / "ttest:welch" reducers on every
# input, so gds data get the same engines, contrasts and options as in-memory
# blocks. Each test below failed on the code before P2.

skip_if_not_installed("fmrigds")

p2_data <- function(S = 14L, P = 6L, seed = 3) {
  set.seed(seed)
  covars <- data.frame(group = factor(rep(c("A", "B"), each = S / 2)),
                       age = stats::rnorm(S, 40, 8))
  rownames(covars) <- sprintf("s%02d", seq_len(S))
  Y <- matrix(stats::rnorm(S * P), S, P) + 0.7 * (covars$group == "B")
  V <- matrix(stats::runif(S * P, 0.2, 1), S, P)
  list(Y = Y, V = V, covars = covars)
}

p2_gds <- function(Y, V, covars, contrasts = NULL) {
  beta <- t(Y)
  var <- t(V)
  dimnames(beta) <- dimnames(var) <- list(paste0("v", seq_len(ncol(Y))), rownames(covars))
  if (!is.null(contrasts)) {
    beta <- array(c(beta, beta + 1), c(dim(beta), 2L),
                  dimnames = c(dimnames(beta), list(contrasts)))
    var <- array(c(var, var), c(dim(var), 2L), dimnames = dimnames(beta))
  }
  fmrigds::as_gds(list(beta = beta, var = var), col_data = covars)
}

p2_block <- function(Y, covars, V = NULL) {
  structure(list(blocks = list(list(Y = Y, V = V, covars = covars))),
            class = "group_data")
}

test_that("P2: classic and welch refuse gds data with several contrasts instead of using the first", {
  d <- p2_data()
  g <- p2_gds(d$Y, d$V, d$covars, contrasts = c("faces", "houses"))
  expect_error(fmri_ttest(g, ~ 1 + group, engine = "classic"),
               "one subject-level contrast.*faces, houses")
  expect_error(fmri_ttest(g, ~ 1 + group, engine = "welch"),
               "one subject-level contrast.*faces, houses")
})

test_that("P2: gds classic engine fits multi-coefficient contrasts exactly (vs lm)", {
  d <- p2_data()
  Y <- d$Y
  Y[3, 2] <- NA
  fit <- fmri_ttest(p2_gds(Y, d$V, d$covars), ~ 1 + group + age,
                    engine = "classic", contrast = c(group = 1, age = 2))
  for (j in seq_len(ncol(Y))) {
    keep <- is.finite(Y[, j])
    m <- stats::lm(y ~ 1 + group + age, data = cbind(d$covars[keep, ], y = Y[keep, j]))
    L <- c(0, -1, 2)  # default sign A - B, so group enters as -groupB
    est <- sum(L * stats::coef(m))
    se <- sqrt(drop(t(L) %*% stats::vcov(m) %*% L))
    expect_equal(unname(fit$beta_contrast[j]), est, tolerance = 1e-10, info = j)
    expect_equal(unname(fit$se_contrast[j]), se, tolerance = 1e-10, info = j)
    expect_equal(unname(fit$df_contrast[j]), m$df.residual, info = j)
    expect_equal(unname(fit$p_contrast[j]),
                 2 * stats::pt(abs(est / se), m$df.residual, lower.tail = FALSE),
                 tolerance = 1e-10, info = j)
  }
})

test_that("P2: welch engine matches t.test per feature, drops non-finite subjects, reports the SE", {
  d <- p2_data()
  Y <- d$Y
  Y[2, 1] <- NA
  Y[c(9, 12), 4] <- c(NaN, Inf)
  for (gd in list(p2_block(Y, d$covars), p2_gds(Y, d$V, d$covars))) {
    fit <- fmri_ttest(gd, ~ 1 + group, engine = "welch")
    for (j in seq_len(ncol(Y))) {
      ok <- is.finite(Y[, j])
      a <- Y[ok & d$covars$group == "A", j]
      b <- Y[ok & d$covars$group == "B", j]
      ref <- stats::t.test(a, b, var.equal = FALSE)  # A - B, the default sign
      expect_equal(unname(fit$beta["group", j]), mean(a) - mean(b), tolerance = 1e-10)
      expect_equal(unname(fit$beta["(Intercept)", j]), (mean(a) + mean(b)) / 2, tolerance = 1e-10)
      expect_equal(unname(fit$se["group", j]), unname(ref$stderr), tolerance = 1e-10)
      expect_equal(unname(fit$t["group", j]), unname(ref$statistic), tolerance = 1e-10)
      expect_equal(unname(fit$df["group", j]), unname(ref$parameter), tolerance = 1e-10)
      expect_equal(unname(fit$p["group", j]), ref$p.value, tolerance = 1e-10)
    }
  }
})

test_that("P2: voxelwise_cov with a group contrast is the signed group coefficient", {
  d <- p2_data()
  C <- matrix(stats::rnorm(length(d$Y)), nrow(d$Y))
  gd <- p2_block(d$Y, d$covars)
  fit <- fmri_ttest(gd, ~ 1 + group, engine = "classic", voxelwise_cov = C,
                    contrast = c(group = 2))
  ref <- fmri_ttest(gd, ~ 1 + group, engine = "classic", voxelwise_cov = C)
  expect_equal(unname(fit$beta_contrast), 2 * unname(ref$beta["group", ]), tolerance = 1e-12)
  expect_equal(unname(fit$se_contrast), 2 * unname(ref$se["group", ]), tolerance = 1e-12)
  expect_equal(unname(fit$t_contrast), unname(ref$t["group", ]), tolerance = 1e-12)

  # The voxelwise coefficient itself can be the contrast, by name or by a
  # full-length weight vector.
  fit_v <- fmri_ttest(gd, ~ 1 + group, engine = "classic", voxelwise_cov = C,
                      contrast = c(voxel_cov = 1))
  fit_u <- fmri_ttest(gd, ~ 1 + group, engine = "classic", voxelwise_cov = C,
                      contrast = c(0, 0, 1))
  expect_equal(unname(fit_v$beta_contrast), unname(ref$beta["voxel_cov", ]), tolerance = 1e-12)
  expect_equal(unname(fit_v$t_contrast), unname(ref$t["voxel_cov", ]), tolerance = 1e-12)
  expect_equal(fit_u$beta_contrast, fit_v$beta_contrast)
})

test_that("P2: paired and mu0 work on gds data (vs t.test)", {
  d <- p2_data()
  g <- p2_gds(d$Y, d$V, d$covars)
  fit_mu <- fmri_ttest(g, ~ 1, engine = "classic", mu0 = 0.3)
  fit_pair <- suppressMessages(fmri_ttest(g, ~ 1, paired = TRUE))
  for (j in seq_len(ncol(d$Y))) {
    ref <- stats::t.test(d$Y[, j], mu = 0.3)
    expect_equal(unname(fit_mu$t[1, j]), unname(ref$statistic), tolerance = 1e-10)
    expect_equal(unname(fit_pair$t[1, j]), unname(stats::t.test(d$Y[, j])$statistic),
                 tolerance = 1e-10)
  }
  expect_identical(fit_pair$engine, "classic")
})

test_that("P2: fmri_meta() accepts a gds built directly with fmrigds", {
  d <- p2_data()
  g <- p2_gds(d$Y, d$V, d$covars)
  g_wrapped <- g
  class(g_wrapped) <- c("group_data_gds", "group_data", class(g))
  fit <- fmri_meta(g, ~ 1 + group, verbose = FALSE)
  ref <- fmri_meta(g_wrapped, ~ 1 + group, verbose = FALSE)
  expect_equal(fit$coefficients, ref$coefficients, tolerance = 1e-12)
  expect_equal(fit$se, ref$se, tolerance = 1e-12)
  tt <- fmri_ttest(g, ~ 1 + group, engine = "meta")
  expect_equal(unname(tt$beta["group", ]), -unname(ref$coefficients[, "groupB"]), tolerance = 1e-12)
})

test_that("P2: engine = 'auto' records the engine that ran", {
  d <- p2_data()
  expect_identical(fmri_ttest(p2_block(d$Y, d$covars), ~ 1)$engine, "classic")
  expect_identical(fmri_ttest(p2_block(d$Y, d$covars, d$V), ~ 1)$engine, "meta")
})

test_that("P2: a rank-deficient design fails fast", {
  d <- p2_data()
  covars <- d$covars
  covars$age2 <- 2 * covars$age
  expect_error(fmri_ttest(p2_block(d$Y, covars), ~ 1 + age + age2, engine = "classic"),
               "rank deficient")
})

# LEGACY-GROUP-DELETE (next minor): tests the deprecated ols_t_cpp()/welch_t_cpp() wrappers
test_that("P2: ols_t_cpp() and welch_t_cpp() are deprecated and unchanged", {
  d <- p2_data()
  X <- stats::model.matrix(~ 1 + group, d$covars)
  expect_warning(ols <- ols_t_cpp(d$Y, X), "deprecated")
  ref <- fmri_ols_fit(d$Y, X)
  expect_equal(unname(ols$beta), unname(ref$beta), tolerance = 1e-12)
  expect_equal(unname(ols$se), unname(ref$se), tolerance = 1e-12)
  expect_equal(ols$df, ref$df)
  expect_warning(w <- welch_t_cpp(d$Y, as.integer(d$covars$group)), "deprecated")
  fit <- fmri_ttest(p2_block(d$Y, d$covars), ~ 1 + group, engine = "welch")
  expect_equal(as.numeric(w$t), unname(fit$t["group", ]), tolerance = 1e-10)
})
