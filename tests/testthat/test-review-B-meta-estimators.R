# Review 2026-10, stream B: meta-analysis estimators against metafor (B5),
# contrasts on fmri_meta fits (B6), dynamic OpenMP schedule (P-B2).

meta_sim <- function(S = 15L, P = 6L, seed = 11) {
  set.seed(seed)
  covars <- data.frame(
    x = stats::rnorm(S),
    g = factor(rep(c("a", "b", "c"), length.out = S))
  )
  V <- matrix(stats::runif(S * P, 0.02, 0.4), S, P)
  tau <- seq(0, 0.6, length.out = P)  # includes a homogeneous voxel
  Y <- matrix(NA_real_, S, P)
  for (j in seq_len(P)) {
    Y[, j] <- 0.3 + 0.5 * covars$x + stats::rnorm(S, sd = sqrt(V[, j] + tau[j]^2))
  }
  list(Y = Y, V = V, covars = covars)
}

# Subject-level effects (subjects x voxels) as in-memory gds group data.
gds_gd <- function(beta, se, covariates = NULL) {
  S <- nrow(beta)
  P <- ncol(beta)
  subj <- paste0("s", seq_len(S))
  B <- t(beta)
  dimnames(B) <- list(paste0("v", seq_len(P)), subj)
  V <- t(se^2)
  dimnames(V) <- dimnames(B)
  cd <- if (is.null(covariates)) data.frame(row.names = subj) else covariates
  rownames(cd) <- subj
  g <- fmrigds::as_gds(list(beta = B, var = V), col_data = cd)
  class(g) <- c("group_data_gds", "group_data", class(g))
  g
}

rma_ref <- function(y, v, covars, mods, method) {
  ctrl <- list(tol = 1e-13, threshold = 1e-13, maxiter = 10000)
  if (is.null(mods)) {
    metafor::rma(yi = y, vi = v, method = method, control = ctrl)
  } else {
    metafor::rma(yi = y, vi = v, mods = mods, data = covars, method = method, control = ctrl)
  }
}

test_that("B5: FE, DL and PM match metafor::rma (intercept-only and meta-regression)", {
  skip_if_not_installed("metafor")
  sim <- meta_sim()
  designs <- list(
    intercept = list(X = matrix(1, nrow(sim$Y), 1), mods = NULL),
    regression = list(X = stats::model.matrix(~ x + g, sim$covars), mods = ~ x + g)
  )
  methods <- c(fe = "FE", dl = "DL", pm = "PM")
  for (dn in names(designs)) {
    d <- designs[[dn]]
    for (m in names(methods)) {
      fit <- fmri_meta_fit(sim$Y, sim$V, d$X, method = m, n_threads = 1)
      for (j in seq_len(ncol(sim$Y))) {
        ref <- rma_ref(sim$Y[, j], sim$V[, j], sim$covars, d$mods, methods[[m]])
        info <- paste(dn, m, "voxel", j)
        expect_equal(as.numeric(fit$beta[, j]), as.numeric(ref$beta), tolerance = 1e-6, info = info)
        expect_equal(as.numeric(fit$se[, j]), as.numeric(ref$se), tolerance = 1e-6, info = info)
        expect_equal(fit$tau2[j], ref$tau2, tolerance = 1e-6, info = info)
        expect_equal(fit$Q_fe[j], ref$QE, tolerance = 1e-6, info = info)
        expect_equal(fit$df[j], ref$k - ref$p, info = info)
      }
    }
  }
})

test_that("B5: method = 'reml' is labelled honestly (Paule-Mandel with a message)", {
  sim <- meta_sim()
  X <- matrix(1, nrow(sim$Y), 1)
  expect_message(
    fit_reml <- fmri_meta_fit(sim$Y, sim$V, X, method = "reml", n_threads = 1),
    "Paule-Mandel"
  )
  fit_pm <- fmri_meta_fit(sim$Y, sim$V, X, method = "pm", n_threads = 1)
  expect_identical(fit_reml$method, "pm")
  expect_equal(fit_reml$tau2, fit_pm$tau2)

  skip_if_not_installed("fmrigds")
  gd <- gds_gd(sim$Y, sqrt(sim$V))
  expect_message(
    meta <- fmri_meta(gd, method = "reml", verbose = FALSE, n_threads = 1),
    "Paule-Mandel"
  )
  expect_identical(meta$method, "pm")
})

test_that("B6: contrast() on an fmri_meta fit uses the coefficient covariance (vs metafor)", {
  skip_if_not_installed("metafor")
  sim <- meta_sim(seed = 12)
  skip_if_not_installed("fmrigds")
  gd <- gds_gd(sim$Y, sqrt(sim$V), covariates = sim$covars)
  for (m in c(fe = "FE", dl = "DL", pm = "PM")) {
    meth <- tolower(m)
    fit <- fmri_meta(gd, formula = ~ x + g, method = meth, verbose = FALSE, n_threads = 1)
    L <- c(`(Intercept)` = 1, x = 2, gb = -1, gc = 0)
    con <- contrast(fit, L)
    for (j in seq_len(ncol(sim$Y))) {
      ref <- rma_ref(sim$Y[, j], sim$V[, j], sim$covars, ~ x + g, m)
      an <- stats::anova(ref, X = matrix(L, nrow = 1))
      expect_equal(con$estimate[j], as.numeric(an$Xb), tolerance = 1e-6, info = paste(meth, j))
      expect_equal(con$se[j], as.numeric(an$se), tolerance = 1e-6, info = paste(meth, j))
      expect_equal(con$z[j], as.numeric(an$zval), tolerance = 1e-6, info = paste(meth, j))
    }
  }
})

test_that("B6: contrast() without stored covariance warns for multi-coefficient contrasts", {
  obj <- structure(
    list(coefficients = matrix(c(0.3, 0.1), nrow = 1, dimnames = list(NULL, c("A", "B"))),
         se = matrix(c(0.05, 0.06), nrow = 1), robust = "none"),
    class = "fmri_meta"
  )
  expect_warning(contrast(obj, c(1, -1)), "covariance")
  # A single-coefficient contrast is exact from the standard errors alone
  expect_silent(single <- contrast(obj, c(0, 2)))
  expect_equal(single$se, 0.12)
})

test_that("P-B2: results do not depend on the number of threads", {
  sim <- meta_sim(P = 40L, seed = 13)
  X <- stats::model.matrix(~ x, sim$covars)
  one <- fmri_meta_fit_cov(sim$Y, sim$V, X, method = "pm", n_threads = 1)
  many <- fmri_meta_fit_cov(sim$Y, sim$V, X, method = "pm", n_threads = 4)
  expect_equal(many$beta, one$beta)
  expect_equal(many$cov_tri, one$cov_tri)
  expect_equal(many$tau2, one$tau2)
})
