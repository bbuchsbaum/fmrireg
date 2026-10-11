# Review 2026-10, final pass: defects found while merging the fix streams.

.fi_fixture <- function(nrun = 2L, n = 120L, V = 8L, seed = 11L, phi = 0.5) {
  set.seed(seed)
  ons <- seq(10, 2 * n - 30, by = 16)
  ev <- do.call(rbind, lapply(seq_len(nrun), function(r) {
    data.frame(onset = ons,
               condition = factor(rep(c("A", "B"), length.out = length(ons))),
               run = r)
  }))
  Y <- matrix(0, n * nrun, V)
  for (r in seq_len(nrun)) for (v in seq_len(V)) {
    Y[(r - 1L) * n + seq_len(n), v] <- 3 * as.numeric(arima.sim(list(ar = phi), n)) + 100
  }
  matrix_frame(Y, TR = 2, run_length = rep(n, nrun), event_table = ev)
}

.fi_fit <- function(dset, backend, noise = "ar1") {
  fmri_lm(
    onset ~ hrf(condition, contrasts = pair_contrast(~condition == "A",
                                                     ~condition == "B",
                                                     name = "AvB")),
    block = ~run, dataset = dset, durations = 0,
    control = fmri_lm_control(estimation = estimation_spec("joint"),
                              noise = noise_spec(noise)),
    compute = compute_spec(voxel_chunks = 1L, backend = backend)
  )
}

test_that("reference backend reports sigma, not sigma^2, for AR fits", {
  # `result$sigma` partially matched `sigma2`, so the reference backend scaled
  # standard errors by a variance. With one chunk it must equal the matrix
  # backend up to AR-estimation detail (the backends iterate phi slightly
  # differently: 0.511 vs 0.512 here), so compare at 1% -- the bug was 3x.
  dset <- .fi_fixture()
  ref <- .fi_fit(dset, "reference")
  mat <- .fi_fit(dset, "matrix")
  expect_equal(as.matrix(coef(ref)), as.matrix(coef(mat)), tolerance = 0.01)
  expect_equal(as.matrix(standard_error(ref)), as.matrix(standard_error(mat)),
               tolerance = 0.01)
  cr <- ref$result$contrasts$data[[1]]
  cm <- mat$result$contrasts$data[[1]]
  expect_equal(cr$se, cm$se, tolerance = 0.01)
})

test_that("soft projection: residual df is tr((I-H)P^2), so sigma^2 is unbiased", {
  # Without the correction runwise kept dfres = n - rank after shrinking the
  # nuisance subspace, biasing sigma^2 low (~8% here) and t-values high.
  set.seed(21)
  n <- 200L; V <- 3000L; k <- 20L
  ons <- seq(10, 2 * n - 30, by = 16)
  ev <- data.frame(onset = ons,
                   condition = factor(rep(c("A", "B"), length.out = length(ons))),
                   run = 1)
  N <- matrix(rnorm(n * k), n, k)
  Y <- matrix(rnorm(n * V), n, V)          # pure noise, unit variance
  dset <- matrix_frame(Y, TR = 2, run_length = n, event_table = ev)
  fit <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = dset,
                 control = fmri_lm_control(
                   estimation = estimation_spec("runwise_meta"),
                   projection = projection_spec("soft_subspace",
                                                nuisance_matrix = N,
                                                lambda = 50)))
  X <- as.matrix(design_matrix(fit$model))
  sp <- soft_projection(N, lambda = 50)
  df_exact <- .soft_projection_residual_df(sp, apply_soft_projection(sp, Y[, 1:2], X)$X)
  expect_lt(df_exact, n - qr(X)$rank - 5)       # the correction is material
  dfres <- unique(as.numeric(fit$result$df$inference))
  expect_equal(dfres, df_exact, tolerance = 1e-8)
  s2 <- as.numeric(fit$result$sigma)^2
  expect_equal(mean(s2), 1, tolerance = 0.02)   # unbiased under the null
})

test_that("fmri_ttest: a contrast on the group coefficient equals that coefficient", {
  # The contrast was computed in the requested sign convention and then had
  # its t and z flipped a second time, leaving t_contrast opposite to
  # beta_contrast whenever the requested sign differed from the engine's.
  set.seed(131)
  S <- 24; P <- 10
  group <- factor(rep(c("A", "B"), each = S / 2))
  Y <- matrix(rnorm(S * P, sd = 0.5), S, P)
  Y[group == "B", ] <- Y[group == "B", ] + 1
  V <- matrix(0.2, S, P)
  gd <- structure(list(blocks = list(list(Y = Y, V = V,
                                          covars = data.frame(group = group),
                                          feature = NULL))),
                  class = "group_data")
  for (engine in c("classic", "welch", "meta")) {
    for (sg in c("AminusB", "BminusA")) {
      fit <- fmri_ttest(gd, formula = ~ 1 + group, engine = engine, sign = sg,
                        contrast = c(group = 1))
      lab <- paste(engine, sg)
      expect_equal(as.numeric(fit$beta_contrast), as.numeric(fit$beta["group", ]),
                   tolerance = 1e-8, info = lab)
      expect_equal(as.numeric(fit$t_contrast), as.numeric(fit$t["group", ]),
                   tolerance = 1e-6, info = lab)
      expect_equal(sign(as.numeric(fit$z_contrast)),
                   sign(as.numeric(fit$beta_contrast)), info = lab)
      truth <- if (sg == "BminusA") 1 else -1
      expect_equal(sign(mean(fit$beta_contrast)), truth, info = lab)
    }
  }
})

test_that("fmri_ttest (GDS input): group contrast equals the group coefficient", {
  skip_if_not_installed("fmrigds")
  set.seed(132)
  S <- 20; R <- 4
  df <- expand.grid(subject = sprintf("s%02d", seq_len(S)), roi = sprintf("r%d", seq_len(R)),
                    stringsAsFactors = FALSE)
  df$group <- factor(ifelse(as.integer(sub("s", "", df$subject)) <= S / 2, "A", "B"))
  df$beta <- rnorm(nrow(df), sd = 0.5) + ifelse(df$group == "B", 1, 0)
  df$se <- 0.4
  df$contrast <- "c1"
  gd <- group_data(df, format = "csv", effect_cols = c(beta = "beta", se = "se"),
                   subject_col = "subject", roi_col = "roi", contrast_col = "contrast",
                   covariate_cols = "group")
  for (engine in c("classic", "meta")) {
    for (sg in c("AminusB", "BminusA")) {
      fit <- fmri_ttest(gd, formula = ~ 1 + group, engine = engine, sign = sg,
                        contrast = c(group = 1))
      lab <- paste("gds", engine, sg)
      expect_equal(as.numeric(fit$beta_contrast), as.numeric(fit$beta["group", ]),
                   tolerance = 1e-8, info = lab)
      expect_equal(sign(as.numeric(fit$z_contrast)), sign(as.numeric(fit$beta_contrast)),
                   info = lab)
      expect_equal(sign(mean(fit$beta_contrast)), if (sg == "BminusA") 1 else -1, info = lab)
    }
  }
})

test_that("runwise multi-run fits store per-voxel sigma, rss and resvar", {
  # pool_runwise_results stacked as.matrix(<V-vector>) columns, so colSums
  # summed over voxels *and* runs: sigma/rss/resvar became one scalar.
  dset <- .fi_fixture(nrun = 2L, V = 6L)
  fit <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = dset,
                 control = fmri_lm_control(estimation = estimation_spec("runwise_meta")))
  r <- fit$result
  Y <- as.matrix(.dset_data_matrix(dset))
  rows <- .dset_run_rows(dset)
  rss_ref <- numeric(ncol(Y)); rdf_ref <- 0
  for (k in seq_along(rows)) {
    tm <- term_matrices(fit$model, k)
    env <- list2env(tm); env[[".y"]] <- rep(0, nrow(tm[[1]]))
    Xk <- model.matrix(get_formula(fit$model), env)
    m <- lm.fit(Xk, Y[rows[[k]], , drop = FALSE])
    rss_ref <- rss_ref + colSums(m$residuals^2)
    rdf_ref <- rdf_ref + (nrow(Xk) - m$rank)
  }
  expect_length(as.numeric(r$sigma), ncol(Y))
  expect_equal(as.numeric(r$rss), rss_ref, tolerance = 1e-8)
  expect_equal(as.numeric(r$sigma), sqrt(rss_ref / rdf_ref), tolerance = 1e-8)
  expect_equal(as.numeric(r$resvar), rss_ref / rdf_ref, tolerance = 1e-8)
})

test_that("post-hoc contrasts on multi-run runwise fits: exact for OLS, refused when runs are whitened", {
  dset <- .fi_fixture(nrun = 2L, V = 6L)
  pc <- pair_contrast(~condition == "A", ~condition == "B", name = "AvB")
  fit_rw <- function(noise) {
    fmri_lm(onset ~ hrf(condition, contrasts = pc), block = ~run, dataset = dset,
            control = fmri_lm_control(estimation = estimation_spec("runwise_meta"),
                                      noise = noise_spec(noise)))
  }
  ols <- fit_rw("iid")
  ph <- fit_contrasts(ols, list(AvB = pc))[[1]]
  infit <- ols$result$contrasts$data[[1]]
  expect_equal(ph$estimate, infit$estimate, tolerance = 1e-8)
  # Pooled-run OLS: global (X'X)^-1 with the pooled per-voxel sigma; the in-fit
  # value pools run-level contrasts by inverse variance. They agree closely
  # (the scalar-sigma bug made this ratio ~3).
  expect_equal(ph$se, infit$se, tolerance = 0.03)

  ar <- fit_rw("ar1")
  expect_error(fit_contrasts(ar, list(AvB = pc)), "multi-run runwise")
  # Joint AR fits remain supported and exact.
  joint <- .fi_fit(dset, "matrix")
  ph_j <- fit_contrasts(joint, list(AvB = pc))[[1]]
  expect_equal(ph_j$se, joint$result$contrasts$data[[1]]$se, tolerance = 1e-8)
})

test_that("meta_betas pools runs by inverse-variance fixed effects (closed form)", {
  # Replaces a coverage test that accepted any return value, including an error.
  set.seed(82)
  V <- 4L; p <- 3L
  mk <- function(b, s) tibble::tibble(
    type = "beta", name = "parameter_estimates",
    data = list(tibble::tibble(estimate = list(b), se = list(s),
                               stat = list(b / s), prob = list(b * 0),
                               sigma = list(rep(1, V)))))
  e1 <- matrix(rnorm(V * p), V, p); s1 <- matrix(runif(V * p, 0.2, 1), V, p)
  e2 <- matrix(rnorm(V * p), V, p); s2 <- matrix(runif(V * p, 0.2, 1), V, p)
  mb <- fmrireg:::meta_betas(list(mk(e1, s1), mk(e2, s2)), colind = 1:p)
  w1 <- 1 / s1^2; w2 <- 1 / s2^2
  est <- (w1 * e1 + w2 * e2) / (w1 + w2)
  se <- sqrt(1 / (w1 + w2))
  got <- mb$data[[1]]
  expect_equal(unname(got$estimate[[1]]), unname(est), tolerance = 1e-12)
  expect_equal(unname(got$se[[1]]), unname(se), tolerance = 1e-12)
  expect_equal(unname(got$stat[[1]]), unname(est / se), tolerance = 1e-12)
  expect_equal(unname(got$prob[[1]]),
               unname(2 * pnorm(abs(est / se), lower.tail = FALSE)), tolerance = 1e-10)
})

test_that("fit_lm_contrasts forms (X'X)^-1 once per model, not once per contrast", {
  set.seed(5)
  n <- 80; X <- cbind(1, matrix(rnorm(n * 3), n, 3)); Y <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- c("(Intercept)", "a", "b", "c")
  fit <- lm(Y ~ X - 1)
  mk_t <- function(w, nm) structure(list(weights = matrix(w, ncol = 1), name = nm),
                                    class = c("contrast", "list"), colind = 2:4)
  mk_f <- function(w, nm) structure(list(weights = w, name = nm),
                                    class = c("Fcontrast", "list"), colind = 2:4)
  cons <- list(mk_t(c(1, -1, 0), "ab"), mk_t(c(0, 1, -1), "bc"),
               mk_f(cbind(c(1, -1, 0), c(0, 1, -1)), "F"))
  calls <- 0L
  real <- fmrireg:::.lm_cov_unscaled
  testthat::local_mocked_bindings(
    .lm_cov_unscaled = function(...) { calls <<- calls + 1L; real(...) },
    .package = "fmrireg")
  res <- fmrireg:::fit_lm_contrasts(fit, cons, NULL, colnames(X))
  expect_equal(calls, 1L)
  # Values are the textbook ones.
  b <- coef(fit); V1 <- vcov(lm(Y[, 1] ~ X - 1))
  w <- c(0, 1, -1, 0)
  expect_equal(as.numeric(res$contrasts$ab$data[[1]]$estimate[1]), sum(w * b[, 1]),
               tolerance = 1e-10)
  expect_equal(as.numeric(res$contrasts$ab$data[[1]]$se[1]),
               sqrt(drop(t(w) %*% V1 %*% w)), tolerance = 1e-10)
  full <- lm(Y[, 1] ~ X - 1); red <- lm(Y[, 1] ~ X[, 1] + I(X[, 2] + X[, 3] + X[, 4]) - 1)
  expect_equal(as.numeric(res$contrasts$F$data[[1]]$stat[1]),
               anova(red, full)$F[2], tolerance = 1e-8)
})
