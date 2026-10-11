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
