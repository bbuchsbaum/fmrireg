test_that("runwise_lm handles rank deficiency gracefully", {
  options(mc.cores = 1)

  set.seed(123)

  sframe <- fmrihrf::sampling_frame(blocklens = c(10, 10), TR = 1)

  # simple event table with one event per run
  etab <- data.frame(onset = c(1, 1),
                     condition = factor(c("A", "A")),
                     run = c(1, 2))

  Y <- matrix(rnorm(sum(fmrihrf::blocklens(sframe)) * 2),
              sum(fmrihrf::blocklens(sframe)), 2)

  dset <- matrix_dataset(Y, TR = 1,
                         run_length = fmrihrf::blocklens(sframe),
                         event_table = etab)

  espec <- event_model(onset ~ hrf(condition), data = etab,
                       block = ~run, sampling_frame = sframe)

  # nuisance columns duplicate the runwise intercept
  nlist <- list(matrix(1, 10, 1), matrix(1, 10, 1))
  bspec <- baseline_model(basis = "poly", degree = 1,
                          sframe = sframe, intercept = "runwise",
                          nuisance_list = nlist)

  fmod <- fmri_model(espec, bspec, dset)

  expect_warning(
    fit <- fmri_lm(onset ~ hrf(condition), block = ~run,
                   dataset = dset, baseline_model = bspec,
                   strategy = "runwise", use_fast_path = TRUE),
    regexp = "rank deficient"
  )

  # cov.unscaled should keep column names
  expect_equal(colnames(fit$result$cov.unscaled), colnames(design_matrix(fmod)))
  expect_false(isTRUE(attr(fit$result$cov.unscaled, "is_full_rank")))
  expect_equal(attr(fit$result$cov.unscaled, "rank"), qr(as.matrix(design_matrix(fmod)))$rank)
  expect_true(length(attr(fit$result$cov.unscaled, "aliased")) > 0)
})

test_that("rank-deficient matrix path preserves columns but marks non-estimable coefficients", {
  set.seed(2026)
  n <- 40
  X <- cbind(intercept = 1, signal = rnorm(n), duplicate = 0)
  X[, "duplicate"] <- X[, "signal"]
  Y <- matrix(rnorm(n * 3), n, 3)

  expect_warning(
    proj <- fmrireg:::.fast_preproject(X),
    regexp = "rank deficient.*Aliased columns: duplicate"
  )
  expect_equal(proj$rank, 2)
  expect_equal(proj$aliased, 3)

  fit <- fmrireg:::solve_glm_core(
    fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  )
  expect_equal(attr(fit$betas, "aliased"), 3)
  expect_true(all(is.finite(fit$betas)))

  bstats <- fmrireg:::beta_stats_matrix(
    fit$betas,
    proj$XtXinv,
    sqrt(fit$sigma2),
    fit$dfres,
    colnames(X)
  )
  estimates <- bstats$data[[1]]$estimate[[1]]
  ses <- bstats$data[[1]]$se[[1]]

  expect_equal(colnames(estimates), colnames(X))
  expect_true(all(is.na(estimates[, "duplicate"])))
  expect_true(all(is.na(ses[, "duplicate"])))
  expect_true(all(is.finite(estimates[, "signal"])))
})

test_that("lm.fit stats avoid chol2inv failure and reject aliased contrasts", {
  set.seed(2027)
  n <- 35
  X <- cbind(intercept = 1, signal = rnorm(n), duplicate = 0)
  X[, "duplicate"] <- X[, "signal"]
  Y <- matrix(rnorm(n * 2), n, 2)
  fit <- lm.fit(X, Y)

  expect_warning(
    bstats <- fmrireg:::beta_stats(fit, colnames(X)),
    regexp = "rank deficient.*Aliased columns: duplicate"
  )
  estimates <- bstats$data[[1]]$estimate[[1]]
  expect_true(all(is.na(estimates[, "duplicate"])))
  expect_true(all(is.finite(estimates[, "signal"])))

  expect_warning(
    tcon <- fmrireg:::fit_contrasts.default(fit, c(0, 0, 1), colind = seq_len(3)),
    regexp = "non-estimable coefficients"
  )
  expect_true(all(is.na(tcon$estimate)))
  expect_true(all(is.na(tcon$se)))
})
