# Variance helpers, update_inference, finalize sandwich, tidy_fitted_hrf branches.

test_that(".fmri_lm_segments / hac_weight / score / residual / update_inference", {
  set.seed(11)
  n <- 24
  p <- 3
  V <- 2
  X <- cbind(1, rnorm(n), rnorm(n))
  runs <- rep(1:2, each = 12)

  segs <- fmrireg:::.fmri_lm_segments(n, runs = runs, censor = c(3L, 15L))
  expect_true(length(segs) >= 2)
  expect_false(any(unlist(segs) %in% c(3L, 15L)))

  segs_l <- fmrireg:::.fmri_lm_segments(
    n, runs = runs, censor = c(rep(FALSE, 5), TRUE, rep(FALSE, n - 6))
  )
  expect_false(6L %in% unlist(segs_l))
  expect_error(
    fmrireg:::.fmri_lm_segments(n, censor = rep(TRUE, 3)),
    "Logical censor"
  )

  expect_equal(fmrireg:::.fmri_lm_hac_weight(0L, 4L, "tukey"), 1)
  expect_equal(fmrireg:::.fmri_lm_hac_weight(1L, 4L, "none"), 1)
  w_t <- fmrireg:::.fmri_lm_hac_weight(1L, 4L, "tukey")
  expect_true(w_t > 0 && w_t < 1)
  w_p <- fmrireg:::.fmri_lm_hac_weight(1L, 4L, "parzen")
  expect_true(w_p > 0 && w_p <= 1)
  w_p2 <- fmrireg:::.fmri_lm_hac_weight(3L, 4L, "parzen")
  expect_true(w_p2 >= 0 && w_p2 <= 1)

  residual <- rnorm(n)
  segs2 <- fmrireg:::.fmri_lm_segments(n, runs = runs)
  Vhat <- fmrireg:::.fmri_lm_score_covariance(
    X, residual, segs2, max_lag = 2L, taper = "tukey", debias = TRUE
  )
  expect_equal(dim(Vhat), c(p, p))
  expect_true(all(eigen(Vhat, symmetric = TRUE)$values >= -1e-8))

  Rcov <- fmrireg:::.fmri_lm_residual_covariance(
    residual, segs2, max_lag = 2L, taper = "parzen"
  )
  expect_equal(dim(Rcov), c(n, n))

  df_h <- fmrireg:::.fmri_lm_hac_df(X, residual, segs2, max_lag = 2L, taper = "tukey")
  expect_true(is.finite(df_h) && df_h >= 1)

  # update_inference with t- and F-contrasts
  estimate <- matrix(rnorm(V * p), V, p)
  payload <- tibble::tibble(
    estimate = list(estimate),
    se = list(matrix(1, V, p)),
    stat = list(matrix(0, V, p)),
    prob = list(matrix(1, V, p))
  )
  result <- list(
    betas = tibble::tibble(
      type = "beta",
      name = "parameter_estimates",
      data = list(payload)
    ),
    contrasts = tibble::tibble(
      type = c("contrast", "Fcontrast"),
      name = c("c1", "f1"),
      conmat = list(matrix(c(0, 1, -1), ncol = 1), diag(p)[2:3, , drop = FALSE]),
      colind = list(1:3, 1:3),
      data = list(
        list(estimate = numeric(V), se = numeric(V), stat = numeric(V), prob = numeric(V)),
        list(estimate = numeric(V), se = numeric(V), stat = numeric(V), prob = numeric(V))
      )
    )
  )
  covs <- lapply(seq_len(V), function(i) diag(p) * 0.25)
  updated <- fmrireg:::.fmri_lm_update_inference(result, covs, df = c(20, 20))
  expect_true(all(is.finite(updated$betas$data[[1]]$se[[1]])))
  expect_true(all(is.finite(updated$contrasts$data[[1]]$stat)))
  expect_true(all(is.finite(updated$contrasts$data[[2]]$prob) |
                    is.na(updated$contrasts$data[[2]]$prob)))
})

test_that(".fmri_lm_finalize_result sandwich path with inference_context", {
  set.seed(12)
  n <- 40
  p <- 2
  V <- 2
  X <- cbind(1, rnorm(n))
  Y <- X %*% matrix(c(1, 0.5), 2, V) + matrix(rnorm(n * V), n, V)
  beta <- t(solve(crossprod(X), crossprod(X, Y)))
  resid <- Y - X %*% t(beta)

  payload <- tibble::tibble(
    estimate = list(beta),
    se = list(matrix(1, V, p)),
    stat = list(matrix(0, V, p)),
    prob = list(matrix(1, V, p)),
    sigma = list(rep(1, V))
  )
  result <- list(
    betas = tibble::tibble(
      type = "beta",
      name = "parameter_estimates",
      data = list(payload),
      df.residual = n - p
    ),
    rdf = n - p,
    sigma = apply(resid, 2, sd),
    cov.unscaled = solve(crossprod(X)),
    inference_context = list(
      X = X,
      residuals = resid,
      runs = list(1:20, 21:40),
      censor = NULL
    )
  )
  cfg <- fmri_lm_control(variance = variance_spec("sandwich"))
  finalized <- fmrireg:::.fmri_lm_finalize_result(result, cfg = cfg, engine = "test")
  expect_equal(finalized$schema_version, 2L)
  expect_equal(finalized$variance_model$covariance_scope, "voxel")
  expect_s3_class(finalized$variance_model, "fmri_lm_variance_model")
  expect_output(print(finalized$variance_model), "fmri_lm_variance_model")

  expect_error(variance_model(list()), "fmri_lm")
})

test_that("tidy_fitted_hrf covers average/exact/regex and error branches", {
  skip_if_not_installed("fmridataset")

  set.seed(101)
  n <- 50L
  Y <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  ev <- data.frame(
    onsets = c(6, 16, 26, 36),
    run = 1,
    condition = factor(c("A", "B", "A", "B"))
  )
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = ev)
  fit <- fmri_lm(
    onsets ~ hrf(condition), block = ~ run, dataset = dset,
    strategy = "runwise", nchunks = 1
  )

  expect_error(tidy_fitted_hrf(list()), "fmri_lm")

  avg <- tidy_fitted_hrf(fit, sample_at = 0:8, average_voxels = TRUE)
  expect_true(all(is.na(avg$voxel)))
  expect_gt(nrow(avg), 0L)

  exact <- tidy_fitted_hrf(
    fit, sample_at = 0:5, term = names(fitted_hrf(fit, sample_at = 0:5))[1],
    term_match = "exact", voxel = 1L
  )
  expect_gt(nrow(exact), 0L)

  rx <- tidy_fitted_hrf(
    fit, sample_at = 0:5, term = "cond", term_match = "regex", voxel = 2L
  )
  expect_gt(nrow(rx), 0L)

  expect_error(
    tidy_fitted_hrf(fit, term = "no-such-term-zzz", term_match = "exact"),
    "No fitted HRF terms"
  )
  expect_error(tidy_fitted_hrf(fit, voxel = 0L), ">= 1")
  expect_error(tidy_fitted_hrf(fit, voxel = c(1L, 2L)), "single numeric")
  expect_error(tidy_fitted_hrf(fit, voxel = 99L), "exceeds prediction")
})

test_that(".fmri_normalize_mask_space drops 4D and checks dims", {
  skip_if_not_installed("neuroim2")
  sp4 <- neuroim2::NeuroSpace(c(4, 4, 2, 10))
  sp3 <- fmrireg:::.fmri_normalize_mask_space(sp4, NULL, "test")
  expect_equal(length(dim(sp3)), 3L)

  expect_error(
    fmrireg:::.fmri_normalize_mask_space(sp3, c(3L, 3L, 3L), "test"),
    "do not match"
  )
})
