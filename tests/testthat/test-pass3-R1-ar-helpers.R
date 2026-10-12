# Pass 3, R1-4: dead and inconsistent AR helpers.

.p3_ar_data <- function(n = 160L, V = 6L, phi = 0.5, seed = 7L) {
  set.seed(seed)
  X <- cbind(1, rnorm(n), sin(seq_len(n) / 7))
  E <- vapply(seq_len(V), function(v) as.numeric(arima.sim(list(ar = phi), n)),
              numeric(n))
  list(X = X, Y = X %*% matrix(c(1, 0.5, 0), 3, V) + E)
}

test_that("fast_rlm_run (a MAD-scale robust fitter with no caller) is gone", {
  ns <- asNamespace("fmrireg")
  expect_false(exists("fast_rlm_run", envir = ns, inherits = FALSE))
})

test_that("AR effective df is not reduced a second time after whitening", {
  # .compute_ar_effective_df_compat shrank n by 1 + 2 sum rho_k although the
  # data had already been whitened; the integrated solver then reported that
  # df as `effective_df` and compute_contrast() used it for p-values.
  ns <- asNamespace("fmrireg")
  expect_false(exists(".compute_ar_effective_df_compat", envir = ns, inherits = FALSE))
  expect_false(exists("compute_ar_effective_df", envir = ns, inherits = FALSE))

  d <- .p3_ar_data()
  n <- nrow(d$X); p <- ncol(d$X)
  cfg <- fmri_lm_control(estimation = estimation_spec("joint"),
                         noise = noise_spec("ar1"))
  fit <- fmrireg:::solve_integrated_glm(d$X, d$Y, cfg)
  expect_gt(abs(unlist(fit$ar_coef)[1]), 0.3)
  expect_equal(fit$effective_df, n - p)
  expect_equal(fit$effective_df,
               fmrireg:::calculate_effective_df(n, p, method = "residual"))
  con <- fmrireg:::compute_contrast(fit, c(0, 1, 0))
  expect_equal(con$df, n - p)
  expect_equal(as.numeric(con$pvalue),
               as.numeric(2 * pt(-abs(con$tstat), n - p)))

  # Whitened then robustly weighted: the one robust effective df, sum(w) - rank.
  cfg_r <- fmri_lm_control(estimation = estimation_spec("joint"),
                           noise = noise_spec("ar1"),
                           robust = robust_spec("huber"))
  fit_r <- fmrireg:::solve_integrated_glm(d$X, d$Y, cfg_r)
  w <- fit_r$robust_weights
  expect_equal(fit_r$effective_df, sum(w) - p)
  expect_equal(fit_r$effective_df, fmrireg:::.robust_effective_df(w, p))
})

test_that("exact_first defaults to FALSE everywhere (documented default)", {
  d <- .p3_ar_data(n = 40L, V = 2L)
  ctx <- fmrireg:::glm_context(X = d$X, Y = d$Y,
                               proj = fmrireg:::.fast_preproject(d$X))
  # No exact_first given: the first observation is not rescaled ...
  w_default <- fmrireg:::whiten_glm_context(ctx, list(struct = "ar1", phi = 0.6))
  w_false <- fmrireg:::whiten_glm_context(
    ctx, list(struct = "ar1", phi = 0.6, exact_first = FALSE))
  w_true <- fmrireg:::whiten_glm_context(
    ctx, list(struct = "ar1", phi = 0.6, exact_first = TRUE))
  expect_equal(w_default$X, w_false$X)
  expect_equal(w_default$X[1, ], d$X[1, ])
  expect_equal(w_true$X[1, ], d$X[1, ] * sqrt(1 - 0.6^2))
  # ... and the fmriAR translation agrees.
  expect_identical(
    fmrireg:::.fmrireg_to_fmriAR_config(list(struct = "ar1"))$exact_first, "none")
  expect_identical(
    fmrireg:::.fmrireg_to_fmriAR_config(list(struct = "ar1", exact_first = TRUE))$exact_first,
    "ar1")
  expect_identical(fmri_lm_control()$ar$exact_first, FALSE)
})

test_that("one AR-order helper; unknown structures are an error, not order 0/NULL", {
  ao <- fmrireg:::.ar_order
  expect_identical(ao(list(struct = "iid")), 0L)
  expect_identical(ao(list(struct = "ar1")), 1L)
  expect_identical(ao(list(struct = "ar2")), 2L)
  expect_identical(ao(list(struct = "arp", p = 5L)), 5L)
  expect_identical(ao(fmri_lm_control(noise = noise_spec("ar2"))), 2L)
  expect_error(ao(list(struct = "ar_bogus")), "Unknown AR structure")
  expect_error(ao(list(struct = "arp", p = 0L)), "positive integer")

  # Former call sites mapped unknown structures to 0 (no whitening) or NULL.
  cfg_bad <- fmri_lm_control(noise = noise_spec("ar1"))
  cfg_bad$ar$struct <- "ar_bogus"
  expect_error(fmrireg:::.rrr_ar_order_from_cfg(cfg_bad), "Unknown AR structure")
  expect_error(fmrireg:::.get_ar_order(cfg = cfg_bad$ar), "Unknown AR structure")
  expect_error(fmrireg:::.fmrireg_to_fmriAR_config(cfg_bad$ar), "Unknown AR structure")

  # No copy of the old switch block survives in the fitting code.
  src <- list.files(test_path("..", "..", "R"), pattern = "[.]R$", full.names = TRUE)
  skip_if(length(src) == 0L, "package sources not available")
  hits <- unlist(lapply(src, function(f) {
    grep("switch\\((cfg\\$ar\\$struct|as.character\\(ar_struct\\)|ar_opts\\$struct),",
         readLines(f, warn = FALSE), value = TRUE)
  }))
  expect_length(hits, 0L)
})
