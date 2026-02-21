test_that("process_run_standard re-estimates AR from raw-domain residuals", {
  n <- 6L
  v <- 2L
  Y <- matrix(seq_len(n * v), nrow = n, ncol = v)
  run_chunk <- list(data = Y, chunk_num = 1L)

  cfg <- list(
    ar = list(struct = "ar1", iter_gls = 2L, exact_first = TRUE, censor = NULL),
    volume_weights = list(enabled = FALSE),
    soft_subspace = list(enabled = FALSE)
  )

  captured <- list()

  with_mocked_bindings(
    {
      out <- fmrireg:::process_run_standard(
        run_chunk = run_chunk,
        model = structure(list(), class = "fmri_model"),
        cfg = cfg,
        phi_fixed = NULL,
        simple_conlist_weights = list(),
        fconlist_weights = list(),
        dataset = NULL
      )
      expect_true(is.list(out))
    },
    term_matrices = function(...) {
      onset <- matrix(seq_len(n), ncol = 1L)
      colnames(onset) <- "onset"
      list(onset = onset)
    },
    get_formula = function(...) .y ~ 1,
    preprocess_run_data = function(X, Y, cfg, dataset = NULL, run_num = NULL) {
      list(X = X, Y = Y, preprocess_info = list())
    },
    resolve_censor = function(...) NULL,
    .fast_preproject = function(X) {
      list(
        XtXinv = diag(ncol(X)),
        dfres = max(1L, nrow(X) - ncol(X)),
        Pinv = MASS::ginv(X),
        rank = ncol(X),
        is_full_rank = TRUE,
        qr = qr(X)
      )
    },
    glm_context = function(X, Y, proj, phi_hat = NULL) {
      list(X = X, Y = Y, proj = proj, phi_hat = phi_hat)
    },
    ar_whiten_transform = function(X, Y, phi_hat, exact_first, censor = NULL) {
      list(X = 2 * X, Y = 3 * Y)
    },
    solve_glm_core = function(ctx, return_fitted = FALSE) {
      p <- ncol(ctx$X)
      vv <- ncol(ctx$Y)
      betas <- matrix(rep(c(0.5, -0.25), length.out = p * vv), nrow = p, ncol = vv)
      fitted <- ctx$X %*% betas
      rss <- colSums((ctx$Y - fitted)^2)
      out <- list(
        betas = betas,
        sigma2 = rss / max(1L, ctx$proj$dfres),
        rss = rss
      )
      if (isTRUE(return_fitted)) {
        out$fitted <- fitted
      }
      out
    },
    .estimate_ar_parameters_routed = function(residuals_vec, ar_order, run_indices = NULL, censor = NULL) {
      captured[[length(captured) + 1L]] <<- as.numeric(residuals_vec)
      0.4
    },
    .package = "fmrireg"
  )

  expect_gte(length(captured), 2L)

  X_raw <- model.matrix(.y ~ 1, data.frame(.y = rep(0, n)))
  betas <- matrix(rep(c(0.5, -0.25), length.out = ncol(X_raw) * v), nrow = ncol(X_raw), ncol = v)

  expected_raw <- rowMeans(Y - X_raw %*% betas)
  expected_white <- rowMeans((3 * Y) - (2 * X_raw) %*% betas)

  expect_equal(unname(captured[[2]]), unname(expected_raw))
  expect_gt(max(abs(captured[[2]] - expected_white)), 1e-6)
})

test_that(".run_lowrank_engine re-estimates AR from raw-domain residuals", {
  n <- 8L
  v <- 3L
  X <- cbind("(Intercept)" = 1, x = seq_len(n))
  set.seed(123)
  Z <- matrix(rnorm(n * v), nrow = n, ncol = v)

  etab <- data.frame(onset = 1, run = 1)
  dataset <- fmridataset::matrix_dataset(Z, TR = 1, run_length = n, event_table = etab)
  fm_dummy <- structure(list(), class = "fmri_model")

  lowrank <- list(time_sketch = list(method = "gaussian", m = n))
  cfg <- list(
    ar = list(struct = "ar1", iter_gls = 2L, exact_first = TRUE, by_cluster = FALSE),
    robust = list()
  )

  captured <- list()

  with_mocked_bindings(
    {
      out <- fmrireg:::.run_lowrank_engine(
        fm = fm_dummy,
        dataset = dataset,
        lowrank = lowrank,
        cfg = cfg
      )
      expect_true(is.list(out))
    },
    design_matrix = function(fm) X,
    term_matrices = function(fm, ...) {
      tm <- list()
      attr(tm, "event_term_indices") <- integer(0)
      attr(tm, "baseline_term_indices") <- integer(0)
      tm
    },
    make_time_sketch = function(Tlen, sk) diag(Tlen),
    ar_whiten_transform = function(X, Y, phi_hat, exact_first, ...) {
      list(X = 2 * X, Y = 3 * Y)
    },
    .estimate_ar_parameters_routed = function(residuals_vec, ar_order, run_indices = NULL, censor = NULL) {
      captured[[length(captured) + 1L]] <<- as.numeric(residuals_vec)
      0.4
    },
    beta_stats_matrix = function(Betas, XtXinv, sigma, dfres, varnames, ...) Betas,
    .package = "fmrireg"
  )

  expect_gte(length(captured), 2L)

  Xw <- 2 * X
  Zw <- 3 * Z
  beta_iter <- solve(crossprod(Xw), crossprod(Xw, Zw))

  expected_raw <- rowMeans(Z - X %*% beta_iter)
  expected_white <- rowMeans(Zw - Xw %*% beta_iter)

  expect_equal(unname(captured[[2]]), unname(expected_raw))
  expect_gt(max(abs(captured[[2]] - expected_white)), 1e-6)
})
