# Review 2026-10, stream D: robust fitting (D2), AR re-estimation (D1) and
# runwise volume weights (D3).
#
# Oracle: the robust fit is a weighted least-squares fit with the returned
# observation weights, so `stats::lm(weights = w)` must reproduce its
# coefficients, standard errors, residual scale and residual df.

.rb_sim <- function(nrun = 1L, n = 120L, V = 5L, phi = 0, seed = 1L,
                    spikes = c(7L, 40L, 41L, 95L), spike_size = 25) {
  set.seed(seed)
  ons <- seq(6, 2 * n - 30, by = 14)
  ev <- do.call(rbind, lapply(seq_len(nrun), function(r) {
    data.frame(onset = ons,
               condition = factor(rep(c("A", "B"), length.out = length(ons))),
               run = r)
  }))
  noise <- matrix(0, n * nrun, V)
  for (r in seq_len(nrun)) {
    for (v in seq_len(V)) {
      e <- if (phi == 0) rnorm(n) else as.numeric(arima.sim(list(ar = phi), n))
      noise[(r - 1L) * n + seq_len(n), v] <- e + 50 + 5 * r
    }
  }
  fr0 <- matrix_frame(noise, TR = 2, run_length = rep(n, nrun), event_table = ev)
  mod <- create_fmri_model(onset ~ hrf(condition), block = ~run,
                           dataset = fr0, durations = 0)
  X <- as.matrix(design_matrix(mod))
  ec <- unlist(attr(mod$event_model$design_matrix, "col_indices"))
  scale <- 1 / max(abs(X[, ec]))
  Y <- noise + X[, ec] %*% (scale * matrix(c(3, -2), 2, V))
  for (r in seq_len(nrun)) {
    rows <- (r - 1L) * n + spikes
    Y[rows, ] <- Y[rows, ] + spike_size
  }
  list(dset = matrix_frame(Y, TR = 2, run_length = rep(n, nrun), event_table = ev),
       Y = Y, X = X, ec = ec, n = n, nrun = nrun)
}

.rb_fit <- function(s, scope = "joint", noise = "iid", robust = "bisquare",
                    chunks = 1L, weights = NULL, ...) {
  fmri_lm(onset ~ hrf(condition), block = ~run, dataset = s$dset, durations = 0,
          control = fmri_lm_control(
            estimation = estimation_spec(scope),
            noise = noise_spec(noise),
            robust = robust_spec(robust, ...),
            weights = weights
          ),
          compute = compute_spec(voxel_chunks = chunks))
}

.rb_beta <- function(fit) as.matrix(fit$result$betas$data[[1]]$estimate[[1]])
.rb_se <- function(fit) as.matrix(fit$result$betas$data[[1]]$se[[1]])

# Per-run design: the run's rows of the joint design, minus the other runs'
# (all-zero) baseline columns.
.rb_run_design <- function(s, r) {
  rows <- (r - 1L) * s$n + seq_len(s$n)
  Xr <- s$X[rows, , drop = FALSE]
  Xr[, colSums(abs(Xr)) > 0, drop = FALSE]
}

# ---- D2: the IRLS engine returns the weights it used -----------------------

test_that("robust_iterative_fitter betas are WLS with the returned weights", {
  s <- .rb_sim()
  X <- s$X[, colSums(abs(s$X)) > 0]
  proj <- fmrireg:::.fast_preproject(X)
  ctx <- fmrireg:::glm_context(X = X, Y = s$Y, proj = proj)
  for (type in c("huber", "bisquare")) {
    cfg <- list(type = type, k_huber = 1.345, c_tukey = 4.685, max_iter = 2L,
                scale_scope = "run")
    rf <- fmrireg:::robust_iterative_fitter(ctx, cfg, X)
    w <- rf$robust_weights_final
    for (v in seq_len(ncol(s$Y))) {
      ref <- lm(s$Y[, v] ~ X - 1, weights = w)
      expect_equal(unname(rf$betas_robust[, v]), unname(coef(ref)),
                   tolerance = 1e-8, info = type)
    }
    expect_equal(rf$dfres, sum(w > 0) - qr(X)$rank, info = type)
  }
})

# ---- D2: runwise and chunkwise agree, and both are WLS ---------------------

test_that("one-run robust fits agree across strategies and match lm(weights=)", {
  s <- .rb_sim(nrun = 1L)
  rw <- .rb_fit(s, scope = "runwise_meta")
  cw <- .rb_fit(s, scope = "joint")
  cw3 <- .rb_fit(s, scope = "joint", chunks = 3L)

  w <- as.numeric(cw$result$robust_weights)
  expect_equal(as.numeric(rw$result$robust_weights), w, tolerance = 1e-12)
  expect_true(any(w == 0))  # bisquare rejects the spikes outright

  expect_equal(.rb_beta(rw), .rb_beta(cw), tolerance = 1e-8)
  expect_equal(.rb_se(rw), .rb_se(cw), tolerance = 1e-8)
  expect_equal(rw$result$sigma, cw$result$sigma, tolerance = 1e-8)
  expect_equal(rw$result$rdf, cw$result$rdf)
  expect_equal(.rb_beta(cw3), .rb_beta(cw), tolerance = 1e-10)
  expect_equal(.rb_se(cw3), .rb_se(cw), tolerance = 1e-10)

  X <- s$X
  for (v in seq_len(ncol(s$Y))) {
    ref <- lm(s$Y[, v] ~ X - 1, weights = w)
    sm <- summary(ref)
    expect_equal(unname(.rb_beta(cw)[v, ]), unname(coef(ref)), tolerance = 1e-8)
    expect_equal(unname(.rb_se(cw)[v, ]), unname(sm$coefficients[, 2]),
                 tolerance = 1e-8)
    expect_equal(cw$result$sigma[v], sm$sigma, tolerance = 1e-8)
    expect_equal(cw$result$rdf, ref$df.residual)
  }
})

test_that("two-run robust fits match WLS: joint on stacked data, runwise by run", {
  s <- .rb_sim(nrun = 2L, seed = 4L)
  cw <- .rb_fit(s, scope = "joint")
  rw <- .rb_fit(s, scope = "runwise_meta")

  # Joint: WLS on the stacked design with run-specific baselines.
  w <- as.numeric(cw$result$robust_weights)
  expect_length(w, nrow(s$Y))
  X <- s$X
  for (v in seq_len(ncol(s$Y))) {
    ref <- lm(s$Y[, v] ~ X - 1, weights = w)
    sm <- summary(ref)
    expect_equal(unname(.rb_beta(cw)[v, ]), unname(coef(ref)), tolerance = 1e-8)
    expect_equal(unname(.rb_se(cw)[v, ]), unname(sm$coefficients[, 2]),
                 tolerance = 1e-8)
    expect_equal(cw$result$rdf, ref$df.residual)
  }

  # Runwise: per-run WLS with the run's weights, then fixed-effect pooling.
  w_runs <- rw$result$robust_weights
  expect_length(w_runs, 2L)
  expect_equal(unlist(w_runs, use.names = FALSE), w, tolerance = 1e-12)
  b_rw <- .rb_beta(rw)
  se_rw <- .rb_se(rw)
  for (v in seq_len(ncol(s$Y))) {
    per_run <- lapply(1:2, function(r) {
      rows <- (r - 1L) * s$n + seq_len(s$n)
      Xr <- .rb_run_design(s, r)
      sm <- summary(lm(s$Y[rows, v] ~ Xr - 1, weights = w_runs[[r]]))
      sm$coefficients[s$ec, 1:2]
    })
    b <- sapply(per_run, function(m) m[, 1])
    se <- sapply(per_run, function(m) m[, 2])
    prec <- 1 / se^2
    expect_equal(unname(b_rw[v, s$ec]), unname(rowSums(b * prec) / rowSums(prec)),
                 tolerance = 1e-8)
    expect_equal(unname(se_rw[v, s$ec]), unname(sqrt(1 / rowSums(prec))),
                 tolerance = 1e-8)
  }

  # The strategies estimate the same effect.
  expect_equal(unname(.rb_beta(rw)[, s$ec]), unname(.rb_beta(cw)[, s$ec]),
               tolerance = 0.05)
})

test_that("zero bisquare weights do not count toward residual df", {
  s <- .rb_sim(nrun = 1L, spike_size = 60)
  cw <- .rb_fit(s, scope = "joint")
  w <- as.numeric(cw$result$robust_weights)
  n_zero <- sum(w == 0)
  expect_gte(n_zero, 4L)
  rank <- qr(s$X)$rank
  expect_equal(cw$result$rdf, length(w) - n_zero - rank)
})

# ---- D1: re-estimating phi uses unwhitened residuals -----------------------

test_that("robust reestimate_phi recovers phi and matches the plain AR fit", {
  s <- .rb_sim(nrun = 1L, n = 300L, V = 40L, phi = 0.5, seed = 21L,
               spike_size = 0)
  base <- .rb_fit(s, scope = "runwise_meta", noise = "ar1", robust = "huber")
  re <- .rb_fit(s, scope = "runwise_meta", noise = "ar1", robust = "huber",
                reestimate_phi = TRUE)
  phi_base <- unlist(base$result$ar_coef)
  phi_re <- unlist(re$result$ar_coef)
  expect_equal(phi_base, 0.5, tolerance = 0.15)
  # Whitened residuals carry no AR(1) structure, so estimating from them
  # collapses phi to ~0. Unwhitened robust residuals recover it.
  expect_equal(phi_re, phi_base, tolerance = 0.05)
  expect_equal(.rb_beta(re)[, s$ec], .rb_beta(base)[, s$ec], tolerance = 0.05)
  expect_equal(.rb_se(re)[, s$ec], .rb_se(base)[, s$ec], tolerance = 0.05)

  # The reported scale belongs to the final (re-whitened) fit.
  w <- as.numeric(re$result$robust_weights)
  Xw <- fmrireg:::ar_whiten_transform(s$X, s$Y, phi_re, FALSE)
  ref <- lm(Xw$Y[, 1] ~ Xw$X - 1, weights = w)
  expect_equal(unname(.rb_beta(re)[1, ]), unname(coef(ref)), tolerance = 1e-8)
  expect_equal(re$result$sigma[1], summary(ref)$sigma, tolerance = 1e-8)
})

# ---- D3: runwise volume weights on multi-run data --------------------------

test_that("runwise volume weights slice by run and match per-run WLS", {
  s <- .rb_sim(nrun = 2L, seed = 8L, spike_size = 0)
  set.seed(30)
  w <- runif(nrow(s$Y), 0.3, 1)
  w[c(3, 4, 50, 121, 122, 200)] <- 0
  rw <- .rb_fit(s, scope = "runwise_meta", robust = "none",
                weights = weights_spec(values = w))
  jt <- .rb_fit(s, scope = "joint", robust = "none",
                weights = weights_spec(values = w))

  X <- s$X
  df_runs <- 0
  b_rw <- .rb_beta(rw)
  se_rw <- .rb_se(rw)
  for (v in seq_len(ncol(s$Y))) {
    # Joint: stacked WLS with run-specific baselines.
    ref <- lm(s$Y[, v] ~ X - 1, weights = w)
    expect_equal(unname(.rb_beta(jt)[v, ]), unname(coef(ref)), tolerance = 1e-8)
    expect_equal(jt$result$rdf, ref$df.residual)

    per_run <- lapply(1:2, function(r) {
      rows <- (r - 1L) * s$n + seq_len(s$n)
      Xr <- .rb_run_design(s, r)
      fit_r <- lm(s$Y[rows, v] ~ Xr - 1, weights = w[rows])
      list(cf = summary(fit_r)$coefficients[s$ec, 1:2], df = fit_r$df.residual)
    })
    if (v == 1L) df_runs <- sum(vapply(per_run, `[[`, numeric(1), "df"))
    b <- sapply(per_run, function(m) m$cf[, 1])
    se <- sapply(per_run, function(m) m$cf[, 2])
    prec <- 1 / se^2
    expect_equal(unname(b_rw[v, s$ec]), unname(rowSums(b * prec) / rowSums(prec)),
                 tolerance = 1e-8)
    expect_equal(unname(se_rw[v, s$ec]), unname(sqrt(1 / rowSums(prec))),
                 tolerance = 1e-8)
  }
  # Zero-weight volumes are dropped from the residual df, as in lm.wfit().
  expect_equal(rw$result$rdf, df_runs)
})
