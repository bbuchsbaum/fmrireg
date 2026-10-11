# Performance pass 2, item P1-1: model-based Satterthwaite df without a
# per-voxel n x n residual covariance or eigen().
#
# Reference: the dense computation as it stood before the change (copied here
# verbatim so it does not depend on package code): build the block-diagonal
# lag-window covariance V of each voxel's residuals, project it onto the PSD
# cone with eigen(), and evaluate tr(RV)^2 / tr(RVRV).

.p11_psd <- function(M) {
  M <- (M + t(M)) / 2
  ev <- eigen(M, symmetric = TRUE)
  tol <- max(abs(ev$values), 1) * .Machine$double.eps^0.5
  ev$values[ev$values < tol] <- 0
  scaled <- sweep(ev$vectors, 2L, sqrt(ev$values), `*`)
  tcrossprod(scaled)
}

.p11_dense_V <- function(residual, segments, max_lag, taper) {
  n <- length(residual)
  V <- matrix(0, n, n)
  for (idx in segments) {
    e <- residual[idx]
    m <- length(e)
    H <- min(as.integer(max_lag), m - 1L)
    gamma <- numeric(H + 1L)
    for (h in 0:H) {
      gamma[h + 1L] <- sum(e[(h + 1L):m] * e[seq_len(m - h)]) / (m - h) *
        fmrireg:::.fmri_lm_hac_weight(h, max_lag, taper)
    }
    distance <- abs(outer(seq_len(m), seq_len(m), "-"))
    block <- matrix(0, m, m)
    inside <- distance <= H
    block[inside] <- gamma[distance[inside] + 1L]
    V[idx, idx] <- block
  }
  V
}

.p11_dense_df <- function(X, E, segments, max_lag, taper) {
  X <- as.matrix(X)
  n <- nrow(X)
  XtXinv <- solve(crossprod(X))
  vapply(seq_len(ncol(E)), function(v) {
    V <- .p11_psd(.p11_dense_V(E[, v], segments, max_lag, taper))
    A <- X %*% XtXinv
    VX <- V %*% X
    B <- crossprod(X, VX)
    tr_RV <- sum(diag(V)) - sum(A * VX)
    tr_HVV <- sum(A * (V %*% VX))
    BG <- B %*% XtXinv
    tr_RVRV <- sum(V * t(V)) - 2 * tr_HVV + sum(BG * t(BG))
    if (!is.finite(tr_RV) || !is.finite(tr_RVRV) || tr_RVRV <= 0) return(1)
    max(tr_RV^2 / tr_RVRV, 1)
  }, numeric(1))
}

# Does the PSD projection change V for this voxel?
.p11_clipped <- function(residual, segments, max_lag, taper) {
  V <- .p11_dense_V(residual, segments, max_lag, taper)
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  any(ev < max(abs(ev), 1) * .Machine$double.eps^0.5)
}

.p11_problem <- function(n = 96L, V = 6L, p = 4L, seed = 1L, scale = 1) {
  set.seed(seed)
  X <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  E <- vapply(seq_len(V), function(v) {
    as.numeric(arima.sim(list(ar = 0.3 + 0.1 * (v %% 4)), n))
  }, numeric(n)) * scale
  # Residuals of the fitted model, as the inference context carries them.
  E <- E - X %*% qr.solve(X, E)
  list(X = X, E = E)
}

test_that("Satterthwaite df equals the dense reference for every taper and lag", {
  pr <- .p11_problem()
  n <- nrow(pr$X)
  layouts <- list(
    one_run = list(runs = list(seq_len(n)), censor = NULL),
    two_runs = list(runs = list(1:50, 51:n), censor = NULL),
    censored = list(runs = list(1:50, 51:n), censor = c(10L, 11L, 70L))
  )
  n_clipped <- 0L
  for (lay in names(layouts)) {
    runs <- layouts[[lay]]$runs
    censor <- layouts[[lay]]$censor
    segs <- fmrireg:::.fmri_lm_segments(n, runs, censor)
    for (taper in c("none", "tukey", "parzen")) {
      for (L in c(0L, 1L, 3L, 8L, 20L)) {
        ctx <- list(X = pr$X, residuals = pr$E, runs = runs, censor = censor)
        got <- fmrireg:::.fmri_lm_variance_from_context(
          ctx, variance_spec("model", max_lag = L, taper = taper,
                             df = "satterthwaite")
        )$df
        ref <- .p11_dense_df(pr$X, pr$E, segs, L, taper)
        expect_equal(got, ref, tolerance = 1e-8,
                     info = paste(lay, taper, L))
        n_clipped <- n_clipped + sum(vapply(seq_len(ncol(pr$E)), function(v) {
          .p11_clipped(pr$E[, v], segs, L, taper)
        }, logical(1)))
      }
    }
  }
  # The comparison above must exercise voxels whose covariance the PSD
  # projection changes (it does for the truncated window at long lags).
  expect_gt(n_clipped, 0L)
})

test_that("Satterthwaite df matches the reference on small-scale data", {
  # With residuals of order 1e-4 the PSD tolerance (absolute, >= 1.5e-8)
  # zeroes the whole covariance; the fast path must notice and agree.
  pr <- .p11_problem(scale = 1e-4, seed = 4L)
  n <- nrow(pr$X)
  segs <- fmrireg:::.fmri_lm_segments(n, list(seq_len(n)))
  for (taper in c("none", "parzen")) {
    ctx <- list(X = pr$X, residuals = pr$E, runs = list(seq_len(n)), censor = NULL)
    got <- fmrireg:::.fmri_lm_variance_from_context(
      ctx, variance_spec("model", max_lag = 4L, taper = taper,
                         df = "satterthwaite")
    )$df
    expect_equal(got, .p11_dense_df(pr$X, pr$E, segs, 4L, taper),
                 tolerance = 1e-8, info = taper)
  }
})

test_that("fmri_lm Satterthwaite df equals the dense reference end to end", {
  set.seed(11)
  n <- 80L
  ev <- data.frame(onset = seq(4, 150, by = 12),
                   condition = factor(rep(c("A", "B"), length.out = 13)),
                   run = 1)
  Y <- vapply(1:5, function(v) as.numeric(arima.sim(list(ar = 0.4), n)) + 10,
              numeric(n))
  dset <- matrix_frame(Y, TR = 2, run_length = n, event_table = ev)
  for (taper in c("none", "parzen")) {
    fit <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = dset,
                   durations = 0,
                   control = fmri_lm_control(
                     estimation = estimation_spec("joint"),
                     variance = variance_spec("model", max_lag = 6,
                                              taper = taper,
                                              df = "satterthwaite")),
                   compute = compute_spec(voxel_chunks = 1))
    X <- as.matrix(design_matrix(fit$model))
    E <- Y - X %*% qr.solve(X, Y)
    ref <- .p11_dense_df(X, E, list(seq_len(n)), 6L, taper)
    expect_equal(as.numeric(fit$result$df$inference), ref, tolerance = 1e-8,
                 info = taper)
  }
})

test_that("Satterthwaite df forms no n x n covariance and calls no eigen()", {
  pr <- .p11_problem(n = 120L, V = 5L)
  n <- nrow(pr$X)
  runs <- list(1:60, 61:120)
  psd_dims <- integer(0)
  projected <- 0L
  orig_psd <- fmrireg:::.fmri_lm_psd
  orig_projected <- fmrireg:::.fmri_lm_projected_traces
  testthat::local_mocked_bindings(
    .fmri_lm_projected_traces = function(...) {
      projected <<- projected + 1L
      orig_projected(...)
    },
    .fmri_lm_psd = function(M) {
      psd_dims <<- c(psd_dims, nrow(M))
      orig_psd(M)
    },
    .fmri_lm_residual_covariance = function(...) {
      stop("dense residual covariance formed")
    },
    .satterthwaite_df = function(...) stop("dense Satterthwaite df called")
  )
  ctx <- list(X = pr$X, residuals = pr$E, runs = runs, censor = NULL)
  # Parzen window at a short lag: every voxel takes the banded fast path,
  # which needs no eigendecomposition at all.
  df <- fmrireg:::.fmri_lm_variance_from_context(
    ctx, variance_spec("model", max_lag = 3L, taper = "parzen",
                       df = "satterthwaite")
  )$df
  expect_length(df, ncol(pr$E))
  expect_length(psd_dims, 0L)
  expect_equal(projected, 0L)
  # Truncated window at a long lag: voxels whose covariance needs the PSD
  # projection are projected one run block at a time, never as n x n.
  df2 <- fmrireg:::.fmri_lm_variance_from_context(
    ctx, variance_spec("model", max_lag = 30L, taper = "none",
                       df = "satterthwaite")
  )$df
  expect_length(df2, ncol(pr$E))
  expect_gt(projected, 0L)
  expect_true(all(psd_dims < n))
})
