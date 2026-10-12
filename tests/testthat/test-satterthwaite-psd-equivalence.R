# Model-based Satterthwaite df with the PSD projection of the residual
# covariance: the fast banded path must equal the dense reference to 1e-10
# relative, including voxels whose covariance the projection changes.
#
# Reference (copied, so it does not depend on package code): build the
# block-diagonal lag-window covariance V of each voxel's residuals as an n x n
# matrix, project it onto the PSD cone with a dense eigen() (eigenvalues below
# max(max|lambda|, 1) * sqrt(eps) are set to 0), and evaluate
# tr(RV)^2 / tr(RVRV).

.sp_hac_weight <- function(h, H, taper) {
  if (h == 0L || identical(taper, "none")) return(1)
  x <- h / (H + 1)
  if (identical(taper, "tukey")) return(0.5 * (1 + cos(pi * x)))
  if (x <= 0.5) 1 - 6 * x^2 + 6 * x^3 else 2 * (1 - x)^3
}

.sp_dense_V <- function(residual, segments, max_lag, taper) {
  n <- length(residual)
  V <- matrix(0, n, n)
  for (idx in segments) {
    e <- residual[idx]
    m <- length(e)
    H <- min(as.integer(max_lag), m - 1L)
    gamma <- numeric(H + 1L)
    for (h in 0:H) {
      gamma[h + 1L] <- sum(e[(h + 1L):m] * e[seq_len(m - h)]) / (m - h) *
        .sp_hac_weight(h, max_lag, taper)
    }
    distance <- abs(outer(seq_len(m), seq_len(m), "-"))
    block <- matrix(0, m, m)
    inside <- distance <= H
    block[inside] <- gamma[distance[inside] + 1L]
    V[idx, idx] <- block
  }
  V
}

.sp_psd <- function(M, tol_mult = 1) {
  M <- (M + t(M)) / 2
  ev <- eigen(M, symmetric = TRUE)
  tol <- max(abs(ev$values), 1) * .Machine$double.eps^0.5 * tol_mult
  ev$values[ev$values < tol] <- 0
  scaled <- sweep(ev$vectors, 2L, sqrt(ev$values), `*`)
  tcrossprod(scaled)
}

.sp_dense_df <- function(X, E, segments, max_lag, taper, tol_mult = 1) {
  X <- as.matrix(X)
  XtXinv <- solve(crossprod(X))
  A <- X %*% XtXinv
  vapply(seq_len(ncol(E)), function(v) {
    V <- .sp_psd(.sp_dense_V(E[, v], segments, max_lag, taper), tol_mult)
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

# Eigenvalues of the (unprojected) reference covariance and its tolerance.
.sp_spectrum <- function(residual, segments, max_lag, taper) {
  ev <- eigen(.sp_dense_V(residual, segments, max_lag, taper),
              symmetric = TRUE, only.values = TRUE)$values
  list(values = ev, tol = max(abs(ev), 1) * .Machine$double.eps^0.5)
}

.sp_n_clipped <- function(residual, segments, max_lag, taper) {
  s <- .sp_spectrum(residual, segments, max_lag, taper)
  sum(s$values < s$tol)
}

.sp_fast_df <- function(X, E, runs, censor, max_lag, taper) {
  ctx <- list(X = X, residuals = E, runs = runs, censor = censor)
  fmrireg:::.fmri_lm_variance_from_context(
    ctx, variance_spec("model", max_lag = max_lag, taper = taper,
                       df = "satterthwaite")
  )$df
}

.sp_max_rel <- function(a, b) max(abs(a - b) / abs(b))

.sp_problem <- function(n = 120L, p = 5L, seed = 1L,
                        ar = c(0.1, 0.3, 0.5, 0.7, 0.8, 0.9, -0.4, 0.6)) {
  set.seed(seed)
  X <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  E <- vapply(ar, function(a) as.numeric(arima.sim(list(ar = a), n)),
              numeric(n))
  E <- E - X %*% qr.solve(X, E)
  list(X = X, E = E)
}

test_that("Satterthwaite df equals the dense PSD reference for every taper, lag and layout", {
  pr <- .sp_problem()
  n <- nrow(pr$X)
  layouts <- list(
    one_run = list(runs = list(seq_len(n)), censor = NULL),
    two_runs = list(runs = list(1:60, 61:n), censor = NULL),
    # Odd-length segments and an isolated single row (row 11).
    censored = list(runs = list(1:60, 61:n),
                    censor = c(10L, 12L, 70L, 100:103))
  )
  clipped <- integer(0)
  for (lay in names(layouts)) {
    runs <- layouts[[lay]]$runs
    censor <- layouts[[lay]]$censor
    segs <- fmrireg:::.fmri_lm_segments(n, runs, censor)
    # variance_spec() offers exactly these lag windows (no Bartlett).
    for (taper in c("none", "parzen", "tukey")) {
      for (L in c(0L, 1L, 3L, 8L, 20L, 40L)) {
        got <- .sp_fast_df(pr$X, pr$E, runs, censor, L, taper)
        ref <- .sp_dense_df(pr$X, pr$E, segs, L, taper)
        expect_lt(.sp_max_rel(got, ref), 1e-10,
                  label = paste("max rel diff", lay, taper, L))
        clipped <- c(clipped, vapply(seq_len(ncol(pr$E)), function(v) {
          .sp_n_clipped(pr$E[, v], segs, L, taper)
        }, integer(1)))
      }
    }
  }
  # The grid must exercise voxels with no, one/few, and many clipped
  # eigenvalues.
  expect_true(any(clipped == 0L))
  expect_true(any(clipped >= 1L & clipped <= 3L))
  expect_true(any(clipped >= 20L))
})

test_that("Satterthwaite df matches the reference with many clipped eigenvalues", {
  # Strongly autocorrelated residuals with the truncated window at a short lag
  # make a large part of the spectrum negative.
  set.seed(7)
  n <- 150L
  X <- cbind(1, matrix(rnorm(n * 3L), n, 3L))
  E <- vapply(c(0.85, 0.9, 0.95), function(a) {
    as.numeric(arima.sim(list(ar = a), n))
  }, numeric(n))
  E <- E - X %*% qr.solve(X, E)
  for (lay in list(list(runs = list(seq_len(n)), censor = NULL),
                   list(runs = list(1:75, 76:n), censor = c(30L, 31L)))) {
    segs <- fmrireg:::.fmri_lm_segments(n, lay$runs, lay$censor)
    for (L in c(3L, 8L)) {
      nclip <- vapply(1:3, function(v) .sp_n_clipped(E[, v], segs, L, "none"),
                      integer(1))
      got <- .sp_fast_df(X, E, lay$runs, lay$censor, L, "none")
      ref <- .sp_dense_df(X, E, segs, L, "none")
      expect_lt(.sp_max_rel(got, ref), 1e-10,
                label = paste("max rel diff, clipped:", toString(nclip)))
    }
  }
})

test_that("Satterthwaite df matches the reference for tiny-scale residuals", {
  pr <- .sp_problem(seed = 3L)
  n <- nrow(pr$X)
  runs <- list(1:60, 61:n)
  segs <- fmrireg:::.fmri_lm_segments(n, runs)
  # 1e-6: the absolute tolerance (>= 1.5e-8) exceeds every eigenvalue, so V
  # projects to 0 and df = 1. 1e-4: the tolerance cuts through the spectrum.
  for (scale in c(1e-6, 1e-4, 3e-4)) {
    E <- pr$E * scale
    for (taper in c("none", "parzen")) {
      for (L in c(3L, 20L)) {
        got <- .sp_fast_df(pr$X, E, runs, NULL, L, taper)
        ref <- .sp_dense_df(pr$X, E, segs, L, taper)
        expect_lt(.sp_max_rel(got, ref), 1e-10,
                  label = paste("max rel diff", scale, taper, L))
      }
    }
  }
  expect_equal(.sp_fast_df(pr$X, pr$E * 1e-6, runs, NULL, 8L, "none"),
               rep(1, ncol(pr$E)))
})

test_that("clip decisions match the reference next to the tolerance", {
  # Scale one voxel's residuals so that a chosen eigenvalue of its covariance
  # sits just below or just above tol = sqrt(eps) (the covariance is scaled so
  # that max|lambda| < 1, which fixes tol and makes the clipped eigenvalue a
  # visible part of the traces). The gap to tol is a small multiple of the
  # eigenvalue rounding level eps * max|lambda|.
  pr <- .sp_problem(n = 100L, seed = 5L, ar = c(0.5, 0.8))
  n <- nrow(pr$X)
  runs <- list(1:50, 51:n)
  segs <- fmrireg:::.fmri_lm_segments(n, runs)
  tol <- .Machine$double.eps^0.5
  checked <- 0L
  for (v in 1:2) {
    for (j in c(1L, 2L, 5L)) {
      for (side in c(-1, 1)) {
        e <- pr$E[, v]
        s <- .sp_spectrum(e, segs, 8L, "parzen")
        pos <- sort(s$values[s$values > 0])
        target <- pos[j]
        for (iter in 1:3) {
          # Aim at tol + side * 200 ulp-scale units, refining for rounding.
          spread <- max(abs(s$values))
          goal <- tol + side * 200 * .Machine$double.eps * spread
          e <- e * sqrt(goal / target)
          s <- .sp_spectrum(e, segs, 8L, "parzen")
          target <- sort(s$values)[which.min(abs(sort(s$values) - goal))]
        }
        spread <- max(abs(s$values))
        gap <- (target - tol) / (.Machine$double.eps * spread)
        expect_true(spread < 1)
        expect_true(sign(gap) == side && abs(gap) > 20 && abs(gap) < 2000,
                    label = paste("placement gap", gap))
        E <- cbind(e)
        got <- .sp_fast_df(pr$X, E, runs, NULL, 8L, "parzen")
        ref <- .sp_dense_df(pr$X, E, segs, 8L, "parzen")
        expect_lt(.sp_max_rel(got, ref), 1e-10,
                  label = paste("max rel diff v", v, "j", j, "side", side))
        # The decision matters: flipping it (moving tol across the
        # eigenvalue) changes the df far beyond the 1e-10 tolerance.
        flipped <- .sp_dense_df(pr$X, E, segs, 8L, "parzen",
                                tol_mult = 1 + side * 1e-4)
        expect_gt(.sp_max_rel(flipped, ref), 1e-8)
        checked <- checked + 1L
      }
    }
  }
  expect_equal(checked, 12L)
})

test_that("fmri_lm Satterthwaite df equals the dense PSD reference end to end", {
  set.seed(11)
  n <- 90L
  ev <- data.frame(onset = seq(4, 170, by = 12),
                   condition = factor(rep(c("A", "B"), length.out = 14)),
                   run = 1)
  Y <- vapply(c(0.2, 0.5, 0.8, 0.9), function(a) {
    as.numeric(arima.sim(list(ar = a), n)) + 10
  }, numeric(n))
  dset <- matrix_frame(Y, TR = 2, run_length = n, event_table = ev)
  for (cfg in list(c("none", 20), c("none", 8), c("parzen", 30),
                   c("tukey", 40))) {
    taper <- cfg[1]
    L <- as.integer(cfg[2])
    fit <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = dset,
                   durations = 0,
                   control = fmri_lm_control(
                     estimation = estimation_spec("joint"),
                     variance = variance_spec("model", max_lag = L,
                                              taper = taper,
                                              df = "satterthwaite")),
                   compute = compute_spec(voxel_chunks = 1))
    X <- as.matrix(design_matrix(fit$model))
    E <- Y - X %*% qr.solve(X, Y)
    ref <- .sp_dense_df(X, E, list(seq_len(n)), L, taper)
    expect_lt(.sp_max_rel(as.numeric(fit$result$df$inference), ref), 1e-10,
              label = paste("max rel diff", taper, L))
  }
})
