# Performance pass 2, P2-1: estimate_betas(method = "mixed") decomposes the
# shared design once and runs a 1-D REML search per voxel, instead of a full
# fmrilss::mixed_solve() (two n x n eigen/SVD decompositions) per voxel.
# Reference: fmrilss::mixed_solve() per voxel.

.p2_rel_err <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  max(abs(a - b)) / max(max(abs(b)), .Machine$double.xmin)
}

# Per-voxel reference: mixed_solve with the same fixed design the mixed path
# builds (fixed effects of interest, then the non-redundant baseline columns).
.p2_mixed_reference <- function(Y, Z, X_fixed) {
  lapply(seq_len(ncol(Y)), function(v) {
    fmrilss::mixed_solve(Y = Y[, v], Z = Z, X = X_fixed)
  })
}

.p2_sim_design <- function(n, m, nuis = c("constant", "poly", "none"),
                           n_fixed = 0L, seed = 1L) {
  nuis <- match.arg(nuis)
  set.seed(seed)
  tt <- seq_len(n)
  # Smooth, overlapping trial regressors (like HRF-convolved trial onsets).
  onsets <- sort(sample(seq(1, n - 2), m))
  Z <- vapply(onsets, function(o) {
    x <- pmax(tt - o, 0)
    x^2 * exp(-x / 1.5) / 2
  }, numeric(n))
  Xi <- if (n_fixed > 0L) matrix(rnorm(n * n_fixed), n, n_fixed) else matrix(0, n, 0)
  Base <- switch(nuis,
    constant = matrix(1, n, 1),
    poly = cbind(1, stats::poly(tt, 2)),
    none = NULL
  )
  list(Z = Z, Xi = Xi, Base = Base, n = n, m = m)
}

.p2_sim_data <- function(des, V, signal_sd, noise_sd, scale = 1, seed = 2L) {
  set.seed(seed)
  n <- des$n
  Y <- des$Z %*% matrix(rnorm(des$m * V, sd = signal_sd), des$m, V) +
    noise_sd * matrix(rnorm(n * V), n, V)
  if (ncol(des$Xi) > 0) {
    Y <- Y + des$Xi %*% matrix(rnorm(ncol(des$Xi) * V), ncol(des$Xi), V)
  }
  if (!is.null(des$Base)) {
    Y <- Y + des$Base %*% matrix(rnorm(ncol(des$Base) * V, sd = 3), ncol(des$Base), V)
  }
  scale * Y
}

test_that("P2-1: batched mixed solver matches fmrilss::mixed_solve per voxel", {
  skip_if_not_installed("fmrilss")
  cases <- list(
    list(nuis = "poly", n = 120, m = 12, n_fixed = 1L, signal = 1, noise = 0.5, scale = 1),
    list(nuis = "constant", n = 90, m = 10, n_fixed = 0L, signal = 1, noise = 1, scale = 1e-4),
    list(nuis = "none", n = 80, m = 8, n_fixed = 2L, signal = 2, noise = 0.3, scale = 1e4),
    list(nuis = "poly", n = 100, m = 15, n_fixed = 0L, signal = 0.05, noise = 1, scale = 1),
    # n <= m + p: fmrilss takes its eigen (not Cholesky/SVD) branch.
    list(nuis = "constant", n = 30, m = 28, n_fixed = 1L, signal = 1, noise = 0.5, scale = 1)
  )
  for (cs in cases) {
    des <- .p2_sim_design(cs$n, cs$m, cs$nuis, cs$n_fixed, seed = cs$n)
    Y <- .p2_sim_data(des, V = 6L, cs$signal, cs$noise, cs$scale, seed = cs$m)
    X_fixed <- fmrireg:::.mixed_fixed_design(des$Xi, des$Base)
    batch <- fmrireg:::.mixed_solve_batch(Y, X = X_fixed, Z = des$Z)
    ref <- .p2_mixed_reference(Y, des$Z, X_fixed)
    for (v in seq_len(ncol(Y))) {
      lab <- paste(cs$nuis, cs$n, cs$m, "voxel", v)
      expect_lt(.p2_rel_err(batch$u[, v], ref[[v]]$u), 1e-6, label = lab)
      expect_lt(.p2_rel_err(batch$beta[, v], ref[[v]]$beta), 1e-6, label = lab)
      expect_lt(.p2_rel_err(batch$Vu[v], ref[[v]]$Vu), 1e-6, label = lab)
      expect_lt(.p2_rel_err(batch$Ve[v], ref[[v]]$Ve), 1e-6, label = lab)
    }
  }
})

test_that("P2-1: bound and flat-likelihood voxels match fmrilss::mixed_solve", {
  skip_if_not_installed("fmrilss")
  des <- .p2_sim_design(100, 10, "constant", 0L, seed = 11L)
  X_fixed <- fmrireg:::.mixed_fixed_design(des$Xi, des$Base)
  set.seed(2)
  y_signal_only <- 4 + des$Z %*% rnorm(des$m)          # Ve -> 0: lower bound
  # Pure noise / negligible trial signal: the REML surface is flat in lambda
  # and L-BFGS-B stops on the plateau wherever fmrilss's arithmetic leads it.
  y_noise_only <- 4 + rnorm(des$n)
  y_weak <- 4 + 1e-3 * des$Z %*% rep(1, des$m) + rnorm(des$n)
  Y <- cbind(y_signal_only, y_noise_only, y_weak)
  batch <- fmrireg:::.mixed_solve_batch(Y, X = X_fixed, Z = des$Z)
  ref <- .p2_mixed_reference(Y, des$Z, X_fixed)
  lam_ref <- vapply(ref, function(r) r$Ve / r$Vu, numeric(1))
  # The fixture really does put a voxel on the lower bound and voxels on the
  # plateau (lambda far above every eigenvalue theta).
  expect_lte(lam_ref[1], 1e-9 * (1 + 1e-6))
  expect_true(all(lam_ref[2:3] > 1e4))
  for (v in seq_len(ncol(Y))) {
    expect_lt(.p2_rel_err(batch$lambda[v], lam_ref[v]), 1e-6)
    expect_lt(.p2_rel_err(batch$u[, v], ref[[v]]$u), 1e-6)
    expect_lt(.p2_rel_err(batch$beta[, v], ref[[v]]$beta), 1e-6)
    expect_lt(.p2_rel_err(batch$Vu[v], ref[[v]]$Vu), 1e-6)
  }
})

.p2_beta_fixture <- function(basis = "poly", degree = 2L, V = 5L, seed = 7L,
                             with_fixed = FALSE) {
  set.seed(seed)
  n <- 140L
  etab <- data.frame(onset = seq(6, 126, by = 8), run = 1L)
  etab$cond <- factor(rep(c("A", "B"), length.out = nrow(etab)))
  etab$rt <- seq(-1, 1, length.out = nrow(etab))
  dset0 <- matrix_frame(matrix(0, n, V), TR = 1, run_length = n, event_table = etab)
  bmod <- baseline_model(basis = basis, degree = degree, sframe = frame_sframe(dset0))
  fixed <- if (with_fixed) onset ~ hrf(rt) else NULL
  bdes <- fmrireg:::gen_beta_design(fixed = fixed, ran = onset ~ trialwise(),
                                    block = ~ run, bmod = bmod, dset = dset0)
  Z <- as.matrix(bdes$dmat_ran)
  Y <- Z %*% matrix(rnorm(ncol(Z) * V), ncol(Z), V) +
    as.matrix(bdes$dmat_base) %*% matrix(rnorm(ncol(bdes$dmat_base) * V, sd = 2),
                                         ncol(bdes$dmat_base), V) +
    matrix(rnorm(n * V, sd = 0.7), n, V)
  Y[, V] <- 3 + rnorm(n)   # a noise-only voxel
  dset <- matrix_frame(Y, TR = 1, run_length = n, event_table = etab)
  list(dset = dset, bmod = bmod, fixed = fixed, Y = Y, bdes = bdes)
}

test_that("P2-1: estimate_betas(method = 'mixed') equals per-voxel mixed_solve", {
  skip_if_not_installed("fmrilss")
  for (wf in c(FALSE, TRUE)) {
    fx <- .p2_beta_fixture(with_fixed = wf)
    out <- estimate_betas(fx$dset, fixed = fx$fixed, ran = onset ~ trialwise(),
                          block = ~ run, method = "mixed", basemod = fx$bmod,
                          progress = FALSE)
    Z <- as.matrix(out$design_ran)
    Xi <- if (wf) as.matrix(out$design_fixed) else matrix(0, nrow(Z), 0)
    X_fixed <- fmrireg:::.mixed_fixed_design(Xi, as.matrix(out$design_base))
    ref <- .p2_mixed_reference(fx$Y, Z, X_fixed)
    for (v in seq_len(ncol(fx$Y))) {
      expect_lt(.p2_rel_err(out$betas_ran[, v], ref[[v]]$u), 1e-6)
      if (wf) {
        expect_lt(.p2_rel_err(out$betas_fixed[, v],
                              ref[[v]]$beta[seq_len(ncol(Xi))]), 1e-6)
      }
    }
  }
})

test_that("P2-1 guard: the mixed design is decomposed once, not per voxel", {
  skip_if_not_installed("fmrilss")
  fx <- .p2_beta_fixture(V = 6L)
  n_solve <- 0L
  n_spectral <- 0L
  real_solve <- fmrilss::mixed_solve
  testthat::local_mocked_bindings(
    mixed_solve = function(...) {
      n_solve <<- n_solve + 1L
      real_solve(...)
    },
    .package = "fmrilss"
  )
  real_spectral <- get(".mixed_spectral_design", envir = asNamespace("fmrireg"))
  testthat::local_mocked_bindings(
    .mixed_spectral_design = function(...) {
      n_spectral <<- n_spectral + 1L
      real_spectral(...)
    }
  )
  out <- estimate_betas(fx$dset, ran = onset ~ trialwise(), block = ~ run,
                        method = "mixed", basemod = fx$bmod, progress = FALSE)
  expect_equal(ncol(out$betas_ran), 6L)
  expect_identical(n_solve, 0L)
  expect_identical(n_spectral, 1L)
})

test_that("P2-1: a rank-deficient fixed design gives NA betas with one warning", {
  skip_if_not_installed("fmrilss")
  des <- .p2_sim_design(60, 6, "none", 0L, seed = 3L)
  Xf <- cbind(1, 1)  # exactly collinear fixed effects of interest
  X <- cbind(des$Z, matrix(1, 60, 2))
  Y <- matrix(rnorm(60 * 3), 60, 3)
  expect_warning(
    out <- fmrireg:::mixed_betas_matrix(X, Y, ran_ind = 1:6, fixed_ind = 7:8),
    "X not full rank"
  )
  expect_equal(dim(out), c(8L, 3L))
  expect_true(all(is.na(out)))
  # Same verdict as the per-voxel path.
  expect_warning(one <- fmrireg:::mixed_betas(X, Y[, 1], 1:6, 7:8), "not full rank")
  expect_true(all(is.na(one)))
})
