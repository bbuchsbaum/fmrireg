# Review 2026-10, stream C: estimate_betas() OLS / mixed / LSS paths.
#
# C1: the OLS and mixed paths removed the baseline from Y but not from the
# event regressors, so the betas were biased whenever the baseline was not
# orthogonal to the events (always, for a constant baseline and HRF regressors).
# P-C1: the LSS fixed-effect betas and the OLS path were per-voxel loops.

.beta_review_fixture <- function(n = 120L, V = 4L, basis = "poly", degree = 2L,
                                 seed = 2026L, noise_sd = 0) {
  set.seed(seed)
  etab <- data.frame(
    onset = seq(6, 96, by = 10),
    cond = factor(rep(c("A", "B"), 5)),
    run = 1L
  )
  Y0 <- matrix(0, n, V)
  dset0 <- matrix_frame(Y0, TR = 1, run_length = n, event_table = etab)
  bmod <- baseline_model(basis = basis, degree = degree,
                         sframe = frame_sframe(dset0))
  bdes <- fmrireg:::gen_beta_design(
    fixed = NULL, ran = onset ~ hrf(cond), block = ~ run,
    bmod = bmod, dset = dset0
  )
  X <- as.matrix(bdes$dmat_ran)
  Base <- as.matrix(bdes$dmat_base)
  B <- rbind(seq(2, by = 0.5, length.out = V), seq(1, by = -0.25, length.out = V))
  G <- matrix(rnorm(ncol(Base) * V, sd = 3), ncol(Base), V)
  Y <- X %*% B + Base %*% G + noise_sd * matrix(rnorm(n * V), n, V)
  dset <- matrix_frame(Y, TR = 1, run_length = n, event_table = etab)
  list(dset = dset, bmod = bmod, X = X, Base = Base, B = B, Y = Y, etab = etab)
}

test_that("C1: estimate_betas(method = 'ols') recovers noiseless betas (constant baseline)", {
  fx <- .beta_review_fixture(basis = "constant", degree = 1L)
  out <- estimate_betas(fx$dset, ran = onset ~ hrf(cond), block = ~ run,
                        method = "ols", progress = FALSE)
  expect_equal(unname(out$betas_ran), unname(fx$B), tolerance = 1e-8)

  # Same as stats::lm with the baseline as covariates.
  ref <- coef(lm(fx$Y ~ 0 + fx$X + fx$Base))[seq_len(ncol(fx$X)), , drop = FALSE]
  expect_equal(unname(out$betas_ran), unname(ref), tolerance = 1e-8)
})

test_that("C1: glm_ols recovers noiseless betas with a polynomial drift baseline", {
  fx <- .beta_review_fixture(basis = "poly", degree = 2L)
  sf <- frame_sframe(fx$dset)
  emod <- event_model(onset ~ hrf(cond), data = fx$etab, block = ~ run,
                      sampling_frame = sf)
  out <- glm_ols(fx$dset, emod, fmrihrf::HRF_SPMG1, basemod = fx$bmod,
                 block = ~ run, progress = FALSE)
  expect_equal(unname(out$betas_ran), unname(fx$B), tolerance = 1e-8)
})

test_that("C1: OLS path with fixed effects equals lm on the joint design", {
  fx <- .beta_review_fixture(basis = "poly", degree = 1L, noise_sd = 0.5)
  etab <- fx$etab
  etab$rt <- seq(-1, 1, length.out = nrow(etab))
  dset <- matrix_frame(fx$Y, TR = 1, run_length = nrow(fx$Y), event_table = etab)
  out <- estimate_betas(dset, fixed = onset ~ hrf(rt), ran = onset ~ hrf(cond),
                        block = ~ run, method = "ols", basemod = fx$bmod,
                        progress = FALSE)
  Xf <- as.matrix(out$design_fixed)
  ref <- coef(lm(fx$Y ~ 0 + fx$X + Xf + fx$Base))
  k <- ncol(fx$X)
  expect_equal(unname(out$betas_ran), unname(ref[seq_len(k), , drop = FALSE]),
               tolerance = 1e-8)
  expect_equal(unname(out$betas_fixed),
               unname(ref[k + seq_len(ncol(Xf)), , drop = FALSE]),
               tolerance = 1e-8)
})

test_that("C1: mixed path treats the baseline as fixed effects (REML/BLUP reference)", {
  skip_if_not_installed("fmrilss")
  fx <- .beta_review_fixture(basis = "poly", degree = 2L, noise_sd = 0.3, V = 3L)
  out <- estimate_betas(fx$dset, ran = onset ~ hrf(cond), block = ~ run,
                        method = "mixed", basemod = fx$bmod, progress = FALSE)

  for (v in seq_len(ncol(fx$Y))) {
    y <- fx$Y[, v]
    ref <- fmrilss::mixed_solve(Y = y, Z = fx$X, X = fx$Base)
    expect_equal(unname(out$betas_ran[, v]), as.numeric(ref$u), tolerance = 1e-6)

    # Independent check: Henderson's mixed-model equations at the REML
    # variance ratio, with the baseline as fixed effects.
    lambda <- ref$Ve / ref$Vu
    W <- cbind(fx$Base, fx$X)
    q <- ncol(fx$Base)
    lhs <- crossprod(W)
    idx <- q + seq_len(ncol(fx$X))
    lhs[idx, idx] <- lhs[idx, idx] + diag(lambda, ncol(fx$X))
    sol <- solve(lhs, crossprod(W, y))
    expect_equal(unname(out$betas_ran[, v]), unname(sol[idx]), tolerance = 1e-6)
  }
})

test_that("P-C1: LSS fixed-effect betas equal per-voxel lm.fit results", {
  skip_if_not_installed("fmrilss")
  fx <- .beta_review_fixture(basis = "poly", degree = 1L, noise_sd = 0.5, V = 5L)
  etab <- fx$etab
  etab$rt <- seq(-1, 1, length.out = nrow(etab))
  dset <- matrix_frame(fx$Y, TR = 1, run_length = nrow(fx$Y), event_table = etab)
  out <- estimate_betas(dset, fixed = onset ~ hrf(rt), ran = onset ~ trialwise(),
                        block = ~ run, method = "lss", basemod = fx$bmod,
                        progress = FALSE)
  Xbf <- cbind(as.matrix(out$design_base), as.matrix(out$design_fixed))
  nb <- ncol(as.matrix(out$design_base))
  ref <- vapply(seq_len(ncol(fx$Y)), function(v) {
    cf <- coef(lm.fit(Xbf, fx$Y[, v]))
    cf[(nb + 1):length(cf)]
  }, numeric(ncol(Xbf) - nb))
  expect_equal(unname(out$betas_fixed), unname(matrix(ref, ncol = ncol(fx$Y))),
               tolerance = 1e-10)
})
