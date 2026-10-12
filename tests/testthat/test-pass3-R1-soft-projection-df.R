# Pass 3, R1-1: residual df after soft-subspace projection beyond runwise iid.
#
# Runwise iid used tr((I - H) P^2) since pass 2. The AR runwise path still
# used n - rank after projecting and whitening, volume weights + projection
# ignored the zero-weight volumes, and the robust, chunkwise (joint),
# reference-backend and MA paths either ignored the projection in the final
# fit or had no tractable df -- they now refuse the combination.

.p3_soft_fixture <- function(n = 150L, V = 3000L, k = 15L, phi = 0.3,
                             seed = 31L) {
  set.seed(seed)
  ons <- seq(10, 2 * n - 30, by = 18)
  ev <- data.frame(onset = ons,
                   condition = factor(rep(c("A", "B"), length.out = length(ons))),
                   run = 1)
  N <- matrix(rnorm(n * k), n, k) %*% diag(seq(0.5, 6, length.out = k))
  E <- matrix(rnorm(n * V), n, V)                # innovations, sd 1
  if (phi != 0) for (t in 2:n) E[t, ] <- phi * E[t - 1L, ] + E[t, ]
  list(dset = matrix_frame(E, TR = 2, run_length = n, event_table = ev),
       N = N, n = n)
}

.p3_soft_control <- function(N, noise = noise_spec("iid"), lambda = 20,
                             scope = "runwise_meta", ...) {
  fmri_lm_control(
    estimation = estimation_spec(scope),
    noise = noise,
    projection = projection_spec("soft_subspace", nuisance_matrix = N,
                                 lambda = lambda),
    ...
  )
}

.p3_task_t <- function(fit) {
  tt <- fit$result$betas$data[[1]]$stat[[1]]
  tt[, seq_along(fit$result$event_indices), drop = FALSE]
}

test_that("AR df helper equals the dense trace tr((I-H) W P Sigma P W')", {
  set.seed(5)
  n <- 90L; k <- 8L
  N <- matrix(rnorm(n * k), n, k)
  X <- cbind(1, rnorm(n), sin(seq_len(n) / 4))
  for (lam in list(0.5, 30, "auto")) {
    sp <- soft_projection(N, lambda = lam)
    P <- sp$P_lambda(diag(n))
    for (ef in c(FALSE, TRUE)) {
      W <- fmrireg:::ar_whiten_transform(matrix(0, n, 0), diag(n), c(0.4, 0.15), ef)$Y
      Xw <- W %*% P %*% X
      H <- Xw %*% solve(crossprod(Xw), t(Xw))
      Sigma <- solve(crossprod(W))
      dense <- sum(diag((diag(n) - H) %*% W %*% P %*% Sigma %*% P %*% t(W)))
      expect_equal(fmrireg:::.soft_projection_ar_residual_df(sp, Xw, W), dense,
                   tolerance = 1e-8)
    }
  }
  # No projection -> n - rank; no whitening -> the iid formula.
  sp0 <- soft_projection(N, lambda = 1e12)
  W <- fmrireg:::ar_whiten_transform(matrix(0, n, 0), diag(n), 0.5, FALSE)$Y
  expect_equal(fmrireg:::.soft_projection_ar_residual_df(sp0, W %*% X, W), n - 3,
               tolerance = 1e-6)
  sp <- soft_projection(N, lambda = 2)
  Xp <- sp$P_lambda(X)
  expect_equal(fmrireg:::.soft_projection_ar_residual_df(sp, Xp, diag(n)),
               fmrireg:::.soft_projection_residual_df(sp, Xp), tolerance = 1e-10)
})

test_that("volume-weighted df helper equals tr((I-H) P A P)", {
  set.seed(6)
  n <- 70L; k <- 6L
  N <- matrix(rnorm(n * k), n, k)
  X <- cbind(1, rnorm(n))
  active <- rep(TRUE, n); active[c(3, 20, 21, 50)] <- FALSE
  D <- diag(as.numeric(active))
  sp <- soft_projection(N, lambda = 1)
  Xp <- sp$P_lambda(D %*% X)
  P <- sp$P_lambda(diag(n))
  H <- Xp %*% MASS::ginv(crossprod(Xp)) %*% t(Xp)
  dense <- sum(diag((diag(n) - H) %*% P %*% D %*% P))
  expect_equal(fmrireg:::.soft_projection_residual_df(sp, Xp, active = active),
               dense, tolerance = 1e-8)
})

test_that("ACVF correction for the projected operator reduces to fmriAR's design correction", {
  set.seed(7)
  n <- 80L
  X <- cbind(1, rnorm(n), cos(seq_len(n) / 6))
  N <- matrix(rnorm(n * 4), n, 4)
  sp_none <- soft_projection(N, lambda = 1e14)        # P = I
  for (cens in list(NULL, c(10L, 11L, 40L))) {
    ours <- fmrireg:::.soft_projection_acvf_correction(sp_none, X, 6L, censor = cens)
    ref <- fmriAR::acvf_bias_matrix(X, censor = cens, max_lag = 6L)[[1]]
    expect_equal(ours, unname(ref), tolerance = 1e-8, ignore_attr = TRUE)
  }
  # With a real projection: E[residual ACVF] = A gamma, checked against the
  # dense expectation of the residual covariance G Sigma G'.
  sp <- soft_projection(N, lambda = 3)
  G <- sp$P_lambda(diag(n))
  Q <- qr.Q(qr(sp$P_lambda(X)))
  G <- G - Q %*% crossprod(Q, G)
  gam <- 0.5^(0:(n - 1)) / (1 - 0.25)                 # AR(1), phi = 0.5
  Sigma <- toeplitz(gam)
  C <- G %*% Sigma %*% t(G)
  C <- C - outer(rowMeans(C), rep(1, n)) - outer(rep(1, n), colMeans(C)) + mean(C)
  emp <- vapply(0:6, function(h) mean(C[cbind((h + 1):n, 1:(n - h))]), numeric(1))
  A <- fmrireg:::.soft_projection_acvf_correction(sp, sp$P_lambda(X), 6L)
  # A ignores the autocovariance beyond lag 6 (as fmriAR does), which is
  # 0.5^7 / 0.75 ~ 0.01 here.
  expect_equal(as.numeric(A %*% gam[1:7]), emp, tolerance = 0.02)
})

test_that("runwise AR + soft projection: phi and sigma^2 unbiased, t calibrated on pure AR noise", {
  fx <- .p3_soft_fixture()
  fit <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = fx$dset,
                 control = .p3_soft_control(fx$N, noise = noise_spec("ar1")))
  # The residual-ACVF correction for I - H (not (I - H)P) gave 0.271.
  phi <- unlist(fit$result$ar_coef)
  expect_equal(phi, 0.3, tolerance = 0.015 / 0.3)
  s2 <- as.numeric(fit$result$sigma)^2
  # Monte Carlo SE of mean(s2) ~ sqrt(2 / 130 / 3000) ~ 0.002; n - rank gave 0.91.
  expect_equal(mean(s2), 1, tolerance = 0.02)
  tt <- .p3_task_t(fit)
  rate <- mean(abs(tt) > qt(0.975, unique(fit$result$df$inference)))
  # Binomial SE with 6000 tests ~ 0.003; before the fix ~0.08.
  expect_gt(rate, 0.04)
  expect_lt(rate, 0.06)
  # Var(t) ~ 1 (SE ~ 0.02); before ~1.17.
  expect_equal(mean(apply(tt, 2, var)), 1, tolerance = 0.06)

  # The reported df and coefficient covariance are the closed forms at the
  # fitted AR coefficient.
  n <- fx$n
  X <- as.matrix(design_matrix(fit$model))
  sp <- soft_projection(fx$N, lambda = 20)
  P <- sp$P_lambda(diag(n))
  W <- fmrireg:::ar_whiten_transform(matrix(0, n, 0), diag(n), phi, FALSE)$Y
  Sigma <- solve(crossprod(W))
  Xw <- W %*% P %*% X
  Ci <- solve(crossprod(Xw))
  H <- Xw %*% Ci %*% t(Xw)
  M <- W %*% P %*% Sigma %*% P %*% t(W)
  dense_df <- sum(diag((diag(n) - H) %*% M))
  expect_equal(unique(as.numeric(fit$result$df$inference)), dense_df, tolerance = 1e-8)
  expect_lt(dense_df, n - qr(X)$rank - 8)       # the correction is material
  dense_cov <- Ci %*% t(Xw) %*% M %*% Xw %*% Ci
  se <- fit$result$betas$data[[1]]$se[[1]]
  sig <- as.numeric(fit$result$sigma)
  expect_equal(unname(se[1:5, ]), unname(outer(sig[1:5], sqrt(diag(dense_cov)))),
               tolerance = 1e-8)
})

test_that("global AR pooling + soft projection estimates phi from the projected runs", {
  set.seed(34)
  n <- 120L; V <- 1500L; k <- 12L
  ons <- seq(10, 2 * n - 30, by = 18)
  ev <- do.call(rbind, lapply(1:2, function(r) {
    data.frame(onset = ons,
               condition = factor(rep(c("A", "B"), length.out = length(ons))),
               run = r)
  }))
  N <- matrix(rnorm(2 * n * k), 2 * n, k) %*% diag(seq(0.5, 6, length.out = k))
  Y <- matrix(rnorm(2 * n * V), 2 * n, V)
  for (r in 0:1) for (t in 2:n) Y[r * n + t, ] <- 0.3 * Y[r * n + t - 1L, ] + Y[r * n + t, ]
  dset <- matrix_frame(Y, TR = 2, run_length = c(n, n), event_table = ev)
  fit <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = dset,
                 control = .p3_soft_control(N, noise = noise_spec("ar1", pooling = "global")))
  phi <- unique(unlist(fit$result$ar_coef))
  expect_length(phi, 1L)
  expect_equal(phi, 0.3, tolerance = 0.02 / 0.3)
  expect_equal(mean(as.numeric(fit$result$sigma)^2), 1, tolerance = 0.02)
  rate <- mean(abs(.p3_task_t(fit)) > qt(0.975, min(fit$result$df$inference)))
  expect_gt(rate, 0.035)
  expect_lt(rate, 0.065)
})

test_that("runwise iid + volume weights + soft projection: zero-weight volumes leave the df", {
  fx <- .p3_soft_fixture(phi = 0, seed = 32L)
  w <- rep(1, fx$n); w[seq(5, fx$n, by = 9)] <- 0     # 17 censored volumes
  fit <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = fx$dset,
                 control = .p3_soft_control(fx$N, weights = weights_spec(values = w)))
  s2 <- as.numeric(fit$result$sigma)^2
  expect_equal(mean(s2), 1, tolerance = 0.03)
  X <- as.matrix(design_matrix(fit$model))
  sp <- soft_projection(fx$N, lambda = 20)
  Xp <- sp$P_lambda(X * sqrt(w))
  expect_equal(unique(as.numeric(fit$result$df$inference)),
               fmrireg:::.soft_projection_residual_df(sp, Xp, active = w > 0),
               tolerance = 1e-8)
  # Coefficient covariance of the projected fit: C X'PAPX C, not C.
  P <- sp$P_lambda(diag(fx$n))
  Ci <- solve(crossprod(Xp))
  dense_cov <- Ci %*% t(Xp) %*% P %*% diag(as.numeric(w > 0)) %*% P %*% Xp %*% Ci
  se <- fit$result$betas$data[[1]]$se[[1]]
  sig <- as.numeric(fit$result$sigma)
  expect_equal(unname(se[1:5, ]), unname(outer(sig[1:5], sqrt(diag(dense_cov)))),
               tolerance = 1e-8)
  rate <- mean(abs(.p3_task_t(fit)) > qt(0.975, unique(fit$result$df$inference)))
  expect_gt(rate, 0.035)
  expect_lt(rate, 0.065)
})

test_that("soft projection is refused where no exact df exists or the projection is not applied", {
  fx <- .p3_soft_fixture(n = 80L, V = 10L, phi = 0, seed = 33L)
  run <- function(control, ...) {
    fmri_lm(onset ~ hrf(condition), block = ~run, dataset = fx$dset,
            control = control, ...)
  }
  msg <- "soft-subspace projection"
  # Robust weights: the weighted df convention has no counterpart once P
  # mixes down-weighted volumes into the others.
  expect_error(run(.p3_soft_control(fx$N, robust = robust_spec("huber"))), msg)
  expect_error(run(.p3_soft_control(fx$N, noise = noise_spec("ar1"),
                                    robust = robust_spec("huber"))), msg)
  # Joint (chunkwise) estimation never applied the projection to the fit.
  expect_error(run(.p3_soft_control(fx$N, scope = "joint")), msg)
  expect_error(run(.p3_soft_control(fx$N, scope = "joint",
                                    noise = noise_spec("ar1"))), msg)
  # The reference backend projected Y but not X.
  expect_error(run(.p3_soft_control(fx$N),
                   compute = compute_spec(backend = "reference")), msg)
  # MA noise and AR with volume weights: no derived df.
  expect_error(run(.p3_soft_control(fx$N, noise = noise_spec("iid", q = 1L))), msg)
  expect_error(run(.p3_soft_control(fx$N, noise = noise_spec("ar1"),
                                    weights = weights_spec(values = rep(1, 80)))),
               msg)
})
