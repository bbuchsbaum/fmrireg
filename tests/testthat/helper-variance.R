# Substrate for first-level variance validation.
#
# Everything here is a self-contained reference implementation: it does not call
# fmrireg, so it can serve as an independent oracle for the package's variance
# arithmetic. Noise covariances are returned with unit lag-0 so a caller can
# scale them freely.

# ---- noise covariance constructors -------------------------------------------

.vt_ar1 <- function(n_t, rho) rho^abs(outer(seq_len(n_t), seq_len(n_t), "-"))

.vt_arma11 <- function(n_t, phi, theta) {
  h <- seq_len(n_t) - 1L
  g0 <- (1 + 2 * phi * theta + theta^2) / (1 - phi^2)
  g1 <- (1 + phi * theta) * (phi + theta) / (1 - phi^2)
  stats::toeplitz(c(g0, g1 * phi^(h[-1] - 1)) / g0)
}

# 1/f^alpha noise built from its power spectrum, with a white floor
.vt_1f <- function(n_t, alpha = 1) {
  f <- seq(0, 0.5, length.out = n_t)
  S <- 1 / (pmax(f, 1 / (2 * n_t))^alpha) + 0.5
  g <- vapply(seq_len(n_t) - 1L, function(h) sum(S * cos(2 * pi * f * h)), numeric(1))
  stats::toeplitz(g / g[1])
}

# short-memory AR plus a damped oscillation (physio-like); PSD-projected
.vt_osc <- function(n_t) {
  h <- seq_len(n_t) - 1L
  g <- 0.4^h + 0.35 * cos(2 * pi * h / 12) * exp(-h / 40)
  S <- stats::toeplitz(g / g[1])
  e <- eigen(S, symmetric = TRUE)
  e$values <- pmax(e$values, 1e-8)
  S <- e$vectors %*% diag(e$values) %*% t(e$vectors)
  S / S[1, 1]
}

.vt_regimes <- function(n_t) {
  list(
    `AR(1) rho=0.4`    = .vt_ar1(n_t, 0.4),
    `ARMA(1,1)`        = .vt_arma11(n_t, 0.5, 0.4),
    `1/f long memory`  = .vt_1f(n_t),
    `AR+oscillatory`   = .vt_osc(n_t)
  )
}

# ---- design construction -----------------------------------------------------

.vt_hrf <- function(t) {
  ifelse(t < 0, 0, stats::dgamma(t, 6, scale = 1) - (1 / 6) * stats::dgamma(t, 16, scale = 1))
}

.vt_convolve <- function(onsets, n_t, TR = 2) {
  s <- numeric(n_t)
  idx <- pmin(n_t, pmax(1L, round(onsets / TR) + 1L))
  for (i in idx) s[i] <- s[i] + 1
  k <- .vt_hrf(seq(0, 30, by = TR))
  r <- as.numeric(stats::filter(c(s, rep(0, length(k))), k, sides = 1)[seq_len(n_t)])
  r[is.na(r)] <- 0
  r
}

.vt_dct <- function(n_t, cutoff_s, TR = 2) {
  n <- max(1L, floor(2 * n_t * TR / cutoff_s))
  vapply(seq_len(n), function(k) {
    sqrt(2 / n_t) * cos(pi * (2 * seq_len(n_t) - 1) * k / (2 * n_t))
  }, numeric(n_t))
}

# type: "er" (rapid randomized event-related) or "block"
.vt_design <- function(n_t, type = c("er", "block"), n_cond = 2L,
                       cutoff_s = 128, TR = 2, seed = 1L) {
  type <- match.arg(type)
  set.seed(seed)
  if (type == "er") {
    onsets <- sort(stats::runif(20 * n_cond, 0, n_t * TR - 20))
    grp <- rep(seq_len(n_cond), length.out = length(onsets))
    Xt <- vapply(seq_len(n_cond),
                 function(g) .vt_convolve(onsets[grp == g], n_t, TR), numeric(n_t))
  } else {
    starts <- seq(0, n_t * TR - 60, by = 60)
    grp <- rep(seq_len(n_cond), length.out = length(starts))
    Xt <- vapply(seq_len(n_cond), function(g) {
      on <- unlist(lapply(starts[grp == g], function(s) seq(s, s + 28, by = 2)))
      .vt_convolve(on, n_t, TR)
    }, numeric(n_t))
  }
  cbind(Xt, 1, .vt_dct(n_t, cutoff_s, TR))
}

# ---- reduced-space variance oracle -------------------------------------------

# Cov(C beta) = sum_h gamma_h * (L T_h L'), computed per run with no lag product
# spanning a run boundary. gamma_by_run is a named list keyed by run label.
.vt_reduced_var <- function(L, gamma_by_run, run_id) {
  q <- nrow(L)
  acc <- matrix(0, q, q)
  for (r in sort(unique(run_id))) {
    Lr <- L[, run_id == r, drop = FALSE]
    n <- ncol(Lr)
    g <- gamma_by_run[[as.character(r)]]
    H <- min(length(g) - 1L, n - 1L)
    acc <- acc + g[1] * tcrossprod(Lr)
    if (H >= 1L) {
      for (h in seq_len(H)) {
        A <- Lr[, seq_len(n - h), drop = FALSE] %*% t(Lr[, seq(h + 1L, n), drop = FALSE])
        acc <- acc + g[h + 1L] * (A + t(A))
      }
    }
  }
  acc
}

# mean of the entries of Rm at each lag 0..H
.vt_lag_mean <- function(Rm, H, d) {
  vapply(0:H, function(h) mean(Rm[d == h]), numeric(1))
}

# Residual-forming bias map: E[acvf(resid)] = A %*% gamma_true, A from the design
# only. Solve A for a debiased autocovariance. Exact when the true covariance is
# supported within lag H.
.vt_bias_map <- function(X, H) {
  n_t <- nrow(X)
  M <- diag(n_t) - X %*% solve(crossprod(X), t(X))
  d <- abs(outer(seq_len(n_t), seq_len(n_t), "-"))
  vapply(0:H, function(h2) {
    Th <- matrix(0, n_t, n_t)
    Th[d == h2] <- 1
    .vt_lag_mean(M %*% Th %*% M, H, d)
  }, numeric(H + 1L))
}

# ---- analytic references -----------------------------------------------------

# exact SE of a contrast under OLS with known Sigma:  sqrt(c' P Sigma P' c)
.vt_true_se <- function(X, Sigma, cvec) {
  P <- solve(crossprod(X), t(X))
  L <- matrix(cvec, nrow = 1) %*% P
  sqrt(drop(L %*% Sigma %*% t(L)))
}

# classical OLS SE assuming iid errors with variance s2
.vt_iid_se <- function(X, s2, cvec) {
  sqrt(drop(s2 * matrix(cvec, nrow = 1) %*% solve(crossprod(X)) %*% matrix(cvec, ncol = 1)))
}
