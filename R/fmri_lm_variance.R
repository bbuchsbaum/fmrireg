# Variance estimation and centralized inference for fmri_lm result schema v2.

.fmri_lm_segments <- function(n, runs = NULL, censor = NULL) {
  if (is.null(runs)) {
    runs <- list(seq_len(n))
  } else if (!is.list(runs)) {
    runs <- split(seq_len(n), runs)
  }
  censored <- rep(FALSE, n)
  if (!is.null(censor)) {
    if (is.logical(censor)) {
      if (length(censor) != n) stop("Logical censor mask must match observations.", call. = FALSE)
      censored <- censor
    } else {
      ind <- as.integer(censor)
      ind <- ind[is.finite(ind) & ind >= 1L & ind <= n]
      censored[ind] <- TRUE
    }
  }
  out <- list()
  for (idx in runs) {
    idx <- as.integer(idx)
    idx <- idx[idx >= 1L & idx <= n & !censored[idx]]
    if (!length(idx)) next
    groups <- cumsum(c(TRUE, diff(idx) != 1L))
    pieces <- split(idx, groups)
    out <- c(out, pieces[vapply(pieces, length, integer(1)) > 0L])
  }
  out
}

.fmri_lm_hac_weight <- function(h, H, taper) {
  if (h == 0L || identical(taper, "none")) return(1)
  x <- h / (H + 1)
  if (identical(taper, "tukey")) {
    return(0.5 * (1 + cos(pi * x)))
  }
  # Parzen lag window.
  if (x <= 0.5) 1 - 6 * x^2 + 6 * x^3 else 2 * (1 - x)^3
}

.fmri_lm_auto_max_lag <- function(X, runs, censor = NULL,
                                  condition_limit = 50, cap = 60L) {
  n <- nrow(X)
  segments <- .fmri_lm_segments(n, runs, censor)
  longest <- max(vapply(segments, length, integer(1)), 1L)
  cap <- min(as.integer(cap), longest - 1L)
  if (cap <= 0L) return(0L)
  run_id <- integer(n)
  if (is.list(runs)) {
    for (i in seq_along(runs)) run_id[runs[[i]]] <- i
  } else if (is.null(runs)) {
    run_id[] <- 1L
  } else {
    run_id <- as.integer(factor(runs))
  }
  best <- 0L
  for (h in seq_len(cap)) {
    mats <- tryCatch(
      fmriAR::acvf_bias_matrix(X, runs = run_id, censor = censor, max_lag = h),
      error = function(e) NULL
    )
    if (is.null(mats) || !length(mats)) break
    cond <- vapply(mats, function(A) tryCatch(kappa(A), error = function(e) Inf), numeric(1))
    if (all(is.finite(cond)) && max(cond) <= condition_limit) best <- h else break
  }
  best
}

.fmri_lm_psd <- function(M) {
  M <- (M + t(M)) / 2
  ev <- eigen(M, symmetric = TRUE)
  tol <- max(abs(ev$values), 1) * .Machine$double.eps^0.5
  ev$values[ev$values < tol] <- 0
  scaled <- sweep(ev$vectors, 2L, sqrt(ev$values), `*`)
  tcrossprod(scaled)
}

#' Bread of the sandwich estimator, computed once per design
#' @keywords internal
#' @noRd
.fmri_lm_bread <- function(X) {
  XtX <- crossprod(X)
  tryCatch(solve(XtX), error = function(e) qr.solve(XtX))
}

.fmri_lm_score_covariance <- function(X, residual, segments,
                                      max_lag = 0L, taper = "none",
                                      debias = TRUE, bread = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(bread)) bread <- .fmri_lm_bread(X)
  score <- X * as.numeric(residual)
  meat <- matrix(0, p, p)
  for (idx in segments) {
    S <- score[idx, , drop = FALSE]
    meat <- meat + crossprod(S)
    H <- min(as.integer(max_lag), nrow(S) - 1L)
    if (H <= 0L) next
    for (h in seq_len(H)) {
      gamma <- crossprod(S[(h + 1L):nrow(S), , drop = FALSE],
                         S[seq_len(nrow(S) - h), , drop = FALSE])
      weight <- .fmri_lm_hac_weight(h, max_lag, taper)
      meat <- meat + weight * (gamma + t(gamma))
    }
  }
  meat <- .fmri_lm_psd(meat)
  V <- bread %*% meat %*% bread
  if (isTRUE(debias) && n > p) V <- V * n / (n - p)
  .fmri_lm_psd(V)
}

.fmri_lm_residual_covariance <- function(residual, segments, max_lag,
                                         taper = "none", distances = NULL) {
  n <- length(residual)
  V <- matrix(0, n, n)
  for (k in seq_along(segments)) {
    idx <- segments[[k]]
    e <- residual[idx]
    m <- length(e)
    H <- min(as.integer(max_lag), m - 1L)
    gamma <- numeric(H + 1L)
    for (h in 0:H) {
      gamma[h + 1L] <- sum(e[(h + 1L):m] * e[seq_len(m - h)]) / (m - h)
      gamma[h + 1L] <- gamma[h + 1L] * .fmri_lm_hac_weight(h, max_lag, taper)
    }
    distance <- if (!is.null(distances)) {
      distances[[k]]
    } else {
      abs(outer(seq_len(m), seq_len(m), "-"))
    }
    block <- matrix(0, m, m)
    inside <- distance <= H
    block[inside] <- gamma[distance[inside] + 1L]
    V[idx, idx] <- block
  }
  .fmri_lm_psd(V)
}

#' Lag-window autocovariances of each voxel's residuals, per segment
#'
#' Row `offsets[s] + h + 1` of `gamma` holds, for every voxel, the tapered
#' lag-`h` autocovariance of segment `s`: the entries of the block-Toeplitz
#' matrix `.fmri_lm_residual_covariance()` builds.
#' @keywords internal
#' @noRd
.fmri_lm_segment_gammas <- function(E, segments, max_lag, taper) {
  lags <- lapply(segments, function(idx) {
    0:min(as.integer(max_lag), length(idx) - 1L)
  })
  Gamma <- matrix(0, sum(lengths(lags)), ncol(E))
  row <- 0L
  for (s in seq_along(segments)) {
    idx <- segments[[s]]
    m <- length(idx)
    Es <- E[idx, , drop = FALSE]
    for (h in lags[[s]]) {
      row <- row + 1L
      prod <- Es[(h + 1L):m, , drop = FALSE] * Es[seq_len(m - h), , drop = FALSE]
      Gamma[row, ] <- colSums(prod) / (m - h) *
        .fmri_lm_hac_weight(h, max_lag, taper)
    }
  }
  list(gamma = Gamma, lags = lags, offsets = cumsum(c(0L, lengths(lags))))
}

#' Symmetric lag-h shift of the rows of a matrix
#'
#' `T_0 = I`; for `h > 0`, `(T_h M)[i, ] = M[i + h, ] + M[i - h, ]`, terms
#' outside the segment omitted. A segment's banded Toeplitz covariance is
#' `sum_h gamma_h T_h`.
#' @keywords internal
#' @noRd
.fmri_lm_lag_shift <- function(M, h) {
  if (h == 0L) return(M)
  m <- nrow(M)
  out <- matrix(0, m, ncol(M))
  if (h < m) {
    out[seq_len(m - h), ] <- M[(h + 1L):m, , drop = FALSE]
    out[(h + 1L):m, ] <- out[(h + 1L):m, , drop = FALSE] +
      M[seq_len(m - h), , drop = FALSE]
  }
  out
}

#' Model-based Satterthwaite df for every voxel, from lag sums
#'
#' Equals, voxel by voxel, `calculate_effective_df(method = "satterthwaite")`
#' applied to `.fmri_lm_residual_covariance()` (the dense formula used before
#' performance pass 2), without forming an n x n matrix.
#'
#' A voxel's residual covariance is block diagonal over segments, each block
#' the banded Toeplitz matrix `V_s = sum_h g_{s,h} T_h`. Every trace in
#' `tr(RV)^2 / tr(RVRV)` is linear or quadratic in the stacked lag
#' autocovariances `g`: `tr(RV) = l'g` and `tr(RVRV) = g'Mg`, where `l` and
#' `M` (K x K, `K = sum_s (L_s + 1)`) depend only on the design and the
#' segments. They are built once, in O(n p L^2); each voxel then costs O(K^2)
#' (one matrix product over all voxels) instead of an O(n^3) eigen().
#'
#' The dense formula projects V onto the PSD cone (`.fmri_lm_psd()`), and that
#' projection does change V for some voxels: the truncated window
#' (`taper = "none"`) at moderate or long lags, or residuals so small that the
#' projection's absolute tolerance exceeds their variance. Each voxel is
#' therefore certified first (`.fmri_lm_toeplitz_needs_projection()`); voxels
#' that cannot be certified are projected exactly as before, but one segment
#' block at a time (`.fmri_lm_projected_traces()`), so that path costs
#' O(sum_s m_s^3) rather than O(n^3) and still never forms an n x n matrix.
#' @keywords internal
#' @noRd
.fmri_lm_satterthwaite_banded <- function(X, E, segments, max_lag, taper,
                                          XtXinv = NULL) {
  X <- as.matrix(X)
  E <- as.matrix(E)
  nvox <- ncol(E)
  # calculate_effective_df() returns 1 when n <= p.
  if (nrow(X) <= ncol(X)) return(rep(1, nvox))
  if (is.null(XtXinv)) XtXinv <- .satterthwaite_xtxinv(X)
  A <- X %*% XtXinv
  p <- ncol(X)

  gl <- .fmri_lm_segment_gammas(E, segments, max_lag, taper)
  Gamma <- gl$gamma
  lags <- gl$lags
  offsets <- gl$offsets
  K <- nrow(Gamma)

  # Design-only terms. With H = A X' and G = (X'X)^-1:
  #   tr(V)    = sum_s m_s g_{s,0}
  #   tr(HV)   = sum g_{s,h} sum(A_s * T_h X_s)
  #   tr(VV)   = sum_s [m_s g_{s,0}^2 + sum_h 2 (m_s - h) g_{s,h}^2]
  #   tr(HVV)  = sum g_{s,h} g_{s,k} sum((T_h A_s) * (T_k X_s))
  #   tr(HVHV) = sum g g' tr(C_{s,h} C_{s',k}),  C_{s,h} = X_s' T_h X_s G
  l <- numeric(K)
  M <- matrix(0, K, K)
  Cvec <- matrix(0, p * p, K)
  Ctvec <- matrix(0, p * p, K)
  for (s in seq_along(segments)) {
    idx <- segments[[s]]
    m <- length(idx)
    hs <- lags[[s]]
    rows <- offsets[s] + seq_along(hs)
    Xs <- X[idx, , drop = FALSE]
    As <- A[idx, , drop = FALSE]
    TX <- matrix(vapply(hs, function(h) as.numeric(.fmri_lm_lag_shift(Xs, h)),
                        numeric(m * p)), ncol = length(hs))
    TA <- matrix(vapply(hs, function(h) as.numeric(.fmri_lm_lag_shift(As, h)),
                        numeric(m * p)), ncol = length(hs))
    Q <- crossprod(TA, TX)
    Q <- (Q + t(Q)) / 2
    trTT <- ifelse(hs == 0L, m, 2 * (m - hs))
    M[rows, rows] <- diag(trTT, nrow = length(hs)) - 2 * Q
    l[rows] <- ifelse(hs == 0L, m, 0) - as.numeric(crossprod(as.numeric(As), TX))
    for (j in seq_along(hs)) {
      C <- crossprod(Xs, matrix(TX[, j], m, p)) %*% XtXinv
      Cvec[, rows[j]] <- as.numeric(C)
      Ctvec[, rows[j]] <- as.numeric(t(C))
    }
  }
  P <- crossprod(Cvec, Ctvec)
  M <- M + (P + t(P)) / 2

  tr_RV <- as.numeric(crossprod(l, Gamma))
  tr_RVRV <- colSums(Gamma * (M %*% Gamma))

  for (v in which(.fmri_lm_toeplitz_needs_projection(Gamma, lags, offsets))) {
    tr <- .fmri_lm_projected_traces(X, A, XtXinv, Gamma[, v], segments,
                                    lags, offsets)
    tr_RV[v] <- tr[[1L]]
    tr_RVRV[v] <- tr[[2L]]
  }

  ok <- is.finite(tr_RV) & is.finite(tr_RVRV) & tr_RVRV > 0
  out <- rep(1, nvox)
  out[ok] <- pmax(tr_RV[ok]^2 / tr_RVRV[ok], 1)
  out
}

#' Voxels whose PSD projection cannot be certified to be the identity
#'
#' `.fmri_lm_psd()` zeroes eigenvalues below `max(max|lambda|, 1) sqrt(eps)`.
#' Every eigenvalue of a finite section of a symmetric Toeplitz matrix lies in
#' the range of its symbol `f(w) = g_0 + 2 sum_h g_h cos(h w)`. `f` is
#' evaluated on a grid of `[0, pi]` and lowered by its Lipschitz constant
#' (`2 sum h |g_h|`) times the grid half-spacing, a rigorous lower bound;
#' `|g_0| + 2 sum |g_h|` bounds `|lambda|`, hence the tolerance, from above. A
#' voxel is certified when every segment's lower bound exceeds the upper
#' tolerance; otherwise (`TRUE`) it needs the exact projection.
#' @keywords internal
#' @noRd
.fmri_lm_toeplitz_needs_projection <- function(Gamma, lags, offsets) {
  nvox <- ncol(Gamma)
  lower <- rep(Inf, nvox)
  upper <- rep(1, nvox)
  for (s in seq_along(lags)) {
    hs <- lags[[s]]
    G <- Gamma[offsets[s] + seq_along(hs), , drop = FALSE]
    G_lag <- G[-1L, , drop = FALSE]
    upper <- pmax(upper, abs(G[1L, ]) + 2 * colSums(abs(G_lag)))
    if (length(hs) == 1L) {
      lower <- pmin(lower, G[1L, ])
      next
    }
    N <- max(64L, 8L * length(hs))
    w <- pi * (0:N) / N
    Cos <- cos(outer(w, hs)) * rep(c(1, rep(2, length(hs) - 1L)), each = N + 1L)
    f_min <- matrixStats::colMins(Cos %*% G)
    lipschitz <- 2 * colSums(abs(G_lag) * hs[-1L])
    lower <- pmin(lower, f_min - lipschitz * pi / (2 * N))
  }
  !(lower > upper * .Machine$double.eps^0.5)
}

#' tr(RV) and tr(RVRV) with V projected onto the PSD cone, block by block
#'
#' The dense computation (`.fmri_lm_residual_covariance()`, then
#' `.satterthwaite_df()`) for one voxel, using that the eigenvalues of a
#' block-diagonal matrix are those of its blocks: only segment-sized blocks
#' are decomposed, with the tolerance taken over all blocks as before.
#' @keywords internal
#' @noRd
.fmri_lm_projected_traces <- function(X, A, XtXinv, gamma, segments, lags,
                                      offsets) {
  eig <- vector("list", length(segments))
  max_abs <- 0
  for (s in seq_along(segments)) {
    m <- length(segments[[s]])
    g <- gamma[offsets[s] + seq_along(lags[[s]])]
    block <- stats::toeplitz(c(g, rep(0, m - length(g))))
    eig[[s]] <- eigen(block, symmetric = TRUE)
    max_abs <- max(max_abs, abs(eig[[s]]$values))
  }
  tol <- max(max_abs, 1) * .Machine$double.eps^0.5
  tr_V <- tr_VV <- tr_HV <- tr_HVV <- 0
  B <- matrix(0, ncol(X), ncol(X))
  for (s in seq_along(segments)) {
    idx <- segments[[s]]
    ev <- eig[[s]]
    ev$values[ev$values < tol] <- 0
    Vs <- tcrossprod(sweep(ev$vectors, 2L, sqrt(ev$values), `*`))
    Xs <- X[idx, , drop = FALSE]
    As <- A[idx, , drop = FALSE]
    VX <- Vs %*% Xs
    tr_V <- tr_V + sum(diag(Vs))
    tr_VV <- tr_VV + sum(Vs * t(Vs))
    tr_HV <- tr_HV + sum(As * VX)
    tr_HVV <- tr_HVV + sum(As * (Vs %*% VX))
    B <- B + crossprod(Xs, VX)
  }
  BG <- B %*% XtXinv
  list(tr_V - tr_HV, tr_VV - 2 * tr_HVV + sum(BG * t(BG)))
}

.fmri_lm_hac_df <- function(X, residual, segments, max_lag,
                            taper = "none", nominal = NULL) {
  X <- as.matrix(X)
  n_valid <- sum(vapply(segments, length, integer(1)))
  if (is.null(nominal)) {
    nominal <- n_valid - qr(X[unlist(segments), , drop = FALSE])$rank
  }
  if (n_valid <= 1L || max_lag <= 0L) return(max(nominal, 1))

  # A lag-window variance estimate averages fewer effectively independent
  # spectral ordinates as its bandwidth grows. Matching its variance to a
  # scaled chi-square gives the equivalent df below; squared window weights
  # are the spectral-window energy and the finite-segment factors ensure that
  # lags never borrow information across run or censor boundaries.
  energy <- 1
  for (h in seq_len(as.integer(max_lag))) {
    available <- sum(vapply(segments, function(idx) max(length(idx) - h, 0L),
                            integer(1)))
    if (available <= 0L) next
    w <- .fmri_lm_hac_weight(h, max_lag, taper)
    energy <- energy + 2 * w^2 * available / n_valid
  }
  max(min(n_valid / energy, nominal), 1)
}

#' Covariance and df for the joint-fit inference context
#'
#' Everything that depends only on the design (the sandwich bread, the rank,
#' the valid-segment rank and the lag-distance matrices) is computed once,
#' outside the voxel loop. The per-voxel sandwich covariance is computed only
#' when it is used, i.e. when `spec$method` is not `"model"`; model-based
#' Satterthwaite inference needs only the df.
#' @keywords internal
#' @noRd
.fmri_lm_variance_from_context <- function(context, spec) {
  X <- as.matrix(context$X)
  E <- as.matrix(context$residuals)
  if (nrow(X) != nrow(E)) stop("Inference context has incompatible X/residual rows.", call. = FALSE)
  segments <- .fmri_lm_segments(nrow(X), context$runs, context$censor)
  max_lag <- if (identical(spec$max_lag, "auto")) {
    .fmri_lm_auto_max_lag(X, context$runs, context$censor)
  } else {
    as.integer(spec$max_lag)
  }
  if (identical(spec$method, "sandwich")) max_lag <- 0L
  model_based <- identical(spec$method, "model")
  rank_X <- qr(X)$rank

  covariance <- NULL
  if (!model_based) {
    bread <- .fmri_lm_bread(X)
    covariance <- lapply(seq_len(ncol(E)), function(v) {
      .fmri_lm_score_covariance(
        X, E[, v], segments, max_lag = max_lag,
        taper = spec$taper, debias = spec$debias, bread = bread
      )
    })
  }
  df <- rep(nrow(X) - rank_X, ncol(E))
  if (identical(spec$df, "satterthwaite")) {
    if (model_based) {
      # Banded/Toeplitz algebra; no per-voxel n x n covariance or eigen().
      df <- .fmri_lm_satterthwaite_banded(
        X, E, segments, max_lag = max_lag, taper = spec$taper,
        XtXinv = .satterthwaite_xtxinv(X)
      )
    } else {
      nominal <- sum(vapply(segments, length, integer(1))) -
        qr(X[unlist(segments), , drop = FALSE])$rank
      df <- vapply(seq_len(ncol(E)), function(v) {
        .fmri_lm_hac_df(
          X, E[, v], segments, max_lag = max_lag,
          taper = spec$taper, nominal = nominal
        )
      }, numeric(1))
    }
  }
  if (!is.null(context$robust_weights)) {
    # The single robust formula (.robust_effective_df): scale the df by
    # (sum(w) - rank) / (n - rank).
    df <- vapply(df, function(d) {
      .robust_effective_df(context$robust_weights, rank_X, df_base = d)
    }, numeric(1))
  }
  list(covariance = covariance, df = df, max_lag = max_lag,
       segments = segments)
}

.fmri_lm_update_inference <- function(result, covariance, df) {
  payload <- .fmri_lm_beta_payload(result)
  if (is.null(payload) || !length(payload$estimate)) return(result)
  estimate <- as.matrix(payload$estimate[[1L]])
  V <- nrow(estimate)
  p <- ncol(estimate)
  if (length(covariance) != V) stop("One covariance matrix is required per voxel.", call. = FALSE)
  if (length(df) == 1L) df <- rep(df, V)

  se <- stat <- prob <- matrix(NA_real_, V, p)
  for (v in seq_len(V)) {
    cv <- covariance[[v]]
    if (is.null(cv) || any(dim(cv) != c(p, p))) next
    se[v, ] <- sqrt(pmax(diag(cv), 0))
    valid <- is.finite(se[v, ]) & se[v, ] > .Machine$double.eps^0.5
    stat[v, valid] <- estimate[v, valid] / se[v, valid]
    prob[v, valid] <- 2 * stats::pt(-abs(stat[v, valid]), df = df[v])
  }
  payload$se[[1L]] <- se
  payload$stat[[1L]] <- stat
  payload$prob[[1L]] <- prob
  result$betas$data[[1L]] <- payload

  if (!is.null(result$contrasts) && nrow(result$contrasts)) {
    for (i in seq_len(nrow(result$contrasts))) {
      row <- result$contrasts[i, , drop = FALSE]
      type <- row$type[[1L]]
      data <- row$data[[1L]]
      colind <- row$colind[[1L]]
      C <- row$conmat[[1L]]
      if (identical(type, "contrast")) {
        full <- numeric(p)
        full[colind] <- as.numeric(C)
        est <- drop(estimate %*% full)
        cse <- vapply(covariance, function(cv) {
          sqrt(pmax(drop(crossprod(full, cv %*% full)), 0))
        }, numeric(1))
        cstat <- ifelse(cse > .Machine$double.eps^0.5, est / cse, NA_real_)
        data$estimate <- est
        data$se <- cse
        data$stat <- cstat
        data$prob <- 2 * stats::pt(-abs(cstat), df = df)
      } else if (identical(type, "Fcontrast")) {
        C <- .orient_fcontrast(C, colind)
        full <- matrix(0, nrow(C), p)
        full[, colind] <- C
        rank_c <- qr(full)$rank
        fstat <- rep(NA_real_, V)
        for (v in seq_len(V)) {
          u <- full %*% estimate[v, ]
          M <- full %*% covariance[[v]] %*% t(full)
          fstat[v] <- tryCatch(drop(crossprod(u, solve(M, u))) / rank_c,
                              error = function(e) NA_real_)
        }
        data$estimate <- fstat
        data$se <- rep(1, V)
        data$stat <- fstat
        data$prob <- stats::pf(fstat, rank_c, df, lower.tail = FALSE)
      }
      result$contrasts$data[[i]] <- data
    }
  }
  result
}
