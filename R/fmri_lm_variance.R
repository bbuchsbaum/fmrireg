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

.fmri_lm_score_covariance <- function(X, residual, segments,
                                      max_lag = 0L, taper = "none",
                                      debias = TRUE) {
  n <- nrow(X)
  p <- ncol(X)
  bread <- tryCatch(solve(crossprod(X)), error = function(e) qr.solve(crossprod(X)))
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
                                         taper = "none") {
  n <- length(residual)
  V <- matrix(0, n, n)
  for (idx in segments) {
    e <- residual[idx]
    m <- length(e)
    H <- min(as.integer(max_lag), m - 1L)
    gamma <- numeric(H + 1L)
    for (h in 0:H) {
      gamma[h + 1L] <- sum(e[(h + 1L):m] * e[seq_len(m - h)]) / (m - h)
      gamma[h + 1L] <- gamma[h + 1L] * .fmri_lm_hac_weight(h, max_lag, taper)
    }
    distance <- abs(row(diag(m)) - col(diag(m)))
    block <- matrix(0, m, m)
    inside <- distance <= H
    block[inside] <- gamma[distance[inside] + 1L]
    V[idx, idx] <- block
  }
  .fmri_lm_psd(V)
}

.fmri_lm_hac_df <- function(X, residual, segments, max_lag,
                            taper = "none") {
  X <- as.matrix(X)
  n_valid <- sum(vapply(segments, length, integer(1)))
  nominal <- n_valid - qr(X[unlist(segments), , drop = FALSE])$rank
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

  covariance <- lapply(seq_len(ncol(E)), function(v) {
    .fmri_lm_score_covariance(
      X, E[, v], segments, max_lag = max_lag,
      taper = spec$taper, debias = spec$debias
    )
  })
  df <- rep(nrow(X) - qr(X)$rank, ncol(E))
  if (identical(spec$df, "satterthwaite")) {
    df <- vapply(seq_len(ncol(E)), function(v) {
      if (identical(spec$method, "model")) {
        resid_cov <- .fmri_lm_residual_covariance(
          E[, v], segments, max_lag = max_lag, taper = spec$taper
        )
        calculate_effective_df(
          nrow(X), ncol(X), method = "satterthwaite",
          X = X, resid_cov = resid_cov
        )
      } else {
        .fmri_lm_hac_df(
          X, E[, v], segments, max_lag = max_lag,
          taper = spec$taper
        )
      }
    }, numeric(1))
  }
  if (!is.null(context$robust_weights)) {
    w <- context$robust_weights[is.finite(context$robust_weights)]
    adjustment <- if (length(w)) sum(w) / length(w) else 0
    df <- pmax(df * adjustment, 1)
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
