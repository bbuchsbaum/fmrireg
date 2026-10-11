# Review 2026-10, stream D: runwise and chunkwise strategies estimate the same
# model. On two-run data the joint (chunkwise) fit equals a hand-coded GLS on
# the stacked data, the runwise fit equals hand-coded per-run GLS pooled by
# fixed effects, and the two strategies agree closely with each other.

.cs_fixture <- function(n = 150L, V = 5L, seed = 12L) {
  set.seed(seed)
  ons <- seq(8, 2 * n - 30, by = 12)
  ev <- do.call(rbind, lapply(1:2, function(r) {
    data.frame(onset = ons,
               condition = factor(rep(c("A", "B"), length.out = length(ons))),
               run = r)
  }))
  noise <- matrix(0, 2 * n, V)
  for (r in 1:2) {
    for (v in seq_len(V)) {
      noise[(r - 1L) * n + seq_len(n), v] <- as.numeric(arima.sim(list(ar = 0.4), n)) +
        40 + 3 * r
    }
  }
  fr0 <- matrix_frame(noise, TR = 2, run_length = rep(n, 2), event_table = ev)
  mod <- create_fmri_model(onset ~ hrf(condition), block = ~run,
                           dataset = fr0, durations = 0)
  X <- as.matrix(design_matrix(mod))
  ec <- unlist(attr(mod$event_model$design_matrix, "col_indices"))
  Y <- noise + X[, ec] %*% ((2 / max(abs(X[, ec]))) * matrix(c(1, -1), 2, V))
  list(dset = matrix_frame(Y, TR = 2, run_length = rep(n, 2), event_table = ev),
       X = X, Y = Y, ec = ec, n = n)
}

.cs_fit <- function(fx, scope, noise, chunks = 1L) {
  fmri_lm(onset ~ hrf(condition), block = ~run, dataset = fx$dset, durations = 0,
          control = fmri_lm_control(estimation = estimation_spec(scope),
                                    noise = noise_spec(noise)),
          compute = compute_spec(voxel_chunks = chunks))
}

.cs_get <- function(fit, what) as.matrix(fit$result$betas$data[[1]][[what]][[1]])

# AR(1) whitening with the conditional first row (exact_first = FALSE).
.cs_whitener <- function(n, phi) {
  W <- diag(n)
  if (phi != 0) W[cbind(2:n, 1:(n - 1))] <- -phi
  W
}

.cs_rows <- function(fx, r) (r - 1L) * fx$n + seq_len(fx$n)

.cs_run_design <- function(fx, r) {
  Xr <- fx$X[.cs_rows(fx, r), , drop = FALSE]
  Xr[, colSums(abs(Xr)) > 0, drop = FALSE]
}

test_that("joint fits equal hand-coded stacked GLS (OLS and AR1)", {
  fx <- .cs_fixture()
  for (noise in c("iid", "ar1")) {
    fit <- .cs_fit(fx, "joint", noise, chunks = 2L)
    phi <- if (noise == "iid") c(0, 0) else unlist(fit$result$ar_coef)
    expect_length(phi, 2L)
    W <- matrix(0, 2 * fx$n, 2 * fx$n)
    for (r in 1:2) W[.cs_rows(fx, r), .cs_rows(fx, r)] <- .cs_whitener(fx$n, phi[r])
    Xw <- W %*% fx$X
    for (v in seq_len(ncol(fx$Y))) {
      ref <- summary(lm(drop(W %*% fx$Y[, v]) ~ Xw - 1))
      expect_equal(unname(.cs_get(fit, "estimate")[v, ]),
                   unname(ref$coefficients[, 1]), tolerance = 1e-8, info = noise)
      expect_equal(unname(.cs_get(fit, "se")[v, ]),
                   unname(ref$coefficients[, 2]), tolerance = 1e-8, info = noise)
      expect_equal(unname(.cs_get(fit, "stat")[v, ]),
                   unname(ref$coefficients[, 3]), tolerance = 1e-8, info = noise)
    }
  }
})

test_that("runwise fits equal hand-coded per-run GLS pooled by fixed effects", {
  fx <- .cs_fixture()
  for (noise in c("iid", "ar1")) {
    fit <- .cs_fit(fx, "runwise_meta", noise)
    phi <- if (noise == "iid") c(0, 0) else unlist(fit$result$ar_coef)
    for (v in seq_len(ncol(fx$Y))) {
      per_run <- lapply(1:2, function(r) {
        W <- .cs_whitener(fx$n, phi[r])
        Xr <- W %*% .cs_run_design(fx, r)
        summary(lm(drop(W %*% fx$Y[.cs_rows(fx, r), v]) ~ Xr - 1))$coefficients[fx$ec, 1:2]
      })
      b <- sapply(per_run, function(m) m[, 1])
      prec <- sapply(per_run, function(m) 1 / m[, 2]^2)
      expect_equal(unname(.cs_get(fit, "estimate")[v, fx$ec]),
                   unname(rowSums(b * prec) / rowSums(prec)),
                   tolerance = 1e-8, info = noise)
      expect_equal(unname(.cs_get(fit, "se")[v, fx$ec]),
                   unname(sqrt(1 / rowSums(prec))), tolerance = 1e-8, info = noise)
    }
  }
})

test_that("runwise and chunkwise agree on two-run data", {
  fx <- .cs_fixture()
  for (noise in c("iid", "ar1")) {
    rw <- .cs_fit(fx, "runwise_meta", noise)
    cw <- .cs_fit(fx, "joint", noise, chunks = 3L)
    if (noise == "ar1") {
      # Both estimate the AR model run by run from the same OLS residuals.
      expect_equal(unlist(rw$result$ar_coef), unlist(cw$result$ar_coef),
                   tolerance = 1e-10)
    }
    # Same effects, standard errors and test statistics up to the difference
    # between pooling run estimates and fitting the runs jointly.
    expect_equal(unname(.cs_get(rw, "estimate")[, fx$ec]),
                 unname(.cs_get(cw, "estimate")[, fx$ec]), tolerance = 0.02,
                 info = noise)
    expect_equal(unname(.cs_get(rw, "se")[, fx$ec]),
                 unname(.cs_get(cw, "se")[, fx$ec]), tolerance = 0.03, info = noise)
    expect_equal(unname(.cs_get(rw, "stat")[, fx$ec]),
                 unname(.cs_get(cw, "stat")[, fx$ec]), tolerance = 0.03, info = noise)
  }
})
