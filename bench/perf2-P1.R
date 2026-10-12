# Performance pass 2, stream P1 (fitting core): old vs new, back to back.
#
# Run from the package root with one thread:
#   OMP_NUM_THREADS=1 Rscript bench/perf2-P1.R
#
# Each "old" implementation is copied below (the code as it stood before
# pass 2), so both run in the same R process on the same data. Every
# comparison also checks that old and new agree.

suppressMessages(devtools::load_all(".", quiet = TRUE))
options(fmrireg.num_threads = 1L)
if (requireNamespace("RcppParallel", quietly = TRUE)) {
  RcppParallel::setThreadOptions(numThreads = 1L)
}
ns <- asNamespace("fmrireg")

timeit <- function(expr, reps = 3L) {
  f <- substitute(expr)
  env <- parent.frame()
  ts <- vapply(seq_len(reps), function(i) {
    gc(FALSE)
    system.time(eval(f, env))[["elapsed"]]
  }, numeric(1))
  median(ts)
}
rows <- list()
record <- function(item, case, old, new, agree) {
  rows[[length(rows) + 1L]] <<- data.frame(
    item = item, case = case, old_s = round(old, 3), new_s = round(new, 3),
    speedup = round(old / new, 1), max_rel_diff = signif(agree, 2)
  )
  print(rows[[length(rows)]], row.names = FALSE)
}
rel <- function(a, b) max(abs(a - b) / pmax(1, abs(b)), na.rm = TRUE)

# ---------------------------------------------------------------------------
# P1-1 Satterthwaite df: dense n x n covariance + eigen() per voxel (old)
# ---------------------------------------------------------------------------
old_psd <- function(M) {
  M <- (M + t(M)) / 2
  ev <- eigen(M, symmetric = TRUE)
  tol <- max(abs(ev$values), 1) * .Machine$double.eps^0.5
  ev$values[ev$values < tol] <- 0
  tcrossprod(sweep(ev$vectors, 2L, sqrt(ev$values), `*`))
}
old_satt_df <- function(X, E, segments, max_lag, taper) {
  XtXinv <- ns$.satterthwaite_xtxinv(X)
  distances <- lapply(segments, function(idx) {
    m <- length(idx); abs(outer(seq_len(m), seq_len(m), "-"))
  })
  vapply(seq_len(ncol(E)), function(v) {
    n <- nrow(E); V <- matrix(0, n, n)
    for (k in seq_along(segments)) {
      idx <- segments[[k]]; e <- E[idx, v]; m <- length(e)
      H <- min(as.integer(max_lag), m - 1L)
      g <- vapply(0:H, function(h) {
        sum(e[(h + 1L):m] * e[seq_len(m - h)]) / (m - h) *
          ns$.fmri_lm_hac_weight(h, max_lag, taper)
      }, numeric(1))
      D <- distances[[k]]; B <- matrix(0, m, m); inside <- D <= H
      B[inside] <- g[D[inside] + 1L]; V[idx, idx] <- B
    }
    V <- old_psd(V)
    ns$calculate_effective_df(nrow(X), ncol(X), method = "satterthwaite",
                              X = X, resid_cov = V, XtXinv = XtXinv)
  }, numeric(1))
}
set.seed(1)
n <- 400L; p <- 12L; V <- 60L
X <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
E <- vapply(seq_len(V), function(v) as.numeric(arima.sim(list(ar = 0.4), n)), numeric(n))
E <- E - X %*% qr.solve(X, E)
runs <- list(1:200, 201:400)
segs <- ns$.fmri_lm_segments(n, runs)
for (cfg in list(c("parzen", 8), c("none", 8), c("none", 30))) {
  taper <- cfg[1]; L <- as.integer(cfg[2])
  spec <- variance_spec("model", max_lag = L, taper = taper, df = "satterthwaite")
  ctx <- list(X = X, residuals = E, runs = runs, censor = NULL)
  ref <- old_satt_df(X, E, segs, L, taper)
  got <- ns$.fmri_lm_variance_from_context(ctx, spec)$df
  t_old <- timeit(old_satt_df(X, E, segs, L, taper), reps = 1L)
  t_new <- timeit(ns$.fmri_lm_variance_from_context(ctx, spec))
  nproj <- sum(ns$.fmri_lm_toeplitz_needs_projection(
    ns$.fmri_lm_segment_gammas(E, segs, L, taper)$gamma,
    ns$.fmri_lm_segment_gammas(E, segs, L, taper)$lags,
    ns$.fmri_lm_segment_gammas(E, segs, L, taper)$offsets))
  record("P1-1", sprintf("n=%d p=%d V=%d %s lag=%d (%d/%d voxels need PSD)",
                         n, p, V, taper, L, nproj, V),
         t_old, t_new, rel(got, ref))
}

# ---------------------------------------------------------------------------
# P1-2 voxelwise AR statistics: per-voxel .fast_t/.fast_F calls (old)
# ---------------------------------------------------------------------------
old_beta_stats <- function(Betas, XtXinv_list, sigma, dfres) {
  V <- ncol(Betas); p <- nrow(Betas)
  se <- tt <- pr <- matrix(NA_real_, V, p)
  for (v in seq_len(V)) {
    M <- XtXinv_list[[v]]
    if (is.null(M) || !is.finite(sigma[v])) next
    s <- sqrt(diag(M)) * sigma[v]
    t_v <- ifelse(abs(s) < .Machine$double.eps^0.5, 0, Betas[, v] / s)
    se[v, ] <- s; tt[v, ] <- t_v; pr[v, ] <- 2 * pt(-abs(t_v), dfres)
  }
  list(se = se, stat = tt, prob = pr)
}
old_contrasts <- function(Betas, sigma2, XtXinv_list, conlist, fconlist, dfres) {
  p <- nrow(Betas); V <- ncol(Betas); out <- list()
  for (nm in names(conlist)) {
    l <- matrix(0, 1, p); l[, attr(conlist[[nm]], "colind")] <- conlist[[nm]]
    out[[nm]] <- vapply(seq_len(V), function(v) {
      ns$.fast_t_contrast(Betas[, v, drop = FALSE], sigma2[v], XtXinv_list[[v]],
                          l, dfres, NULL, 1L)$stat
    }, numeric(1))
  }
  for (nm in names(fconlist)) {
    L <- fconlist[[nm]]; ci <- attr(L, "colind")
    Cw <- ns$.orient_fcontrast(L, ci); full <- matrix(0, nrow(Cw), p); full[, ci] <- Cw
    out[[nm]] <- vapply(seq_len(V), function(v) {
      ns$.fast_F_contrast(Betas[, v, drop = FALSE], sigma2[v], XtXinv_list[[v]],
                          full, dfres, NULL, 1L)$stat
    }, numeric(1))
  }
  out
}
set.seed(2)
V <- 5000L; p <- 10L
Betas <- matrix(rnorm(p * V), p, V)
base <- solve(crossprod(matrix(rnorm(200 * p), 200, p)))
XtXinv_list <- lapply(seq_len(V), function(v) base * (1 + v / V))
sigma <- runif(V, 0.5, 2)
tcon <- list(a = structure(1, colind = 2L), ab = structure(c(1, -1), colind = 2:3),
             c = structure(c(1, 1, -2), colind = 4:6))
fcon <- list(f2 = structure(diag(2), colind = 2:3), f3 = structure(diag(3), colind = 4:6))
ref_b <- old_beta_stats(Betas, XtXinv_list, sigma, 190)
got_b <- ns$beta_stats_matrix_voxelwise(Betas, XtXinv_list, sigma, 190,
                                        paste0("x", 1:p), ar_order = 1L)$data[[1]]
record("P1-2", sprintf("beta stats, V=%d p=%d", V, p),
       timeit(old_beta_stats(Betas, XtXinv_list, sigma, 190)),
       timeit(ns$beta_stats_matrix_voxelwise(Betas, XtXinv_list, sigma, 190,
                                             paste0("x", 1:p), ar_order = 1L)),
       rel(got_b$stat[[1]], ref_b$stat))
ref_c <- old_contrasts(Betas, sigma^2, XtXinv_list, tcon, fcon, 190)
got_c <- ns$fit_lm_contrasts_voxelwise(Betas, sigma^2, XtXinv_list, tcon, fcon, 190,
                                       ar_order = 1L)
record("P1-2", sprintf("3 t + 2 F contrasts, V=%d", V),
       timeit(old_contrasts(Betas, sigma^2, XtXinv_list, tcon, fcon, 190), reps = 1L),
       timeit(ns$fit_lm_contrasts_voxelwise(Betas, sigma^2, XtXinv_list, tcon, fcon,
                                            190, ar_order = 1L)),
       max(vapply(names(ref_c), function(nm) rel(got_c[[nm]]$data[[1]]$stat, ref_c[[nm]]),
                  numeric(1))))
# Data shipped to future workers by the parallel voxel path: the old closure
# serialised runwise_lm_voxelwise()'s frame, including the whole run matrix,
# with every worker's chunk; the new path ships one series per voxel.
Y_run <- matrix(rnorm(300 * 2000), 300, 2000)
old_fun <- local({ Y_run <- Y_run; function(v) Y_run[, v] })
cat(sprintf("P1-2 future payload per worker chunk: old >= %.1f MB (closure env), new %.3f MB per voxel\n",
            length(serialize(old_fun, NULL)) / 2^20,
            length(serialize(Y_run[, 1], NULL)) / 2^20))

# ---------------------------------------------------------------------------
# P1-3 reference backend: AR estimated per chunk (old) vs once (new)
# ---------------------------------------------------------------------------
set.seed(3)
n_run <- 150L; V <- 2000L; nchunks <- 10L
ev <- data.frame(onset = rep(seq(5, 2 * n_run - 30, by = 15), 2),
                 cond = factor(rep(c("A", "B"), length.out = 2 * length(seq(5, 2 * n_run - 30, by = 15)))),
                 run = rep(1:2, each = length(seq(5, 2 * n_run - 30, by = 15))))
Y <- matrix(as.numeric(arima.sim(list(ar = 0.4), 2 * n_run * V)), 2 * n_run, V) + 10
dset <- matrix_frame(Y, TR = 2, run_length = c(n_run, n_run), event_table = ev)
for (cs in list(list(name = "AR1", noise = noise_spec("ar1"), robust = robust_spec("none")),
                list(name = "AR1 + huber", noise = noise_spec("ar1"),
                     robust = robust_spec("huber", scale_scope = "run")))) {
  ctrl <- fmri_lm_control(estimation = estimation_spec("joint"), noise = cs$noise,
                          robust = cs$robust)
  # The chunk loop of chunkwise_lm_slow(): old = solve_integrated_glm() per
  # chunk (re-estimating the noise model), new = estimate once, fixed per chunk.
  model <- fmri_lm(onset ~ hrf(cond), block = ~run, dataset = dset, durations = 0,
                   control = ctrl,
                   compute = compute_spec(voxel_chunks = 1L, backend = "reference"))$model
  modmat <- as.matrix(design_matrix(model))
  cfg <- ctrl
  run_idx <- ns$.dset_run_rows(dset)
  chunks <- ns$.dset_chunks(dset, nchunks)
  old_path <- function() {
    lapply(chunks, function(ch) {
      ns$solve_integrated_glm(modmat, as.matrix(ns$.chunk_data(ch)), cfg, run_idx)$betas
    })
  }
  new_path <- function() {
    fixed <- ns$.integrated_noise_estimate(
      modmat, cfg, run_idx, ncol(dset),
      function(cols) as.matrix(ns$.dset_read(dset, cols = cols)))
    lapply(chunks, function(ch) {
      ns$solve_integrated_glm(modmat, as.matrix(ns$.chunk_data(ch)), cfg, run_idx,
                              fixed = fixed)$betas
    })
  }
  one <- ns$solve_integrated_glm(modmat, Y, cfg, run_idx)$betas
  b_old <- do.call(cbind, old_path())
  b_new <- do.call(cbind, new_path())
  cat(sprintf("P1-3 %s: max |beta(nchunks=%d) - beta(nchunks=1)|: old %.2e, new %.2e\n",
              cs$name, nchunks, max(abs(b_old - one)), max(abs(b_new - one))))
  record("P1-3", sprintf("%s, T=%d V=%d nchunks=%d", cs$name, 2 * n_run, V, nchunks),
         timeit(old_path()), timeit(new_path()), rel(b_new, one))
}

# ---------------------------------------------------------------------------
# P1-4 robust IRLS: weighted response and sweep() copies per iteration (old)
# ---------------------------------------------------------------------------
src <- readLines("tests/testthat/test-perf2-P1-robust-irls.R")
start <- grep("^\\.p14_old_fitter <- function", src)
end <- grep("^}", src)
end <- end[end > start][1]
eval(parse(text = src[start:end]))
environment(.p14_old_fitter) <- ns
set.seed(4)
n <- 400L; V <- 5000L
X <- cbind(1, matrix(rnorm(n * 9), n, 9))
Y <- X %*% matrix(rnorm(10 * V), 10, V) + matrix(rt(n * V, df = 3), n, V)
# Volume-wide artefacts (the robust weights are per volume): 5% of volumes
# carry a global spike, so the weights move across iterations.
spikes <- sample(n, n %/% 20)
Y[spikes, ] <- Y[spikes, ] + rnorm(length(spikes), 6, 2)
ctx <- ns$glm_context(X = X, Y = Y, proj = ns$.fast_preproject(X))
for (scope in c("run", "voxel")) {
  cfg <- list(type = "huber", k_huber = 1.345, c_tukey = 4.685, max_iter = 5L,
              tol = 0, scale_scope = scope)
  ref <- .p14_old_fitter(ctx, cfg, X)
  got <- ns$robust_iterative_fitter(ctx, cfg, X)
  record("P1-4", sprintf("huber %s scope, T=%d V=%d, 5 iterations", scope, n, V),
         timeit(.p14_old_fitter(ctx, cfg, X)), timeit(ns$robust_iterative_fitter(ctx, cfg, X)),
         rel(got$betas_robust, ref$betas_robust))
}

# ---------------------------------------------------------------------------
# P1-5 whitening: per-voxel legacy calls; whiten_glm_context with fixed phi
# ---------------------------------------------------------------------------
old_whiten_transform <- function(X, Y, phi, exact_first = FALSE, censor = NULL,
                                 run_indices = NULL) {
  if (anyNA(X) || anyNA(Y)) stop("NA values")
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  phi_list <- if (is.list(phi)) phi else list(phi)
  theta_list <- replicate(length(phi_list), numeric(0), simplify = FALSE)
  pooling_mode <- if (length(phi_list) > 1L) "run" else "global"
  phi_arg <- if (identical(pooling_mode, "global")) phi_list[[1]] else phi_list
  runs <- ns$.build_run_labels(nrow(X), run_indices)
  plan <- fmriAR::compat$plan_from_phi(phi = phi_arg, theta = theta_list,
                                       runs = runs, pooling = pooling_mode,
                                       exact_first = exact_first)
  if (is.null(plan$theta) || length(plan$theta) == 0L) plan$theta <- theta_list
  result <- fmriAR::whiten_apply(plan = plan, X = X, Y = Y, runs = runs,
                                 parallel = TRUE, censor = censor)
  list(X = result$X, Y = result$Y)
}
set.seed(5)
n <- 300L
X <- cbind(1, matrix(rnorm(n * 9), n, 9))
ys <- lapply(1:2000, function(i) matrix(rnorm(n)))
ref <- old_whiten_transform(X, ys[[1]], 0.4)
got <- ns$ar_whiten_transform(X, ys[[1]], 0.4)
record("P1-5", "ar_whiten_transform, 2000 single-voxel calls",
       timeit(for (y in ys) old_whiten_transform(X, y, 0.4)),
       timeit(for (y in ys) ns$ar_whiten_transform(X, y, 0.4)),
       rel(got$Y, ref$Y))
old_whiten_ctx <- function(ctx, phi) {
  coef <- ctx$proj$Pinv %*% ctx$Y
  residuals <- ctx$Y - ctx$X %*% coef
  ns$whiten_glm_context(ctx, list(cor_struct = "ar1", phi = phi))
}
Yb <- matrix(rnorm(n * 20000), n, 20000)
ctx <- ns$glm_context(X = X, Y = Yb, proj = ns$.fast_preproject(X))
record("P1-5", "whiten_glm_context fixed phi, T=300 V=20000",
       timeit(old_whiten_ctx(ctx, 0.4)),
       timeit(ns$whiten_glm_context(ctx, list(cor_struct = "ar1", phi = 0.4))),
       0)

cat("\n")
print(do.call(rbind, rows), row.names = FALSE)
