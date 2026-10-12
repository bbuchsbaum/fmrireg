# Model-based Satterthwaite df with the PSD projection: dense per-block
# eigen() (old) vs banded eigensolver + low-rank trace correction (new).
#
# Run from the package root with one thread:
#   OMP_NUM_THREADS=1 Rscript bench/satterthwaite-psd.R
#
# "old" is .fmri_lm_satterthwaite_banded() with the projected path replaced
# by the per-voxel dense per-block helper (.fmri_lm_projected_traces_dense(),
# the code as it stood before this change); everything else is shared, so
# the difference is the projected path alone. Reported: ms per voxel over all
# voxels, ms per projected voxel, and the max relative df difference.

# load_all() compiles with -O0 (pkgbuild debug flags); build the shared
# library with the installed-package flags first so the C++ kernel is timed
# as users run it.
pkgbuild::compile_dll(".", debug = FALSE, force = TRUE, quiet = TRUE)
suppressMessages(devtools::load_all(".", quiet = TRUE, compile = FALSE))
options(fmrireg.num_threads = 1L)
ns <- asNamespace("fmrireg")

old_env <- new.env(parent = ns)
old_env$.fmri_lm_projected_traces <- function(X, A, XtXinv, Gamma, segments,
                                              lags, offsets, tr_RV, tr_RVRV) {
  for (v in seq_len(ncol(Gamma))) {
    tr <- ns$.fmri_lm_projected_traces_dense(X, A, XtXinv, Gamma[, v],
                                             segments, lags, offsets)
    tr_RV[v] <- tr[[1L]]
    tr_RVRV[v] <- tr[[2L]]
  }
  list(tr_RV = tr_RV, tr_RVRV = tr_RVRV)
}
old_df <- ns$.fmri_lm_satterthwaite_banded
environment(old_df) <- old_env
new_df <- ns$.fmri_lm_satterthwaite_banded

elapsed <- function(expr) {
  gc(FALSE)
  system.time(expr)[["elapsed"]]
}

set.seed(1)
n <- 400L
p <- 12L
X <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
Vmax <- 2000L
E_all <- vapply(seq_len(Vmax), function(v) {
  as.numeric(arima.sim(list(ar = 0.4), n))
}, numeric(n))
E_all <- E_all - X %*% qr.solve(X, E_all)
XtXinv <- ns$.satterthwaite_xtxinv(X)

rows <- list()
for (layout in c("one run", "two runs")) {
  runs <- if (layout == "one run") list(seq_len(n)) else list(1:200, 201:400)
  segs <- ns$.fmri_lm_segments(n, runs)
  for (V in c(60L, Vmax)) {
    E <- E_all[, seq_len(V), drop = FALSE]
    for (cfg in list(c("none", 8), c("none", 30), c("none", 60),
                     c("parzen", 8))) {
      taper <- cfg[1]
      L <- as.integer(cfg[2])
      gl <- ns$.fmri_lm_segment_gammas(E, segs, L, taper)
      nproj <- sum(ns$.fmri_lm_toeplitz_needs_projection(gl$gamma, gl$lags,
                                                         gl$offsets))
      t_old <- elapsed(ref <- old_df(X, E, segs, L, taper, XtXinv = XtXinv))
      t_new <- elapsed(got <- new_df(X, E, segs, L, taper, XtXinv = XtXinv))
      row <- data.frame(
        layout = layout, V = V, taper = taper, lag = L, projected = nproj,
        old_ms_vox = round(1000 * t_old / V, 3),
        new_ms_vox = round(1000 * t_new / V, 3),
        old_ms_proj = if (nproj) round(1000 * t_old / nproj, 2) else NA,
        new_ms_proj = if (nproj) round(1000 * t_new / nproj, 2) else NA,
        speedup = round(t_old / t_new, 1),
        max_rel_diff = signif(max(abs(got - ref) / abs(ref)), 2)
      )
      print(row, row.names = FALSE)
      rows[[length(rows) + 1L]] <- row
    }
  }
}
cat("\n")
print(do.call(rbind, rows), row.names = FALSE)
