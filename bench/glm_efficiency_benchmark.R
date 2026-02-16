#!/usr/bin/env Rscript

# Reproducible GLM efficiency benchmark for PRD workstreams.
# Measures wall time for runwise/chunkwise AR+robust fitting and compares
# default fixed-order AR behavior against legacy fallback mode.

suppressPackageStartupMessages({
  library(pkgload)
  library(fmridataset)
  library(fmrihrf)
})

load_fmrireg_for_bench <- function() {
  # Prefer benchmarking the local working tree to avoid stale installed-package timings.
  if (file.exists("DESCRIPTION") && dir.exists("R")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(fmrireg)
  }
}

make_dataset <- function(n_runs = 3L, run_len = 160L, n_vox = 1200L, seed = 123L) {
  set.seed(seed)
  n <- n_runs * run_len
  onsets <- rep(seq(12, run_len - 20, by = 16), n_runs)
  event_table <- data.frame(
    onset = onsets,
    condition = factor(rep(c("A", "B"), length.out = length(onsets))),
    run = rep(seq_len(n_runs), each = length(seq(12, run_len - 20, by = 16)))
  )

  Y <- matrix(rnorm(n * n_vox), n, n_vox)
  dset <- matrix_dataset(
    Y,
    TR = 1,
    run_length = rep(run_len, n_runs),
    event_table = event_table
  )

  sframe <- sampling_frame(rep(run_len, n_runs), TR = 1)
  base <- baseline_model(basis = "poly", degree = 2, sframe = sframe)

  list(dataset = dset, baseline = base)
}

bench_fit <- function(strategy, dataset, baseline, legacy_fallback = FALSE, nchunks = 6L) {
  old_opt <- getOption("fmrireg.ar.fixed_order_legacy_fallback")
  options(fmrireg.ar.fixed_order_legacy_fallback = legacy_fallback)
  on.exit(options(fmrireg.ar.fixed_order_legacy_fallback = old_opt), add = TRUE)

  gc()
  tm <- system.time({
    fmri_lm(
      onset ~ hrf(condition),
      block = ~ run,
      dataset = dataset,
      baseline_model = baseline,
      strategy = strategy,
      nchunks = nchunks,
      use_fast_path = TRUE,
      ar_options = list(struct = "ar1", iter_gls = 3),
      robust = "huber",
      robust_options = list(max_iter = 5, scale_scope = "run"),
      progress = FALSE
    )
  })

  data.frame(
    strategy = strategy,
    legacy_fallback = legacy_fallback,
    elapsed_sec = unname(tm[["elapsed"]]),
    user_sec = unname(tm[["user.self"]]),
    system_sec = unname(tm[["sys.self"]]),
    stringsAsFactors = FALSE
  )
}

run_benchmark <- function(n_runs = 3L, run_len = 160L, n_vox = 1200L, seed = 123L,
                          nchunks = 6L, n_reps = 5L, warmup = 1L,
                          out_csv = "bench/glm_efficiency_benchmark.csv") {
  load_fmrireg_for_bench()
  dat <- make_dataset(n_runs = n_runs, run_len = run_len, n_vox = n_vox, seed = seed)

  cases <- list(
    list(strategy = "runwise", legacy_fallback = FALSE),
    list(strategy = "runwise", legacy_fallback = TRUE),
    list(strategy = "chunkwise", legacy_fallback = FALSE),
    list(strategy = "chunkwise", legacy_fallback = TRUE)
  )

  raw_results <- do.call(
    rbind,
    lapply(cases, function(cs) {
      reps <- lapply(seq_len(n_reps + warmup), function(rep_i) {
        out <- bench_fit(
          strategy = cs$strategy,
          dataset = dat$dataset,
          baseline = dat$baseline,
          legacy_fallback = cs$legacy_fallback,
          nchunks = nchunks
        )
        out$rep <- rep_i
        out
      })
      do.call(rbind, reps)
    })
  )

  kept <- subset(raw_results, rep > warmup)
  med <- aggregate(cbind(elapsed_sec, user_sec, system_sec) ~ strategy + legacy_fallback,
                   data = kept, FUN = median)
  mean_vals <- aggregate(cbind(elapsed_sec, user_sec, system_sec) ~ strategy + legacy_fallback,
                         data = kept, FUN = mean)
  sd_vals <- aggregate(cbind(elapsed_sec, user_sec, system_sec) ~ strategy + legacy_fallback,
                       data = kept, FUN = sd)

  results <- med
  names(results)[names(results) %in% c("elapsed_sec", "user_sec", "system_sec")] <-
    c("elapsed_sec_median", "user_sec_median", "system_sec_median")
  results$elapsed_sec_mean <- mean_vals$elapsed_sec
  results$elapsed_sec_sd <- sd_vals$elapsed_sec
  results$n_reps <- n_reps

  dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
  write.csv(results, out_csv, row.names = FALSE)
  write.csv(raw_results, sub("[.]csv$", "_raw.csv", out_csv), row.names = FALSE)

  cat("GLM efficiency benchmark complete\n")
  print(results)
  cat("\nSaved:", out_csv, "\n")
  cat("Saved:", sub("[.]csv$", "_raw.csv", out_csv), "\n")
}

if (sys.nframe() == 0L) {
  run_benchmark()
}
