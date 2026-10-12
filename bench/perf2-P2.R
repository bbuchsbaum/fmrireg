# Performance pass 2, stream P2 (trial betas and outputs): old vs new, back to
# back in one R process, single-threaded.
#
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript bench/perf2-P2.R
#
# Run from the package root of a checkout that contains the P2 commits. The
# pre-P2 write_results() code is taken from git (`OLD_SHA`, the commit before
# stream P2); the old mixed path is the per-voxel loop the package used
# (mixed_betas() -> fmrilss::mixed_solve() per voxel), which still exists.

OLD_SHA <- Sys.getenv("P2_OLD_SHA", "206600c")

suppressMessages(devtools::load_all(".", quiet = TRUE))
options(fmrireg.num_threads = 1)
ns <- asNamespace("fmrireg")

timeit <- function(expr, reps = 1L) {
  expr <- substitute(expr)
  env <- parent.frame()
  best <- Inf
  for (r in seq_len(reps)) {
    gc(FALSE)
    t0 <- proc.time()[["elapsed"]]
    eval(expr, env)
    best <- min(best, proc.time()[["elapsed"]] - t0)
  }
  best
}
results <- list()
add <- function(item, case, old, new, check) {
  results[[length(results) + 1L]] <<- data.frame(
    item = item, case = case, old_s = round(old, 3), new_s = round(new, 3),
    speedup = round(old / new, 1), check = check, stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# P2-1: mixed-model trial betas
# ---------------------------------------------------------------------------
mixed_case <- function(n_run, run_len, m_per_run, V, seed = 1) {
  set.seed(seed)
  n <- n_run * run_len
  tt <- seq_len(n)
  onsets <- sort(sample(seq(1, n - 10), m_per_run * n_run))
  Z <- vapply(onsets, function(o) { x <- pmax(tt - o, 0); x^2 * exp(-x / 1.5) / 2 }, numeric(n))
  run <- rep(seq_len(n_run), each = run_len)
  Base <- do.call(cbind, lapply(seq_len(n_run), function(r) {
    idx <- run == r
    B <- cbind(1, stats::poly(seq_len(run_len), 2))
    out <- matrix(0, n, ncol(B)); out[idx, ] <- B; out
  }))
  Xi <- matrix(rnorm(n), n, 1)
  X <- cbind(Z, Xi)
  Y <- Z %*% matrix(rnorm(ncol(Z) * V, sd = 0.7), ncol(Z), V) +
    Base %*% matrix(rnorm(ncol(Base) * V, sd = 3), ncol(Base), V) +
    matrix(rnorm(n * V), n, V)
  list(X = X, Y = Y, Base = Base, ran = seq_len(ncol(Z)), fixed = ncol(Z) + 1L)
}

old_mixed <- function(cs) {
  vapply(seq_len(ncol(cs$Y)), function(v) {
    mixed_betas(cs$X, cs$Y[, v], ran_ind = cs$ran, fixed_ind = cs$fixed,
                nuisance = cs$Base)
  }, numeric(length(cs$ran) + length(cs$fixed)))
}
new_mixed <- function(cs) {
  mixed_betas_matrix(cs$X, cs$Y, ran_ind = cs$ran, fixed_ind = cs$fixed,
                     nuisance = cs$Base)
}

for (spec in list(list(n_run = 1, run_len = 300, m = 40, V = 200),
                  list(n_run = 3, run_len = 300, m = 30, V = 30))) {
  cs <- mixed_case(spec$n_run, spec$run_len, spec$m, spec$V)
  old_out <- NULL; new_out <- NULL
  t_old <- timeit(old_out <- old_mixed(cs))
  t_new <- timeit(new_out <- new_mixed(cs))
  rel <- max(abs(old_out - new_out)) / max(abs(old_out))
  add("P2-1 mixed betas",
      sprintf("n=%d, trials=%d, V=%d", nrow(cs$X), length(cs$ran), ncol(cs$Y)),
      t_old, t_new, sprintf("max rel diff %.1e", rel))
}
# New path only, at a whole-brain-ish voxel count.
cs <- mixed_case(1, 300, 40, 20000, seed = 2)
t_big <- timeit(new_mixed(cs))
add("P2-1 mixed betas", "n=300, trials=40, V=20000 (new only)", NA, t_big, "-")

# ---------------------------------------------------------------------------
# Volumetric fit shared by P2-2 / P2-3
# ---------------------------------------------------------------------------
set.seed(31)
dims <- c(24, 24, 16); nt <- 120
scans <- lapply(1:2, function(r) neuroim2::NeuroVec(array(rnorm(prod(dims) * nt), c(dims, nt)),
                                                    neuroim2::NeuroSpace(dim = c(dims, nt))))
mask <- neuroim2::LogicalNeuroVol(array(runif(prod(dims)) > 0.3, dims), neuroim2::NeuroSpace(dim = dims))
etab <- data.frame(onset = rep(seq(5, 105, by = 10), 2),
                   condition = factor(rep(c("A", "B", "C", "D"), length.out = 22)),
                   run = rep(1:2, each = 11))
dset <- neurovec_frame(scans = scans, mask = mask, TR = 1.5, event_table = etab)
con <- contrast_set(pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B"),
                    pair_contrast(~ condition == "C", ~ condition == "D", name = "C_vs_D"),
                    oneway_contrast(~ condition, name = "main"))
fit <- fmri_lm(onset ~ hrf(condition, contrasts = con), block = ~ run, dataset = dset,
               progress = FALSE)

# ---------------------------------------------------------------------------
# P2-2: coef_images (old = one coef_image() per coefficient, which is what the
# old coef_images() did)
# ---------------------------------------------------------------------------
old_ci <- function() lapply(coef_names(fit), function(nm) coef_image(fit, coef = nm))
t_old <- timeit(a <- old_ci(), reps = 3)
t_new <- timeit(b <- coef_images(fit), reps = 3)
add("P2-2 coef_images", sprintf("%d coefs, %d voxels", length(b), sum(mask)), t_old, t_new,
    sprintf("identical: %s", identical(unname(a), unname(b))))

# ---------------------------------------------------------------------------
# P2-3: write_results(format = c("h5", "nifti"))
# ---------------------------------------------------------------------------
if (requireNamespace("fmristore", quietly = TRUE)) {
  old_file <- tempfile(fileext = ".R")
  system2("git", c("show", paste0(OLD_SHA, ":R/bids_export.R")), stdout = old_file)
  old_env <- new.env(parent = ns)
  sys.source(old_file, envir = old_env)
  write_old <- get("write_results.fmri_lm", envir = old_env)
  write_new <- get("write_results.fmri_lm", envir = ns)
  for (strategy in c("by_stat", "by_contrast")) {
    run_one <- function(fun) {
      d <- tempfile("p2bench-"); dir.create(d)
      on.exit(unlink(d, recursive = TRUE))
      suppressWarnings(fun(fit, path = d, subject = "01", task = "t", space = "MNI",
                           format = c("h5", "nifti"), strategy = strategy))
      invisible(NULL)
    }
    t_old <- timeit(run_one(write_old), reps = 2)
    t_new <- timeit(run_one(write_new), reps = 2)
    add("P2-3 write_results h5+nifti", strategy, t_old, t_new,
        "outputs identical (see test-perf2-P2-write-results.R)")
  }
}

out <- do.call(rbind, results)
print(out, row.names = FALSE)
