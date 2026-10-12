# Performance pass 2, stream P3: old vs new timings, back to back, one thread.
#
# Run from the package root:  Rscript bench/perf2-P3.R
# The "old" code is the pre-P3 revision (OLD_REV), taken from git:
#   - src/meta_kernels.cpp is compiled with Rcpp::sourceCpp;
#   - the old R implementations are sourced into an environment whose parent
#     is the fmrireg namespace (so they call the same helpers as the new code).
# Every timing also checks that old and new give the same numbers.

OLD_REV <- Sys.getenv("P3_OLD_REV", "206600c")
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
options(fmrireg.num_threads = 1)
suppressMessages(devtools::load_all(".", quiet = TRUE))
ns <- asNamespace("fmrireg")

git_show <- function(path) {
  out <- tempfile(fileext = paste0(".", tools::file_ext(path)))
  writeLines(system2("git", c("show", paste0(OLD_REV, ":", path)), stdout = TRUE), out)
  out
}
old_env <- new.env(parent = ns)
for (f in c("R/fmri_meta_gds.R", "R/fmri_lm_engine_rrr.R", "R/soft_subspace_projection.R")) {
  sys.source(git_show(f), envir = old_env)
}
# Both kernels are compiled the same way (sourceCpp, default optimisation);
# load_all() builds the package DLL without optimisation, so the package's
# own meta_fit_cov_cpp is not a fair "new" timing.
old_cpp <- new.env()
Rcpp::sourceCpp(git_show("src/meta_kernels.cpp"), env = old_cpp, rebuild = TRUE)
new_cpp <- new.env()
Rcpp::sourceCpp("src/meta_kernels.cpp", env = new_cpp, rebuild = TRUE)

timeit <- function(expr, reps = 3L) {
  expr <- substitute(expr)
  env <- parent.frame()
  t <- vapply(seq_len(reps), function(i) {
    gc(FALSE)
    system.time(eval(expr, env))[["elapsed"]]
  }, numeric(1))
  stats::median(t)
}
results <- list()
record <- function(item, case, old, new, maxdiff) {
  results[[length(results) + 1L]] <<- data.frame(
    item = item, case = case, old_s = round(old, 3), new_s = round(new, 3),
    speedup = round(old / new, 1), max_abs_diff = signif(maxdiff, 3)
  )
  cat(sprintf("%-5s %-34s old %8.3fs  new %8.3fs  x%5.1f  maxdiff %.2g\n",
              item, case, old, new, old / new, maxdiff))
}
maxdiff <- function(a, b) {
  num <- function(x) unlist(rapply(list(x), function(z) if (is.numeric(z)) as.numeric(z) else NULL,
                                   how = "unlist"))
  a <- num(a); b <- num(b)
  ok <- is.finite(a) & is.finite(b)
  stopifnot(identical(is.na(a), is.na(b)))
  max(c(0, abs(a[ok] - b[ok])))
}

# ---- P3-1: meta kernels ------------------------------------------------------
set.seed(1)
S <- 30L; P <- 20000L
cov <- data.frame(x = rnorm(S), g = factor(rep(c("a", "b", "c"), length.out = S)))
V <- matrix(runif(S * P, 0.02, 0.5), S, P)
Y <- 0.2 + 0.3 * cov$x + matrix(rnorm(S * P), S, P) * sqrt(V + 0.1)
Y[sample(length(Y), 2000)] <- NA
for (xn in c("intercept", "regression")) {
  X <- if (xn == "intercept") matrix(1, S, 1) else model.matrix(~ x + g, cov)
  for (m in c("fe", "dl", "pm")) {
    o <- old_cpp$meta_fit_cov_cpp(Y, V, X, m, "none", 1.345, 2L, 1L)
    n <- new_cpp$meta_fit_cov_cpp(Y, V, X, m, "none", 1.345, 2L, 1L)
    record("P3-1", paste(m, xn, "P=20000"),
           timeit(old_cpp$meta_fit_cov_cpp(Y, V, X, m, "none", 1.345, 2L, 1L)),
           timeit(new_cpp$meta_fit_cov_cpp(Y, V, X, m, "none", 1.345, 2L, 1L)),
           maxdiff(o, n))
  }
}

# ---- P3-2 / P3-3: GDS meta path ---------------------------------------------
set.seed(2)
S <- 20L; P <- 3000L
df <- expand.grid(subject = sprintf("s%02d", seq_len(S)), roi = sprintf("R%05d", seq_len(P)),
                  contrast = "c1", stringsAsFactors = FALSE)
sid <- as.integer(sub("s", "", df$subject))
df$group <- factor(ifelse(sid <= S / 2, "A", "B"))
df$se <- runif(nrow(df), 0.1, 0.4)
df$beta <- 0.2 + 0.3 * (df$group == "B") + rnorm(nrow(df), sd = sqrt(df$se^2 + 0.05))
df$df <- 30
df$t <- df$beta / df$se
gd <- ns$.group_data_gds(df, format = "csv",
  effect_cols = c(beta = "beta", se = "se", t = "t", df = "df"),
  subject_col = "subject", roi_col = "roi", contrast_col = "contrast",
  covariate_cols = "group")
strip <- function(o) { o$data <- NULL; unclass(o) }
run_meta <- function(impl, ...) strip(impl(gd, verbose = FALSE, ...))
for (cs in list(
  list(case = "pm ~group + contrast", args = list(formula = ~ group, method = "pm",
                                                  contrasts = c(groupB = 1))),
  list(case = "fe ~1 (no contrast)", args = list(formula = ~ 1, method = "fe"))
)) {
  o <- do.call(run_meta, c(list(old_env$.fmri_meta_group_data_gds_impl), cs$args))
  n <- do.call(run_meta, c(list(ns$.fmri_meta_group_data_gds_impl), cs$args))
  stopifnot(isTRUE(all.equal(o, n, tolerance = 1e-12)))
  record("P3-2", paste(cs$case, "P=3000"),
         timeit(do.call(run_meta, c(list(old_env$.fmri_meta_group_data_gds_impl), cs$args))),
         timeit(do.call(run_meta, c(list(ns$.fmri_meta_group_data_gds_impl), cs$args))),
         maxdiff(o[c("coefficients", "se", "tau2", "Q", "cov", "contrasts")],
                 n[c("coefficients", "se", "tau2", "Q", "cov", "contrasts")]))
}
for (cm in c("stouffer", "lancaster")) {
  a <- list(formula = ~ 1, combine = cm, weights = "equal")
  o <- do.call(run_meta, c(list(old_env$.fmri_meta_group_data_gds_impl), a))
  n <- do.call(run_meta, c(list(ns$.fmri_meta_group_data_gds_impl), a))
  stopifnot(isTRUE(all.equal(o, n, tolerance = 1e-12)))
  record("P3-3", paste(cm, "combine P=3000"),
         timeit(do.call(run_meta, c(list(old_env$.fmri_meta_group_data_gds_impl), a))),
         timeit(do.call(run_meta, c(list(ns$.fmri_meta_group_data_gds_impl), a))),
         maxdiff(o$coefficients, n$coefficients))
}

# ---- P3-4: rrr_gls on latent data -------------------------------------------
set.seed(3)
n_time <- 400L; n_comp <- 20L; n_vox <- 60000L
ev <- data.frame(onset = seq(6, 760, by = 14), run = 1L)
ev$condition <- factor(rep(c("a", "b", "c"), length.out = nrow(ev)))
basis <- matrix(rnorm(n_time * n_comp), n_time, n_comp)
loadings <- matrix(rnorm(n_vox * n_comp), n_vox, n_comp)
lvec <- fmristore::LatentNeuroVec(basis = basis, loadings = loadings,
  space = neuroim2::NeuroSpace(c(n_vox, 1, 1, n_time)),
  mask = rep(TRUE, n_vox), offset = rep(0, n_vox))
lds <- latent_frame(lvec, TR = 2, run_length = n_time, event_table = ev)
for (a in list(list(rank = 2L), list(rank = 2L, se_mode = "bootstrap", bootstrap_n = 20L,
                                      bootstrap_seed = 1L))) {
  fit_new <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = lds,
                     engine = "rrr_gls", engine_args = a)
  model <- fit_new$model
  cfg <- attr(fit_new, "config")
  o <- old_env$.fit_rrr_gls_engine(model, lds, a, cfg)
  n <- ns$.fit_rrr_gls_engine(model, lds, a, cfg)
  record("P3-4", paste(if (is.null(a$se_mode)) "conditional" else "bootstrap(20)",
                       "T=400 k=20 V=60000"),
         timeit(old_env$.fit_rrr_gls_engine(model, lds, a, cfg)),
         timeit(ns$.fit_rrr_gls_engine(model, lds, a, cfg)),
         max(maxdiff(o$betas_fixed, n$betas_fixed),
             maxdiff(o$result$betas$data[[1]]$se[[1]], n$result$betas$data[[1]]$se[[1]]),
             maxdiff(o$result$betas$data[[1]]$stat[[1]], n$result$betas$data[[1]]$stat[[1]])))
}

# ---- P3-5: soft projection GCV ----------------------------------------------
set.seed(4)
n <- 300L; p <- 400L; V <- 20000L
N <- matrix(rnorm(n * p), n, p)
Yd <- matrix(rnorm(n * V), n, V)
o <- old_env$soft_projection(N, lambda = "gcv", Y = Yd)
nn <- ns$soft_projection(N, lambda = "gcv", Y = Yd)
record("P3-5", "soft_projection(gcv) n=300 p=400 V=20000",
       timeit(old_env$soft_projection(N, lambda = "gcv", Y = Yd)),
       timeit(ns$soft_projection(N, lambda = "gcv", Y = Yd)),
       abs(o$lambda - nn$lambda))
sv <- svd(N, nu = min(n, p), nv = 0)
record("P3-5", ".select_lambda_gcv only",
       timeit(old_env$.select_lambda_gcv(N, Yd, sv)),
       timeit(ns$.select_lambda_gcv(N, Yd, sv)),
       abs(old_env$.select_lambda_gcv(N, Yd, sv) - ns$.select_lambda_gcv(N, Yd, sv)))

cat("\n")
print(do.call(rbind, results), row.names = FALSE)
