# Thread count for the OpenMP kernels (meta-analysis, t-tests, spatial FDR):
# options(fmrireg.num_threads) if set, else the environment variable
# FMRIREG_NUM_THREADS if set and positive, else 0 (the OpenMP default).
# Read at call time, so either can be changed after the package is loaded.
.fmrireg_num_threads <- function() {
  opt <- getOption("fmrireg.num_threads")
  if (!is.null(opt)) {
    n <- suppressWarnings(as.integer(opt)[1L])
    return(if (length(n) != 1L || is.na(n) || n < 0L) 0L else n)
  }
  env <- suppressWarnings(as.integer(Sys.getenv("FMRIREG_NUM_THREADS", unset = "")))
  if (length(env) == 1L && !is.na(env) && env > 0L) env else 0L
}

.onLoad <- function(libname, pkgname) {
  # Expose fmrireg's spatial FDR to fmrigds pipelines as the post-hoc
  # "fmrireg:fdr:spatial". fmrigds's own reducers/post-hocs are never replaced.
  try(.register_fmrireg_posthoc(), silent = TRUE)

  if (exists(".register_builtin_engines", mode = "function")) {
    try(.register_builtin_engines(), silent = TRUE)
  }
}
