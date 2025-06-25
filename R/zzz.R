#' @importFrom RcppParallel setThreadOptions
.onLoad <- function(libname, pkgname) {
  env_threads <- Sys.getenv("FMRIREG_NUM_THREADS", unset = NA)
  opt_threads <- getOption("fmrireg.num_threads", default = NA)
  val <- NA
  if (!is.na(opt_threads)) {
    val <- opt_threads
  } else if (!is.na(env_threads) && nzchar(env_threads)) {
    val <- as.numeric(env_threads)
  }
  if (!is.na(val) && val > 0) {
    try({
      RcppParallel::setThreadOptions(numThreads = as.integer(val))
    }, silent = TRUE)
  }
}
