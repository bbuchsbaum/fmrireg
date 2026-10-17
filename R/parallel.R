## Parallel execution
##
## fmrireg has two layers of parallelism and one switch for each:
## - R-level loops (voxel chunks, voxelwise AR fits, trial-wise mixed models,
##   run_jobs() fan-out) run through .fmrireg_lapply(), which uses the user's
##   future::plan() when parallel. future is already a dependency of
##   neuroim2 and fmridataset, so this costs nothing extra.
## - C++ kernels use OpenMP threads, set by .fmrireg_num_threads()
##   (options(fmrireg.num_threads) or FMRIREG_NUM_THREADS).

#' Map over a list sequentially or through the active future plan
#'
#' The single R-level parallel entry point. With `parallel = FALSE` this is
#' `lapply()`; with `parallel = TRUE` it is `future.apply::future_lapply()`
#' under the user's `future::plan()` (sequential unless they set one), with
#' parallel-safe RNG streams. `progress = TRUE` shows a cli progress bar when
#' sequential, and a progressr progressor when parallel and progressr is
#' installed.
#' @keywords internal
#' @noRd
.fmrireg_lapply <- function(X, FUN, parallel = FALSE, progress = FALSE,
                            label = "Processing") {
  if (!isTRUE(parallel)) {
    if (!isTRUE(progress)) {
      return(lapply(X, FUN))
    }
    pb <- cli::cli_progress_bar(label, total = length(X), clear = FALSE)
    on.exit(cli::cli_progress_done(id = pb), add = TRUE)
    return(lapply(X, function(x) {
      out <- FUN(x)
      cli::cli_progress_update(id = pb)
      out
    }))
  }

  if (isTRUE(progress) && requireNamespace("progressr", quietly = TRUE)) {
    return(progressr::with_progress({
      p <- progressr::progressor(steps = length(X))
      future.apply::future_lapply(X, function(x) {
        out <- FUN(x)
        p()
        out
      }, future.seed = TRUE)
    }))
  }
  future.apply::future_lapply(X, FUN, future.seed = TRUE)
}
