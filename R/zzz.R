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
  
  # Register S3 methods from fmridesign for Fcontrasts
  # This ensures the methods are available when fmrireg is loaded
  # The methods are not exported, so we get them from the namespace
  if (requireNamespace("fmridesign", quietly = TRUE)) {
    registerS3method("Fcontrasts", "event_model", 
                     utils::getFromNamespace("Fcontrasts.event_model", "fmridesign"))
    registerS3method("Fcontrasts", "event_term", 
                     utils::getFromNamespace("Fcontrasts.event_term", "fmridesign"))
    registerS3method("Fcontrasts", "convolved_term", 
                     utils::getFromNamespace("Fcontrasts.convolved_term", "fmridesign"))

    # Register correlation_map.event_model from fmridesign
    # This method doesn't depend on internal functions so it can be safely registered
    registerS3method("correlation_map", "event_model",
                     utils::getFromNamespace("correlation_map.event_model", "fmridesign"))

    # Register our columns() method for event_model objects with the fmridesign generic
    registerS3method("columns", "event_model", columns.event_model,
                     envir = asNamespace("fmridesign"))
    # Note: correlation_map.baseline_model is defined locally in correlation_map_methods.R
    # because it needs access to fmrireg's internal .correlation_map_common function
  }

  # fmrigds bridge registrations happen via fmrigds .onLoad when available
}
