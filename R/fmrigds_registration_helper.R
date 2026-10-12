# Session state for the fmrigds bindings: registered once per session.
.fmrigds_reg_state <- new.env(parent = emptyenv())
.fmrigds_reg_state$done <- FALSE

# Register fmrireg's fmrigds reducers/posthoc on first use. Idempotent: later
# calls return immediately. A failure is reported as a warning (fmrireg's GDS
# paths fit subject-level data directly and do not depend on the reducers)
# and leaves the state unset so the next call retries.
.ensure_fmrigds_registered <- function() {
  if (isTRUE(.fmrigds_reg_state$done)) return(invisible(TRUE))
  if (!requireNamespace("fmrigds", quietly = TRUE)) return(invisible(FALSE))
  ok <- tryCatch(
    isTRUE(.register_fmrireg_bindings()),
    error = function(e) {
      warning("fmrireg could not register its fmrigds reducers: ",
              conditionMessage(e), call. = FALSE)
      FALSE
    }
  )
  .fmrigds_reg_state$done <- ok
  invisible(ok)
}
