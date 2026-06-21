fmrireg_cli_bootstrap <- function(command = "fmrireg") {
  status <- tryCatch(
    {
      args <- commandArgs(trailingOnly = TRUE)
      root <- .find_fmrireg_source_root(command)

      if (!is.null(root) && requireNamespace("pkgload", quietly = TRUE)) {
        pkgload::load_all(root, export_all = FALSE, helpers = FALSE, quiet = TRUE)
        return(getExportedValue("fmrireg", "fmrireg_cli")(args))
      }

      if (requireNamespace("fmrireg", quietly = TRUE)) {
        return(fmrireg::fmrireg_cli(args))
      }

      stop(
        "Package 'fmrireg' is not installed, and source-checkout fallback ",
        "requires the 'pkgload' package."
      )
    },
    error = function(err) {
      message(command, ": ", conditionMessage(err))
      2L
    }
  )

  quit(save = "no", status = as.integer(status))
}

.find_fmrireg_source_root <- function(command) {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0L) {
    return(NULL)
  }

  script_path <- normalizePath(
    sub("^--file=", "", file_arg[[1]]),
    winslash = "/",
    mustWork = FALSE
  )
  root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)

  # Only treat `root` as a development source checkout when R/ actually holds
  # `.R` sources. An *installed* package also has DESCRIPTION + an R/ directory,
  # but R/ there contains only the compiled lazy-load DB (no `.R` files), so
  # keying on DESCRIPTION + R/ alone would mis-fire and trigger
  # pkgload::load_all() on the installed tree, which errors. When there are no
  # `.R` sources, return NULL so the caller falls back to the installed namespace.
  r_dir <- file.path(root, "R")
  is_source_checkout <-
    file.exists(file.path(root, "DESCRIPTION")) &&
    dir.exists(r_dir) &&
    length(list.files(r_dir, pattern = "\\.[Rr]$")) > 0L

  if (is_source_checkout) {
    return(root)
  }

  NULL
}
