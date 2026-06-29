###############################################################################
## export_jobs.R
##
## Serialize a batch of jobs to disk for execution by ANY external array
## scheduler. No scheduler-specific code is emitted: `export_jobs()` writes a
## serialized manifest, a list of ids, and a small backend-agnostic `run_one.R`
## that runs one job by (1-based) index. The user wires this into SLURM/SGE/PBS
## /a local loop by passing the array index as the first argument.
##
## See `docs/plans/model-templates-multisubject-fanout.md`.
###############################################################################

#' Export jobs for external / array execution
#'
#' Writes a serialized job manifest plus a minimal, backend-agnostic
#' \code{run_one.R} runner to \code{dir}. Deliberately emits no scheduler code:
#' drive it with whatever array system you have, passing the 1-based job index
#' as the first argument, e.g. \code{Rscript run_one.R $SLURM_ARRAY_TASK_ID}.
#'
#' @param jobs A single [fmri_job] or a list of them (e.g. from [instantiate()]).
#' @param dir Output directory (created if needed).
#' @param overwrite Overwrite an existing manifest. Default \code{FALSE}.
#' @param setup Character vector of R lines run at the top of \code{run_one.R}
#'   before jobs are executed (load packages, set threads, etc.). Default loads
#'   \pkg{fmrireg}.
#' @param results_dir Default output subdirectory written by \code{run_one.R}.
#' @return Invisibly, a list describing what was written (\code{dir},
#'   \code{manifest}, \code{runner}, \code{n}, \code{ids}).
#' @seealso [read_jobs()], [run_jobs()], [instantiate()]
#' @export
#' @examples
#' \dontrun{
#' jobs <- instantiate(tmpl, manifest)
#' export_jobs(jobs, "study/jobs")
#' # then, per array task:  Rscript study/jobs/run_one.R $SLURM_ARRAY_TASK_ID
#' }
export_jobs <- function(jobs, dir, overwrite = FALSE,
                        setup = "library(fmrireg)",
                        results_dir = "results") {
  if (inherits(jobs, "fmri_job")) jobs <- list(jobs)
  assert_that(is.list(jobs) && length(jobs) > 0,
              msg = "'jobs' must be an fmri_job or a non-empty list of fmri_jobs")
  assert_that(all(vapply(jobs, inherits, logical(1), "fmri_job")),
              msg = "every element of 'jobs' must be an 'fmri_job'")
  assert_that(is.character(setup), msg = "'setup' must be a character vector")

  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  manifest_path <- file.path(dir, "manifest.rds")
  if (file.exists(manifest_path) && !isTRUE(overwrite)) {
    stop(sprintf("manifest already exists at %s (use overwrite = TRUE)", manifest_path),
         call. = FALSE)
  }

  ids <- vapply(jobs, `[[`, character(1), "id")
  assert_that(!anyDuplicated(ids), msg = "job ids must be unique")

  saveRDS(jobs, manifest_path)
  writeLines(ids, file.path(dir, "job_ids.txt"))
  runner_path <- file.path(dir, "run_one.R")
  writeLines(.run_one_script(setup = setup, results_dir = results_dir), runner_path)

  invisible(list(dir = dir, manifest = manifest_path, runner = runner_path,
                 n = length(jobs), ids = ids))
}

#' Read jobs exported by export_jobs()
#'
#' @param dir Directory containing \code{manifest.rds}.
#' @return The list of [fmri_job]s.
#' @seealso [export_jobs()]
#' @export
read_jobs <- function(dir) {
  manifest_path <- if (file.exists(file.path(dir, "manifest.rds"))) {
    file.path(dir, "manifest.rds")
  } else dir
  assert_that(file.exists(manifest_path), msg = "no manifest.rds found")
  readRDS(manifest_path)
}

#' @keywords internal
#' @noRd
.run_one_script <- function(setup = "library(fmrireg)", results_dir = "results") {
  c(
    "#!/usr/bin/env Rscript",
    "## Backend-agnostic runner: execute ONE fmri_job by 1-based index.",
    "## Usage: Rscript run_one.R <index> [results_dir]",
    "## Drive with any array scheduler, e.g.:",
    "##   SLURM: Rscript run_one.R $SLURM_ARRAY_TASK_ID",
    "##   SGE:   Rscript run_one.R $SGE_TASK_ID",
    "##   local: for i in $(seq 1 N); do Rscript run_one.R $i; done",
    "",
    "args <- commandArgs(trailingOnly = TRUE)",
    "if (length(args) < 1L) stop('usage: Rscript run_one.R <index> [results_dir]')",
    "idx <- as.integer(args[[1L]])",
    sprintf("results_dir <- if (length(args) >= 2L) args[[2L]] else %s",
            deparse(results_dir)),
    "",
    "## --- setup ---",
    setup,
    "",
    "## locate this script's directory to find manifest.rds",
    "argv <- commandArgs(FALSE)",
    "this <- sub('^--file=', '', grep('^--file=', argv, value = TRUE))",
    "job_dir <- if (length(this)) dirname(normalizePath(this)) else getwd()",
    "## a relative results_dir is resolved against the manifest's directory",
    "if (!grepl('^(/|[A-Za-z]:)', results_dir)) results_dir <- file.path(job_dir, results_dir)",
    "jobs <- readRDS(file.path(job_dir, 'manifest.rds'))",
    "if (is.na(idx) || idx < 1L || idx > length(jobs)) stop(sprintf('index must be an integer in 1..%d', length(jobs)))",
    "job <- jobs[[idx]]",
    "",
    "value <- fmrireg::run_job(job)",
    "if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)",
    "## defense-in-depth: keep the output filename inside results_dir",
    "safe_id <- gsub('[^A-Za-z0-9._-]', '_', job$id)",
    "out <- file.path(results_dir, paste0(safe_id, '.rds'))",
    "saveRDS(value, out)",
    "message(sprintf('[run_one] job %s -> %s', job$id, out))"
  )
}
