###############################################################################
## collect_results.R
##
## The inverse of the fan-out: gather the per-subject maps that
## reduce_write_results() wrote to disk and hand them to the group level as a
## group_data object ready for fmri_meta(). This owns the BIDS naming convention
## end-to-end (write_results writes it; collect_results reads it back), closing
## the fan-out -> fan-in loop in one call.
##
## See `docs/plans/model-templates-multisubject-fanout.md`.
###############################################################################

#' Collect per-subject result maps into a group_data object
#'
#' Globs the per-subject statistical maps written by [reduce_write_results()]
#' under \code{dir}, pairs them by subject, and builds a [group_data()] object
#' for [fmri_meta()]. When both beta and se maps are present the result carries
#' variance (suitable for inverse-variance meta-analysis); a beta-only result is
#' returned otherwise.
#'
#' @param dir Directory (searched recursively) containing the written maps.
#' @param space Optional BIDS space label to filter on (e.g.
#'   \code{"MNI152NLin2009cAsym"}); \code{NULL} takes all.
#' @param beta_desc,se_desc The BIDS \code{desc-} labels of the effect and
#'   standard-error maps (defaults \code{"beta"} / \code{"se"}).
#' @param format Group-data backend format (currently \code{"nifti"}).
#' @return A \code{group_data} object.
#' @note Pairing by BIDS subject requires \pkg{fmrigds} with BIDS-aware subject
#'   keying. Multi-file \code{"h5"} aggregation is not yet supported upstream.
#' @seealso [reduce_write_results()], [run_jobs()], [fmri_meta()]
#' @export
#' @examples
#' \dontrun{
#' tmpl <- fmri_template(onset ~ hrf(trial_type, contrasts = con), ~ run,
#'                       reducer = reduce_write_results(format = "nifti",
#'                                                      stats = c("beta", "se"),
#'                                                      path = "study/glm"))
#' run_jobs(instantiate(tmpl, manifest))
#' gd <- collect_results("study/glm", space = "MNI152NLin2009cAsym")
#' fm <- fmri_meta(gd, ~ 1, method = "fe")
#' }
collect_results <- function(dir, space = NULL,
                            beta_desc = "beta", se_desc = "se",
                            format = c("nifti")) {
  format <- match.arg(format)
  assert_that(is.character(dir), length(dir) == 1, dir.exists(dir),
              msg = "'dir' must be an existing directory")

  pat <- function(d) {
    sp <- if (!is.null(space)) sprintf("space-%s_.*", space) else ""
    sprintf("%sdesc-%s_bold[.]nii([.]gz)?$", sp, d)
  }
  beta <- sort(list.files(dir, pattern = pat(beta_desc), full.names = TRUE,
                          recursive = TRUE))
  se <- sort(list.files(dir, pattern = pat(se_desc), full.names = TRUE,
                        recursive = TRUE))

  assert_that(length(beta) > 0,
              msg = sprintf("no 'desc-%s' maps found under %s", beta_desc, dir))

  data <- if (length(se) == length(beta) && length(se) > 0) {
    list(beta = beta, se = se)
  } else {
    list(beta = beta)
  }
  group_data(data, format = format)
}
