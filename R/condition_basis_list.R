#' Convert an `event_term` to a per-condition basis list
#'
#' A lightweight wrapper around [`fmrireg::convolve()`] that post-processes the
#' resulting design matrix into a named list of T x d matrices - one per
#' experimental condition ("base condition tag").  This keeps **all** of the
#' heavy lifting inside **fmrireg** while exposing a minimal, pipe-friendly API
#' that can be used anywhere a condition -> basis split is required (e.g. for
#' CFALS).
#'
#' @param x            An [`event_term`] object.
#' @param hrf          An [`HRF`] object to apply.
#' @param sampling_frame A [`sampling_frame`] object defining the temporal grid.
#' @param ...          Further arguments passed on to
#'                     [`fmrireg::convolve()`] (e.g. `drop.empty = FALSE`).
#' @param output       Either "matrix" (default) for the ordinary design matrix
#'                     or "condition_list" for the split-by-condition list.
#'
#' @return A numeric *matrix* or a named *list* of matrices, depending on
#'         `output`.
#' @export
condition_basis_list <- function(x, hrf, sampling_frame, ...,
                                 output = c("condition_list", "matrix")) {
  stopifnot(inherits(x, "event_term"),
            inherits(hrf, "HRF"),
            inherits(sampling_frame, "sampling_frame"))

  output <- match.arg(output)

  # 1. Convolve in the usual way ------------------------------------------------
  dm <- convolve(x, hrf = hrf, sampling_frame = sampling_frame, ...)
  if (output == "matrix") return(dm)

  # 2. Derive condition tags & column groups -----------------------------------
  base_tags <- fmrireg::conditions(x, expand_basis = FALSE, drop.empty = FALSE)
  if (length(base_tags) == 0L || ncol(dm) == 0L) return(list())

  nb <- fmrihrf::nbasis(hrf)
  term_tag <- attr(x, "term_tag")
  if (is.null(term_tag) && nzchar(x$varname)) {
    term_tag <- fmrireg::sanitize(x$varname, allow_dot = FALSE)
  }

  cols_by_cond_matrix <- vapply(base_tags, function(ct) {
    make_column_names(term_tag, ct, nb)
  }, character(nb))
  
  # Convert matrix to list by column, with proper names
  cols_by_cond <- lapply(seq_len(ncol(cols_by_cond_matrix)), function(i) {
    cols_by_cond_matrix[, i]
  })
  names(cols_by_cond) <- colnames(cols_by_cond_matrix)

  # 3. Split into list, silently dropping incomplete sets -----------------------
  valid <- vapply(cols_by_cond, function(nms) all(nms %in% colnames(dm)), logical(1))
  if (!any(valid)) return(list())
  cols_by_cond <- cols_by_cond[valid]

  lapply(cols_by_cond, function(nms) dm[, nms, drop = FALSE])
}
