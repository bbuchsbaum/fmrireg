#' Tidy fitted HRF curves from an fmri_lm fit
#'
#' Convert [`fitted_hrf()`] output into a long tibble suitable for plotting.
#'
#' @param x An `fmri_lm` object.
#' @param sample_at Numeric vector of time points passed to [`fitted_hrf()`].
#' @param term Optional term selector. When `NULL`, all event terms are returned.
#' @param term_match Matching mode for `term`: `"contains"` (default), `"exact"`,
#'   or `"regex"`.
#' @param voxel Integer voxel index to extract from each term's prediction matrix.
#' @param average_voxels Logical; if `TRUE`, average across voxels instead of
#'   selecting one voxel.
#' @param ... Additional arguments passed to [`fitted_hrf()`].
#'
#' @return A tibble with columns `term`, `time`, `condition`, `estimate`, and
#'   `voxel`.
#' @export
tidy_fitted_hrf <- function(
    x,
    sample_at = seq(0, 24, by = 1),
    term = NULL,
    term_match = c("contains", "exact", "regex"),
    voxel = 1L,
    average_voxels = FALSE,
    ...) {
  if (!inherits(x, "fmri_lm")) {
    stop("`x` must be an 'fmri_lm' object.", call. = FALSE)
  }
  term_match <- match.arg(term_match)
  average_voxels <- isTRUE(average_voxels)

  hrf_list <- fitted_hrf(x, sample_at = sample_at, ...)
  if (length(hrf_list) == 0) {
    return(tibble::tibble(
      term = character(),
      time = numeric(),
      condition = character(),
      estimate = numeric(),
      voxel = integer()
    ))
  }

  term_names <- names(hrf_list)
  if (is.null(term_names) || any(!nzchar(term_names))) {
    term_names <- paste0("term_", seq_along(hrf_list))
    names(hrf_list) <- term_names
  }

  keep <- rep(TRUE, length(term_names))
  if (!is.null(term)) {
    keep <- switch(
      term_match,
      contains = grepl(term, term_names, fixed = TRUE),
      exact = term_names %in% term,
      regex = grepl(term, term_names, perl = TRUE)
    )
    if (!any(keep)) {
      stop("No fitted HRF terms matched `term`.", call. = FALSE)
    }
  }

  hrf_list <- hrf_list[keep]
  term_names <- term_names[keep]

  if (!average_voxels) {
    if (!is.numeric(voxel) || length(voxel) != 1L || is.na(voxel)) {
      stop("`voxel` must be a single numeric index when `average_voxels = FALSE`.", call. = FALSE)
    }
    voxel <- as.integer(voxel)
    if (voxel < 1L) {
      stop("`voxel` must be >= 1.", call. = FALSE)
    }
  } else {
    voxel <- NA_integer_
  }

  rows <- lapply(seq_along(hrf_list), function(i) {
    nm <- term_names[[i]]
    obj <- hrf_list[[i]]

    design_tbl <- tibble::as_tibble(obj$design)
    pred <- as.matrix(obj$pred)
    if (nrow(pred) != nrow(design_tbl)) {
      stop(
        sprintf("Prediction/design row mismatch for term '%s'.", nm),
        call. = FALSE
      )
    }

    if (average_voxels) {
      est <- rowMeans(pred)
      voxel_out <- rep(NA_integer_, nrow(pred))
    } else {
      if (voxel > ncol(pred)) {
        stop(
          sprintf("Requested voxel index %d exceeds prediction columns (%d) for term '%s'.", voxel, ncol(pred), nm),
          call. = FALSE
        )
      }
      est <- pred[, voxel]
      voxel_out <- rep(voxel, nrow(pred))
    }

    time_vals <- if ("time" %in% names(design_tbl)) {
      as.numeric(design_tbl$time)
    } else {
      as.numeric(sample_at)[seq_len(nrow(design_tbl))]
    }

    cond_vals <- if ("condition" %in% names(design_tbl)) {
      as.character(design_tbl$condition)
    } else {
      other_cols <- setdiff(names(design_tbl), "time")
      if (length(other_cols) > 0) {
        apply(as.data.frame(design_tbl[, other_cols, drop = FALSE]), 1, function(row) {
          paste(row, collapse = ":")
        })
      } else {
        rep(nm, nrow(design_tbl))
      }
    }

    tibble::tibble(
      term = nm,
      time = time_vals,
      condition = cond_vals,
      estimate = as.numeric(est),
      voxel = voxel_out
    )
  })

  dplyr::bind_rows(rows)
}
