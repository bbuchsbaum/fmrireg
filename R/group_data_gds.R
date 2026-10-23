#' Group-level data via fmrigds (canonical)
#'
#' Entry point compatible with existing \code{group_data()} API. This always
#' constructs a gds-backed object via fmrigds and returns it with class
#' \code{c("group_data_gds","group_data", ...)} for S3 compatibility.
#'
#' @param data See legacy \code{group_data()}.
#' @param format One of "auto","h5","nifti","csv","fmrilm".
#' @param ... Passed to the underlying constructor.
#' @return A gds-backed group-level data object (class includes 'group_data_gds').
#' @noRd
# There is no in-memory route from fitted fmri_lm objects to a gds object yet;
# per-subject maps written by write_results() are the supported input.
.group_data_fmrilm_unsupported_msg <- paste(
  "group_data(format = \"fmrilm\") is not supported.",
  "Save each subject's fit with write_results() (e.g. format = \"nifti\" with",
  "beta and se maps) and load the files with",
  "group_data(list(beta = beta_paths, se = se_paths), format = \"nifti\")",
  "or, for HDF5 output, group_data(h5_paths, format = \"h5\")."
)

.is_fmri_lm_list <- function(x) {
  is.list(x) && !is.data.frame(x) && length(x) > 0L &&
    all(vapply(x, inherits, logical(1), what = "fmri_lm"))
}

.group_data_gds <- function(data, format = c("auto","h5","nifti","csv","fmrilm"), ...) {
  format <- match.arg(format)
  if (identical(format, "fmrilm") ||
      (identical(format, "auto") && .is_fmri_lm_list(data))) {
    stop(.group_data_fmrilm_unsupported_msg, call. = FALSE)
  }
  dots <- list(...)
  fmrireg_mask <- dots$mask %||% NULL
  legacy_sample_labels <- NULL
  legacy_feature_group <- NULL
  # Map legacy names to fmrigds adapter ids
  fformat <- switch(format,
                   csv = "tabular",
                   format)

  src <- data

  # For in-memory tabular data, prefer direct passing; fall back to temp CSV if needed
  if (fformat == "tabular" && is.data.frame(data)) {
    # Work on a local copy we can normalize for fmrigds
    data_local <- data
    # Map ROI/parcel column to fmrigds' sample axis when provided in legacy arg
    if (!is.null(dots$roi_col) && is.null(dots$sample_col)) {
      dots$sample_col <- dots$roi_col
      dots$roi_col <- NULL
    }
    # A table with no ROI/sample column is a single ROI, "overall", as the
    # legacy csv reader treated it.
    if (is.null(dots$sample_col) && !("sample" %in% names(data_local))) {
      data_local$sample <- "overall"
      dots$sample_col <- "sample"
    }
    # Capture sample labels (e.g., ROI names) to preserve as rownames later
    sample_labels <- NULL
    if (!is.null(dots$sample_col) && dots$sample_col %in% names(data_local)) {
      sc <- dots$sample_col
      sample_labels <- as.character(data_local[[sc]][!duplicated(data_local[[sc]])])
      legacy_sample_labels <- sample_labels
    }
    if (is.null(dots$row_data) && !is.null(sample_labels)) {
      row_df <- data.frame(sample = sample_labels, row.names = sample_labels, check.names = FALSE)
      if (!is.null(dots$feature_group)) {
        feature_group <- dots$feature_group
        legacy_feature_group <- feature_group
        if (length(feature_group) != length(sample_labels)) {
          stop("feature_group must have length equal to number of samples", call. = FALSE)
        }
        row_df$feature_group <- feature_group
        dots$feature_group <- NULL
      }
      dots$row_data <- row_df
    } else if (!is.null(dots$row_data) && !is.null(sample_labels) && !"sample" %in% names(dots$row_data)) {
      dots$row_data$sample <- sample_labels
      if (is.null(rownames(dots$row_data))) {
        rownames(dots$row_data) <- sample_labels
      }
    }
    if (is.null(legacy_sample_labels) && !is.null(dots$row_data)) {
      if ("sample" %in% names(dots$row_data)) {
        legacy_sample_labels <- as.character(dots$row_data$sample)
      } else if (!is.null(rownames(dots$row_data))) {
        legacy_sample_labels <- rownames(dots$row_data)
      }
    }
    if (is.null(legacy_feature_group) && !is.null(dots$row_data) && "feature_group" %in% names(dots$row_data)) {
      legacy_feature_group <- dots$row_data$feature_group
    }
    # Build col_data from provided covariate_cols if present
    if (!is.null(dots$covariate_cols)) {
      covar_cols <- dots$covariate_cols
      subj_col <- dots$subject_col %||% "subject"
      if (all(covar_cols %in% names(data)) && subj_col %in% names(data)) {
        subj_rows <- !duplicated(data[[subj_col]])
        col_df <- data[subj_rows, c(subj_col, covar_cols), drop = FALSE]
        rownames(col_df) <- col_df[[subj_col]]
        dots$col_data <- col_df[, covar_cols, drop = FALSE]
      }
      dots$covariate_cols <- NULL
    }

    # Ensure a contrast column exists for tabular adapter (use single default level)
    if (is.null(dots$contrast_col)) {
      if (!("contrast" %in% names(data_local))) {
        data_local$contrast <- "contrast1"
      }
      dots$contrast_col <- "contrast"
    }

    # Try direct data.frame first
    gd_try <- try(do.call(fmrigds::gds, c(list(source = data_local, format = fformat), dots)), silent = TRUE)
    if (!inherits(gd_try, "try-error")) {
      return(.annotate_group_data_gds(gd_try, legacy_sample_labels, legacy_feature_group))
    }
    # Try list(data=df) variant
    gd_try2 <- try(do.call(fmrigds::gds, c(list(source = list(data = data_local), format = fformat), dots)), silent = TRUE)
    if (!inherits(gd_try2, "try-error")) {
      return(.annotate_group_data_gds(gd_try2, legacy_sample_labels, legacy_feature_group))
    }
    # Fallback: temp CSV path
    tmpfile <- tempfile(fileext = ".csv")
    utils::write.csv(data_local, tmpfile, row.names = FALSE)
    src <- tmpfile
  }

  # NIfTI sources: list(beta=, se= | var=) or list(t=, df=), also with the
  # legacy *_paths names. fmrigds's NIfTI adapter reads all of them; names are
  # passed exactly (its source list is read with `$`, which would otherwise
  # partially match "beta" to "beta_paths").
  if (fformat == "nifti" && is.list(data) && !inherits(data, "group_data_nifti")) {
    src <- .group_data_nifti_source(
      beta = data[["beta"]] %||% data[["beta_paths"]],
      se = data[["se"]] %||% data[["se_paths"]],
      var = data[["var"]] %||% data[["var_paths"]],
      t = data[["t"]] %||% data[["t_paths"]],
      df = data[["df"]]
    )
  }

  # Translate the legacy public argument name to fmrigds metadata so formulas
  # such as ~ 1 + group have a real subject-level design matrix.
  if (fformat == "nifti" && !is.null(dots$covariates) && is.null(dots$col_data)) {
    dots$col_data <- as.data.frame(dots$covariates)
    if (!is.null(dots$subjects) &&
        nrow(dots$col_data) == length(dots$subjects)) {
      rownames(dots$col_data) <- as.character(dots$subjects)
    }
    dots$covariates <- NULL
  }

  # LEGACY-GROUP-DELETE (next minor): this conversion of legacy group_data_nifti objects.
  # For nifti when passed a legacy group_data_nifti object, adapt its fields
  if (fformat == "nifti" && inherits(data, "group_data_nifti")) {
    src <- .group_data_nifti_source(
      beta = data$beta_paths, se = data$se_paths, var = data$var_paths,
      t = data$t_paths, df = data$df
    )
  }

  gd <- do.call(fmrigds::gds, c(list(source = src, format = fformat), dots))
  gd <- .annotate_group_data_gds(gd, legacy_sample_labels, legacy_feature_group)
  if (!is.null(fmrireg_mask)) attr(gd, "fmrireg_mask") <- fmrireg_mask
  gd
}

.annotate_group_data_gds <- function(gd, sample_labels = NULL, feature_group = NULL) {
  if (!is.null(sample_labels)) {
    attr(gd, "fmrireg_sample_labels") <- as.character(sample_labels)
  }
  if (!is.null(feature_group)) {
    attr(gd, "fmrireg_feature_group") <- feature_group
  }
  class(gd) <- c("group_data_gds", "group_data", class(gd))
  gd
}

# Basic S3 bridges to keep helpers working on gds-backed objects
.n_subjects_group_data_gds <- function(x) {
  x <- if (inherits(x, "gds_plan")) fmrigds::compute(x) else x
  length(fmrigds::subjects(x))
}

.get_subjects_group_data_gds <- function(x) {
  x <- if (inherits(x, "gds_plan")) fmrigds::compute(x) else x
  fmrigds::subjects(x)
}

.get_covariates_group_data_gds <- function(x) {
  x <- if (inherits(x, "gds_plan")) fmrigds::compute(x) else x
  fmrigds::col_data(x)
}

.print_group_data_gds <- function(x, ...) {
  cat("group_data (backed by gds)\n")
  if (inherits(x, "gds_plan")) {
    try(fmrigds::explain(x), silent = TRUE)
  } else {
    assays <- .gds_safe_assay_names(x)
    if (length(assays) > 0) {
      cat(" assays: ", paste(assays, collapse = ", "), "\n", sep = "")
    }
    ns <- try(length(fmrigds::subjects(x)), silent = TRUE)
    if (!inherits(ns, "try-error")) {
      cat(" subjects:", ns, "\n")
    }
  }
  invisible(x)
}

.group_data_nifti_source <- function(beta = NULL, se = NULL, var = NULL, t = NULL, df = NULL) {
  if (!is.null(t)) {
    return(list(t = t, df = df))
  }
  src <- list(beta = beta, se = se, var = var)
  src[!vapply(src, is.null, logical(1))]
}
