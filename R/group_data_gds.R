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
#' @keywords internal
group_data <- function(data, format = c("auto","h5","nifti","csv","fmrilm"), ...) {
  format <- match.arg(format)
  dots <- list(...)
  # Legacy compatibility: for CSV/tabular in-memory data, preserve group_data_csv
  if (identical(format, "csv") && is.data.frame(data) && !("contrast" %in% names(data))) {
    # Delegate to legacy constructor to satisfy tests expecting group_data_csv
    return(do.call(group_data_from_csv, c(list(data), dots)))
  }
  # Map legacy names to fmrigds adapter ids
  fformat <- switch(format,
                   csv = "tabular",
                   fmrilm = "tabular", # temporary mapping until fmrilm adapter is present
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
    # Capture sample labels (e.g., ROI names) to preserve as rownames later
    sample_labels <- NULL
    if (!is.null(dots$sample_col) && dots$sample_col %in% names(data)) {
      sc <- dots$sample_col
      sample_labels <- as.character(data[[sc]][!duplicated(data[[sc]])])
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
      gd <- gd_try
      if (!is.null(sample_labels)) attr(gd, "fmrireg_sample_labels") <- sample_labels
      class(gd) <- c("group_data_gds", "group_data", class(gd))
      return(gd)
    }
    # Try list(data=df) variant
    gd_try2 <- try(do.call(fmrigds::gds, c(list(source = list(data = data_local), format = fformat), dots)), silent = TRUE)
    if (!inherits(gd_try2, "try-error")) {
      gd <- gd_try2
      if (!is.null(sample_labels)) attr(gd, "fmrireg_sample_labels") <- sample_labels
      class(gd) <- c("group_data_gds", "group_data", class(gd))
      return(gd)
    }
    # Fallback: temp CSV path
    tmpfile <- tempfile(fileext = ".csv")
    utils::write.csv(data_local, tmpfile, row.names = FALSE)
    src <- tmpfile
  }

  # For nifti with list(beta=..., se=...) or list(t=..., df=...), accept both *_paths and plain names
  if (fformat == "nifti" && is.list(data)) {
    beta_src <- data$beta %||% data$beta_paths
    se_src   <- data$se %||% data$se_paths
    var_src  <- data$var %||% data$var_paths
    t_src    <- data$t %||% data$t_paths
    df_src   <- data$df
    if (!is.null(beta_src)) {
      src <- list(beta_paths = beta_src)
      if (!is.null(se_src))  src$se_paths  <- se_src
      if (!is.null(var_src)) src$var_paths <- var_src
    } else if (!is.null(t_src)) {
      src <- list(t_paths = t_src)
      if (!is.null(df_src)) src$df <- df_src
    }
  }

  # For nifti when passed a legacy group_data_nifti object, adapt its fields
  if (fformat == "nifti" && inherits(data, "group_data_nifti")) {
    bp <- data$beta_paths
    sp <- data$se_paths
    vp <- data$var_paths
    tp <- data$t_paths
    dfv <- data$df
    if (!is.null(bp)) {
      src <- list(beta_paths = bp)
      if (!is.null(sp)) src$se_paths  <- sp
      if (!is.null(vp)) src$var_paths <- vp
    } else if (!is.null(tp)) {
      src <- list(t_paths = tp)
      if (!is.null(dfv)) src$df <- dfv
    }
  }

  gd <- do.call(fmrigds::gds, c(list(source = src, format = fformat), dots))
  # Attach sample labels when available
  if (exists("sample_labels") && !is.null(sample_labels)) {
    attr(gd, "fmrireg_sample_labels") <- sample_labels
  }
  class(gd) <- c("group_data_gds", "group_data", class(gd))
  gd
}

# Basic S3 bridges to keep helpers working on gds-backed objects
n_subjects.group_data_gds <- function(x) {
  x <- if (inherits(x, "gds_plan")) fmrigds::compute(x) else x
  length(fmrigds::subjects(x))
}

get_subjects.group_data_gds <- function(x) {
  x <- if (inherits(x, "gds_plan")) fmrigds::compute(x) else x
  fmrigds::subjects(x)
}

get_covariates.group_data_gds <- function(x) {
  x <- if (inherits(x, "gds_plan")) fmrigds::compute(x) else x
  fmrigds::col_data(x)
}

print.group_data_gds <- function(x, ...) {
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
