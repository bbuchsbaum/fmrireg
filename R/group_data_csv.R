# Group Data Support for CSV/Tabular Data
# Handles ROI summaries and pre-extracted statistics for group analysis

#' Create Group Dataset from CSV File or Data Frame
#'
#' Creates a group dataset from tabular data containing pre-extracted statistics
#' such as ROI means, effect sizes, and standard errors. This format is useful
#' for ROI-based analyses or when working with summary statistics.
#'
#' @param data Either a path to a CSV file or a data frame containing the data
#' @param effect_cols Named vector or list specifying column names for effect statistics.
#'   E.g., c(beta = "mean_activation", se = "std_error") or c(t = "t_stat", df = "df")
#' @param subject_col Character string specifying the column containing subject IDs
#' @param roi_col Character string specifying the column containing ROI names (optional)
#' @param contrast_col Character string specifying the column containing contrast names (optional)
#' @param covariate_cols Character vector of column names to use as covariates (optional)
#' @param wide_format Logical. If TRUE, expects wide format with ROIs as columns (default: FALSE)
#'
#' @return A group_data_csv object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Long format: one row per subject-ROI combination
#' gd <- group_data_from_csv(
#'   "roi_statistics.csv",
#'   effect_cols = c(beta = "mean_beta", se = "se_beta"),
#'   subject_col = "participant_id",
#'   roi_col = "roi_name",
#'   covariate_cols = c("age", "sex", "group")
#' )
#'
#' # Wide format: one row per subject, ROIs as columns
#' gd <- group_data_from_csv(
#'   "subject_summary.csv",
#'   effect_cols = c(beta = "roi_"),  # Prefix for ROI columns
#'   subject_col = "subject",
#'   wide_format = TRUE
#' )
#'
#' # From data frame with multiple contrasts
#' df <- read.csv("contrast_results.csv")
#' gd <- group_data_from_csv(
#'   df,
#'   effect_cols = c(beta = "estimate", se = "std_error", t = "t_value"),
#'   subject_col = "subject_id",
#'   contrast_col = "contrast_name",
#'   roi_col = "region"
#' )
#' }
group_data_from_csv <- function(data,
                                effect_cols,
                                subject_col = "subject",
                                roi_col = NULL,
                                contrast_col = NULL,
                                covariate_cols = NULL,
                                wide_format = FALSE) {
  
  # Read data if it's a file path
  if (is.character(data) && length(data) == 1) {
    if (!file.exists(data)) {
      stop("File does not exist: ", data, call. = FALSE)
    }
    df <- read.csv(data, stringsAsFactors = FALSE)
  } else if (is.data.frame(data)) {
    df <- data
  } else {
    stop("'data' must be a file path or data frame", call. = FALSE)
  }
  
  # Validate required columns exist
  if (!subject_col %in% names(df)) {
    stop("Subject column '", subject_col, "' not found in data", call. = FALSE)
  }
  
  # Validate effect columns
  effect_cols <- validate_effect_cols(effect_cols, df)
  
  # Process wide format if specified
  if (wide_format) {
    df <- wide_to_long_format(df, effect_cols, subject_col)
    roi_col <- "roi"  # Created by wide_to_long_format
  }
  
  # Validate optional columns
  if (!is.null(roi_col) && !roi_col %in% names(df)) {
    stop("ROI column '", roi_col, "' not found in data", call. = FALSE)
  }
  
  if (!is.null(contrast_col) && !contrast_col %in% names(df)) {
    stop("Contrast column '", contrast_col, "' not found in data", call. = FALSE)
  }
  
  if (!is.null(covariate_cols)) {
    missing_cols <- setdiff(covariate_cols, names(df))
    if (length(missing_cols) > 0) {
      stop("Covariate columns not found: ", paste(missing_cols, collapse = ", "), 
           call. = FALSE)
    }
  }
  
  # Extract unique subjects
  subjects <- unique(df[[subject_col]])
  n_subjects <- length(subjects)
  
  # Extract unique ROIs if present
  rois <- NULL
  n_rois <- 1
  if (!is.null(roi_col)) {
    rois <- unique(df[[roi_col]])
    n_rois <- length(rois)
  }
  
  # Extract unique contrasts if present
  contrasts <- NULL
  if (!is.null(contrast_col)) {
    contrasts <- unique(df[[contrast_col]])
  }
  
  # Build covariates data frame
  covariates <- NULL
  if (!is.null(covariate_cols)) {
    # Get unique subject rows for covariates
    subject_rows <- !duplicated(df[[subject_col]])
    covariates <- df[subject_rows, covariate_cols, drop = FALSE]
    rownames(covariates) <- df[[subject_col]][subject_rows]
    
    # Ensure order matches subjects vector
    covariates <- covariates[subjects, , drop = FALSE]
  }
  
  # Validate data completeness
  expected_rows <- n_subjects
  if (!is.null(roi_col)) expected_rows <- expected_rows * n_rois
  if (!is.null(contrast_col)) expected_rows <- expected_rows * length(contrasts)
  
  if (nrow(df) != expected_rows) {
    warning("Data may be incomplete. Expected ", expected_rows, " rows but found ", nrow(df))
  }
  
  # Create group_data_csv object
  structure(
    list(
      data = df,
      effect_cols = effect_cols,
      subject_col = subject_col,
      roi_col = roi_col,
      contrast_col = contrast_col,
      covariate_cols = covariate_cols,
      subjects = subjects,
      rois = rois,
      contrasts = contrasts,
      covariates = covariates,
      format = "csv",
      n_subjects = n_subjects,
      n_rois = n_rois,
      n_contrasts = length(contrasts)
    ),
    class = c("group_data_csv", "group_data")
  )
}

#' Validate Effect Column Specification
#'
#' @param effect_cols Named vector or list of effect column names
#' @param df Data frame to validate against
#' @return Validated effect_cols list
#' @keywords internal
#' @noRd
validate_effect_cols <- function(effect_cols, df) {
  # Convert to list if needed
  if (!is.list(effect_cols)) {
    effect_cols <- as.list(effect_cols)
  }
  
  # Check that names are provided
  if (is.null(names(effect_cols))) {
    stop("'effect_cols' must be a named vector or list (e.g., c(beta = 'column_name'))", 
         call. = FALSE)
  }
  
  # Validate column names exist
  for (effect_type in names(effect_cols)) {
    col_name <- effect_cols[[effect_type]]
    
    # Handle column prefixes for wide format
    if (!col_name %in% names(df)) {
      # Check if it's a prefix
      matching_cols <- grep(paste0("^", col_name), names(df), value = TRUE)
      if (length(matching_cols) == 0) {
        stop("Effect column '", col_name, "' (for ", effect_type, ") not found in data", 
             call. = FALSE)
      }
    }
  }
  
  # Validate effect type combinations
  has_beta <- "beta" %in% names(effect_cols)
  has_se <- "se" %in% names(effect_cols)
  has_var <- "var" %in% names(effect_cols)
  has_t <- "t" %in% names(effect_cols)
  
  if (has_beta && !has_se && !has_var) {
    stop("When 'beta' is provided, must also provide 'se' or 'var'", call. = FALSE)
  }
  
  if (has_t && !"df" %in% names(effect_cols)) {
    warning("T-statistics provided without degrees of freedom. Meta-analysis options will be limited.")
  }
  
  return(effect_cols)
}

#' Convert Wide Format to Long Format
#'
#' @param df Data frame in wide format
#' @param effect_cols Effect column specification
#' @param subject_col Subject column name
#' @return Data frame in long format
#' @keywords internal
wide_to_long_format <- function(df, effect_cols, subject_col) {
  # This is a simplified implementation
  # In practice, might want to use tidyr::pivot_longer or similar
  
  stop("Wide format conversion not yet implemented. Please reshape your data to long format.", 
       call. = FALSE)
}

#' Extract Data for Meta-Analysis from CSV
#'
#' @param gd A group_data_csv object
#' @param roi Optional ROI name to extract
#' @param contrast Optional contrast name to extract
#' @return List with effect sizes and variances
#' @examples
#' gd <- fmrireg:::.demo_group_data_csv()
#' extract_csv_data(gd, roi = "ROI1")
#' @export
extract_csv_data <- function(gd, roi = NULL, contrast = NULL) {
  if (!inherits(gd, "group_data_csv")) {
    stop("Input must be a group_data_csv object", call. = FALSE)
  }
  
  df <- gd$data
  
  # Filter by ROI if specified
  if (!is.null(roi)) {
    if (is.null(gd$roi_col)) {
      stop("Cannot filter by ROI: no ROI column specified", call. = FALSE)
    }
    df <- df[df[[gd$roi_col]] == roi, ]
    if (nrow(df) == 0) {
      stop("No data found for ROI: ", roi, call. = FALSE)
    }
  }
  
  # Filter by contrast if specified
  if (!is.null(contrast)) {
    if (is.null(gd$contrast_col)) {
      stop("Cannot filter by contrast: no contrast column specified", call. = FALSE)
    }
    df <- df[df[[gd$contrast_col]] == contrast, ]
    if (nrow(df) == 0) {
      stop("No data found for contrast: ", contrast, call. = FALSE)
    }
  }
  
  # Extract effect statistics
  result <- list()
  
  # Extract beta/effect sizes
  if ("beta" %in% names(gd$effect_cols)) {
    col_name <- gd$effect_cols[["beta"]]
    result$beta <- df[[col_name]]
  }
  
  # Extract standard errors
  if ("se" %in% names(gd$effect_cols)) {
    col_name <- gd$effect_cols[["se"]]
    result$se <- df[[col_name]]
    result$var <- result$se^2
  }
  
  # Extract variances
  if ("var" %in% names(gd$effect_cols)) {
    col_name <- gd$effect_cols[["var"]]
    result$var <- df[[col_name]]
    if (is.null(result$se)) {
      result$se <- sqrt(result$var)
    }
  }
  
  # Extract t-statistics
  if ("t" %in% names(gd$effect_cols)) {
    col_name <- gd$effect_cols[["t"]]
    result$t <- df[[col_name]]
  }
  
  # Extract degrees of freedom
  if ("df" %in% names(gd$effect_cols)) {
    col_name <- gd$effect_cols[["df"]]
    result$df <- df[[col_name]]
  }
  
  # Add subject information
  result$subjects <- df[[gd$subject_col]]
  
  return(result)
}

#' Get Available ROIs
#'
#' @param gd A group_data_csv object
#' @return Character vector of ROI names
#' @examples
#' gd <- fmrireg:::.demo_group_data_csv()
#' get_rois(gd)
#' @export
get_rois <- function(gd) {
  if (!inherits(gd, "group_data_csv")) {
    stop("Input must be a group_data_csv object", call. = FALSE)
  }
  return(gd$rois)
}

#' Get Available Contrasts
#'
#' @param gd A group_data_csv object
#' @return Character vector of contrast names
#' @examples
#' gd <- fmrireg:::.demo_group_data_csv()
#' get_contrasts(gd)
#' @export
get_contrasts <- function(gd) {
  if (!inherits(gd, "group_data_csv")) {
    stop("Input must be a group_data_csv object", call. = FALSE)
  }
  return(gd$contrasts)
}

#' Validate group_data_csv Object
#'
#' @param x A group_data_csv object
#' @keywords internal
#' @noRd
validate_group_data_csv <- function(x) {
  # Check for required fields
  required_fields <- c("data", "effect_cols", "subject_col", "subjects", "format")
  missing_fields <- setdiff(required_fields, names(x))
  if (length(missing_fields) > 0) {
    stop("Missing required fields for group_data_csv: ",
         paste(missing_fields, collapse = ", "), call. = FALSE)
  }
  
  # Validate data is a data frame
  if (!is.data.frame(x$data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }
  
  # Validate effect_cols structure
  if (!is.list(x$effect_cols) || is.null(names(x$effect_cols))) {
    stop("'effect_cols' must be a named list", call. = FALSE)
  }
  
  invisible(TRUE)
}
