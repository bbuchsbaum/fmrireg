#' Create fMRI linear model configuration
#'
#' @description
#' Wrapper function to create configuration for fMRI linear models.
#' This provides a simpler interface to fmri_lm_control.
#'
#' @param robust Logical or character. FALSE for no robust fitting,
#'   "huber" or "bisquare" for robust methods.
#' @param ar_options List with AR options including:
#'   - cor_struct: "none", "ar1", "ar2", etc.
#'   - iter: Number of iterations
#' @param weights Optional weights vector
#' @param method Solving method ("auto", "ols", "ar", "robust", "ar_robust")
#'
#' @return An fmri_lm_config object
#' @keywords internal
#' @noRd
fmri_lm_config <- function(robust = FALSE, 
                          ar_options = NULL,
                          weights = NULL,
                          method = "auto") {
  
  # Convert simplified robust spec
  if (is.logical(robust)) {
    robust_opts <- list(type = ifelse(robust, "huber", FALSE))
  } else if (is.character(robust)) {
    robust_opts <- list(type = robust)
  } else {
    robust_opts <- robust
  }
  
  # Convert simplified AR spec
  if (is.null(ar_options)) {
    ar_opts <- list(struct = "iid")
  } else {
    # Map cor_struct to struct
    if (!is.null(ar_options$cor_struct)) {
      ar_opts <- ar_options
      ar_opts$struct <- switch(ar_options$cor_struct,
        "none" = "iid",
        "ar1" = "ar1",
        "ar2" = "ar2",
        "ar3" = "arp",
        "ar4" = "arp",
        "iid"
      )
      
      # Set p for arp
      if (ar_opts$struct == "arp" && is.null(ar_opts$p)) {
        ar_opts$p <- as.numeric(gsub("ar", "", ar_options$cor_struct))
      }
      
      # Map iter to iter_gls
      if (!is.null(ar_options$iter)) {
        ar_opts$iter_gls <- ar_options$iter
      }
    } else {
      ar_opts <- ar_options
    }
  }
  
  # Create config
  cfg <- fmri_lm_control(
    robust_options = robust_opts,
    ar_options = ar_opts
  )
  
  # Add extra fields
  cfg$weights <- weights
  cfg$method <- method
  
  cfg
}