# Helper constructors for documentation examples
# These functions create tiny objects that can be reused across roxygen examples
# to keep the snippets short and avoid heavyweight setup.

#' Internal: demo event data for examples
#' @keywords internal
#' @noRd
.demo_event_data <- function() {
  data.frame(
    onsets = c(0, 4, 0, 4),
    condition = factor(rep(c("A", "B"), each = 2)),
    run = rep(1:2, each = 2)
  )
}

#' Internal: demo sampling frame with two short runs
#' @keywords internal
#' @noRd
.demo_sampling_frame <- function() {
  fmridesign::sampling_frame(blocklens = c(4, 4), TR = 2)
}

#' Internal: demo event model for examples
#' @keywords internal
#' @noRd
.demo_event_model <- function() {
  fmridesign::event_model(
    onsets ~ hrf(condition),
    data = .demo_event_data(),
    block = ~run,
    sampling_frame = .demo_sampling_frame()
  )
}

#' Internal: demo matrix_dataset with two voxels
#' @keywords internal
#' @noRd
.demo_matrix_dataset <- function() {
  signals <- matrix(seq_len(16), nrow = 8, ncol = 2)
  fmridataset::matrix_dataset(
    signals,
    TR = 2,
    run_length = c(4, 4),
    event_table = .demo_event_data()
  )
}

#' Internal: demo fmri_model object
#' @keywords internal
#' @noRd
.demo_fmri_model <- function() {
  create_fmri_model(
    formula = onsets ~ hrf(condition),
    block = ~run,
    dataset = .demo_matrix_dataset()
  )
}

#' Internal: demo fmri_lm fit
#' @keywords internal
#' @noRd
.demo_fmri_lm <- function() {
  base <- baseline_model(basis = "poly", degree = 1, sframe = .demo_sampling_frame())
  fmri_lm(
    formula = onsets ~ hrf(condition),
    block = ~run,
    dataset = .demo_matrix_dataset(),
    baseline_model = base,
    progress = FALSE
  )
}

#' Internal: demo group_data_csv object
#' @keywords internal
#' @noRd
.demo_group_data_csv <- function() {
  df <- data.frame(
    subject = rep(paste0("s", 1:3), each = 2),
    roi = rep(c("ROI1", "ROI2"), times = 3),
    contrast = "A_vs_B",
    age = rep(c(30, 32, 34), each = 2),
    beta = c(0.2, 0.1, 0.3, 0.15, 0.25, 0.2),
    se = rep(0.1, 6)
  )
  group_data_from_csv(
    df,
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    contrast_col = "contrast",
    covariate_cols = "age"
  )
}

#' Internal: demo fmri_meta fit using CSV group data
#' @keywords internal
#' @noRd
.demo_fmri_meta <- function() {
  fmri_meta(
    data = .demo_group_data_csv(),
    formula = ~ 1,
    method = "fe",
    robust = "none",
    verbose = FALSE
  )
}

#' Internal: demo fmri_ttest fit from small matrices
#' @keywords internal
#' @noRd
.demo_fmri_ttest_fit <- function() {
  coeff_names <- c("conditionA", "conditionB")
  structure(
    list(
      beta = matrix(c(0.2, -0.1), nrow = 2, ncol = 1,
                    dimnames = list(coeff_names, "voxel1")),
      se = matrix(0.1, nrow = 2, ncol = 1,
                  dimnames = list(coeff_names, "voxel1")),
      t = matrix(c(2, -1), nrow = 2, ncol = 1,
                 dimnames = list(coeff_names, "voxel1")),
      z = matrix(c(2.1, -1.05), nrow = 2, ncol = 1,
                 dimnames = list(coeff_names, "voxel1")),
      p = matrix(c(0.04, 0.32), nrow = 2, ncol = 1,
                 dimnames = list(coeff_names, "voxel1")),
      df = 8,
      engine = "classic",
      formula = ~1,
      n_subjects = 5,
      n_features = 1
    ),
    class = "fmri_ttest_fit"
  )
}
