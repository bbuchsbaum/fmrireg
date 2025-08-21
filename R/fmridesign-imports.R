# Re-exports from fmridesign for backward compatibility

#' @importFrom fmridesign event_model baseline_model event_factor event_variable event_matrix event_basis condition_basis_list
#' @importFrom fmridesign covariate nuisance hrf design_matrix conditions cells convolve sanitize onsets durations
#' @importFrom fmridesign elements blockids event_terms baseline_terms term_matrices is_categorical is_continuous
#' @importFrom fmridesign longnames shortnames columns term_indices split_by_block event_table event_term construct
#' @importFrom fmridesign convolve_design sub_basis basis_suffix feature_suffix trialwise
#' @importFrom fmridesign Ident Poly BSpline Scale ScaleWithin Standardized RobustScale baseline block
#' @importFrom fmridesign contrast pair_contrast column_contrast unit_contrast oneway_contrast poly_contrast
#' @importFrom fmridesign interaction_contrast one_against_all_contrast contrast_set contrast_weights plot_contrasts design_map
NULL

# Re-export key functions from fmridesign
#' @export
fmridesign::event_model

#' @export
fmridesign::baseline_model

#' @export
fmridesign::event_factor

#' @export
fmridesign::event_variable

#' @export
fmridesign::event_matrix

#' @export
fmridesign::event_basis

#' @export
fmridesign::condition_basis_list

#' @export
fmridesign::covariate

#' @export
fmridesign::nuisance

#' @export
fmridesign::hrf

#' @export
fmridesign::design_matrix

#' @export
fmridesign::conditions

#' @export
fmridesign::cells

#' @export
fmridesign::convolve

#' @export
fmridesign::sanitize

#' @export
fmridesign::onsets

#' @export
fmridesign::durations

#' @export
fmridesign::elements

#' @export
fmridesign::blockids

#' @export
fmridesign::event_terms

#' @export
fmridesign::baseline_terms

#' @export
fmridesign::term_matrices

#' @export
fmridesign::is_categorical

#' @export
fmridesign::is_continuous

#' @export
fmridesign::longnames

#' @export
fmridesign::shortnames

#' @export
fmridesign::columns

#' @export
fmridesign::term_indices

#' @export
fmridesign::split_by_block

#' @export
fmridesign::event_table

#' @export
fmridesign::event_term

#' @export
fmridesign::construct

#' @export
fmridesign::convolve_design

#' @export
fmridesign::sub_basis

#' @export
fmridesign::basis_suffix

#' @export
fmridesign::feature_suffix

#' @export
fmridesign::trialwise

# Re-export transformation functions
#' @export
fmridesign::Ident

#' @export
fmridesign::Poly

#' @export
fmridesign::BSpline

#' @export
fmridesign::Scale

#' @export
fmridesign::ScaleWithin

#' @export
fmridesign::Standardized

#' @export
fmridesign::RobustScale

#' @export
fmridesign::baseline

#' @export
fmridesign::block

# Re-export contrast functions
#' @export
fmridesign::contrast

#' @export
fmridesign::pair_contrast

#' @export
fmridesign::column_contrast

#' @export
fmridesign::unit_contrast

#' @export
fmridesign::oneway_contrast

#' @export
fmridesign::poly_contrast

#' @export
fmridesign::interaction_contrast

#' @export
fmridesign::one_against_all_contrast

#' Create a contrast set
#' 
#' @param ... contrast specifications
#' @export
contrast_set <- function(...) {
  fmridesign::contrast_set(...)
}

#' @export
fmridesign::contrast_weights

#' @export
fmridesign::plot_contrasts

#' @export
fmridesign::design_map