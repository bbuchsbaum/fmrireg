# Shared frame fixtures and accessors for the test suite.
#
# fmrireg fits models to fmridataset::fmri_frame objects. Tests build them
# with fmridataset's constructors (as_fmri_frame() and
# read_nifti_bold(), re-exported by fmrireg) and read them back through the public
# fmridataset API; the three accessors below just keep that one line long.

# Event table stored on a frame, as a data frame (empty when absent).
frame_events <- function(x) {
  fmrireg:::.dset_event_table(x)
}

# Dense observation-by-feature matrix of the active assay.
frame_data <- function(x) {
  fmridataset::collect_assay(x)
}

# fmrihrf sampling frame implied by the frame's run_id / TR columns.
frame_sframe <- function(x) {
  fmridataset::as_sampling_frame(x)
}

# A latent frame from list(basis, loadings, mask), the form only the deprecated
# latent_frame() accepts; built here without the deprecation warning.
latent_list_frame <- function(x, ..., events = NULL) {
  .latent_frame_impl(x, ..., event_table = events)
}
