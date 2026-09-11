# Shared frame fixtures and accessors for the test suite.
#
# fmrireg fits models to fmridataset::fmri_frame objects. Tests build them
# through the package's own constructors (matrix_frame(), neurovec_frame(),
# latent_frame(), nifti_frame()) and read them back through the public
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
