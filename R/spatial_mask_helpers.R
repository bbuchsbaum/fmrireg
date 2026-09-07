#' Resolve Dataset Mask and Spatial Metadata
#'
#' Recovers the neuroim2 spatial reference of an `fmri_frame` whose feature
#' space is a `volume_space` (directly, or as the parent of a `basis_space`).
#' The reconstruction goes through `fmridataset::reconstruct_space()`, so the
#' mask and `NeuroSpace` agree exactly with what fmridataset itself would use
#' for `spatial_map()`.
#'
#' @return A list with `mask` (a `LogicalNeuroVol`), `mask_array` (3D logical
#'   array), and `space` (a 3D `NeuroSpace`).
#' @keywords internal
#' @noRd
.fmri_dataset_mask_space <- function(dataset, context = "spatial reconstruction") {
  if (!.is_fmri_frame(dataset)) {
    stop("Expected an 'fmri_frame' object for ", context, ".", call. = FALSE)
  }

  vspace <- .dset_volume_space(dataset)
  if (is.null(vspace)) {
    stop(
      "Spatial metadata is required for ", context,
      ", but the dataset's feature space is not a volume_space ",
      "(build the frame with neurovec_frame(), nifti_frame(), or a volume_space).",
      call. = FALSE
    )
  }

  n_feat <- fmridataset::n_features(vspace)
  vol <- fmridataset::reconstruct_space(vspace, rep(1, n_feat))
  space <- neuroim2::space(vol)
  mask_array <- !is.na(as.array(vol))
  mask_dim <- dim(mask_array)
  if (length(mask_dim) != 3L) {
    stop("Dataset mask must be a 3D array for ", context, ".", call. = FALSE)
  }
  space <- .fmri_normalize_mask_space(space, mask_dim, context)
  mask <- neuroim2::LogicalNeuroVol(mask_array, space)
  list(mask = mask, mask_array = mask_array, space = space)
}

#' @keywords internal
#' @noRd
.fmri_try_space <- function(x) {
  if (is.null(x) || is.character(x)) return(NULL)
  tryCatch(neuroim2::space(x), error = function(e) NULL)
}

#' @keywords internal
#' @noRd
.fmri_normalize_mask_space <- function(space, mask_dim = NULL, context) {
  space_dim <- dim(space)
  if (length(space_dim) > 3L) {
    space <- tryCatch(
      neuroim2::drop_dim(space, seq.int(4L, length(space_dim))),
      error = function(e) NULL
    )
    if (is.null(space)) {
      stop("Could not reduce 4D spatial metadata to a 3D mask space for ", context, ".", call. = FALSE)
    }
    space_dim <- dim(space)
  }

  if (!is.null(mask_dim) && !identical(as.integer(space_dim[1:3]), as.integer(mask_dim))) {
    stop(
      "Mask dimensions (", paste(mask_dim, collapse = " x "),
      ") do not match spatial metadata dimensions (",
      paste(space_dim[1:3], collapse = " x "), ") for ", context, ".",
      call. = FALSE
    )
  }

  space
}
