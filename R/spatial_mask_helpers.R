#' Resolve Dataset Mask and Spatial Metadata
#' @keywords internal
#' @noRd
.fmri_dataset_mask_space <- function(dataset, context = "spatial reconstruction") {
  if (!inherits(dataset, "fmri_dataset")) {
    stop("Expected an 'fmri_dataset' object for ", context, ".", call. = FALSE)
  }

  mask <- fmridataset::get_mask(dataset)
  mask_raw <- as.array(mask)
  space <- .fmri_mask_space_from_sources(mask, dataset)
  if (is.null(space)) {
    stop(
      "Spatial metadata is required for ", context,
      ", but the dataset mask does not carry a neuroim2::NeuroSpace and ",
      "no recoverable backend spatial metadata was found.",
      call. = FALSE
    )
  }

  mask_dim <- dim(mask_raw)
  if (is.null(mask_dim)) {
    mask_dim <- dim(.fmri_normalize_mask_space(space, NULL, context))[1:3]
  }
  if (length(mask_dim) != 3L) {
    stop("Dataset mask must be a 3D array for ", context, ".", call. = FALSE)
  }

  mask_array <- array(as.logical(mask_raw), dim = mask_dim)
  space <- .fmri_normalize_mask_space(space, dim(mask_array), context)
  list(mask = mask, mask_array = mask_array, space = space)
}

#' @keywords internal
#' @noRd
.fmri_mask_space_from_sources <- function(mask, dataset) {
  space <- .fmri_try_space(mask)
  if (!is.null(space)) return(space)

  backend <- dataset$backend
  if (!is.null(backend)) {
    metadata <- tryCatch(fmridataset::backend_get_metadata(backend), error = function(e) NULL)
    if (!is.null(metadata$space)) return(metadata$space)

    for (field in c("mask", "mask_source", "data")) {
      obj <- backend[[field]]
      space <- .fmri_try_space(obj)
      if (!is.null(space)) return(space)
    }

    source <- backend$source
    if (is.list(source) && length(source) > 0L) {
      space <- .fmri_try_space(source[[1]])
      if (!is.null(space)) return(space)
    } else {
      space <- .fmri_try_space(source)
      if (!is.null(space)) return(space)
    }
  }

  for (field in c("mask", "mask_source", "data")) {
    obj <- dataset[[field]]
    space <- .fmri_try_space(obj)
    if (!is.null(space)) return(space)
  }

  NULL
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
