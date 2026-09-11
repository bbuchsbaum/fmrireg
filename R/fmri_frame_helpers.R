###############################################################################
## fmri_frame_helpers.R
##
## fmrireg consumes `fmridataset::fmri_frame` objects (observations by
## features, with the acquisition timing carried as observation metadata and
## the experimental events carried as a keyed event table). This file holds
##
##   * exported convenience constructors that build a frame from the inputs
##     fmrireg users typically hold (a matrix, NeuroVec objects, NIfTI files,
##     or a fmristore::LatentNeuroVec), and
##   * the small internal accessor layer the fitting code uses so that the
##     frame representation is touched in one place.
##
## Every accessor calls the public fmridataset API; nothing here reaches into
## fmridataset internals.
###############################################################################

#' @importFrom fmridataset fmri_frame collect_assay temporal_schema
#'   as_sampling_frame
NULL

# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------

#' Build an fmri_frame from an in-memory matrix
#'
#' fmrireg fits models to [fmridataset::fmri_frame()] objects. This helper
#' builds one from a time-by-feature matrix plus the run structure and event
#' table that the design machinery needs. Observation IDs are derived
#' deterministically as `run-<r>-vol-<index>`; feature IDs default to the
#' matrix column names when they are unique and to `feature-<j>` otherwise.
#'
#' @param datamat A numeric matrix with one row per acquired volume and one
#'   column per feature (voxel, vertex, parcel, or component).
#' @param TR Repetition time in seconds; one value, or one per run.
#' @param run_length Integer vector giving the number of volumes in each run.
#'   Must sum to `nrow(datamat)`.
#' @param event_table Optional data frame of events (onsets, conditions, and
#'   the block variable). It is stored on the frame as a keyed
#'   [fmridataset::event_table()]; an `event_id` column is added when absent.
#' @param censor Optional censoring indicator with one entry per volume: a
#'   logical vector, a 0/1 vector, or a vector of 1-based volume indices.
#'   Stored as the logical `censor` observation column that
#'   [fmridataset::temporal_schema()] recognises. [fmri_lm()] uses it without
#'   further configuration; it feeds AR estimation and whitening, and does not
#'   drop the censored volumes from the regression.
#' @param feature_ids Optional stable feature IDs (one per column).
#' @param assay Name of the assay holding `datamat`.
#' @return An `fmri_frame` whose feature space is an
#'   [fmridataset::index_space()].
#' @seealso [neurovec_frame()], [nifti_frame()], [latent_frame()],
#'   [fmridataset::fmri_frame()]
#' @export
#' @examples
#' Y <- matrix(rnorm(80 * 3), 80, 3)
#' events <- data.frame(onset = c(5, 25, 45, 65),
#'                      condition = factor(c("A", "B", "A", "B")),
#'                      run = c(1, 1, 2, 2))
#' frame <- matrix_frame(Y, TR = 2, run_length = c(40, 40), event_table = events)
#' fmridataset::temporal_schema(frame)$run_lengths
matrix_frame <- function(datamat, TR, run_length, event_table = NULL,
                         censor = NULL, feature_ids = NULL, assay = "signal") {
  if (is.vector(datamat) && !is.list(datamat)) {
    datamat <- as.matrix(datamat)
  }
  assert_that(is.matrix(datamat), msg = "'datamat' must be a matrix")
  run_length <- .frame_check_run_length(run_length, nrow(datamat))
  if (is.null(feature_ids)) {
    cn <- colnames(datamat)
    feature_ids <- if (!is.null(cn) && !anyDuplicated(cn) && all(nzchar(cn))) {
      cn
    } else {
      sprintf("feature-%d", seq_len(ncol(datamat)))
    }
  }
  assert_that(length(feature_ids) == ncol(datamat),
              msg = "'feature_ids' must have one entry per column of 'datamat'")
  # Column names, when usable, became the feature IDs just above; row names
  # never describe the generated volume IDs. Either way the axis IDs are
  # authoritative, and leaving dimnames on the assay would assert a second
  # naming of the same axes that fmri_frame() rejects as ambiguous.
  dimnames(datamat) <- NULL
  assays <- list(datamat)
  names(assays) <- assay
  fmridataset::fmri_frame(
    assays = assays,
    observations = .frame_observations(run_length, TR, censor),
    space = fmridataset::index_space(ncol(datamat), ids = as.character(feature_ids)),
    tables = .frame_tables(event_table)
  )
}

#' Build an fmri_frame from in-memory NeuroVec objects
#'
#' Extracts the masked time series of one or more `neuroim2::NeuroVec` runs
#' and builds an [fmridataset::fmri_frame()] whose feature space is the
#' [fmridataset::volume_space()] described by `mask`. Voxel time series are
#' read with `neuroim2::series()`, so a `SparseNeuroVec` does not need to be
#' densified.
#'
#' @param scans A `NeuroVec`, or a list of them, one per run.
#' @param mask A `neuroim2::LogicalNeuroVol` (or a `NeuroVol` whose non-zero
#'   voxels form the mask) on the same 3D grid as `scans`.
#' @param TR Repetition time in seconds; one value, or one per run.
#' @param run_length Volumes per run; defaults to the fourth dimension of each
#'   scan.
#' @param event_table,censor As in [matrix_frame()].
#' @param assay Name of the assay.
#' @return An `fmri_frame` over a `volume_space`.
#' @seealso [matrix_frame()], [nifti_frame()]
#' @export
neurovec_frame <- function(scans, mask, TR, run_length = NULL, event_table = NULL,
                           censor = NULL, assay = "bold") {
  if (inherits(scans, "NeuroVec")) scans <- list(scans)
  assert_that(is.list(scans) && length(scans) > 0 &&
                all(vapply(scans, function(x) inherits(x, "NeuroVec"), logical(1))),
              msg = "'scans' must be a NeuroVec or a list of NeuroVec objects")
  assert_that(inherits(mask, "NeuroVol"), msg = "'mask' must be a neuroim2 NeuroVol")
  mask_dim <- dim(mask)
  for (s in scans) {
    assert_that(all(dim(s)[1:3] == mask_dim),
                msg = "every scan must share the spatial dimensions of 'mask'")
  }
  scan_lengths <- vapply(scans, function(x) as.integer(dim(x)[4]), integer(1))
  if (is.null(run_length)) run_length <- scan_lengths
  run_length <- .frame_check_run_length(run_length, sum(scan_lengths))

  voxel_ind <- which(as.vector(mask) != 0)
  datamat <- do.call(rbind, lapply(scans, function(s) neuroim2::series(s, voxel_ind)))
  # series() labels its columns by voxel index; the frame's feature IDs come
  # from the volume space, so the labels are a second naming of the same axis.
  dimnames(datamat) <- NULL
  assays <- list(datamat)
  names(assays) <- assay
  fmridataset::fmri_frame(
    assays = assays,
    observations = .frame_observations(run_length, TR, censor),
    space = .volume_space_from_mask(mask),
    tables = .frame_tables(event_table)
  )
}

#' Build a lazy fmri_frame from NIfTI files
#'
#' Wraps [fmridataset::nifti_array_source()], so headers and the mask are
#' read at construction but no BOLD volumes are read until the fit needs them.
#' For fMRIPrep derivatives see [fmridataset::read_bids_bold()].
#'
#' @param scans Character vector of NIfTI paths, one per run, all on the same
#'   grid.
#' @param mask A NIfTI mask path, a `neuroim2::LogicalNeuroVol`, or a
#'   [fmridataset::volume_space()].
#' @param TR Repetition time in seconds; one value, or one per run.
#' @param run_length Volumes per run; defaults to the fourth header dimension
#'   of each file.
#' @param event_table,censor As in [matrix_frame()].
#' @param base_path Directory that relative `scans` and `mask` paths are
#'   resolved against.
#' @param assay Name of the assay.
#' @return An `fmri_frame` over the NIfTI source's `volume_space`.
#' @seealso [matrix_frame()], [neurovec_frame()],
#'   [fmridataset::nifti_array_source()]
#' @export
nifti_frame <- function(scans, mask, TR, run_length = NULL, event_table = NULL,
                        censor = NULL, base_path = ".", assay = "bold") {
  assert_that(is.character(scans) && length(scans) > 0,
              msg = "'scans' must be a character vector of NIfTI paths")
  paths <- .frame_resolve_paths(scans, base_path)
  if (is.character(mask)) {
    mask <- .frame_resolve_paths(mask, base_path)
  } else if (inherits(mask, "NeuroVol")) {
    mask <- .volume_space_from_mask(mask)
  }
  source <- fmridataset::nifti_array_source(paths, mask)
  space <- fmridataset::nifti_source_space(source)
  n_obs <- fmridataset::source_shape(source)[[1L]]
  if (is.null(run_length)) {
    run_length <- vapply(paths, function(p) {
      hdr <- neuroim2::read_header(p)
      dims <- methods::slot(hdr, "dims")
      if (length(dims) >= 4L) as.integer(dims[4L]) else 1L
    }, integer(1))
  }
  run_length <- .frame_check_run_length(run_length, n_obs)
  assays <- list(source)
  names(assays) <- assay
  fmridataset::fmri_frame(
    assays = assays,
    observations = .frame_observations(run_length, TR, censor),
    space = space,
    tables = .frame_tables(event_table)
  )
}

#' Build an fmri_frame from a latent (basis) decomposition
#'
#' A `fmristore::LatentNeuroVec` stores temporal scores (time by component)
#' and spatial loadings (voxel by component). This helper places the scores in
#' the assay and the loadings in a [fmridataset::basis_space()] whose parent
#' is the mask's `volume_space`, so fmrireg's latent engines can fit in
#' component space and reconstruct voxel maps through
#' [fmridataset::basis_synthesis()].
#'
#' @param x A `fmristore::LatentNeuroVec`, or a list with elements `basis`
#'   (time by component scores), `loadings` (masked-voxel by component), and
#'   `mask` (a `neuroim2::LogicalNeuroVol`).
#' @param TR Repetition time in seconds; one value, or one per run.
#' @param run_length Volumes per run; defaults to a single run spanning all
#'   score rows.
#' @param event_table,censor As in [matrix_frame()].
#' @param assay Name of the assay.
#' @return An `fmri_frame` over a synthesis-only `basis_space`.
#' @seealso [matrix_frame()], [fmridataset::basis_space_from_decoder()]
#' @export
latent_frame <- function(x, TR, run_length = NULL, event_table = NULL,
                         censor = NULL, assay = "scores") {
  if (methods::is(x, "LatentNeuroVec")) {
    parts <- list(
      basis = methods::slot(x, "basis"),
      loadings = methods::slot(x, "loadings"),
      mask = methods::slot(x, "mask")
    )
  } else if (is.list(x) && all(c("basis", "loadings", "mask") %in% names(x))) {
    parts <- x
  } else {
    stop("'x' must be a fmristore::LatentNeuroVec or a list with 'basis', ",
         "'loadings', and 'mask'", call. = FALSE)
  }
  scores <- as.matrix(parts$basis)
  loadings <- parts$loadings
  if (!inherits(loadings, "Matrix")) loadings <- as.matrix(loadings)
  assert_that(ncol(scores) == ncol(loadings),
              msg = "'basis' and 'loadings' must have the same number of components")
  mask <- parts$mask
  assert_that(inherits(mask, "NeuroVol"), msg = "the latent mask must be a neuroim2 NeuroVol")
  if (is.null(run_length)) run_length <- nrow(scores)
  run_length <- .frame_check_run_length(run_length, nrow(scores))

  parent <- .volume_space_from_mask(mask)
  assert_that(fmridataset::n_features(parent) == nrow(loadings),
              msg = "'loadings' must have one row per voxel in the mask")
  component_ids <- .frame_component_ids(scores, loadings)
  space <- fmridataset::basis_space_from_decoder(
    parent, component_ids = component_ids, decoder = loadings, encoder = "none"
  )
  # The axis IDs are authoritative. A decomposition routinely labels its own
  # columns (prcomp writes PC1..PCr) and those labels either became the
  # component IDs just above or were not usable; either way carrying them on
  # the assay as well would assert a second, unchecked naming of the same axis,
  # which fmri_frame() rejects as ambiguous.
  dimnames(scores) <- NULL
  assays <- list(scores)
  names(assays) <- assay
  fmridataset::fmri_frame(
    assays = assays,
    observations = .frame_observations(run_length, TR, censor),
    space = space,
    tables = .frame_tables(event_table)
  )
}

# ---------------------------------------------------------------------------
# Construction helpers
# ---------------------------------------------------------------------------

# Component IDs for a latent decomposition. A decomposition usually names its
# own components (prcomp writes PC1..PCr); those names are preferred, as
# matrix_frame() prefers column names, and generated only when they are absent,
# duplicated, or empty. Scores are checked before loadings so that the assay's
# own labelling wins when the two disagree.
#' @keywords internal
#' @noRd
.frame_component_ids <- function(scores, loadings) {
  usable <- function(nm, n) {
    !is.null(nm) && length(nm) == n && !anyNA(nm) && all(nzchar(nm)) && !anyDuplicated(nm)
  }
  n <- ncol(loadings)
  for (nm in list(colnames(scores), colnames(loadings))) {
    if (usable(nm, n)) {
      return(as.character(nm))
    }
  }
  sprintf("component-%d", seq_len(n))
}

#' @keywords internal
#' @noRd
.frame_check_run_length <- function(run_length, n_obs) {
  assert_that(is.numeric(run_length) && length(run_length) > 0 && all(run_length > 0) &&
                all(run_length == as.integer(run_length)),
              msg = "'run_length' must be a vector of positive integers")
  run_length <- as.integer(run_length)
  if (sum(run_length) != n_obs) {
    stop(sprintf("sum(run_length) = %d does not match the number of volumes (%d)",
                 sum(run_length), n_obs), call. = FALSE)
  }
  run_length
}

#' @keywords internal
#' @noRd
.frame_observations <- function(run_length, TR, censor = NULL) {
  n_runs <- length(run_length)
  n_obs <- sum(run_length)
  assert_that(is.numeric(TR) && all(is.finite(TR)) && all(TR > 0),
              msg = "'TR' must be positive and finite")
  if (length(TR) == 1L) {
    TR <- rep(TR, n_runs)
  }
  assert_that(length(TR) == n_runs, msg = "'TR' must have length one or one entry per run")
  run <- rep(seq_len(n_runs), run_length)
  index <- unlist(lapply(run_length, seq_len), use.names = FALSE)
  obs <- data.frame(
    .obs_id = sprintf("run-%d-vol-%05d", run, index),
    run_id = sprintf("run-%d", run),
    TR = rep(as.numeric(TR), run_length),
    stringsAsFactors = FALSE
  )
  if (!is.null(censor)) {
    obs$censor <- .frame_censor_logical(censor, n_obs)
  }
  obs
}

#' @keywords internal
#' @noRd
.frame_censor_logical <- function(censor, n_obs) {
  if (is.logical(censor)) {
    assert_that(length(censor) == n_obs && !anyNA(censor),
                msg = "a logical 'censor' must have one non-missing entry per volume")
    return(censor)
  }
  assert_that(is.numeric(censor) && !anyNA(censor), msg = "'censor' must be logical or numeric")
  if (length(censor) == n_obs && all(censor %in% c(0, 1))) {
    return(censor == 1)
  }
  assert_that(all(censor >= 1 & censor <= n_obs & censor == as.integer(censor)),
              msg = "a numeric 'censor' must be a 0/1 vector or 1-based volume indices")
  out <- rep(FALSE, n_obs)
  out[as.integer(censor)] <- TRUE
  out
}

#' @keywords internal
#' @noRd
.frame_tables <- function(event_table) {
  if (is.null(event_table)) return(list())
  assert_that(is.data.frame(event_table), msg = "'event_table' must be a data frame")
  if (nrow(event_table) == 0L) return(list())
  if (!"event_id" %in% names(event_table)) {
    event_table <- cbind(
      data.frame(event_id = sprintf("event-%d", seq_len(nrow(event_table))),
                 stringsAsFactors = FALSE),
      as.data.frame(event_table)
    )
  }
  list(events = fmridataset::event_table(event_table, key = "event_id"))
}

#' @keywords internal
#' @noRd
.frame_resolve_paths <- function(paths, base_path = ".") {
  ifelse(grepl("^(/|~|[A-Za-z]:)", paths), paths, file.path(base_path, paths))
}

#' Volume space for a neuroim2 mask volume
#' @keywords internal
#' @noRd
.volume_space_from_mask <- function(mask) {
  sp <- neuroim2::space(mask)
  fmridataset::volume_space(
    dim = as.integer(dim(mask)[1:3]),
    affine = unname(neuroim2::trans(sp)),
    support = as.logical(as.vector(mask) != 0)
  )
}

# ---------------------------------------------------------------------------
# Internal accessors used by the fitting code
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.is_fmri_frame <- function(x) inherits(x, "fmri_frame")

#' @keywords internal
#' @noRd
.dset_root <- function(x) {
  while (inherits(x, "fmri_view")) x <- x$base
  x
}

#' Sampling frame implied by the frame's observation metadata
#' @keywords internal
#' @noRd
.dset_sampling_frame <- function(x) fmridataset::as_sampling_frame(x)

#' Run lengths (in acquisition order) of a contiguous frame
#' @keywords internal
#' @noRd
.dset_run_lengths <- function(x) {
  unname(as.integer(fmridataset::temporal_schema(x)$run_lengths))
}

#' Integer block id (1-based, order of first appearance) per observation
#' @keywords internal
#' @noRd
.dset_block_ids <- function(x) fmridataset::temporal_schema(x)$block_ids

#' Event table stored on the frame, as a data frame (empty when absent)
#' @keywords internal
#' @noRd
.dset_event_table <- function(x) {
  tables <- .dset_root(x)$tables
  events <- if (is.list(tables)) tables[["events"]] else NULL
  if (is.null(events)) return(data.frame())
  fmridataset::event_data(events)
}

#' Logical censor vector from the frame's temporal schema, or NULL
#' @keywords internal
#' @noRd
.dset_censor <- function(x) fmridataset::temporal_schema(x)$censor

#' Dense observation-by-feature matrix of the active assay
#' @keywords internal
#' @noRd
.dset_data_matrix <- function(x, ...) fmridataset::collect_assay(x, ...)

#' Feature space (or the parent volume space of a basis space)
#' @keywords internal
#' @noRd
.dset_space <- function(x) fmridataset::space(x)

#' @keywords internal
#' @noRd
.dset_is_latent <- function(x) {
  .is_fmri_frame(x) && inherits(fmridataset::space(x), "basis_space")
}

#' Voxel-by-component loadings of a latent frame
#' @keywords internal
#' @noRd
.dset_loadings <- function(x) {
  sp <- fmridataset::space(x)
  if (!inherits(sp, "basis_space")) {
    stop("the dataset's feature space is not a basis_space; no latent loadings are available",
         call. = FALSE)
  }
  loadings <- fmridataset::basis_synthesis(sp)
  if (is.null(loadings)) {
    stop("the dataset's basis_space carries no synthesis operator (loadings)", call. = FALSE)
  }
  loadings
}

#' The volume_space underlying a frame, or NULL for non-volumetric features
#' @keywords internal
#' @noRd
.dset_volume_space <- function(x) {
  sp <- fmridataset::space(x)
  if (inherits(sp, "basis_space")) sp <- fmridataset::parent_space(sp)
  if (inherits(sp, "volume_space")) sp else NULL
}

#' A chunk of data: the shape fmrireg's chunkwise and runwise engines consume
#' @keywords internal
#' @noRd
.data_chunk <- function(mat, voxel_ind, row_ind, chunk_num) {
  structure(
    list(data = mat, voxel_ind = voxel_ind, row_ind = row_ind, chunk_num = chunk_num),
    class = c("data_chunk", "list")
  )
}

#' One data chunk per run, in acquisition order
#'
#' Each chunk carries the run's observation-by-feature matrix, the global row
#' indices it occupies, the feature positions, and the run number. Reads are
#' performed one run at a time through a lazy frame view.
#' @keywords internal
#' @noRd
.dset_run_chunks <- function(x) {
  schema <- fmridataset::temporal_schema(x)
  if (!isTRUE(schema$contiguous)) {
    stop("runwise processing requires a frame whose runs occupy contiguous ",
         "observation ranges", call. = FALSE)
  }
  rows <- split(seq_len(nrow(x)), schema$block_ids)
  nvox <- ncol(x)
  lapply(seq_along(rows), function(i) {
    ind <- rows[[i]]
    mat <- fmridataset::collect_assay(x[ind, ])
    .data_chunk(mat, voxel_ind = seq_len(nvox), row_ind = ind, chunk_num = i)
  })
}

#' Feature-wise data chunks covering all observations
#' @keywords internal
#' @noRd
.dset_chunks <- function(x, nchunks = 1) {
  nvox <- ncol(x)
  nobs <- nrow(x)
  assert_that(is.numeric(nchunks) && length(nchunks) == 1L && nchunks >= 1,
              msg = "'nchunks' must be a positive number")
  if (nchunks > nvox) {
    warning("requested number of chunks is greater than number of voxels in mask",
            call. = FALSE)
    nchunks <- nvox
  }
  if (nchunks == 1) {
    return(list(.data_chunk(fmridataset::collect_assay(x), voxel_ind = seq_len(nvox),
                            row_ind = seq_len(nobs), chunk_num = 1L)))
  }
  sidx <- split(seq_len(nvox), sort(rep(seq_len(nchunks), length.out = nvox)))
  lapply(seq_along(sidx), function(i) {
    cols <- sidx[[i]]
    mat <- fmridataset::collect_assay(x[, cols])
    .data_chunk(mat, voxel_ind = cols, row_ind = seq_len(nobs), chunk_num = i)
  })
}
