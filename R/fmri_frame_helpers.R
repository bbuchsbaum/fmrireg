###############################################################################
## fmri_frame_helpers.R
##
## fmrireg consumes `fmridataset::fmri_frame` objects (observations by
## features, with the acquisition timing carried as observation metadata and
## the experimental events carried as a keyed event table). This file holds
##
##   * the deprecated frame constructors (matrix_frame() and friends), which
##     forward to fmridataset::as_fmri_frame() and read_nifti_bold(), and
##   * the small internal accessor layer the fitting code uses so that the
##     frame representation is touched in one place.
##
## Every accessor calls the public fmridataset API; nothing here reaches into
## fmridataset internals.
###############################################################################

#' @importFrom fmridataset fmri_frame collect_assay
NULL

#' @importFrom fmridataset as_fmri_frame
#' @export
fmridataset::as_fmri_frame

#' @importFrom fmridataset read_nifti_bold
#' @export
fmridataset::read_nifti_bold

# ---------------------------------------------------------------------------
# Constructors (deprecated)
#
# fmridataset now builds run-structured frames itself: as_fmri_frame() for a
# matrix, NeuroVec volumes, or a LatentNeuroVec, and read_nifti_bold() for
# NIfTI files. The four helpers below forward to those and are kept for one
# release. Saved jobs that name them still realize silently; see
# .dataset_constructors.
# ---------------------------------------------------------------------------

#' Deprecated frame constructors
#'
#' `matrix_frame()`, `neurovec_frame()`, `nifti_frame()`, and `latent_frame()`
#' are deprecated. fmridataset now builds the same run-structured frames:
#'
#' | Deprecated | Use instead |
#' |---|---|
#' | `matrix_frame(datamat, ...)` | [fmridataset::as_fmri_frame()]`(datamat, ...)` |
#' | `neurovec_frame(scans, mask, ...)` | [fmridataset::as_fmri_frame()]`(scans, mask = mask, ...)` |
#' | `nifti_frame(scans, mask, ...)` | [fmridataset::read_nifti_bold()]`(scans, mask, ...)` |
#' | `latent_frame(x, ...)` | [fmridataset::as_fmri_frame()]`(x, ...)` |
#'
#' The argument `event_table` is called `events` in the replacements. Each
#' helper warns once per call and forwards to its replacement, so it returns
#' what the replacement returns: observation IDs such as `run-1_vol-0001`, a
#' `volume_index` column, and, for volumes, a mask that must match the data's
#' grid (dimensions and affine).
#'
#' @param datamat A numeric matrix (time by feature), or a vector for one
#'   feature.
#' @param scans For `neurovec_frame()`, a `NeuroVec` or a list of them, one per
#'   run; for `nifti_frame()`, NIfTI file paths, one per run.
#' @param mask A `neuroim2::NeuroVol` mask, a [fmridataset::volume_space()],
#'   or (for `nifti_frame()`) a NIfTI mask path.
#' @param x A `LatentNeuroVec`, or a list with elements `basis` (time by
#'   component scores), `loadings` (masked voxel by component), and `mask`.
#' @param TR Repetition time in seconds: one value, or one per run.
#' @param run_length Number of volumes in each run.
#' @param event_table Optional events: a data frame or an
#'   [fmridataset::event_table()].
#' @param censor Optional censoring indicator: logical or 0/1 per volume, or
#'   1-based volume indices.
#' @param feature_ids Optional feature IDs, one per column of `datamat`.
#' @param base_path Directory that relative `scans` and `mask` paths are
#'   resolved against.
#' @param assay Name of the assay.
#' @return An `fmri_frame`.
#' @name deprecated-frame-constructors
#' @keywords internal
NULL

.deprecate_frame_constructor <- function(old, new) {
  .Deprecated(new, package = "fmrireg", old = old)
}

#' @rdname deprecated-frame-constructors
#' @export
matrix_frame <- function(datamat, TR, run_length, event_table = NULL,
                         censor = NULL, feature_ids = NULL, assay = "signal") {
  .deprecate_frame_constructor("matrix_frame", "fmridataset::as_fmri_frame")
  .matrix_frame_impl(datamat, TR, run_length, event_table, censor, feature_ids, assay)
}

#' @rdname deprecated-frame-constructors
#' @export
neurovec_frame <- function(scans, mask, TR, run_length = NULL, event_table = NULL,
                           censor = NULL, assay = "bold") {
  .deprecate_frame_constructor("neurovec_frame", "fmridataset::as_fmri_frame")
  .neurovec_frame_impl(scans, mask, TR, run_length, event_table, censor, assay)
}

#' @rdname deprecated-frame-constructors
#' @export
nifti_frame <- function(scans, mask, TR, run_length = NULL, event_table = NULL,
                        censor = NULL, base_path = ".", assay = "bold") {
  .deprecate_frame_constructor("nifti_frame", "fmridataset::read_nifti_bold")
  .nifti_frame_impl(scans, mask, TR, run_length, event_table, censor, base_path, assay)
}

#' @rdname deprecated-frame-constructors
#' @export
latent_frame <- function(x, TR, run_length = NULL, event_table = NULL,
                         censor = NULL, assay = "scores") {
  .deprecate_frame_constructor("latent_frame", "fmridataset::as_fmri_frame")
  .latent_frame_impl(x, TR, run_length, event_table, censor, assay)
}

# The forwarding itself, shared by the deprecated helpers and by saved jobs
# that name them (which realize without warnings).

#' @keywords internal
#' @noRd
.matrix_frame_impl <- function(datamat, TR, run_length, event_table = NULL,
                               censor = NULL, feature_ids = NULL, assay = "signal") {
  if (is.vector(datamat) && !is.list(datamat)) datamat <- as.matrix(datamat)
  fmridataset::as_fmri_frame(
    datamat,
    TR = TR, run_length = run_length, events = event_table, censor = censor,
    feature_ids = feature_ids, assay = assay
  )
}

#' @keywords internal
#' @noRd
.neurovec_frame_impl <- function(scans, mask, TR, run_length = NULL, event_table = NULL,
                                 censor = NULL, assay = "bold") {
  fmridataset::as_fmri_frame(
    scans,
    mask = mask, TR = TR, run_length = run_length, events = event_table,
    censor = censor, assay = assay
  )
}

#' @keywords internal
#' @noRd
.nifti_frame_impl <- function(scans, mask, TR, run_length = NULL, event_table = NULL,
                              censor = NULL, base_path = ".", assay = "bold") {
  assert_that(is.character(scans) && length(scans) > 0,
              msg = "'scans' must be a character vector of NIfTI paths")
  if (is.character(mask)) mask <- .frame_resolve_paths(mask, base_path)
  fmridataset::read_nifti_bold(
    .frame_resolve_paths(scans, base_path), mask,
    TR = TR, run_length = run_length, events = event_table, censor = censor,
    assay = assay
  )
}

#' @keywords internal
#' @noRd
.latent_frame_impl <- function(x, TR, run_length = NULL, event_table = NULL,
                               censor = NULL, assay = "scores") {
  if (methods::is(x, "LatentNeuroVec")) {
    return(fmridataset::as_fmri_frame(
      x,
      TR = TR, run_length = run_length, events = event_table, censor = censor,
      assay = assay
    ))
  }
  # The list form (basis, loadings, mask) has no fmridataset counterpart: build
  # the basis space here and hand the scores to as_fmri_frame().
  if (!(is.list(x) && all(c("basis", "loadings", "mask") %in% names(x)))) {
    stop("'x' must be a fmristore::LatentNeuroVec or a list with 'basis', ",
         "'loadings', and 'mask'", call. = FALSE)
  }
  scores <- as.matrix(x$basis)
  loadings <- x$loadings
  if (!inherits(loadings, "Matrix")) loadings <- as.matrix(loadings)
  assert_that(ncol(scores) == ncol(loadings),
              msg = "'basis' and 'loadings' must have the same number of components")
  assert_that(inherits(x$mask, "NeuroVol"), msg = "the latent mask must be a neuroim2 NeuroVol")
  parent <- .volume_space_from_mask(x$mask)
  assert_that(fmridataset::n_features(parent) == nrow(loadings),
              msg = "'loadings' must have one row per voxel in the mask")
  component_ids <- .frame_component_ids(scores, loadings)
  dimnames(loadings) <- NULL
  space <- fmridataset::basis_space_from_decoder(
    parent, component_ids = component_ids, decoder = loadings, encoder = "none"
  )
  fmridataset::as_fmri_frame(
    unname(scores),
    TR = TR, run_length = run_length %||% nrow(scores), events = event_table,
    censor = censor, space = space, assay = assay
  )
}

# ---------------------------------------------------------------------------
# Construction helpers
# ---------------------------------------------------------------------------

# Component IDs for a latent decomposition. A decomposition usually names its
# own components (prcomp writes PC1..PCr); those names are preferred, and
# generated only when they are absent, duplicated, or empty. Scores are checked
# before loadings so that the assay's own labelling wins when the two disagree.
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

#' Read a block of the active assay
#'
#' The single data-access point of the fitting engines: every read of voxel
#' data goes through here, one (rows x cols) block at a time.
#' @param x An `fmri_frame`.
#' @param rows,cols Integer row / feature positions, or `NULL` for all.
#' @keywords internal
#' @noRd
.dset_read <- function(x, rows = NULL, cols = NULL) {
  view <- if (is.null(rows) && is.null(cols)) {
    x
  } else if (is.null(rows)) {
    x[, cols]
  } else if (is.null(cols)) {
    x[rows, ]
  } else {
    x[rows, cols]
  }
  fmridataset::collect_assay(view)
}

#' A chunk of data: the shape fmrireg's chunkwise and runwise engines consume
#'
#' A chunk is either eager (it carries `mat`) or a lazy descriptor that holds
#' a source and the row / feature positions and reads them only when its data
#' is requested. The source is an environment holding the frame, shared by all
#' chunks of one split, so serialising the chunk list (e.g. to future workers)
#' carries one reference to the frame rather than one copy per chunk.
#' `chunk$data` reads a lazy chunk on every access, so engines read it once
#' with `.chunk_data()` and keep the matrix for as long as they need it.
#' @keywords internal
#' @noRd
.data_chunk <- function(mat, voxel_ind, row_ind, chunk_num, source = NULL,
                        rows = NULL, cols = NULL) {
  structure(
    list(data = mat, voxel_ind = voxel_ind, row_ind = row_ind,
         chunk_num = chunk_num, source = source, rows = rows, cols = cols),
    class = c("data_chunk", "list")
  )
}

#' @keywords internal
#' @noRd
.chunk_source <- function(x) {
  src <- new.env(parent = emptyenv())
  src$frame <- x
  src
}

#' Data of a (possibly lazy) chunk
#' @keywords internal
#' @noRd
.chunk_data <- function(chunk) {
  mat <- .subset2(chunk, "data")
  if (!is.null(mat)) return(mat)
  src <- .subset2(chunk, "source")
  if (is.null(src)) return(NULL)
  .dset_read(src$frame, .subset2(chunk, "rows"), .subset2(chunk, "cols"))
}

#' @keywords internal
#' @noRd
#' @rawNamespace S3method("$",data_chunk)
`$.data_chunk` <- function(x, name) {
  if (identical(name, "data")) return(.chunk_data(x))
  .subset2(x, name)
}

#' Global row indices of each run, from the temporal schema (no data read)
#' @keywords internal
#' @noRd
.dset_run_rows <- function(x) {
  schema <- fmridataset::temporal_schema(x)
  if (!isTRUE(schema$contiguous)) {
    stop("runwise processing requires a frame whose runs occupy contiguous ",
         "observation ranges", call. = FALSE)
  }
  unname(split(seq_len(nrow(x)), schema$block_ids))
}

#' One data chunk per run, in acquisition order
#'
#' Each chunk carries the global row indices the run occupies, the feature
#' positions and the run number. The run's observation-by-feature matrix is
#' read lazily (`.chunk_data()`), one run at a time.
#' @keywords internal
#' @noRd
.dset_run_chunks <- function(x) {
  rows <- .dset_run_rows(x)
  nvox <- ncol(x)
  src <- .chunk_source(x)
  lapply(seq_along(rows), function(i) {
    ind <- rows[[i]]
    .data_chunk(NULL, voxel_ind = seq_len(nvox), row_ind = ind, chunk_num = i,
                source = src, rows = ind, cols = NULL)
  })
}

#' Feature-wise data chunks covering all observations
#'
#' Lazy descriptors: each chunk's data is read only when an engine asks for it,
#' so at most one chunk is resident at a time.
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
    return(list(.data_chunk(NULL, voxel_ind = seq_len(nvox),
                            row_ind = seq_len(nobs), chunk_num = 1L,
                            source = .chunk_source(x), rows = NULL, cols = NULL)))
  }
  sidx <- split(seq_len(nvox), sort(rep(seq_len(nchunks), length.out = nvox)))
  src <- .chunk_source(x)
  lapply(seq_along(sidx), function(i) {
    cols <- sidx[[i]]
    .data_chunk(NULL, voxel_ind = cols, row_ind = seq_len(nobs), chunk_num = i,
                source = src, rows = NULL, cols = cols)
  })
}
