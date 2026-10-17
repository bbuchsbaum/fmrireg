# BIDS Export Functions using LabeledVolumeSet
# Functions for saving fMRI analysis results in BIDS-compliant format
#' @importFrom jsonlite write_json
NULL

#' Write Results from fMRI Linear Model
#'
#' Exports statistical maps from an fmri_lm object with BIDS-compliant naming
#' and JSON metadata sidecars. Supported image outputs are HDF5
#' LabeledVolumeSet files and NIfTI files; an fmrigds representation can also
#' be requested.
#'
#' @param x An fmri_lm object containing fitted model results
#' @param path Output directory path. If NULL, uses current working directory
#' @param subject Subject identifier (e.g., "01", "1001"). Required.
#' @param task Task identifier (e.g., "nback", "rest"). Required for BIDS compliance.
#' @param space Spatial reference (e.g., "MNI152NLin2009cAsym"). Optional but recommended.
#' @param desc BIDS \code{desc-} label of the analysis (default: "GLM"),
#'   sanitised to letters and digits. Beta files are labelled
#'   \code{desc-<desc>}; statistic maps \code{desc-<desc><stat>} (e.g.
#'   \code{desc-run1tstat}), except that the default \code{"GLM"} is omitted
#'   (\code{desc-tstat}), the naming [collect_results()] reads by default.
#' @param format Output format(s). One or more of `"h5"` (BIDS-style HDF5
#'   outputs), `"nifti"` (BIDS-style NIfTI outputs, written as `.nii.gz`), and
#'   `"gds"` (fmrigds-compatible assays plus an `.rds` plan). Pass a character
#'   vector to write multiple formats at once (e.g. `c("h5", "nifti")`).
#'   Defaults to `"h5"` only; `"h5"` requires the suggested package fmristore,
#'   `"nifti"` does not. NIfTI maps are stored as 32-bit floats
#'   (datatype FLOAT32), so their values are the fitted values rounded to
#'   single precision.
#' @param strategy Storage strategy: "by_stat" (group contrasts by statistic) or "by_contrast" (separate files)
#' @param save_betas Logical. Save raw regressor betas (default: TRUE)
#' @param contrasts Character vector of contrast names to save. NULL saves all contrasts
#' @param contrast_match How to match `contrasts`: `"auto"` first matches exact
#'   contrast names and treats unmatched selectors as regular expressions;
#'   `"exact"` uses literal names only; `"regex"` treats all selectors as regular
#'   expressions.
#' @param contrast_stats Character vector of contrast statistics to save (default: c("beta", "tstat", "pval", "se"))
#' @param overwrite Logical. Overwrite existing files (default: FALSE)
#' @param validate_inputs Logical. Validate fmrilm object structure (default: TRUE)
#' @param ... Additional arguments passed to internal functions
#'
#' @return Invisible list of file paths created
#' @export
#'
#' @examples
#' \dontrun{
#' # Save all results using default settings
#' write_results(fitted_model, subject = "01", task = "nback")
#'
#' # Save only specific contrasts and statistics
#' write_results(fitted_model,
#'               subject = "01", task = "nback", space = "MNI152NLin2009cAsym",
#'               contrasts = c("FacesVsPlaces", "GoVsNoGo"),
#'               contrast_stats = c("beta", "tstat"))
#' }
#' @method write_results fmri_lm
write_results.fmri_lm <- function(x,
                                  path = NULL,
                                  subject = NULL,
                                  task = NULL,
                                  space = NULL,
                                  desc = "GLM",
                                  format = c("h5", "nifti", "gds"),
                                  strategy = c("by_stat", "by_contrast"),
                                  save_betas = TRUE,
                                  contrasts = NULL,
                                  contrast_match = c("auto", "exact", "regex"),
                                  contrast_stats = c("beta", "tstat", "pval", "se"),
                                  overwrite = FALSE,
                                  validate_inputs = TRUE,
                                  ...) {

  strategy <- match.arg(strategy)
  contrast_match <- match.arg(contrast_match)
  # The full choice vector is the formal default purely so it appears in the
  # usage/help; an unspecified `format` still defaults to "h5" alone. Because
  # `several.ok = TRUE` disables match.arg's "arg identical to choices ->
  # first element" shortcut, detect the default explicitly with missing().
  if (missing(format)) {
    format <- "h5"
  } else {
    format <- match.arg(format, c("h5", "nifti", "gds"), several.ok = TRUE)
  }
  if (length(format) == 0) format <- "h5"
  produce_gds <- "gds" %in% format
  # fmristore (Suggests) is only needed for HDF5 output; check before staging.
  if ("h5" %in% format) .require_fmristore()

  # Input validation
  if (validate_inputs) {
    .validate_fmrilm_object(x)
  }

  # Set default path
  if (is.null(path)) {
    path <- getwd()
  }

  # Create BIDS entities list
  entities <- .create_bids_entities(subject, task, space)

  # Validate required entities
  .validate_required_entities(entities)

  # One plan drives both the overwrite check and the writes.
  plan <- .bids_output_plan(
    fmrilm_obj = x,
    entities = entities,
    desc = desc,
    format = format,
    strategy = strategy,
    save_betas = save_betas,
    contrasts = contrasts,
    contrast_match = contrast_match,
    contrast_stats = contrast_stats
  )
  predicted_files <- .predict_output_files(path, plan, entities, desc, format)

  if (!overwrite) {
    existing_files <- predicted_files[file.exists(predicted_files)]
    if (length(existing_files) > 0) {
      stop(
        "Output files already exist: ",
        paste(existing_files, collapse = ", "),
        ". Set overwrite = TRUE to replace them.",
        call. = FALSE
      )
    }
  }

  # Create output directory if needed
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }

  # Use atomic write pattern with temporary directory
  temp_dir <- .create_temp_write_dir(path)
  finalized <- NULL

  # A contrast that cannot be written must not destroy outputs that were
  # written successfully. Betas are staged first and are frequently the
  # expensive deliverable, so a contrast failure is recorded here and re-raised
  # only after everything that did succeed has been moved into place. This
  # applies only when writing into a clean destination; see the handler below
  # for why `overwrite = TRUE` keeps the all-or-nothing behaviour.
  contrast_error <- NULL

  created_files <- tryCatch({
    created_files <- list()

    created_files <- .merge_created_files(
      created_files,
      .write_bids_plan(plan$betas, temp_dir, format)
    )

    # Everything staged up to this point is complete and safe to finalize.
    staged_before_contrasts <- .list_staged_files(temp_dir)

    created_files <- .merge_created_files(
      created_files,
      tryCatch(
        .write_bids_plan(plan$contrasts, temp_dir, format),
        error = function(e) {
          # When overwriting, `path` already holds a self-consistent set of
          # outputs from an earlier run. Keeping only the new betas would leave
          # betas and contrasts from different fits side by side, which is a
          # worse failure than losing the write, so let the error propagate and
          # take the full rollback below.
          if (isTRUE(overwrite)) stop(e)

          # A half-written file must never be promoted, so discard whatever the
          # failed stage staged while keeping the completed outputs.
          .discard_staged_since(temp_dir, staged_before_contrasts)
          contrast_error <<- conditionMessage(e)
          list()
        }
      )
    )

    if (produce_gds) {
      if (!requireNamespace("fmrigds", quietly = TRUE)) {
        stop("fmrigds package is required for GDS export. Please install fmrigds.", call. = FALSE)
      }
      gds_files <- .save_gds_outputs(x, temp_dir, entities, desc)
      created_files <- c(created_files, list(gds = gds_files))
    }

    # Atomic move: if all writes succeeded, move temp_dir contents to final location
    finalized <- .finalize_atomic_write(temp_dir, path, created_files, overwrite = overwrite)
    created_files <- finalized$files

    if (produce_gds && "gds" %in% names(created_files)) {
      .write_gds_plan_file(
        gds_h5_final = created_files$gds$h5,
        gds_plan_final = created_files$gds$plan,
        overwrite = overwrite,
        backup_root = file.path(temp_dir, ".overwrite_backup")
      )
    }

    .cleanup_temp_write_dir(temp_dir)

    created_files

  }, error = function(e) {
    if (!is.null(finalized)) {
      .rollback_finalized_write(finalized$transaction)
    }
    .cleanup_temp_write_dir(temp_dir)
    stop("Failed to write BIDS results: ", e$message, call. = FALSE)
  })

  # Raised here, not at the point of failure, so that the outputs that did
  # succeed are already committed to `path`. The call still fails loudly: the
  # requested set of results is incomplete.
  if (!is.null(contrast_error)) {
    stop(
      "Failed to write BIDS contrast results: ", contrast_error,
      .describe_written_outputs(created_files),
      call. = FALSE
    )
  }

  invisible(created_files)
}

#' Files currently staged in an atomic-write directory
#' @keywords internal
#' @noRd
.list_staged_files <- function(temp_dir) {
  list.files(temp_dir, full.names = TRUE, recursive = TRUE)
}

#' Remove files staged since a snapshot, keeping the completed ones
#'
#' Uses `file.remove()` rather than `unlink()`: `unlink()` glob-expands its
#' argument, and `desc` reaches the staged filenames unsanitized, so a `desc`
#' containing `*`, `?` or `[` would let a stale path match completed outputs.
#'
#' Signals an error if any file survives. The caller's invariant is that a
#' partially written file is never promoted to the output directory, and
#' `.finalize_atomic_write()` moves whatever is left in the staging directory.
#' If the invariant cannot be met -- a file held open by a failed HDF5 write on
#' a platform that refuses to unlink open files, say -- the write must fail
#' whole rather than publish a truncated file.
#'
#' @keywords internal
#' @noRd
.discard_staged_since <- function(temp_dir, keep) {
  stale <- setdiff(.list_staged_files(temp_dir), keep)
  if (length(stale) == 0) {
    return(invisible(character(0)))
  }
  removed <- suppressWarnings(file.remove(stale))
  if (!all(removed)) {
    stop(
      "Could not discard partially written file(s): ",
      paste(stale[!removed], collapse = ", "),
      call. = FALSE
    )
  }
  invisible(stale)
}

#' Describe outputs that survived a partial write, for error messages
#' @keywords internal
#' @noRd
.describe_written_outputs <- function(created_files) {
  paths <- unique(as.character(unlist(created_files, use.names = FALSE)))
  paths <- paths[!is.na(paths) & nzchar(paths)]
  if (length(paths) == 0) {
    return(" No output files were written.")
  }
  paste0(
    " The other requested outputs were written and retained: ",
    paste(basename(paths), collapse = ", "), "."
  )
}

#' Plan the image outputs of write_results()
#'
#' One entry per map: the BIDS filename stem (no extension), the slot it is
#' reported under in the returned file list, a function computing the image
#' (a 4-D NeuroVec, its volume labels and the mask) and a function building
#' its JSON sidecar. The same plan predicts the files for the overwrite check
#' and drives the single writer, `.write_bids_plan()`, so prediction and
#' writing cannot drift apart and each image and sidecar is produced once,
#' whatever the number of formats.
#'
#' @return `list(betas = <entries>, contrasts = <entries>)`.
#' @keywords internal
#' @noRd
.bids_output_plan <- function(fmrilm_obj, entities, desc, format, strategy,
                              save_betas, contrasts, contrast_match, contrast_stats) {
  plan <- list(betas = list(), contrasts = list())
  if (length(.image_output_formats(format)) == 0) {
    return(plan)
  }

  # Spatial metadata is looked up once, when the first image is computed.
  spatial_cache <- NULL
  spatial <- function() {
    if (is.null(spatial_cache)) {
      spatial_cache <<- .fmri_dataset_mask_space(fmrilm_obj$dataset, "BIDS export")
    }
    spatial_cache
  }

  if (isTRUE(save_betas)) {
    plan$betas <- list(list(
      slot = "betas",
      stem = .generate_bids_filename(entities, desc = desc, suffix = "betas", extension = NULL),
      context = "raw regressor beta",
      image = function() .compute_beta_volumes(fmrilm_obj, spatial()),
      json = function(labels, output_formats) {
        .betas_json_metadata(fmrilm_obj, entities, labels, output_formats)
      }
    ))
  }

  all_contrasts <- fmrilm_obj$result$contrasts
  if (is.null(all_contrasts) || nrow(all_contrasts) == 0) {
    return(plan)
  }
  if (!is.null(contrasts)) {
    all_contrasts <- .select_contrasts(all_contrasts, contrasts, contrast_match)
    if (nrow(all_contrasts) == 0) {
      warning("None of the specified contrasts found in model")
      return(plan)
    }
  }
  if (strategy == "by_contrast") {
    .check_contrast_label_collisions(all_contrasts$name)
  }

  # t- and F-contrasts are written to separate, family-labelled files.
  for (fam in .resolve_contrast_families(all_contrasts, contrast_stats)) {
    entries <- if (strategy == "by_stat") {
      .bids_by_stat_entries(fmrilm_obj, fam$contrasts, fam$stats, entities, desc, spatial)
    } else {
      .bids_by_contrast_entries(fmrilm_obj, fam$contrasts, fam$stats, entities, desc, spatial)
    }
    plan$contrasts <- c(plan$contrasts, entries)
  }

  plan
}

#' Plan entries for strategy "by_stat": one map per statistic, one volume
#' per contrast
#' @keywords internal
#' @noRd
.bids_by_stat_entries <- function(fmrilm_obj, fam_contrasts, fam_stats, entities, desc, spatial) {
  entries <- list()
  for (stat in fam_stats) {
    has_stat <- any(vapply(fam_contrasts$data, function(entry) {
      .map_stat_name(stat) %in% names(entry)
    }, logical(1)))
    if (!has_stat) next

    entries[[length(entries) + 1L]] <- local({
      stat <- stat
      list(
        slot = stat,
        stem = .generate_bids_filename(entities, desc = .glm_stat_desc(desc, stat),
                                       suffix = "bold", extension = NULL),
        context = paste0("contrast ", stat),
        image = function() {
          sp <- spatial()
          vols <- .compute_statistical_volumes(fam_contrasts, stat, dim(sp$mask_array),
                                               sp$mask_array, sp$space)
          if (is.null(vols)) return(NULL)
          list(neurovec = vols$neurovec, labels = vols$contrast_names, mask = sp$mask)
        },
        json = function(labels, output_formats) {
          .contrasts_json_metadata(fmrilm_obj, entities, stat, labels, fam_contrasts,
                                   output_formats)
        }
      )
    })
  }
  entries
}

#' Plan entries for strategy "by_contrast": one map per contrast, one volume
#' per statistic
#' @keywords internal
#' @noRd
.bids_by_contrast_entries <- function(fmrilm_obj, fam_contrasts, fam_stats, entities, desc, spatial) {
  entries <- list()
  for (i in seq_len(nrow(fam_contrasts))) {
    contrast_data <- fam_contrasts$data[[i]]
    available_stats <- fam_stats[vapply(fam_stats, function(stat) {
      .map_stat_name(stat) %in% names(contrast_data)
    }, logical(1))]
    if (length(available_stats) == 0) next

    entries[[length(entries) + 1L]] <- local({
      contrast_name <- fam_contrasts$name[i]
      contrast_row <- fam_contrasts[i, ]
      stats_i <- available_stats
      data_i <- contrast_data
      stat_desc <- if (length(stats_i) == 1) .glm_stat_desc(desc, stats_i[[1]]) else desc
      list(
        slot = contrast_name,
        stem = .generate_bids_filename(entities, desc = stat_desc, contrast = contrast_name,
                                       suffix = "bold", extension = NULL),
        context = paste0("contrast ", contrast_name),
        image = function() {
          sp <- spatial()
          columns <- lapply(stats_i, function(stat) as.vector(data_i[[.map_stat_name(stat)]]))
          what <- sprintf("Statistic '%s' for contrast '%s'", stats_i, contrast_name)
          list(neurovec = .pack_mask_volumes(columns, sp$mask_array, sp$space, what),
               labels = stats_i, mask = sp$mask)
        },
        json = function(labels, output_formats) {
          .single_contrast_json_metadata(fmrilm_obj, entities, contrast_name, labels,
                                         contrast_row, output_formats)
        }
      )
    })
  }
  entries
}

#' Predict Output Files for Overwrite Check
#' @keywords internal
#' @noRd
.predict_output_files <- function(path, plan, entities, desc, format) {
  extensions <- c(h5 = ".h5", nifti = ".nii.gz")[.image_output_formats(format)]
  stems <- vapply(c(plan$betas, plan$contrasts), function(entry) entry$stem, character(1))
  predicted_files <- file.path(path, as.vector(outer(stems, c(extensions, ".json"), paste0)))

  if ("gds" %in% format) {
    gds_h5 <- .generate_bids_filename(entities, desc = desc, suffix = "gds", extension = "h5")
    gds_plan <- .generate_bids_filename(entities, desc = desc, suffix = "gds", extension = "rds")
    predicted_files <- c(predicted_files, file.path(path, gds_h5), file.path(path, gds_plan))
  }

  unique(predicted_files)
}

#' Image Output Formats
#' @keywords internal
#' @noRd
.image_output_formats <- function(format) {
  intersect(c("h5", "nifti"), format)
}

#' Write planned maps in every requested image format
#'
#' Each entry's image and JSON sidecar are computed once and written once;
#' only the image write is repeated per format.
#' @return A list keyed by entry slot, each `list(h5 =, json =, nifti =)`.
#' @keywords internal
#' @noRd
.write_bids_plan <- function(entries, path, format) {
  output_formats <- .image_output_formats(format)
  created_files <- list()
  if (length(output_formats) == 0) {
    return(created_files)
  }

  for (entry in entries) {
    image <- entry$image()
    if (is.null(image)) next
    files <- list()

    if ("h5" %in% output_formats) {
      files$h5 <- file.path(path, paste0(entry$stem, ".h5"))
      .require_fmristore()
      tryCatch({
        h5_handle <- fmristore::write_labeled_vec(image$neurovec, image$mask, image$labels,
                                                  file = files$h5)
        h5_handle$close_all()
      }, error = function(e) {
        stop("Failed to write LabeledVolumeSet: ", e$message, call. = FALSE)
      })
    }

    if ("nifti" %in% output_formats) {
      files$nifti <- .write_nifti_volume(image$neurovec,
                                         file.path(path, paste0(entry$stem, ".nii.gz")),
                                         paste(entry$context, "NIfTI"))
    }

    files$json <- file.path(path, paste0(entry$stem, ".json"))
    jsonlite::write_json(entry$json(image$labels, output_formats), files$json,
                         pretty = TRUE, auto_unbox = TRUE)

    # Report files in the established order: the primary image, its sidecar,
    # then a NIfTI copy written alongside HDF5.
    files <- files[intersect(c("h5", "json", "nifti"), names(files))]
    if (is.null(files$h5)) files <- files[c("nifti", "json")]
    created_files[[entry$slot]] <- files
  }

  created_files
}

#' Partition Contrasts into t- and F-Families with Family-Appropriate Stats
#'
#' t-contrasts (\code{type == "contrast"}) and F-contrasts
#' (\code{type == "Fcontrast"}) store their test statistic in the same
#' \code{stat} field, but an F statistic is not a t statistic (different
#' distribution and degrees of freedom). This helper splits a contrast table
#' into the two families and resolves which user-facing statistics apply to
#' each, so they are written to separate, correctly labeled files. For
#' F-contrasts the stored \code{estimate}/\code{se} are mean-square components
#' rather than an effect size and standard error, so only the F statistic
#' (\code{fstat}) and its p-value (\code{fpval}) are exported.
#'
#' @return A list of family specs, each \code{list(contrasts = <tbl>,
#'   stats = <character>)}, omitting families with no rows or no applicable
#'   statistics.
#' @keywords internal
#' @noRd
.resolve_contrast_families <- function(available_contrasts, contrast_stats) {
  families <- list()
  if (is.null(available_contrasts) || nrow(available_contrasts) == 0) {
    return(families)
  }

  ctype <- available_contrasts$type
  if (is.null(ctype)) ctype <- rep("contrast", nrow(available_contrasts))
  is_f <- !is.na(ctype) & ctype == "Fcontrast"

  t_rows <- available_contrasts[!is_f, , drop = FALSE]
  if (nrow(t_rows) > 0) {
    t_stats <- contrast_stats[contrast_stats %in% c("beta", "tstat", "se", "pval")]
    if (length(t_stats) > 0) {
      families[[length(families) + 1L]] <- list(contrasts = t_rows, stats = t_stats)
    }
  }

  f_rows <- available_contrasts[is_f, , drop = FALSE]
  if (nrow(f_rows) > 0) {
    f_stats <- character(0)
    # The primary test statistic for an F-contrast is the F value; treat a
    # request for "tstat" (the generic test statistic) or explicit "fstat"
    # as a request for the F statistic.
    if (any(c("tstat", "fstat") %in% contrast_stats)) f_stats <- c(f_stats, "fstat")
    if (any(c("pval", "fpval") %in% contrast_stats)) f_stats <- c(f_stats, "fpval")
    f_stats <- unique(f_stats)
    if (length(f_stats) > 0) {
      families[[length(families) + 1L]] <- list(contrasts = f_rows, stats = f_stats)
    }
  }

  families
}

#' Merge Nested Created File Lists
#' @keywords internal
#' @noRd
.merge_created_files <- function(existing, additions) {
  if (length(additions) == 0) {
    return(existing)
  }
  for (nm in names(additions)) {
    if (!is.null(existing[[nm]]) && is.list(existing[[nm]]) && is.list(additions[[nm]])) {
      for (child in names(additions[[nm]])) {
        existing[[nm]][[child]] <- additions[[nm]][[child]]
      }
    } else {
      existing[[nm]] <- additions[[nm]]
    }
  }
  existing
}

#' Select Contrasts by Exact Name or Regex
#' @keywords internal
#' @noRd
.select_contrasts <- function(contrast_tbl, contrasts = NULL, contrast_match = c("auto", "exact", "regex")) {
  if (is.null(contrasts)) {
    return(contrast_tbl)
  }

  contrast_match <- match.arg(contrast_match)
  contrast_names <- contrast_tbl$name
  keep <- rep(FALSE, length(contrast_names))

  for (selector in contrasts) {
    exact_keep <- contrast_names == selector

    if (contrast_match == "exact") {
      keep <- keep | exact_keep
      next
    }

    if (contrast_match == "auto" && any(exact_keep)) {
      keep <- keep | exact_keep
      next
    }

    regex_keep <- tryCatch(
      grepl(selector, contrast_names, perl = TRUE),
      error = function(e) {
        stop("Invalid contrast regex '", selector, "': ", e$message, call. = FALSE)
      }
    )
    keep <- keep | regex_keep
  }

  contrast_tbl[keep, , drop = FALSE]
}

#' Build Description Entity for Statistic-Specific Maps
#' @keywords internal
#' @noRd
.stat_desc <- function(desc, stat) {
  .sanitize_label(stat)
}

#' Build the desc- label of an fmri_lm statistic map
#'
#' The analysis label is kept so that results written with different `desc`
#' values into one directory do not overwrite each other: `desc-<desc><stat>`,
#' both sanitised to BIDS alphanumerics (e.g. `desc = "run-A"`, `"tstat"` ->
#' `runAtstat`). The default analysis label `"GLM"` (or `NULL`) is omitted,
#' giving `desc-<stat>`: the naming [collect_results()] reads by default.
#' @keywords internal
#' @noRd
.glm_stat_desc <- function(desc, stat, default_desc = "GLM") {
  stat_label <- .sanitize_label(stat)
  if (is.null(desc) || identical(desc, default_desc)) {
    return(stat_label)
  }
  desc_label <- suppressWarnings(.sanitize_label(desc))
  paste0(desc_label %||% "", stat_label)
}

#' Detect contrast names that map to the same BIDS label
#'
#' `contrast-` labels are sanitised to alphanumerics, so distinct names such
#' as `"A_vs_B"` and `"A-vs-B"` would be written to the same file, the second
#' silently replacing the first. Called before anything is written.
#' @keywords internal
#' @noRd
.check_contrast_label_collisions <- function(contrast_names) {
  contrast_names <- unique(as.character(contrast_names))
  if (length(contrast_names) == 0L) return(invisible(TRUE))
  labels <- vapply(contrast_names, function(nm) {
    suppressWarnings(.sanitize_label(nm)) %||% ""
  }, character(1))
  empty <- contrast_names[labels == ""]
  if (length(empty)) {
    stop("Contrast name(s) ", paste0("'", empty, "'", collapse = ", "),
         " contain no letters or digits and cannot form a BIDS 'contrast-' label. ",
         "Rename the contrast(s).", call. = FALSE)
  }
  dup_labels <- unique(labels[duplicated(labels)])
  if (length(dup_labels)) {
    groups <- vapply(dup_labels, function(lab) {
      paste0("'", contrast_names[labels == lab], "'", collapse = " and ")
    }, character(1))
    stop("Contrast names ", paste(groups, collapse = "; "),
         " map to the same BIDS 'contrast-' label (",
         paste(dup_labels, collapse = ", "),
         ") and would overwrite each other with strategy = 'by_contrast'. ",
         "Rename the contrasts or use strategy = 'by_stat'.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Summarise (possibly voxel-specific) degrees of freedom for a JSON sidecar
#'
#' A scalar when every finite value is the same (to rounding), otherwise
#' `list(Min = , Max = )`; `NULL` when nothing finite is available. Writing
#' the full per-voxel vector into a sidecar bloats it without being usable.
#' @keywords internal
#' @noRd
.summarize_df <- function(df) {
  if (is.null(df)) return(NULL)
  df <- as.numeric(df)
  df <- df[is.finite(df)]
  if (!length(df)) return(NULL)
  rng <- range(df)
  if (isTRUE(all.equal(rng[1], rng[2]))) {
    return(rng[1])
  }
  list(Min = rng[1], Max = rng[2])
}

#' Whether the optional fmristore package (HDF5 writer) is available
#'
#' A thin wrapper so the dependency check can be mocked in tests.
#' @keywords internal
#' @noRd
.has_fmristore <- function() {
  requireNamespace("fmristore", quietly = TRUE)
}

#' Fail with a clear message when fmristore is required but missing
#' @keywords internal
#' @noRd
.require_fmristore <- function() {
  if (!.has_fmristore()) {
    stop("format = 'h5' requires the 'fmristore' package, which is not installed. ",
         "Install fmristore or use format = 'nifti'.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Software versions recorded in JSON sidecars
#' @keywords internal
#' @noRd
.bids_software_versions <- function() {
  out <- list(
    R = paste(R.Version()$major, R.Version()$minor, sep = "."),
    fmrireg = as.character(utils::packageVersion("fmrireg"))
  )
  if (.has_fmristore()) {
    out$fmristore <- as.character(utils::packageVersion("fmristore"))
  }
  out$neuroim2 <- as.character(utils::packageVersion("neuroim2"))
  out
}

#' Validate fmri_lm Object Structure
#' @keywords internal
#' @noRd
.validate_fmrilm_object <- function(x) {
  if (!inherits(x, "fmri_lm")) {
    stop("Input must be an 'fmri_lm' object", call. = FALSE)
  }

  if (is.null(x$result)) {
    stop("fmri_lm object has not been fitted (missing $result component)", call. = FALSE)
  }

  if (is.null(x$result$betas)) {
    stop("fmri_lm object missing beta estimates", call. = FALSE)
  }

  if (is.null(x$model)) {
    stop("fmri_lm object missing model specification", call. = FALSE)
  }

  if (is.null(x$dataset)) {
    stop("fmri_lm object missing dataset information", call. = FALSE)
  }
}

#' Create BIDS Entities List
#' @keywords internal
#' @noRd
.create_bids_entities <- function(subject, task, space) {

  entities <- list()

  # Subject: sanitize and store
  if (!is.null(subject)) {
    entities$subject <- .sanitize_label(subject)
  }

  # Task: sanitize and store
  if (!is.null(task)) {
    entities$task <- .sanitize_label(task)
  }

  # Space: sanitize and store
  if (!is.null(space)) {
    entities$space <- .sanitize_label(space)
  }

  return(entities)
}

#' Sanitize Labels for BIDS Compliance
#' @keywords internal
#' @noRd
.sanitize_label <- function(label) {
  if (is.null(label) || is.na(label)) return(NULL)

  # Convert to character and remove invalid characters
  label <- as.character(label)
  # BIDS labels must be alphanumeric only (no hyphens allowed in labels)
  label <- gsub("[^a-zA-Z0-9]", "", label)

  if (nchar(label) == 0) {
    warning("Label sanitization resulted in empty string")
    return(NULL)  # Return NULL instead of empty string to properly trigger validation
  }

  return(label)
}

#' Validate Required BIDS Entities
#' @keywords internal
#' @noRd
.validate_required_entities <- function(entities) {
  if (is.null(entities$subject)) {
    stop("Subject identifier is required. Please provide 'subject' parameter.",
         call. = FALSE)
  }

  if (is.null(entities$task)) {
    stop("Task identifier is required for BIDS compliance. Please provide 'task' parameter.",
         call. = FALSE)
  }

  # Space is optional but recommended
  if (is.null(entities$space)) {
    warning("Spatial reference not specified. Consider providing 'space' parameter for better BIDS compliance.")
  }
}

#' Generate BIDS-Compliant Filename
#' @keywords internal
#' @noRd
.generate_bids_filename <- function(entities, desc = NULL, contrast = NULL,
                                   stat = NULL, suffix = NULL, extension = "h5") {

  parts <- c()

  # Required/standard entities in BIDS order
  if (!is.null(entities$subject)) {
    parts <- c(parts, paste0("sub-", entities$subject))
  }

  if (!is.null(entities$task)) {
    parts <- c(parts, paste0("task-", entities$task))
  }

  if (!is.null(entities$space)) {
    parts <- c(parts, paste0("space-", entities$space))
  }

  if (!is.null(contrast)) {
    parts <- c(parts, paste0("contrast-", .sanitize_label(contrast)))
  }

  if (!is.null(stat)) {
    parts <- c(parts, paste0("stat-", stat))
  }

  if (!is.null(desc)) {
    # desc is a BIDS label: alphanumerics only (e.g. "my_model" -> "mymodel").
    desc <- .sanitize_label(desc)
    if (!is.null(desc)) parts <- c(parts, paste0("desc-", desc))
  }

  # Add suffix
  if (!is.null(suffix)) {
    parts <- c(parts, suffix)
  }

  # Combine parts
  filename <- paste(parts, collapse = "_")

  # Add extension
  if (!is.null(extension)) {
    filename <- paste0(filename, ".", extension)
  }

  return(filename)
}

#' Save GDS representation of subject results
#' @keywords internal
#' @noRd
.save_gds_outputs <- function(fmrilm_obj, path, entities, desc) {
  # Collect beta-level statistics
  betas_tbl <- fmrilm_obj$result$betas$data[[1]]
  extract_component <- function(col) {
    if (is.null(col)) return(NULL)
    if (is.list(col)) {
      col <- col[[1]]
    }
    col
  }
  beta_mat <- extract_component(betas_tbl$estimate)
  se_mat   <- extract_component(betas_tbl$se)
  stat_mat <- extract_component(betas_tbl$stat)
  prob_mat <- extract_component(betas_tbl$prob)

  if (is.null(beta_mat) || is.null(se_mat)) {
    stop("Unable to extract beta estimates or standard errors for GDS export.", call. = FALSE)
  }

  # Determine dimensions and labels
  n_vox <- nrow(beta_mat)
  n_coef <- ncol(beta_mat)
  contrast_names <- colnames(beta_mat) %||% paste0("beta", seq_len(n_coef))
  subject_label <- entities$subject %||% "subject1"
  subject_dimnames <- subject_label

  sample_labels <- rownames(beta_mat)
  if (is.null(sample_labels)) {
    sample_labels <- sprintf("voxel_%0*d", ceiling(log10(n_vox + 1)), seq_len(n_vox))
  }

  # Assemble assays
  beta_array <- array(beta_mat, dim = c(n_vox, 1, n_coef),
                      dimnames = list(sample_labels, subject_dimnames, contrast_names))
  se_array   <- array(se_mat, dim = c(n_vox, 1, n_coef),
                      dimnames = list(sample_labels, subject_dimnames, contrast_names))
  var_array  <- array(se_mat^2, dim = c(n_vox, 1, n_coef),
                      dimnames = list(sample_labels, subject_dimnames, contrast_names))

  assays <- list(beta = beta_array, var = var_array, se = se_array)

  if (!is.null(stat_mat)) {
    t_array <- array(stat_mat, dim = c(n_vox, 1, n_coef),
                     dimnames = list(sample_labels, subject_dimnames, contrast_names))
    assays$t <- t_array
  }
  if (!is.null(prob_mat)) {
    p_array <- array(prob_mat, dim = c(n_vox, 1, n_coef),
                     dimnames = list(sample_labels, subject_dimnames, contrast_names))
    assays$p <- p_array
  }
  df_resid <- .summarize_df(.extract_degrees_of_freedom(fmrilm_obj))

  # Build spatial description from mask
  spatial <- .fmri_dataset_mask_space(fmrilm_obj$dataset, "GDS export")
  mask_array <- spatial$mask_array
  neuro_space <- spatial$space
  affine <- neuroim2::trans(neuro_space)
  voxel_space <- fmrigds::space_voxel(dim = dim(mask_array), affine = affine, mask_bitmap = mask_array)

  # Assemble column metadata
  col_df <- data.frame(subject = subject_label, stringsAsFactors = FALSE)
  rownames(col_df) <- subject_label
  if (!is.null(entities$task)) col_df$task <- entities$task
  if (!is.null(entities$space)) col_df$space <- entities$space

  # Metadata about the model
  model_formula <- tryCatch({
    form <- get_formula(fmrilm_obj$model)
    if (is.null(form)) NULL else as.character(form)
  }, error = function(e) NULL)

  metadata <- list(
    schema_version = "0.1.0",
    model = list(
      formula = model_formula,
    regressors = contrast_names,
    df_residual = df_resid
  ),
    generated_by = list(
      package = "fmrireg",
      version = as.character(utils::packageVersion("fmrireg"))
    ),
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  )

  gds_obj <- fmrigds::as_gds(
    assays,
    space = voxel_space,
    subjects = subject_label,
    contrasts = contrast_names,
    col_data = col_df,
    metadata = metadata
  )

  gds_h5_filename <- .generate_bids_filename(entities, desc = desc, suffix = "gds", extension = "h5")
  gds_h5_path <- file.path(path, gds_h5_filename)

  write_plan <- fmrigds::write_out(gds_obj, path = gds_h5_path, format = "h5")
  written_gds <- fmrigds::compute(write_plan)
  rm(written_gds, write_plan)
  .release_unused_hdf5_handles()

  plan_filename <- .generate_bids_filename(entities, desc = desc, suffix = "gds", extension = "rds")
  plan_path <- file.path(path, plan_filename)

  list(h5 = gds_h5_path, plan = plan_path)
}

#' Release HDF5 handles left for finalization by an optional writer
#' @keywords internal
#' @noRd
.release_unused_hdf5_handles <- function(collector = NULL) {
  # fmrigds' HDF5 writer closes its file, but hdf5r dataset objects created
  # without explicit bindings can remain pending R finalization. POSIX permits
  # renaming such a file; Windows does not. Finalize unreachable R objects and
  # then trigger HDF5's own garbage collector before the atomic move.
  invisible(gc(verbose = FALSE))
  if (is.null(collector)) {
    if (!requireNamespace("hdf5r", quietly = TRUE)) {
      return(invisible(FALSE))
    }
    collector <- hdf5r::h5garbage_collect
  }
  collector()
  invisible(TRUE)
}


.fmri_lm_bids_inference_metadata <- function(fit) {
  vm <- fit$result$variance_model
  if (!inherits(vm, "fmri_lm_variance_model")) return(NULL)
  df <- as.numeric(vm$df_inference)
  finite_df <- df[is.finite(df)]
  list(
    VarianceMethod = vm$method,
    CovarianceScope = vm$covariance_scope,
    DegreesOfFreedomMethod = vm$metadata$df_method %||% fit$result$df$method,
    NominalDegreesOfFreedom = as.numeric(vm$df_nominal)[1L],
    InferenceDegreesOfFreedom = list(
      Min = if (length(finite_df)) min(finite_df) else NULL,
      Max = if (length(finite_df)) max(finite_df) else NULL
    ),
    EstimationScope = vm$metadata$estimation_scope %||% NULL,
    NoiseStructure = vm$metadata$noise$struct %||% NULL,
    RobustMethod = vm$metadata$robust$type %||% NULL,
    MaxLag = vm$metadata$selected_max_lag %||% NULL,
    Taper = vm$metadata$taper %||% NULL
  )
}

#' JSON Metadata for Betas
#' @keywords internal
#' @noRd
.betas_json_metadata <- function(fmrilm_obj, entities, regressor_names, output_formats = "h5") {
  metadata <- list(
    Description = "Raw regressor beta coefficients from GLM analysis",
    Sources = .extract_source_files(fmrilm_obj),
    SoftwareVersions = .bids_software_versions(),
    ModelInfo = list(
      Type = "General Linear Model",
      Formula = as.character(deparse(fmrilm_obj$model$event_model$model_spec$formula_or_list)),
      NumRegressors = as.integer(length(regressor_names))
    ),
    Inference = .fmri_lm_bids_inference_metadata(fmrilm_obj),
    RegressorOrder = as.character(regressor_names),
    DataInfo = .image_data_info(output_formats, "arbitrary (depends on input data scaling)"),
    GeneratedBy = .bids_generated_by(),
    CreationTime = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  )
  if (!is.null(entities$space)) {
    metadata$SpatialReference <- entities$space
  }
  metadata
}

#' JSON Metadata for Contrasts (by stat)
#' @keywords internal
#' @noRd
.contrasts_json_metadata <- function(fmrilm_obj, entities, stat, contrast_names, contrasts,
                                     output_formats = "h5") {
  contrast_definitions <- list()
  for (i in seq_len(nrow(contrasts))) {
    contrast_name <- contrasts$name[i]
    if (contrast_name %in% contrast_names) {
      contrast_definitions[[contrast_name]] <- list(
        Name = contrast_name,
        Type = .extract_contrast_type(contrasts[i, ])
      )
    }
  }

  metadata <- list(
    Description = paste("Contrast", stat, "statistics from GLM analysis"),
    StatisticType = toupper(stat),
    DegreesOfFreedom = .summarize_df(.extract_degrees_of_freedom(fmrilm_obj)),
    Sources = .extract_source_files(fmrilm_obj),
    SoftwareVersions = .bids_software_versions(),
    ModelInfo = list(
      Type = "General Linear Model",
      Formula = as.character(deparse(get_formula(fmrilm_obj$model)))
    ),
    Inference = .fmri_lm_bids_inference_metadata(fmrilm_obj),
    ContrastOrder = contrast_names,
    ContrastDefinitions = contrast_definitions,
    DataInfo = .image_data_info(output_formats, .get_stat_units(stat)),
    GeneratedBy = .bids_generated_by(),
    CreationTime = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  )
  if (!is.null(entities$space)) {
    metadata$SpatialReference <- entities$space
  }
  metadata
}

#' JSON Metadata for a Single Contrast (by contrast)
#' @keywords internal
#' @noRd
.single_contrast_json_metadata <- function(fmrilm_obj, entities, contrast_name, stat_names,
                                           contrast_row, output_formats = "h5") {
  metadata <- list(
    Description = paste("Statistical maps for contrast", contrast_name, "from GLM analysis"),
    ContrastName = contrast_name,
    DegreesOfFreedom = .summarize_df(.extract_degrees_of_freedom(fmrilm_obj)),
    Sources = .extract_source_files(fmrilm_obj),
    SoftwareVersions = .bids_software_versions(),
    ModelInfo = list(
      Type = "General Linear Model",
      Formula = as.character(deparse(get_formula(fmrilm_obj$model)))
    ),
    Inference = .fmri_lm_bids_inference_metadata(fmrilm_obj),
    ContrastInfo = list(
      Name = contrast_name,
      Type = .extract_contrast_type(contrast_row)
    ),
    StatisticOrder = stat_names,
    DataInfo = .image_data_info(output_formats, "various (see StatisticOrder for mapping)"),
    GeneratedBy = .bids_generated_by(),
    CreationTime = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  )
  if (!is.null(entities$space)) {
    metadata$SpatialReference <- entities$space
  }
  metadata
}

#' GeneratedBy block of the JSON sidecars
#' @keywords internal
#' @noRd
.bids_generated_by <- function() {
  list(
    Name = "fmrireg::write_results",
    Version = as.character(utils::packageVersion("fmrireg")),
    CodeURL = "https://github.com/bbuchsbaum/fmrireg"
  )
}

#' Get Units for Statistic Type
#' @keywords internal
#' @noRd
.get_stat_units <- function(stat) {
  switch(stat,
    "beta" = "arbitrary (depends on input data scaling)",
    "tstat" = "t-statistic",
    "pval" = "probability",
    "fpval" = "probability",
    "se" = "standard error",
    "z" = "z-score",
    "f" = "F-statistic",
    "fstat" = "F-statistic",
    "arbitrary"
  )
}

#' Extract Source Files from fmri_lm Object
#' @keywords internal
.extract_source_files <- function(fmrilm_obj) {
  # Try to extract source file information from dataset
  sources <- character(0)

  tryCatch({
    # Check if dataset has source information
    if (!is.null(fmrilm_obj$dataset) && !is.null(attr(fmrilm_obj$dataset, "source_files"))) {
      sources <- attr(fmrilm_obj$dataset, "source_files")
    } else if (!is.null(fmrilm_obj$dataset) && inherits(fmrilm_obj$dataset, "fmri_frame")) {
      # A file-backed frame knows its NIfTI paths through its array source
      descriptor <- fmridataset::assay(.dset_root(fmrilm_obj$dataset))$source
      if (inherits(descriptor, "nifti_array_source") && !is.null(descriptor$uri)) {
        sources <- paste0("bids::", basename(descriptor$uri))
      } else {
        dataset_path <- attr(fmrilm_obj$dataset, "file_path")
        if (!is.null(dataset_path)) {
          # Convert to relative BIDS URI if possible
          sources <- paste0("bids::", basename(dataset_path))
        }
      }
    }
  }, error = function(e) {
    # If extraction fails, leave empty
    sources <<- character(0)
  })

  return(sources)
}

#' Contrast Type Label for JSON Sidecars
#'
#' Every contrast table built at fit time records `type` ("contrast" or
#' "Fcontrast"); rows without it are t-contrasts, as in
#' `.resolve_contrast_families()`.
#' @keywords internal
#' @noRd
.extract_contrast_type <- function(contrast_row) {
  ctype <- contrast_row[["type"]]
  if (length(ctype) >= 1 && identical(as.character(ctype[[1]]), "Fcontrast")) {
    "F-contrast"
  } else {
    "t-contrast"
  }
}

#' Extract Degrees of Freedom from fmri_lm Object
#' @keywords internal
#' @noRd
.extract_degrees_of_freedom <- function(fmrilm_obj) {
  # Prefer the inference df actually used to compute the reported p-values.
  # This may be voxel-specific and non-integer under Satterthwaite inference.
  df_resid <- NULL

  tryCatch({
    if (!is.null(fmrilm_obj$result$df$inference)) {
      df_resid <- as.numeric(fmrilm_obj$result$df$inference)
    } else if (!is.null(fmrilm_obj$result$df.residual)) {
      df_resid <- as.numeric(fmrilm_obj$result$df.residual)
    } else if (!is.null(fmrilm_obj$result$betas) && "df.residual" %in% names(fmrilm_obj$result$betas)) {
      df_resid <- as.numeric(fmrilm_obj$result$betas$df.residual[1])
    } else if (!is.null(fmrilm_obj$result$contrasts) && "df.residual" %in% names(fmrilm_obj$result$contrasts)) {
      df_resid <- as.numeric(fmrilm_obj$result$contrasts$df.residual[1])
    }
  }, error = function(e) {
    # If extraction fails, leave NULL
    df_resid <<- NULL
  })

  return(df_resid)
}

#' Create Temporary Directory for Atomic Writes
#' @keywords internal
#' @noRd
.create_temp_write_dir <- function(base_path) {
  # Unique staging directory inside the destination (same filesystem, so the
  # final move is a rename). tempfile() does not consume the user's RNG
  # stream, unlike the former sample()-based suffix.
  temp_dir <- tempfile(
    pattern = paste0(".tmp_write_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
    tmpdir = base_path
  )

  # Create temporary directory
  if (!dir.create(temp_dir, recursive = TRUE)) {
    stop("Failed to create temporary directory: ", temp_dir, call. = FALSE)
  }

  return(temp_dir)
}

#' Stage an Existing File for Rollback
#' @keywords internal
#' @noRd
.stash_existing_file <- function(file, backup_root, overwrite = FALSE, relative_name = basename(file)) {
  if (!file.exists(file)) {
    return(NULL)
  }
  if (!overwrite) {
    stop("Output file already exists: ", file, call. = FALSE)
  }

  backup_file <- file.path(backup_root, relative_name)
  if (!dir.exists(dirname(backup_file))) {
    dir.create(dirname(backup_file), recursive = TRUE)
  }

  if (!file.rename(file, backup_file)) {
    if (!file.copy(file, backup_file, overwrite = TRUE)) {
      stop("Failed to back up existing file: ", file, call. = FALSE)
    }
    if (!file.remove(file)) {
      stop("Failed to remove existing file after backup: ", file, call. = FALSE)
    }
  }

  backup_file
}

#' Resolve paths relative to an atomic-write staging directory
#' @keywords internal
#' @noRd
.relative_staged_paths <- function(paths, root) {
  normalized_paths <- chartr("\\", "/", paths)
  normalized_root <- sub("/+$", "", chartr("\\", "/", root))
  root_prefix <- paste0(normalized_root, "/")
  is_descendant <- !is.na(normalized_paths) & startsWith(normalized_paths, root_prefix)

  if (any(!is_descendant)) {
    stop(
      "Refusing to finalize a file outside the staged write directory: ",
      paste(paths[!is_descendant], collapse = ", "),
      call. = FALSE
    )
  }

  substring(normalized_paths, nchar(root_prefix) + 1L)
}

#' Finalize Atomic Write by Moving Files
#' @keywords internal
.finalize_atomic_write <- function(temp_dir, final_path, created_files, overwrite = FALSE) {
  # Move all files from temp_dir to final_path
  temp_files <- list.files(temp_dir, full.names = TRUE, recursive = TRUE)
  rel_paths <- .relative_staged_paths(temp_files, temp_dir)
  final_files <- file.path(final_path, rel_paths)
  existing_targets <- final_files[file.exists(final_files)]

  if (!overwrite && length(existing_targets) > 0) {
    stop(
      "Output files already exist: ",
      paste(existing_targets, collapse = ", "),
      ". Set overwrite = TRUE to replace them.",
      call. = FALSE
    )
  }

  backup_root <- file.path(temp_dir, ".overwrite_backup")
  backup_map <- character(0)
  moved_files <- character(0)
  success <- FALSE

  on.exit({
    if (!success) {
      .rollback_finalized_write(list(
        moved_files = moved_files,
        backup_map = backup_map
      ))
    }
  }, add = TRUE)

  for (i in seq_along(temp_files)) {
    temp_file <- temp_files[[i]]
    # Get relative path within temp directory
    rel_path <- rel_paths[[i]]
    final_file <- final_files[[i]]

    # Ensure destination directory exists
    final_dir <- dirname(final_file)
    if (!dir.exists(final_dir)) {
      dir.create(final_dir, recursive = TRUE)
    }

    backup_file <- .stash_existing_file(
      final_file,
      backup_root = backup_root,
      overwrite = overwrite,
      relative_name = rel_path
    )
    if (!is.null(backup_file)) {
      backup_map[[final_file]] <- backup_file
    }

    # Move file
    if (!file.rename(temp_file, final_file)) {
      # If rename fails, try copy + delete
      if (!file.copy(temp_file, final_file, overwrite = FALSE)) {
        stop("Failed to move file from ", temp_file, " to ", final_file, call. = FALSE)
      }
      file.remove(temp_file)
    }
    moved_files <- c(moved_files, final_file)
  }

  # Update file paths in created_files to point to final locations
  updated_files <- .update_file_paths_in_results(created_files, temp_dir, final_path)
  success <- TRUE

  list(
    files = updated_files,
    transaction = list(
      moved_files = moved_files,
      backup_map = backup_map
    )
  )
}

#' Update File Paths in Results After Atomic Move
#' @keywords internal
#' @noRd
.update_file_paths_in_results <- function(created_files, temp_dir, final_path) {
  # Recursively update file paths in the results structure
  for (i in seq_along(created_files)) {
    if (is.list(created_files[[i]])) {
      created_files[[i]] <- .update_file_paths_in_results(created_files[[i]], temp_dir, final_path)
    } else if (is.character(created_files[[i]]) && length(created_files[[i]]) == 1) {
      normalized_path <- chartr("\\", "/", created_files[[i]])
      normalized_root <- sub("/+$", "", chartr("\\", "/", temp_dir))
      if (startsWith(normalized_path, paste0(normalized_root, "/"))) {
        relative_path <- .relative_staged_paths(created_files[[i]], temp_dir)
        created_files[[i]] <- file.path(final_path, relative_path)
      }
    }
  }
  return(created_files)
}

#' Clean Up Temporary Write Directory
#' @keywords internal
#' @noRd
.cleanup_temp_write_dir <- function(temp_dir) {
  if (dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE)
  }
}

#' Roll Back Finalized Writes
#' @keywords internal
#' @noRd
.rollback_finalized_write <- function(transaction) {
  moved_files <- rev(transaction$moved_files %||% character(0))
  backup_map <- transaction$backup_map %||% character(0)

  for (file in moved_files) {
    if (file.exists(file)) {
      file.remove(file)
    }
  }

  if (length(backup_map) > 0) {
    for (final_file in names(backup_map)) {
      backup_file <- backup_map[[final_file]]
      if (file.exists(backup_file)) {
        if (!dir.exists(dirname(final_file))) {
          dir.create(dirname(final_file), recursive = TRUE)
        }
        if (!file.rename(backup_file, final_file)) {
          file.copy(backup_file, final_file, overwrite = TRUE)
          file.remove(backup_file)
        }
      }
    }
  }
}

#' Write GDS Plan File
#' @keywords internal
#' @noRd
.write_gds_plan_file <- function(gds_h5_final, gds_plan_final, overwrite = FALSE, backup_root = NULL) {
  backup_file <- NULL
  if (!is.null(backup_root)) {
    backup_file <- .stash_existing_file(
      gds_plan_final,
      backup_root = backup_root,
      overwrite = overwrite,
      relative_name = basename(gds_plan_final)
    )
  }

  success <- FALSE
  on.exit({
    if (!success && !is.null(backup_file) && file.exists(backup_file)) {
      if (file.exists(gds_plan_final)) {
        file.remove(gds_plan_final)
      }
      if (!file.rename(backup_file, gds_plan_final)) {
        file.copy(backup_file, gds_plan_final, overwrite = TRUE)
        file.remove(backup_file)
      }
    }
  }, add = TRUE)

  disk_plan <- fmrigds::gds(source = gds_h5_final, format = "h5")
  saveRDS(disk_plan, gds_plan_final)
  success <- TRUE

  invisible(gds_plan_final)
}

#' Write a NIfTI Volume
#'
#' Statistical maps are always stored as 32-bit floats (NIfTI datatype
#' FLOAT32), so values are rounded to single precision on disk.
#' @keywords internal
#' @noRd
.write_nifti_volume <- function(volume, filepath, context) {
  tryCatch({
    if (inherits(volume, "NeuroVec")) {
      neuroim2::write_vec(volume, filepath, data_type = "FLOAT")
    } else {
      neuroim2::write_vol(volume, filepath, data_type = "FLOAT")
    }
  }, error = function(e) {
    stop("Failed to write ", context, ": ", e$message, call. = FALSE)
  })
  filepath
}

#' DataInfo Metadata for Image Outputs
#' @keywords internal
#' @noRd
.image_data_info <- function(output_formats, units) {
  output_formats <- .image_output_formats(output_formats)
  if (length(output_formats) == 0) {
    output_formats <- "h5"
  }

  file_formats <- c(h5 = "HDF5", nifti = "NIfTI")
  storage <- c(h5 = "LabeledVolumeSet", nifti = "NIfTI-1 image")

  list(
    FileFormat = unname(file_formats[output_formats]),
    Storage = unname(storage[output_formats]),
    Units = units
  )
}

#' Compute Beta Volumes from fmri_lm Object
#'
#' @param spatial Output of `.fmri_dataset_mask_space()`; looked up from the
#'   fit's dataset when `NULL`.
#' @return `list(neurovec = <4-D NeuroVec, one volume per regressor>,
#'   labels = <regressor names>, mask = <LogicalNeuroVol>)`.
#' @keywords internal
#' @noRd
.compute_beta_volumes <- function(fmrilm_obj, spatial = NULL) {
  # Input validation
  if (is.null(fmrilm_obj$result$betas$data) ||
      length(fmrilm_obj$result$betas$data) == 0 ||
      is.null(fmrilm_obj$result$betas$data[[1]]$estimate) ||
      length(fmrilm_obj$result$betas$data[[1]]$estimate) == 0) {
    stop("Beta data structure is invalid or empty", call. = FALSE)
  }

  # Extract beta data - use coef() to ensure we get all betas including baseline
  beta_data <- coef(fmrilm_obj, type = "betas", include_baseline = TRUE)

  # Additional validation for beta data structure
  if (!is.matrix(beta_data) || nrow(beta_data) == 0 || ncol(beta_data) == 0) {
    stop("Beta data structure is invalid or empty", call. = FALSE)
  }

  spatial <- spatial %||% .fmri_dataset_mask_space(fmrilm_obj$dataset, "beta volume reconstruction")
  n_mask_voxels <- sum(spatial$mask_array)
  if (nrow(beta_data) != n_mask_voxels) {
    stop(
      "Beta data has ", nrow(beta_data),
      " rows but mask contains ", n_mask_voxels, " in-mask voxels.",
      call. = FALSE
    )
  }

  columns <- lapply(seq_len(ncol(beta_data)), function(i) beta_data[, i])
  list(
    neurovec = .pack_mask_volumes(columns, spatial$mask_array, spatial$space),
    labels = colnames(beta_data) %||% paste0("beta_", seq_len(ncol(beta_data))),
    mask = spatial$mask
  )
}

#' Pack in-mask value vectors into a 4-D NeuroVec
#'
#' Volume `j` holds `columns[[j]]` at the in-mask voxels (in mask order) and
#' zero elsewhere.
#' @param what Optional per-column descriptions used in the length error.
#' @keywords internal
#' @noRd
.pack_mask_volumes <- function(columns, mask_array, space, what = NULL) {
  brain_dims <- dim(mask_array)
  mask_indices <- which(mask_array, arr.ind = TRUE)
  n_mask_voxels <- nrow(mask_indices)
  out <- array(0, dim = c(brain_dims, length(columns)))
  for (j in seq_along(columns)) {
    if (length(columns[[j]]) != n_mask_voxels) {
      stop(
        what[j] %||% paste0("Volume ", j), " has ", length(columns[[j]]),
        " values but mask contains ", n_mask_voxels, " in-mask voxels.",
        call. = FALSE
      )
    }
    out[cbind(mask_indices, j)] <- columns[[j]]
  }
  # add_dim extends the 3-D space to 4-D
  neuroim2::NeuroVec(out, neuroim2::add_dim(space, length(columns)))
}

#' Map User Statistic Names to Internal Field Names
#' @keywords internal
#' @noRd
.map_stat_name <- function(stat) {
  # Map user-friendly names to internal contrast data field names.
  # Note: t- and F-contrasts both store their test statistic in the "stat"
  # field and their p-value in "prob"; the family-specific labels ("tstat"
  # vs "fstat", "pval" vs "fpval") only differ in the BIDS desc token and
  # JSON metadata, never in the underlying field name.
  stat_mapping <- c(
    "beta" = "estimate",
    "tstat" = "stat",
    "fstat" = "stat",
    "se" = "se",
    "pval" = "prob",
    "fpval" = "prob",
    "sigma" = "sigma"
  )

  if (stat %in% names(stat_mapping)) {
    return(stat_mapping[[stat]])
  } else {
    return(stat)  # Return as-is if no mapping found
  }
}

#' Compute Statistical Map Volumes from Contrast Data
#' @keywords internal
#' @noRd
.compute_statistical_volumes <- function(contrast_data, stat, brain_dims, mask, space) {
  # Map user statistic name to internal field name
  internal_stat_name <- .map_stat_name(stat)

  stat_values <- list()
  contrast_names <- character(0)
  for (i in seq_len(nrow(contrast_data))) {
    contrast_entry <- contrast_data$data[[i]]
    if (!(internal_stat_name %in% names(contrast_entry))) next
    # Handle both indexed and named access patterns
    value <- contrast_entry[[internal_stat_name]]
    if (is.list(value) && length(value) > 0) value <- value[[1]]
    if (is.null(value)) next
    stat_values[[length(stat_values) + 1L]] <- value
    contrast_names <- c(contrast_names, contrast_data$name[i])
  }

  if (length(stat_values) == 0) {
    return(NULL)
  }

  if (length(mask) != prod(brain_dims)) {
    stop(
      "Mask length (", length(mask),
      ") does not match product of space dims (", prod(brain_dims), ").",
      call. = FALSE
    )
  }
  # as.logical() on a plain array drops dim; restore it so which(arr.ind=TRUE)
  # returns 3D coordinates rather than linear indices.
  mask_array <- array(as.logical(mask), dim = brain_dims)

  what <- sprintf("Statistic '%s' for contrast '%s'", stat, contrast_names)
  list(
    neurovec = .pack_mask_volumes(stat_values, mask_array, space, what),
    contrast_names = contrast_names
  )
}
