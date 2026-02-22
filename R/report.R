#' Render Analysis Reports
#'
#' Generic report interface.
#'
#' @param x Object to report on.
#' @param ... Additional arguments passed to methods.
#' @return Method-specific report output.
#' @export
report <- function(x, ...) {
  UseMethod("report")
}

#' @export
report.default <- function(x, ...) {
  stop(
    "No report() method for class: ",
    paste(class(x), collapse = ", "),
    call. = FALSE
  )
}

#' Render a PDF Report for an `fmri_lm` Fit
#'
#' Builds a publication-ready PDF report via Quarto (Typst backend). The report
#' combines model metadata, design diagnostics, HRF curves, estimate summaries,
#' contrast maps, and peak coordinate tables.
#'
#' @param x A fitted `fmri_lm` object.
#' @param output_file Output PDF path.
#' @param title Report title.
#' @param author Optional author string.
#' @param sections Sections to include. Any subset of:
#'   `"model"`, `"design"`, `"hrf"`, `"estimates"`, `"contrasts"`, `"diagnostics"`.
#' @param brain_map_stat Statistic used for contrast maps:
#'   `"tstat"`, `"estimate"`, or `"prob"`.
#' @param slice_axis Axis for montage slicing (1, 2, or 3).
#' @param n_slices Number of slices to render in contrast maps.
#' @param threshold Optional hard threshold for contrast map rendering. `NULL`
#'   uses proportional alpha overlays.
#' @param bg_vol Optional background `NeuroVol` for overlay plots.
#' @param atlas Optional atlas object for peak labeling. Expected fields:
#'   `atlas`, `ids`, `labels`.
#' @param cluster_thresh Absolute statistic threshold used for peak table
#'   cluster detection.
#' @param min_cluster_size Minimum cluster size (voxels) in peak table.
#' @param max_peaks Maximum number of peaks to report per contrast.
#' @param open If `TRUE`, open the generated PDF after rendering.
#' @param quiet If `TRUE`, suppress Quarto render logs.
#' @param ... Reserved for future use.
#' @return (Invisibly) the output PDF path.
#' @export
#' @method report fmri_lm
report.fmri_lm <- function(
    x,
    output_file = "fmri_lm_report.pdf",
    title = "fMRI GLM Analysis Report",
    author = NULL,
    sections = c("model", "design", "hrf", "estimates", "contrasts", "diagnostics"),
    brain_map_stat = c("tstat", "estimate", "prob"),
    slice_axis = 3L,
    n_slices = 9L,
    threshold = NULL,
    bg_vol = NULL,
    atlas = NULL,
    cluster_thresh = 3.0,
    min_cluster_size = 10L,
    max_peaks = 15L,
    open = interactive(),
    quiet = TRUE,
    ...) {
  brain_map_stat <- match.arg(brain_map_stat)
  section_choices <- c("model", "design", "hrf", "estimates", "contrasts", "diagnostics")
  sections <- unique(match.arg(sections, choices = section_choices, several.ok = TRUE))

  .check_report_deps()

  output_file <- .normalize_pdf_path(output_file)
  tmp_dir <- tempfile("fmrireg-report-")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  fig_dir <- file.path(tmp_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  report_data <- .prepare_report_data(
    x = x,
    sections = sections,
    brain_map_stat = brain_map_stat,
    slice_axis = as.integer(slice_axis),
    n_slices = as.integer(n_slices),
    threshold = threshold,
    bg_vol = bg_vol,
    atlas = atlas,
    cluster_thresh = cluster_thresh,
    min_cluster_size = as.integer(min_cluster_size),
    max_peaks = as.integer(max_peaks),
    fig_dir = fig_dir
  )

  report_data$meta <- list(
    title = title,
    author = author,
    date = as.character(Sys.Date()),
    package_version = as.character(utils::packageVersion("fmrireg"))
  )

  data_file <- file.path(tmp_dir, "report_data.rds")
  saveRDS(report_data, data_file)

  template <- .locate_report_template()
  qmd_file <- file.path(tmp_dir, "fmri_lm_report.qmd")
  ok <- file.copy(template, qmd_file, overwrite = TRUE)
  if (!isTRUE(ok)) {
    stop("Failed to copy report template: ", template, call. = FALSE)
  }

  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  .run_quarto_render(
    qmd_file = qmd_file,
    data_file = data_file,
    output_file = output_file,
    execute_dir = tmp_dir,
    quiet = isTRUE(quiet)
  )

  if (!file.exists(output_file)) {
    stop("Report rendering did not produce output file: ", output_file, call. = FALSE)
  }

  if (isTRUE(open) && interactive()) {
    utils::browseURL(normalizePath(output_file, mustWork = TRUE))
  }

  invisible(output_file)
}

#' @keywords internal
#' @noRd
.check_report_deps <- function() {
  has_cli <- nzchar(.quarto_cli_path())
  has_pkg <- requireNamespace("quarto", quietly = TRUE)
  pkg_cli <- FALSE
  if (has_pkg) {
    pkg_cli <- isTRUE(tryCatch(quarto::quarto_available(), error = function(e) FALSE))
  }

  if (!has_cli && !pkg_cli) {
    stop(
      "Quarto CLI is not available. Install Quarto from https://quarto.org/",
      call. = FALSE
    )
  }

  if (!requireNamespace("tinytable", quietly = TRUE)) {
    stop(
      "Package 'tinytable' is required for report generation. ",
      "Install it with install.packages('tinytable').",
      call. = FALSE
    )
  }
}

#' @keywords internal
#' @noRd
.quarto_cli_path <- function() {
  bin <- Sys.which("quarto")
  if (nzchar(bin)) {
    return(unname(bin))
  }

  # Common macOS installation path.
  mac_bin <- "/Applications/quarto/bin/quarto"
  if (file.exists(mac_bin)) {
    return(mac_bin)
  }

  ""
}

#' @keywords internal
#' @noRd
.run_quarto_render <- function(qmd_file, data_file, output_file, execute_dir, quiet) {
  use_pkg <- requireNamespace("quarto", quietly = TRUE)

  if (isTRUE(use_pkg)) {
    out <- tryCatch({
      quarto::quarto_render(
        input = qmd_file,
        output_file = output_file,
        execute_params = list(data_file = normalizePath(data_file, mustWork = TRUE)),
        execute_dir = execute_dir,
        quiet = quiet
      )
      TRUE
    }, error = function(e) e)

    if (isTRUE(out)) {
      return(invisible(TRUE))
    }
  }

  .run_quarto_render_cli(
    qmd_file = qmd_file,
    data_file = data_file,
    output_file = output_file,
    execute_dir = execute_dir,
    quiet = quiet
  )
}

#' @keywords internal
#' @noRd
.run_quarto_render_cli <- function(qmd_file, data_file, output_file, execute_dir, quiet) {
  quarto_bin <- .quarto_cli_path()
  if (!nzchar(quarto_bin)) {
    stop("Could not find a 'quarto' executable on PATH.", call. = FALSE)
  }

  qmd_in_dir <- basename(qmd_file)
  out_name <- basename(output_file)
  out_file_tmp <- file.path(execute_dir, out_name)

  old_wd <- setwd(execute_dir)
  on.exit(setwd(old_wd), add = TRUE)

  args <- c(
    "render",
    qmd_in_dir,
    "-P", paste0("data_file:", normalizePath(data_file, mustWork = TRUE)),
    "--output", out_name
  )
  if (isTRUE(quiet)) {
    args <- c(args, "--quiet")
  }

  if (isTRUE(quiet)) {
    cli_out <- suppressWarnings(system2(quarto_bin, args, stdout = TRUE, stderr = TRUE))
    status <- attr(cli_out, "status")
    if (is.null(status)) {
      status <- 0L
    }
  } else {
    status <- suppressWarnings(system2(quarto_bin, args, stdout = "", stderr = ""))
    cli_out <- character(0)
  }

  if (!identical(as.integer(status), 0L)) {
    msg <- paste0("quarto CLI failed with exit status ", status, ".")
    if (length(cli_out)) {
      msg <- paste0(msg, "\n", paste(cli_out, collapse = "\n"))
    }
    stop(msg, call. = FALSE)
  }

  if (!file.exists(out_file_tmp)) {
    stop("quarto CLI did not produce expected output: ", out_file_tmp, call. = FALSE)
  }

  if (normalizePath(out_file_tmp, mustWork = TRUE) != normalizePath(output_file, mustWork = FALSE)) {
    ok <- file.copy(out_file_tmp, output_file, overwrite = TRUE)
    if (!isTRUE(ok)) {
      stop("Failed to copy report output to target location: ", output_file, call. = FALSE)
    }
  }

  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.normalize_pdf_path <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (!nzchar(ext)) {
    path <- paste0(path, ".pdf")
  } else if (ext != "pdf") {
    path <- paste0(path, ".pdf")
  }
  normalizePath(path, mustWork = FALSE)
}

#' @keywords internal
#' @noRd
.locate_report_template <- function() {
  template <- system.file("report", "fmri_lm_report.qmd", package = "fmrireg")
  if (nzchar(template) && file.exists(template)) {
    return(template)
  }

  local_template <- file.path("inst", "report", "fmri_lm_report.qmd")
  if (file.exists(local_template)) {
    return(local_template)
  }

  stop("Could not locate report template 'inst/report/fmri_lm_report.qmd'.", call. = FALSE)
}

#' @keywords internal
#' @noRd
.prepare_report_data <- function(
    x,
    sections,
    brain_map_stat,
    slice_axis,
    n_slices,
    threshold,
    bg_vol,
    atlas,
    cluster_thresh,
    min_cluster_size,
    max_peaks,
    fig_dir) {
  spatial <- .is_spatial_dataset(x$dataset)

  out <- list(
    sections = sections,
    is_spatial = spatial,
    model = NULL,
    design = NULL,
    hrf = NULL,
    estimates = NULL,
    contrasts = NULL,
    diagnostics = NULL
  )

  if ("model" %in% sections) {
    out$model <- .report_model_info(x)
  }
  if ("design" %in% sections) {
    out$design <- .report_design_data(x)
  }
  if ("hrf" %in% sections) {
    out$hrf <- .report_hrf_data(x)
  }
  if ("estimates" %in% sections) {
    out$estimates <- .report_estimates(x)
  }
  if ("contrasts" %in% sections) {
    out$contrasts <- .report_contrasts(
      x = x,
      brain_map_stat = brain_map_stat,
      slice_axis = slice_axis,
      n_slices = n_slices,
      threshold = threshold,
      bg_vol = bg_vol,
      atlas = atlas,
      cluster_thresh = cluster_thresh,
      min_cluster_size = min_cluster_size,
      max_peaks = max_peaks,
      fig_dir = fig_dir
    )
  }
  if ("diagnostics" %in% sections) {
    out$diagnostics <- .report_diagnostics(x)
  }

  out
}

#' @keywords internal
#' @noRd
.is_spatial_dataset <- function(dataset) {
  mask <- tryCatch(fmridataset::get_mask(dataset), error = function(e) NULL)
  if (is.null(mask)) {
    return(FALSE)
  }

  if (inherits(mask, c("NeuroVol", "LogicalNeuroVol", "ROIVol"))) {
    return(TRUE)
  }

  d <- tryCatch(dim(mask), error = function(e) NULL)
  if (is.null(d)) {
    return(FALSE)
  }

  length(d) >= 3 && all(d[1:3] > 1)
}

#' @keywords internal
#' @noRd
.report_model_info <- function(x) {
  formula_chr <- tryCatch({
    form <- get_formula(x$model)
    if (is.null(form)) "Unknown" else paste(deparse(form), collapse = "")
  }, error = function(e) "Unknown")

  dataset_class <- paste(class(x$dataset), collapse = ", ")
  spatial <- .is_spatial_dataset(x$dataset)

  etab <- .extract_event_table(x$dataset)
  n_events <- if (is.data.frame(etab)) nrow(etab) else NA_integer_

  n_vox <- tryCatch(nrow(.as_matrix(coef(x, type = "betas"))), error = function(e) NA_integer_)
  n_event_reg <- tryCatch(length(x$result$event_indices), error = function(e) NA_integer_)
  n_baseline_reg <- tryCatch(length(x$result$baseline_indices), error = function(e) NA_integer_)

  df_resid <- .extract_residual_df(x)
  conds <- tryCatch(conditions(x$model$event_model), error = function(e) character(0))
  cond_str <- if (length(conds)) paste(conds, collapse = ", ") else "Unknown"

  ar_desc <- tryCatch({
    ar <- ar_parameters(x, scope = "average")
    if (is.null(ar)) "None"
    else paste(signif(as.numeric(ar), 3), collapse = ", ")
  }, error = function(e) "None")

  kv <- data.frame(
    Parameter = c(
      "Formula",
      "Strategy",
      "Dataset class",
      "Spatial dataset",
      "Voxels",
      "Events",
      "Event regressors",
      "Baseline regressors",
      "Residual df",
      "AR parameters",
      "Conditions"
    ),
    Value = c(
      formula_chr,
      as.character(x$strategy %||% "unknown"),
      dataset_class,
      if (isTRUE(spatial)) "yes" else "no",
      .fmt_int(n_vox),
      .fmt_int(n_events),
      .fmt_int(n_event_reg),
      .fmt_int(n_baseline_reg),
      .fmt_num(df_resid, 1),
      ar_desc,
      cond_str
    ),
    stringsAsFactors = FALSE
  )

  list(table = kv)
}

#' @keywords internal
#' @noRd
.report_design_data <- function(x) {
  dm <- tryCatch(design_matrix(x), error = function(e) NULL)
  if (is.null(dm)) {
    return(list(available = FALSE, matrix_long = data.frame(), corr_long = data.frame()))
  }

  dm <- .as_matrix(dm)
  if (is.null(colnames(dm))) {
    colnames(dm) <- paste0("X", seq_len(ncol(dm)))
  }

  max_rows <- 400L
  max_cols <- 60L

  row_idx <- if (nrow(dm) > max_rows) unique(round(seq(1, nrow(dm), length.out = max_rows))) else seq_len(nrow(dm))
  col_idx <- if (ncol(dm) > max_cols) seq_len(max_cols) else seq_len(ncol(dm))
  dm_sub <- dm[row_idx, col_idx, drop = FALSE]

  matrix_long <- as.data.frame(as.table(dm_sub), stringsAsFactors = FALSE)
  names(matrix_long) <- c("time", "regressor", "value")
  time_chr <- as.character(matrix_long$time)
  time_num <- suppressWarnings(as.integer(time_chr))
  if (all(is.na(time_num))) {
    time_num <- match(time_chr, unique(time_chr))
  }
  matrix_long$time <- time_num
  matrix_long$value <- as.numeric(matrix_long$value)

  corr <- suppressWarnings(stats::cor(dm_sub, use = "pairwise.complete.obs"))
  if (!is.matrix(corr)) {
    corr <- matrix(corr, nrow = 1L, ncol = 1L)
    dimnames(corr) <- list(colnames(dm_sub)[1], colnames(dm_sub)[1])
  }
  corr[upper.tri(corr)] <- NA_real_
  corr_long <- as.data.frame(as.table(corr), stringsAsFactors = FALSE)
  names(corr_long) <- c("x", "y", "corr")
  corr_long$corr <- as.numeric(corr_long$corr)
  corr_long <- corr_long[is.finite(corr_long$corr), , drop = FALSE]

  list(
    available = TRUE,
    n_time = nrow(dm),
    n_regressors = ncol(dm),
    matrix_long = matrix_long,
    corr_long = corr_long
  )
}

#' @keywords internal
#' @noRd
.report_hrf_data <- function(x) {
  hrf_obj <- tryCatch(fitted_hrf(x), error = function(e) NULL)
  if (is.null(hrf_obj) || !length(hrf_obj)) {
    return(list(available = FALSE, data = data.frame()))
  }

  out <- list()
  for (term_name in names(hrf_obj)) {
    term <- hrf_obj[[term_name]]
    if (!is.list(term) || is.null(term$pred)) {
      next
    }

    pred <- .as_matrix(term$pred)
    design_df <- tryCatch(as.data.frame(term$design), error = function(e) NULL)
    if (!nrow(pred)) {
      next
    }

    time_vals <- seq_len(nrow(pred)) - 1
    if (!is.null(design_df)) {
      numeric_cols <- names(design_df)[vapply(design_df, is.numeric, logical(1))]
      if (length(numeric_cols)) {
        time_vals <- design_df[[numeric_cols[1]]]
      }
    }

    for (j in seq_len(ncol(pred))) {
      cname <- colnames(pred)[j]
      if (is.null(cname) || !nzchar(cname)) {
        cname <- paste0("curve_", j)
      }
      cname <- sub("^.*[._]", "", cname)

      out[[length(out) + 1L]] <- data.frame(
        term = term_name,
        time = as.numeric(time_vals),
        condition = as.character(cname),
        amplitude = as.numeric(pred[, j]),
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(out)) {
    return(list(available = FALSE, data = data.frame()))
  }

  hrf_df <- do.call(rbind, out)
  hrf_df <- hrf_df[is.finite(hrf_df$amplitude), , drop = FALSE]
  list(available = nrow(hrf_df) > 0, data = hrf_df)
}

#' @keywords internal
#' @noRd
.report_estimates <- function(x) {
  betas <- tryCatch(.as_matrix(coef(x, type = "betas")), error = function(e) NULL)
  tstats <- tryCatch(.as_matrix(stats(x, type = "estimates")), error = function(e) NULL)
  pvals <- tryCatch(.as_matrix(p_values(x, type = "estimates")), error = function(e) NULL)

  if (is.null(betas) || !ncol(betas)) {
    return(list(available = FALSE, table = data.frame()))
  }

  coef_nms <- colnames(betas)
  if (is.null(coef_nms) || length(coef_nms) != ncol(betas)) {
    coef_nms <- tryCatch(coef_names(x, type = "estimates"), error = function(e) character(0))
    if (length(coef_nms) != ncol(betas)) {
      coef_nms <- paste0("coef_", seq_len(ncol(betas)))
    }
  }

  if (is.null(tstats) || ncol(tstats) != ncol(betas)) {
    tstats <- matrix(NA_real_, nrow = nrow(betas), ncol = ncol(betas))
  }
  if (is.null(pvals) || ncol(pvals) != ncol(betas)) {
    pvals <- matrix(NA_real_, nrow = nrow(betas), ncol = ncol(betas))
  }

  tab <- data.frame(
    Coefficient = coef_nms,
    Mean_Beta = apply(betas, 2L, .safe_mean),
    SD_Beta = apply(betas, 2L, .safe_sd),
    Median_Beta = apply(betas, 2L, .safe_median),
    Mean_T = apply(tstats, 2L, .safe_mean),
    Max_Abs_T = apply(tstats, 2L, .safe_max_abs),
    Pct_P_lt_0_05 = apply(pvals, 2L, function(v) .safe_pct(v < 0.05)),
    stringsAsFactors = FALSE
  )

  list(available = TRUE, table = tab)
}

#' @keywords internal
#' @noRd
.report_contrasts <- function(
    x,
    brain_map_stat,
    slice_axis,
    n_slices,
    threshold,
    bg_vol,
    atlas,
    cluster_thresh,
    min_cluster_size,
    max_peaks,
    fig_dir) {
  ctab <- tryCatch(x$result$contrasts, error = function(e) NULL)
  if (is.null(ctab) || !is.data.frame(ctab) || !nrow(ctab)) {
    return(list(available = FALSE, items = list()))
  }

  is_spatial <- .is_spatial_dataset(x$dataset)
  stat_field <- switch(
    brain_map_stat,
    tstat = "stat",
    estimate = "estimate",
    prob = "prob"
  )

  items <- vector("list", nrow(ctab))
  for (i in seq_len(nrow(ctab))) {
    name <- as.character(ctab$name[i])
    type <- as.character(ctab$type[i])
    cdata <- ctab$data[[i]]

    est <- .extract_contrast_field(cdata, "estimate")
    stat <- .extract_contrast_field(cdata, "stat")
    prob <- .extract_contrast_field(cdata, "prob")

    summary_tbl <- data.frame(
      Metric = c("Mean estimate", "SD estimate", "Mean statistic", "Max |stat|", "% p < 0.05"),
      Value = c(
        .safe_mean(est),
        .safe_sd(est),
        .safe_mean(stat),
        .safe_max_abs(stat),
        .safe_pct(prob < 0.05)
      ),
      stringsAsFactors = FALSE
    )

    map_file <- ""
    peak_table <- .empty_peak_table()
    selected_slices <- integer(0)
    vol_for_peaks <- NULL

    if (isTRUE(is_spatial)) {
      map_vals <- .extract_contrast_field(cdata, stat_field)
      if (!length(map_vals) && stat_field == "stat") {
        map_vals <- .extract_contrast_field(cdata, "fstat")
      }
      vol <- .vector_to_volume(x$dataset, map_vals)
      vol_for_peaks <- .vector_to_volume(x$dataset, stat)

      if (!is.null(vol)) {
        selected_slices <- .pick_peak_slices(vol, along = slice_axis, n = n_slices)
        map_file <- .render_contrast_map(
          vol = vol,
          bg_vol = bg_vol,
          name = name,
          slice_axis = slice_axis,
          n_slices = n_slices,
          threshold = threshold,
          fig_dir = fig_dir
        )
      }
      if (is.null(vol_for_peaks)) {
        vol_for_peaks <- vol
      }
    }

    if (!is.null(vol_for_peaks)) {
      peak_table <- .build_peak_table(
        vol = vol_for_peaks,
        atlas = atlas,
        thresh = cluster_thresh,
        min_size = min_cluster_size,
        max_peaks = max_peaks
      )
    }

    items[[i]] <- list(
      name = name,
      type = type,
      summary = summary_tbl,
      map_file = map_file,
      slices = selected_slices,
      peaks = peak_table
    )
  }

  list(
    available = length(items) > 0,
    items = items
  )
}

#' @keywords internal
#' @noRd
.report_diagnostics <- function(x) {
  ar_avg <- tryCatch(ar_parameters(x, scope = "average"), error = function(e) NULL)
  sigma <- tryCatch(as.numeric(x$result$sigma), error = function(e) numeric(0))

  ar_tbl <- if (!is.null(ar_avg) && length(ar_avg)) {
    data.frame(
      Parameter = paste0("AR", seq_along(ar_avg)),
      Value = as.numeric(ar_avg),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(Parameter = "AR", Value = NA_real_, stringsAsFactors = FALSE)
  }

  sigma_tbl <- if (length(sigma)) {
    sigma <- sigma[is.finite(sigma)]
    data.frame(
      Metric = c("Mean", "SD", "Median", "Min", "Max"),
      Value = c(
        .safe_mean(sigma),
        .safe_sd(sigma),
        .safe_median(sigma),
        .safe_min(sigma),
        .safe_max(sigma)
      ),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(Metric = "No sigma estimates", Value = NA_real_, stringsAsFactors = FALSE)
  }

  list(ar = ar_tbl, sigma = sigma_tbl)
}

#' @keywords internal
#' @noRd
.render_contrast_map <- function(
    vol,
    bg_vol,
    name,
    slice_axis,
    n_slices,
    threshold,
    fig_dir) {
  file_name <- sprintf("contrast_%s.png", .safe_file_component(name))
  out <- file.path(fig_dir, file_name)
  slices <- .pick_peak_slices(vol, along = slice_axis, n = n_slices)

  grDevices::png(out, width = 1400, height = 600, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)

  if (!is.null(bg_vol)) {
    supports_alpha_mode <- "ov_alpha_mode" %in% names(formals(neuroim2::plot_overlay))
    if (supports_alpha_mode) {
      ov_mode <- if (is.null(threshold)) "proportional" else "binary"
      ov_thresh <- if (is.null(threshold)) 0 else as.numeric(threshold)
      print(
        neuroim2::plot_overlay(
          bgvol = bg_vol,
          overlay = vol,
          zlevels = slices,
          along = slice_axis,
          ov_alpha_mode = ov_mode,
          ov_thresh = ov_thresh,
          ncol = min(n_slices, 5L),
          title = name
        )
      )
    } else {
      ov <- vol
      if (!is.null(threshold)) {
        arr <- as.array(ov)
        arr[abs(arr) < threshold] <- NA_real_
        ov <- neuroim2::NeuroVol(arr, neuroim2::space(ov))
      }
      print(
        neuroim2::plot_overlay(
          bgvol = bg_vol,
          overlay = ov,
          zlevels = slices,
          along = slice_axis,
          ov_thresh = 0,
          ncol = min(n_slices, 5L),
          title = name
        )
      )
    }
  } else {
    vv <- vol
    if (!is.null(threshold)) {
      arr <- as.array(vv)
      arr[abs(arr) < threshold] <- NA_real_
      vv <- neuroim2::NeuroVol(arr, neuroim2::space(vv))
    }
    print(
      neuroim2::plot_montage(
        vv,
        zlevels = slices,
        along = slice_axis,
        ncol = min(n_slices, 5L),
        title = name
      )
    )
  }

  out
}

#' @keywords internal
#' @noRd
.pick_peak_slices <- function(vol, along = 3L, n = 9L) {
  arr <- as.array(vol)
  if (length(dim(arr)) < 3) {
    return(seq_len(max(1L, min(n, length(arr)))))
  }
  along <- as.integer(along)
  n <- as.integer(max(1L, n))
  if (along < 1L || along > 3L) {
    along <- 3L
  }

  sums <- apply(abs(arr), along, sum, na.rm = TRUE)
  top_idx <- order(sums, decreasing = TRUE)
  top_idx <- top_idx[sums[top_idx] > 0]

  selected <- sort(head(top_idx, n))
  if (!length(selected) || length(selected) < n) {
    fallback <- unique(round(seq(1L, dim(arr)[along], length.out = n)))
    selected <- sort(unique(c(selected, fallback)))
  }

  as.integer(head(selected, n))
}

#' @keywords internal
#' @noRd
.build_peak_table <- function(vol, atlas = NULL, thresh = 3.0, min_size = 10L, max_peaks = 15L) {
  arr <- as.array(vol)
  if (length(dim(arr)) != 3L) {
    return(.empty_peak_table())
  }

  mask <- is.finite(arr) & (abs(arr) >= thresh)
  if (!any(mask)) {
    return(.empty_peak_table())
  }

  comps <- .connected_components_3d(mask)
  if (!length(comps)) {
    return(.empty_peak_table())
  }

  rows <- list()
  for (cid in seq_along(comps)) {
    vox <- comps[[cid]]
    if (length(vox) < min_size) {
      next
    }

    vals <- arr[vox]
    peak_idx <- vox[which.max(abs(vals))]
    ijk <- arrayInd(peak_idx, dim(arr))[1, ]
    xyz <- .index_to_xyz(vol, peak_idx, ijk)
    label <- .atlas_label(atlas, ijk)

    rows[[length(rows) + 1L]] <- data.frame(
      Cluster = cid,
      X = xyz[1],
      Y = xyz[2],
      Z = xyz[3],
      Peak_Stat = as.numeric(arr[peak_idx]),
      Size = length(vox),
      Label = label,
      stringsAsFactors = FALSE
    )
  }

  if (!length(rows)) {
    return(.empty_peak_table())
  }

  out <- do.call(rbind, rows)
  out <- out[order(abs(out$Peak_Stat), decreasing = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  head(out, max_peaks)
}

#' @keywords internal
#' @noRd
.connected_components_3d <- function(mask) {
  dims <- dim(mask)
  visited <- array(FALSE, dim = dims)
  offsets <- as.matrix(expand.grid(-1:1, -1:1, -1:1))
  offsets <- offsets[rowSums(abs(offsets)) > 0, , drop = FALSE]

  components <- list()
  starts <- which(mask)

  for (start in starts) {
    if (visited[start] || !mask[start]) {
      next
    }

    queue <- c(start)
    head <- 1L
    comp <- integer(0)

    while (head <= length(queue)) {
      lin <- queue[head]
      head <- head + 1L

      if (visited[lin]) {
        next
      }
      visited[lin] <- TRUE
      if (!mask[lin]) {
        next
      }
      comp <- c(comp, lin)

      ijk <- arrayInd(lin, .dim = dims)[1, ]
      neigh <- sweep(offsets, 2L, ijk, "+")
      keep <- neigh[, 1] >= 1 & neigh[, 1] <= dims[1] &
        neigh[, 2] >= 1 & neigh[, 2] <= dims[2] &
        neigh[, 3] >= 1 & neigh[, 3] <= dims[3]

      if (!any(keep)) {
        next
      }

      neigh <- neigh[keep, , drop = FALSE]
      nlin <- neigh[, 1] +
        (neigh[, 2] - 1L) * dims[1] +
        (neigh[, 3] - 1L) * dims[1] * dims[2]
      nlin <- nlin[mask[nlin] & !visited[nlin]]

      if (length(nlin)) {
        queue <- c(queue, nlin)
      }
    }

    if (length(comp)) {
      components[[length(components) + 1L]] <- comp
    }
  }

  components
}

#' @keywords internal
#' @noRd
.extract_event_table <- function(dataset) {
  et <- tryCatch(fmridataset::event_table(dataset), error = function(e) NULL)
  if (is.data.frame(et)) {
    return(et)
  }

  et <- tryCatch(dataset$event_table, error = function(e) NULL)
  if (is.data.frame(et)) {
    return(et)
  }

  NULL
}

#' @keywords internal
#' @noRd
.extract_residual_df <- function(x) {
  df <- tryCatch(x$result$rdf, error = function(e) NA_real_)
  if (length(df) > 1L) {
    df <- df[1]
  }
  suppressWarnings(as.numeric(df))
}

#' @keywords internal
#' @noRd
.extract_contrast_field <- function(cdata, field) {
  if (is.null(cdata) || !is.list(cdata) || !(field %in% names(cdata))) {
    return(numeric(0))
  }
  val <- cdata[[field]]
  if (is.list(val) && length(val)) {
    val <- val[[1]]
  }
  suppressWarnings(as.numeric(val))
}

#' @keywords internal
#' @noRd
.vector_to_volume <- function(dataset, values) {
  if (!length(values)) {
    return(NULL)
  }

  tryCatch({
    mask <- fmridataset::get_mask(dataset)
    sp <- neuroim2::space(mask)
    mask_idx <- which(as.logical(mask))
    arr <- array(NA_real_, dim = dim(sp))
    arr[mask_idx] <- values
    neuroim2::NeuroVol(arr, sp)
  }, error = function(e) NULL)
}

#' @keywords internal
#' @noRd
.index_to_xyz <- function(vol, lin_idx, ijk_fallback) {
  xyz <- tryCatch(neuroim2::index_to_coord(vol, lin_idx), error = function(e) NULL)
  if (!is.null(xyz) && is.matrix(xyz) && ncol(xyz) >= 3) {
    return(as.numeric(xyz[1, 1:3]))
  }
  as.numeric(ijk_fallback)
}

#' @keywords internal
#' @noRd
.atlas_label <- function(atlas, ijk) {
  if (is.null(atlas) || !is.list(atlas)) {
    return(NA_character_)
  }
  if (is.null(atlas$atlas)) {
    return(NA_character_)
  }

  atlas_arr <- tryCatch({
    if (inherits(atlas$atlas, c("NeuroVol", "LogicalNeuroVol"))) {
      as.array(atlas$atlas)
    } else {
      as.array(atlas$atlas)
    }
  }, error = function(e) NULL)

  if (is.null(atlas_arr) || length(dim(atlas_arr)) != 3L) {
    return(NA_character_)
  }

  d <- dim(atlas_arr)
  if (ijk[1] < 1 || ijk[1] > d[1] || ijk[2] < 1 || ijk[2] > d[2] || ijk[3] < 1 || ijk[3] > d[3]) {
    return(NA_character_)
  }

  region_id <- atlas_arr[ijk[1], ijk[2], ijk[3]]
  if (is.na(region_id) || region_id == 0) {
    return(NA_character_)
  }

  ids <- atlas$ids
  labels <- atlas$labels
  if (is.null(ids) || is.null(labels)) {
    return(as.character(region_id))
  }

  idx <- match(region_id, ids)
  if (is.na(idx)) {
    return(as.character(region_id))
  }
  as.character(labels[idx])
}

#' @keywords internal
#' @noRd
.empty_peak_table <- function() {
  data.frame(
    Cluster = integer(0),
    X = numeric(0),
    Y = numeric(0),
    Z = numeric(0),
    Peak_Stat = numeric(0),
    Size = integer(0),
    Label = character(0),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
.safe_file_component <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) "contrast" else x
}

#' @keywords internal
#' @noRd
.as_matrix <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.matrix(x)) {
    return(x)
  }
  if (is.data.frame(x)) {
    return(as.matrix(x))
  }
  if (is.vector(x)) {
    return(matrix(as.numeric(x), ncol = 1L))
  }
  as.matrix(x)
}

#' @keywords internal
#' @noRd
.safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  mean(x)
}

#' @keywords internal
#' @noRd
.safe_sd <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) <= 1L) {
    return(NA_real_)
  }
  stats::sd(x)
}

#' @keywords internal
#' @noRd
.safe_median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  stats::median(x)
}

#' @keywords internal
#' @noRd
.safe_max_abs <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  max(abs(x))
}

#' @keywords internal
#' @noRd
.safe_min <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  min(x)
}

#' @keywords internal
#' @noRd
.safe_max <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  max(x)
}

#' @keywords internal
#' @noRd
.safe_pct <- function(x) {
  x <- as.logical(x)
  x <- x[!is.na(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  100 * mean(x)
}

#' @keywords internal
#' @noRd
.fmt_num <- function(x, digits = 3L) {
  if (length(x) == 0L || !is.finite(x)) {
    return("NA")
  }
  format(round(as.numeric(x), digits = digits), trim = TRUE, scientific = FALSE)
}

#' @keywords internal
#' @noRd
.fmt_int <- function(x) {
  if (length(x) == 0L || !is.finite(x)) {
    return("NA")
  }
  format(as.integer(round(x)), big.mark = ",", scientific = FALSE, trim = TRUE)
}
