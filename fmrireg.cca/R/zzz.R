# nocov start
.onLoad <- function(libname, pkgname) {
  # Register basis and engine with fmrireg when this extension loads
  try({
    fmrireg::register_basis("cca2", cca_basis)
  }, silent = TRUE)

  try({
    fmrireg::register_engine(
      name = "cca",
      fit = cca_fit,
      preflight = cca_preflight,
      capabilities = list(output = "fmri_lm", needs_volumes = TRUE)
    )
  }, silent = TRUE)
}
# nocov end

# Local helper to avoid depending on host internals
`%||%` <- function(x, y) if (!is.null(x)) x else y

#' @keywords internal
cca_preflight <- function(model, dataset, args, cfg) {
  if (!inherits(dataset, "fmri_dataset")) {
    stop("cca engine requires an 'fmri_dataset'", call. = FALSE)
  }

  # Parse options to determine required capabilities
  if (is.null(args) || !is.list(args)) args <- list()
  mode <- tryCatch(match.arg(args$mode %||% "scale-only", c("scale-only", "full")), error = function(e) "scale-only")

  # For full adaptive mode, require volumetric backend (NeuroVec)
  if (identical(mode, "full")) {
    data_obj <- try(fmridataset::get_data(dataset), silent = TRUE)
    if (inherits(data_obj, "try-error") || !inherits(data_obj, "NeuroVec")) {
      stop("engine='cca' with mode='full' requires a volumetric NeuroVec-backed dataset.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' @keywords internal
cca_fit <- function(model, dataset, args, cfg) {
  # Options with defaults
  if (is.null(args) || !is.list(args)) args <- list()
  opts <- do.call(cca_options, args[intersect(names(args), names(formals(cca_options)))])
  if (is.null(opts)) opts <- cca_options()

  # Try volumetric adaptive path (center-only simple variant); fallback to passthrough
  data_obj <- try(fmridataset::get_data(dataset), silent = TRUE)
  if (inherits(data_obj, "try-error") || !inherits(data_obj, "NeuroVec")) {
    # Fallback: raw Y passthrough
    Y <- as.matrix(fmridataset::get_data_matrix(dataset))
    return(fmrireg::fit_glm_on_transformed_series(
      model, Y, cfg = cfg, dataset = dataset,
      strategy = "engine", engine = "cca"
    ))
  }

  # Volumetric path
  arr4d <- as.array(data_obj)  # x,y,z,t
  if (length(dim(arr4d)) != 4L) {
    Y <- as.matrix(fmridataset::get_data_matrix(dataset))
    return(fmrireg::fit_glm_on_transformed_series(model, Y, cfg = cfg, dataset = dataset,
      strategy = "engine", engine = "cca"))
  }

  # Spacing in mm per axis
  spacing_mm <- try(neuroim2::spacing(data_obj), silent = TRUE)
  if (inherits(spacing_mm, "try-error")) {
    msk <- try(fmridataset::get_mask(dataset), silent = TRUE)
    spacing_mm <- if (!inherits(msk, "try-error")) try(neuroim2::spacing(msk), silent = TRUE) else spacing_mm
  }
  if (inherits(spacing_mm, "try-error") || length(spacing_mm) < 3L) spacing_mm <- c(1,1,1)

  # Mask indices
  mask <- try(fmridataset::get_mask(dataset), silent = TRUE)
  mask_arr <- if (!inherits(mask, "try-error") && inherits(mask, "NeuroVol")) as.array(mask) else mask
  mask_idx <- if (is.null(mask_arr) || length(dim(mask_arr)) != 3L) seq_len(prod(dim(arr4d)[1:3])) else which(as.vector(mask_arr) > 0)

  # Build volumes list for Rcpp (one cube per TR)
  Tlen <- dim(arr4d)[4]
  # Build a list of 3D numeric arrays (one per TR)
  vols <- lapply(seq_len(Tlen), function(t) arr4d[,,,t])

  if (identical(opts$mode, "scale-only")) {
    # Center-only: just smooth (G0) and extract S
    Vlen <- length(mask_idx)
    S <- matrix(0, nrow = Tlen, ncol = Vlen)
    for (t in seq_len(Tlen)) {
      resp <- try(friman_base_responses_3d(vols[[t]], spacing_mm = spacing_mm, fwhm_mm = opts$fwhm), silent = TRUE)
      G0 <- if (inherits(resp, "try-error")) vols[[t]] else resp[["G0"]]
      S[t, ] <- G0[mask_idx]
    }
    # Delegate GLM with AR/robust using host helper on transformed series
    return(fmrireg::fit_glm_with_config(
      model, S, cfg = cfg, dataset = dataset,
      strategy = "engine", engine = "cca"
    ))
  }

  # Full adaptive path: compute Ytilde (T×2) from cca2 basis with PCA+alpha
  X <- fmrireg::design_matrix(model)
  Tlen2 <- nrow(X)
  TR <- dataset$TR %||% 1
  tt <- seq(0, by = TR, length.out = Tlen2)
  basis <- cca_basis(alpha = args$alpha %||% NULL,
                     span = max(tt), TR = TR,
                     library_fun = args$library_fun %||% NULL,
                     pgrid = args$pgrid %||% NULL)
  vals <- fmrihrf::evaluate(basis, tt)
  if (is.null(dim(vals)) || ncol(vals) < 2L) {
    # fallback to two SPMG2 columns
    vals <- fmrihrf::evaluate(fmrihrf::HRF_SPMG2, tt)
    if (is.null(dim(vals)) || ncol(vals) < 2L) stop("Failed to build cca2 temporal basis (2 cols)")
  }
  Ytilde <- cbind(vals[,1], vals[,2])

  # Choose spatial implementation: 2d (fast) or 3d
  spatial_dim <- opts$spatial %||% "3d"
  if (identical(spatial_dim, "2d")) {
    # Prepare mask_z (slice index) for each masked voxel
    nx <- dim(arr4d)[1]; ny <- dim(arr4d)[2]
    mask_z <- as.integer((as.integer(mask_idx) - 1L) %/% (nx * ny))

    # Stream Step-1/2 weights in 2D
    wres <- friman_stream_weights_2d(vols, spacing_mm2d = spacing_mm[1:2], fwhm_mm = opts$fwhm,
                                     Ytilde = Ytilde, mask_lin = as.integer(mask_idx), mask_z = as.integer(mask_z),
                                     shrink = opts$shrink, simplex = TRUE)
    w_dir <- wres$w_dir; w_step2 <- wres$w_step2

    # Decide OLS (suffstats) vs full integrated GLM
    ar_struct <- (cfg$ar$struct %||% cfg$ar_options$cor_struct %||% "iid")
    robust_on <- {
      rb <- cfg$robust; if (is.null(rb)) FALSE else if (is.list(rb)) isTRUE(rb$type != FALSE) else if (is.logical(rb)) isTRUE(rb) else if (is.character(rb)) rb %in% c("huber","bisquare") else FALSE
    }
    ar_on <- !(ar_struct %in% c("iid", "none"))

    if (!ar_on && !robust_on) {
      acc <- friman_pass2_xts_sts_2d(vols, spacing_mm2d = spacing_mm[1:2], fwhm_mm = opts$fwhm,
                                     X = X, w_dir = w_dir, w_step2 = w_step2,
                                     mask_lin = as.integer(mask_idx), mask_z = as.integer(mask_z))
      XtX <- crossprod(X); df <- nrow(X) - qr(X)$rank
      return(fmrireg::fit_glm_from_suffstats(model, XtX = XtX, XtS = acc$XtS, StS = acc$StS, df = df,
                                             cfg = cfg, dataset = dataset, strategy = "engine", engine = "cca"))
    } else {
      S <- friman_apply_series_2d(vols, spacing_mm2d = spacing_mm[1:2], fwhm_mm = opts$fwhm,
                                  w_dir = w_dir, w_step2 = w_step2,
                                  mask_lin = as.integer(mask_idx), mask_z = as.integer(mask_z))
      return(fmrireg::fit_glm_with_config(model, S, cfg = cfg, dataset = dataset, strategy = "engine", engine = "cca"))
    }
  } else {
    # 3D path
    wres <- friman_stream_weights_3d(vols, spacing_mm = spacing_mm, fwhm_mm = opts$fwhm,
                                     Ytilde = Ytilde, mask_linear_idx = as.integer(mask_idx),
                                     shrink = opts$shrink, simplex = TRUE)
    w_dir <- wres$w_dir; w_step2 <- wres$w_step2
    S <- friman_apply_series_3d(vols, spacing_mm = spacing_mm, fwhm_mm = opts$fwhm,
                                w_dir = w_dir, w_step2 = w_step2, mask_linear_idx = as.integer(mask_idx))
    return(fmrireg::fit_glm_with_config(model, S, cfg = cfg, dataset = dataset, strategy = "engine", engine = "cca"))
  }
}
