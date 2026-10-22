# Implementation for fmri_meta.group_data_gds (wrapper lives in fmri_meta.R).
#
# Subject-level effects (beta with se or var) are reduced by the fmrigds
# reducer "fmrireg:meta_reg" (registered in fmrigds_registration.R), fit-time
# contrasts come from fmrigds::coef_contrast(), and the realised result is
# wrapped into an `fmri_meta` object per subject-level contrast. t-only data
# with `combine =` use combine_t_statistics().
.fmri_meta_group_data_gds_impl <- function(data,
                          formula = ~ 1,
                          method = c("pm", "fe", "dl", "reml"),
                          robust = c("none", "huber", "t"),
                          weights = c("ivw", "equal", "custom"),
                          weights_custom = NULL,
                          combine = NULL,
                          contrasts = NULL,
                          return_cov = NULL,
                          chunk_size = 10000,
                          n_threads = .fmrireg_num_threads(),
                          verbose = TRUE) {

  method  <- .meta_resolve_method(match.arg(method))
  robust  <- match.arg(robust)
  weights <- match.arg(weights)
  if (!is.null(return_cov)) return_cov <- match.arg(return_cov, c("tri"))

  if (weights == "custom" && is.null(weights_custom)) {
    stop("weights='custom' requires 'weights_custom'", call. = FALSE)
  }
  if (!is.null(contrasts)) {
    if (!is.null(combine)) {
      stop("Fit-time coefficient contrasts are unavailable for t-only combination models.",
           call. = FALSE)
    }
    if (identical(robust, "t")) {
      stop("GDS-backed fit-time contrasts do not support robust = 't'.",
           call. = FALSE)
    }
  }

  gcomp <- if (inherits(data, "gds")) data else fmrigds::compute(data)
  an <- .gds_safe_assay_names(gcomp)
  sample_labels <- .fmri_ttest_sample_labels(data, NULL)
  X_data <- tryCatch(fmrigds::model_matrix(data, formula), error = function(e) NULL)

  if (!is.null(combine) && "t" %in% an) {
    fits <- .fmri_meta_gds_combine(gcomp, formula, combine, weights, weights_custom)
    coef_terms <- colnames(X_data) %||% "(Intercept)"
  } else {
    if (!("beta" %in% an && ("se" %in% an || "var" %in% an))) {
      stop("fmri_meta() needs subject-level 'beta' with 'se' or 'var' assays ",
           "(or 't' and 'df' with combine =).", call. = FALSE)
    }
    .register_fmrireg_meta_reg()
    res <- fmrigds::compute(fmrigds::reduce(
      gcomp,
      method = .fmrireg_meta_reg_name,
      formula = formula,
      weights = if (identical(weights, "ivw")) "1/var" else weights,
      options = list(method = method, robust = robust,
                     custom_weights = weights_custom, n_threads = n_threads)
    ))
    coef_terms <- sub("^coef:", "", grep("^coef:", names(fmrigds::assays(res)), value = TRUE))
    con <- NULL
    if (!is.null(contrasts)) {
      con <- .fmri_meta_gds_contrast_matrix(contrasts, coef_terms)
      res <- fmrigds::coef_contrast(res, con$C, stat = "z", name = con$names)
    }
    fits <- .fmri_meta_gds_fits(res, coef_terms, con, with_cov = is.null(combine))
  }

  wrap <- function(fit) {
    dn <- list(sample_labels, coef_terms)
    out <- list(
      coefficients = matrix(fit$coef, ncol = length(coef_terms), dimnames = dn),
      se           = matrix(fit$se, ncol = length(coef_terms), dimnames = dn),
      tau2         = fit$tau2,
      I2           = fit$I2,
      Q            = fit$Q,
      Q_df         = fit$Q_df,
      model        = list(X = X_data, formula = formula),
      method       = if (!is.null(combine)) paste0("combine:", tolower(combine)) else method,
      robust       = robust,
      weights      = weights,
      data         = data,
      formula      = formula,
      n_voxels     = nrow(fit$coef),
      n_subjects   = tryCatch(length(fmrigds::subjects(data)), error = function(e) NA_integer_),
      roi_names    = sample_labels
    )
    if (!is.null(fit$cov_tri)) {
      out$cov <- list(type = "tri", tri = fit$cov_tri, coef_names = coef_terms)
    }
    out$contrasts <- fit$contrasts
    class(out) <- "fmri_meta"
    out
  }

  out <- lapply(fits, wrap)
  if (length(out) == 1L) return(out[[1L]])
  structure(out, class = "fmri_meta_list")
}

# Per subject-level contrast k: coefficient, SE (features x terms), the
# heterogeneity vectors, the packed covariance and fit-time contrasts.
.fmri_meta_gds_fits <- function(res, terms, con, with_cov) {
  slice <- function(nm, k) as.numeric(fmrigds::assay(res, nm)[, 1L, k])
  stack <- function(prefix, k) {
    do.call(cbind, lapply(terms, function(t) slice(paste0(prefix, t), k)))
  }
  labels <- res$contrasts %||% "1"
  fits <- lapply(seq_along(labels), function(k) {
    fit <- list(
      coef = stack("coef:", k),
      se = stack("se_coef:", k),
      tau2 = slice("tau2", k),
      I2 = slice("I2", k),
      Q = slice("Q", k),
      Q_df = slice("df_res", k)
    )
    if (with_cov) {
      fit$cov_tri <- fmrigds::coef_cov_tri(res, method = .fmrireg_meta_reg_name,
                                           contrast = k)$cov_tri
    }
    if (!is.null(con)) {
      # A contrast matrix labels the columns; a single vector does not.
      cols <- function(prefix) {
        M <- do.call(cbind, lapply(con$names, function(n) slice(paste0(prefix, n), k)))
        if (con$is_matrix) colnames(M) <- con$names
        M
      }
      fit$contrasts <- list(names = con$names, estimate = cols("contrast:"),
                            se = cols("se_contrast:"), z = cols("z_contrast:"))
    }
    fit
  })
  stats::setNames(fits, labels)
}

# Fit-time contrast weights -> named J x K matrix for fmrigds::coef_contrast().
.fmri_meta_gds_contrast_matrix <- function(contrasts, terms) {
  if (is.matrix(contrasts)) {
    if (ncol(contrasts) != length(terms)) {
      stop("contrasts columns must match predictors", call. = FALSE)
    }
    if (!is.null(colnames(contrasts)) && !identical(colnames(contrasts), terms)) {
      stop("contrast matrix columns must exactly match design-matrix columns",
           call. = FALSE)
    }
    names <- rownames(contrasts) %||% paste0("contrast", seq_len(nrow(contrasts)))
    C <- matrix(as.numeric(contrasts), nrow = nrow(contrasts),
                dimnames = list(names, terms))
    return(list(C = C, names = names, is_matrix = TRUE))
  }
  if (!is.numeric(contrasts)) {
    stop("Unsupported 'contrasts' type; provide a matrix or numeric vector",
         call. = FALSE)
  }
  if (!is.null(names(contrasts))) {
    unknown <- setdiff(names(contrasts), terms)
    if (length(unknown)) {
      stop("Unknown contrast coefficient(s): ", paste(unknown, collapse = ", "),
           call. = FALSE)
    }
    w <- stats::setNames(numeric(length(terms)), terms)
    w[names(contrasts)] <- contrasts
  } else if (length(contrasts) == length(terms)) {
    w <- contrasts
  } else {
    stop("Unnamed contrast vector must have length equal to number of predictors",
         call. = FALSE)
  }
  list(C = matrix(w, nrow = 1L, dimnames = list("contrast1", terms)),
       names = "contrast1", is_matrix = FALSE)
}

# t-only data: combine per-subject t statistics into one z per feature
# (intercept-only model; one fit per subject-level contrast).
.fmri_meta_gds_combine <- function(gcomp, formula, combine, weights, weights_custom) {
  if (!identical(formula, ~ 1)) {
    formula_terms <- stats::terms(formula)
    if (length(attr(formula_terms, "term.labels")) > 0L ||
        !isTRUE(attr(formula_terms, "intercept") == 1L)) {
      stop("combine methods with t-statistics require intercept-only model (~ 1)",
           call. = FALSE)
    }
  }
  # t statistics carry no standard errors, so "ivw" (the default) cannot be
  # honoured; combine with equal weights and say so, as group_data_nifti did.
  if (identical(weights, "ivw")) {
    warning("IVW requested for combine but no SE available; using equal weights",
            call. = FALSE)
    weights <- "equal"
  }
  if (identical(tolower(combine), "fisher") && !identical(weights, "equal")) {
    stop("Fisher combination supports equal weights only; use Lancaster for weighted p-value combination.",
         call. = FALSE)
  }
  t_assay <- fmrigds::assay(gcomp, "t")
  df_assay <- .gds_safe_assay(gcomp, "df")
  if (is.null(df_assay)) {
    stop("t-statistics provided without 'df'; cannot combine", call. = FALSE)
  }
  if (!identical(dim(t_assay), dim(df_assay))) {
    stop("t-statistic and df assays must have identical subject-by-feature dimensions.",
         call. = FALSE)
  }
  labels <- dimnames(t_assay)[[3L]] %||% as.character(seq_len(dim(t_assay)[3L]))
  fits <- lapply(seq_len(dim(t_assay)[3L]), function(k) {
    t_mat <- t(matrix(t_assay[, , k], nrow = dim(t_assay)[1L]))
    df_mat <- t(matrix(df_assay[, , k], nrow = dim(df_assay)[1L]))
    custom <- if (identical(weights, "custom")) {
      .fmri_ttest_validate_weights_custom(weights_custom, nrow(t_mat), ncol(t_mat))
    }
    z <- vapply(seq_len(ncol(t_mat)), function(j) {
      feature_weights <- if (is.null(custom)) {
        NULL
      } else if (is.vector(custom)) {
        custom
      } else {
        custom[, j]
      }
      combine_t_statistics(
        matrix(t_mat[, j], ncol = 1L),
        df = df_mat[, j],
        method = tolower(combine),
        weights = weights,
        weights_custom = feature_weights
      )[[1L]]
    }, numeric(1))
    na <- rep(NA_real_, length(z))
    list(coef = matrix(z, ncol = 1L), se = matrix(1, length(z), 1L),
         tau2 = na, I2 = na, Q = na, Q_df = na)
  })
  stats::setNames(fits, labels)
}

#' @export
print.fmri_meta_list <- function(x, ...) {
  cat("fmri_meta fits for", length(x), "subject-level contrasts:",
      paste(names(x), collapse = ", "), "\n")
  cat("Select one with x[[\"<contrast>\"]].\n")
  invisible(x)
}
