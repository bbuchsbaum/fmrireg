# Implementation for fmri_meta.group_data_gds (wrapper lives in fmri_meta.R)

# Internal: .meta_fit_core over feature chunks of at most chunk_size columns
# (custom S x P weights are sliced per chunk). Features are independent, so
# the result does not depend on chunk_size.
.gds_meta_fit_chunked <- function(Y, V, X, method, robust, weights, weights_custom,
                                  chunk_size, n_threads) {
  P <- ncol(Y)
  chunk_size <- suppressWarnings(as.integer(chunk_size))
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) chunk_size <- P
  if (P <= chunk_size) {
    return(.meta_fit_core(Y = Y, V = V, X = X, method = method, robust = robust,
                          weights = weights, weights_custom = weights_custom,
                          n_threads = n_threads))
  }
  starts <- seq.int(1L, P, by = chunk_size)
  parts <- lapply(starts, function(st) {
    idx <- st:min(P, st + chunk_size - 1L)
    wc <- if (is.matrix(weights_custom)) weights_custom[, idx, drop = FALSE] else weights_custom
    .meta_fit_core(Y = Y[, idx, drop = FALSE], V = V[, idx, drop = FALSE], X = X,
                   method = method, robust = robust, weights = weights,
                   weights_custom = wc, n_threads = n_threads)
  })
  out <- lapply(stats::setNames(nm = names(parts[[1L]])), function(nm) {
    vals <- lapply(parts, `[[`, nm)
    if (is.matrix(vals[[1L]])) do.call(cbind, vals) else do.call(c, vals)
  })
  out
}
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
                          n_threads = getOption("fmrireg.num_threads", 0),
                          verbose = TRUE) {

  # Ensure our fmrigds registrations are available at first use
  try(.ensure_fmrigds_registered(), silent = TRUE)

  method  <- .meta_resolve_method(match.arg(method))
  robust  <- match.arg(robust)
  weights <- match.arg(weights)
  if (!is.null(return_cov)) return_cov <- match.arg(return_cov, c("tri"))

  if (weights == "custom" && is.null(weights_custom)) {
    stop("weights='custom' requires 'weights_custom'", call. = FALSE)
  }

  # Map fmrireg weight names to fmrigds names
  weights_gds <- switch(weights,
                        ivw = "1/var",
                        equal = "equal",
                        custom = "custom",
                        weights)  # fallback

  # Keep fmrigds weighting semantics; custom weights are passed via options when supported
  gd <- data
  first_contrast_matrix <- function(a) {
    if (length(dim(a)) == 3L) {
      matrix(
        a[, , 1L],
        nrow = dim(a)[1L],
        ncol = dim(a)[2L],
        dimnames = dimnames(a)[1:2]
      )
    } else {
      as.matrix(a)
    }
  }

  reducer <- switch(method,
    fe   = "meta:fe_reg",
    dl   = "meta:re_reg",
    pm   = "meta:re_reg",
    reml = "meta:re_reg"
  )
  if (!is.null(combine)) reducer <- paste0("combine:", tolower(combine))

  opts <- list()
  if (!is.null(return_cov)) opts$return_cov <- return_cov
  if (!identical(robust, "none")) opts$robust <- robust
  if (!is.null(contrasts)) opts$contrasts <- contrasts
  if (method %in% c("dl","pm","reml")) opts$method <- method

  # The group dataset is computed at most once per call; every subject-level
  # fit (coefficients, diagnostics, covariance, fit-time contrasts) reuses
  # one fit of the same Y / V / X.
  gcomp <- try(fmrigds::compute(gd), silent = TRUE)
  gcomp_error <- NULL
  if (inherits(gcomp, "try-error")) {
    gcomp_error <- attr(gcomp, "condition")
    gcomp <- NULL
  }
  require_gcomp <- function() {
    if (is.null(gcomp)) {
      stop(gcomp_error %||% simpleError("could not compute the group dataset"))
    }
    gcomp
  }
  an_try <- .gds_safe_assay_names(gcomp %||% gd)
  X_model <- tryCatch(fmrigds::model_matrix(data, formula), error = function(e) NULL)
  design_matrix <- function(n) {
    if (!is.null(X_model)) X_model else .gds_safe_model_matrix(gd, formula, fallback_n = n)
  }

  manual_fit <- NULL
  compute_manual_fit <- function() {
    if (!is.null(manual_fit)) return(manual_fit)
    g <- require_gcomp()
    B <- fmrigds::assay(g, "beta")
    Varr <- .gds_safe_assay(g, "var")
    SEa  <- .gds_safe_assay(g, "se")
    X <- design_matrix(dim(B)[2])
    Y <- t(first_contrast_matrix(B))
    V <- if (!is.null(Varr)) {
      t(first_contrast_matrix(Varr))
    } else {
      t(first_contrast_matrix(SEa)^2)
    }
    manual_fit <<- .gds_meta_fit_chunked(
      Y, V, X, method = method, robust = robust, weights = weights,
      weights_custom = weights_custom, chunk_size = chunk_size, n_threads = n_threads
    )
    manual_fit
  }
  cov_direct <- NULL

  # Pre-alloc locals for direct path
  beta_local <- se_local <- z_local <- p_local <- NULL
  tau2_local <- I2_local <- Q_local <- df_local <- NULL
  use_manual_fit <- function(rf) {
    cov_direct <<- rf$cov_tri
    beta_local <<- rf$beta
    se_local   <<- rf$se
    z_local    <<- rf$z
    p_local    <<- 2 * stats::pnorm(abs(rf$z), lower.tail = FALSE)
    tau2_local <<- rf$tau2
    I2_local   <<- rf$I2_fe
    Q_local    <<- rf$Q_fe
    df_local   <<- rf$df
  }

  if (!is.null(combine) && "t" %in% an_try) {
    if (!identical(formula, ~ 1)) {
      formula_terms <- stats::terms(formula)
      if (length(attr(formula_terms, "term.labels")) > 0L ||
          !isTRUE(attr(formula_terms, "intercept") == 1L)) {
        stop("combine methods with t-statistics require intercept-only model (~ 1)",
             call. = FALSE)
      }
    }
    if (identical(weights, "ivw")) {
      stop("Inverse-variance weighting is unavailable for t-only combination models; use equal or custom weights.",
           call. = FALSE)
    }
    if (identical(tolower(combine), "fisher") && !identical(weights, "equal")) {
      stop("Fisher combination supports equal weights only; use Lancaster for weighted p-value combination.",
           call. = FALSE)
    }

    t_assay <- .gds_safe_assay(gcomp, "t")
    df_assay <- .gds_safe_assay(gcomp, "df")
    if (is.null(df_assay)) {
      stop("t-statistics provided without 'df'; cannot combine", call. = FALSE)
    }
    t_mat <- t(first_contrast_matrix(t_assay))
    df_mat <- t(first_contrast_matrix(df_assay))
    if (!identical(dim(t_mat), dim(df_mat))) {
      stop("t-statistic and df assays must have identical subject-by-feature dimensions.",
           call. = FALSE)
    }
    custom <- if (identical(weights, "custom")) {
      .fmri_ttest_validate_weights_custom(
        weights_custom, nrow(t_mat), ncol(t_mat)
      )
    } else {
      NULL
    }
    # One vectorised call over all features (subject df may differ per
    # feature; missing subjects are dropped per feature).
    z_combined <- as.numeric(combine_t_statistics(
      t_mat,
      df = df_mat,
      method = tolower(combine),
      weights = weights,
      weights_custom = custom
    ))
    beta_local <- matrix(z_combined, nrow = 1L)
    se_local <- matrix(1, nrow = 1L, ncol = length(z_combined))
    z_local <- beta_local
    p_local <- 2 * stats::pnorm(abs(z_local), lower.tail = FALSE)
    tau2_local <- I2_local <- Q_local <- df_local <-
      rep(NA_real_, length(z_combined))
    res <- NULL
  } else if (length(an_try) > 0 && ("beta" %in% an_try) &&
      (("se" %in% an_try) || ("var" %in% an_try))) {
    use_manual_fit(compute_manual_fit())
    res <- NULL
  } else {
    pl <- fmrigds::as_plan(gd)
    pl <- fmrigds::reduce(pl,
                          method    = reducer,
                          formula   = formula,
                          weights   = weights_gds,
                          options   = opts)
    res_try <- try(fmrigds::compute(pl), silent = TRUE)
    if (inherits(res_try, 'try-error') || isTRUE(is.null(res_try))) {
      # Fallback: manual fit from the subject-level assays
      use_manual_fit(compute_manual_fit())
      res <- NULL
    } else {
      res <- res_try
    }
  }

  # Extract core assays if present
  # Initialize diagnostics to avoid unbound variables in certain paths
  tau2 <- NULL; I2 <- NULL; Q <- NULL; Q_df <- NULL

  # Prefer regression-style param arrays when available
  X_cols <- colnames(X_model)
  term_info <- tryCatch(stats::terms(formula), error = function(e) NULL)
  intercept_formula <- !is.null(term_info) &&
    length(attr(term_info, "term.labels")) == 0 &&
    isTRUE(attr(term_info, "intercept") == 1)
  coef_terms <- X_cols %||% if (intercept_formula) "(Intercept)" else character(0)
  anames <- if (!is.null(res)) names(fmrigds::assays(res)) else character(0)
  if (length(coef_terms) > 0 && all(paste0("coef:", coef_terms) %in% anames)) {
    # Assemble coefficients from param arrays
    get_term <- function(prefix) {
      mats <- lapply(coef_terms, function(tn) .gds_coerce_matrix(.gds_safe_assay(res, paste0(prefix, tn))))
      if (length(mats) == 0 || any(vapply(mats, is.null, logical(1)))) return(NULL)
      M <- do.call(cbind, mats)               # P x K
      t(M)                                    # K x P
    }
    beta <- get_term("coef:")
    se   <- get_term("se_coef:")
    zmat <- get_term("t_coef:")
    pmat <- get_term("p_coef:")
    tau2 <- .gds_safe_assay(res, "tau2")
    I2   <- .gds_safe_assay(res, "I2")
    Q    <- .gds_safe_assay(res, "Q")
    Q_df <- .gds_safe_assay(res, "df")
  } else {
    if (!is.null(beta_local)) {
      beta <- beta_local; se <- se_local; zmat <- z_local; pmat <- p_local
      tau2 <- tau2_local; I2 <- I2_local; Q <- Q_local; Q_df <- df_local
    } else {
      beta <- if (!is.null(res)) .gds_coerce_matrix(.gds_safe_assay(res, "beta")) else NULL
      se   <- if (!is.null(res)) .gds_coerce_matrix(.gds_safe_assay(res, "se"))   else NULL
      zmat <- if (!is.null(res)) .gds_coerce_matrix(.gds_safe_assay(res, "z"))    else NULL
      pmat <- if (!is.null(res)) .gds_coerce_matrix(.gds_safe_assay(res, "p"))    else NULL
      tau2 <- if (!is.null(res)) .gds_safe_assay(res, "tau2") else NULL
      I2   <- if (!is.null(res)) .gds_safe_assay(res, "I2")   else NULL
      Q    <- if (!is.null(res)) .gds_safe_assay(res, "Q")    else NULL
      Q_df <- if (!is.null(res)) {
        .gds_safe_assay(res, "df") %||% .gds_safe_assay(res, "df_res")
      } else NULL
    }
  }
  diag_missing <- is.null(combine) && (
    (method %in% c("dl", "pm", "reml") && is.null(tau2)) ||
      is.null(I2) || is.null(Q) || is.null(Q_df)
  )
  cov_manual <- cov_direct
  if (diag_missing) {
    mf <- compute_manual_fit()
    if (method %in% c("dl","pm","reml") && is.null(tau2)) tau2 <- mf$tau2
    if (is.null(I2) && !is.null(mf$I2_fe)) I2 <- mf$I2_fe
    if (is.null(Q) && !is.null(mf$Q_fe)) Q <- mf$Q_fe
    if (is.null(Q_df) && !is.null(mf$df)) Q_df <- mf$df
    if (is.null(cov_manual) && !is.null(mf$cov_tri)) cov_manual <- mf$cov_tri
  }
  to_num <- function(a) if (is.null(a)) NULL else as.numeric(a)

  # Minimal fmri_meta wrapper
  if (!is.null(beta) && length(coef_terms) == nrow(beta)) {
    rownames(beta) <- coef_terms
    if (!is.null(se) && nrow(se) == length(coef_terms)) rownames(se) <- coef_terms
    if (!is.null(zmat) && nrow(zmat) == length(coef_terms)) rownames(zmat) <- coef_terms
    if (!is.null(pmat) && nrow(pmat) == length(coef_terms)) rownames(pmat) <- coef_terms
  }

  # Orient like legacy: [features x coefficients]
  coef_mat <- if (!is.null(beta)) t(beta) else NULL
  se_mat   <- if (!is.null(se))   t(se)   else NULL
  # Apply dimnames to match legacy (rows = features/ROI, cols = term names)
  # Use the computed dataset for row labels (a gds_plan would be recomputed)
  sample_labels <- .fmri_ttest_sample_labels(
    if (is.null(gcomp) || !is.null(attr(data, "fmrireg_sample_labels", exact = TRUE))) data else gcomp,
    NULL
  )

  if (!is.null(coef_mat)) {
    dimnames(coef_mat) <- list(sample_labels, coef_terms %||% colnames(coef_mat))
  }
  if (!is.null(se_mat)) {
    dimnames(se_mat) <- list(sample_labels, coef_terms %||% colnames(se_mat))
  }

  # Heterogeneity diagnostics: mirror legacy numeric vectors (or NA when unavailable)
  n_features <- if (!is.null(coef_mat)) nrow(coef_mat) else 0
  to_diag <- function(x) {
    if (is.null(x)) {
      if (n_features > 0) rep(NA_real_, n_features) else NULL
    } else {
      vec <- as.numeric(x)
      if (length(vec) == n_features || n_features == 0) vec else rep(vec, length.out = n_features)
    }
  }
  tau2 <- to_diag(tau2)
  I2   <- to_diag(I2)
  Q    <- to_diag(Q)
  Q_df <- to_diag(Q_df)
  if (!is.null(Q) && !is.null(Q_df) && any(is.na(I2))) {
    denom <- pmax(Q, .Machine$double.eps)
    I2_calc <- pmax((Q - Q_df) / denom, 0)
    idx <- is.na(I2)
    if (any(idx)) I2[idx] <- I2_calc[idx]
  }

  out <- list(
    coefficients = coef_mat,
    se           = se_mat,
    tau2         = tau2,
    I2           = I2,
    Q            = Q,
    Q_df         = Q_df,
    model        = list(X = X_model, formula = formula),
    method       = if (!is.null(combine)) paste0("combine:", tolower(combine)) else method,
    robust       = robust,
    weights      = weights,
    data         = data,
    formula      = formula,
    n_voxels     = tryCatch(ncol(beta), error = function(e) NA_integer_),
    n_subjects   = tryCatch(length(fmrigds::subjects(data)), error = function(e) NA_integer_),
    roi_names    = sample_labels
  )

  if (!is.null(out$coefficients)) {
    if (is.null(colnames(out$coefficients)) &&
        length(coef_terms) == ncol(out$coefficients)) {
      colnames(out$coefficients) <- coef_terms
    }
    if (is.null(colnames(out$coefficients)) &&
        intercept_formula && ncol(out$coefficients) == 1L) {
      colnames(out$coefficients) <- "(Intercept)"
    }
    if (!is.null(sample_labels) &&
        is.null(rownames(out$coefficients)) &&
        length(sample_labels) == nrow(out$coefficients)) {
      rownames(out$coefficients) <- sample_labels
    }
  }
  if (!is.null(out$se)) {
    if (is.null(colnames(out$se)) &&
        length(coef_terms) == ncol(out$se)) {
      colnames(out$se) <- coef_terms
    }
    if (is.null(colnames(out$se)) &&
        intercept_formula && ncol(out$se) == 1L) {
      colnames(out$se) <- "(Intercept)"
    }
    if (!is.null(sample_labels) &&
        is.null(rownames(out$se)) &&
        length(sample_labels) == nrow(out$se)) {
      rownames(out$se) <- sample_labels
    }
  }

  design_names <- tryCatch(colnames(out$model$X), error = function(e) NULL)
  if (is.null(design_names) && intercept_formula) {
    n_coef <- if (!is.null(out$coefficients)) {
      ncol(out$coefficients)
    } else if (!is.null(out$se)) {
      ncol(out$se)
    } else {
      1L
    }
    design_names <- rep("(Intercept)", n_coef)
  }
  if (!is.null(out$coefficients) && !is.null(design_names) &&
      is.null(colnames(out$coefficients))) {
    colnames(out$coefficients) <- design_names
  }
  if (!is.null(out$se) && !is.null(design_names) &&
      is.null(colnames(out$se))) {
    colnames(out$se) <- design_names
  }

  # Coefficient covariance (packed upper triangle per feature). It is always
  # attached when available so that contrast() can use exact standard errors.
  if (is.null(combine)) {
    cov_tri <- cov_manual
    if (is.null(cov_tri) && !is.null(res)) {
      # Reducer output: cov:<term_i>:<term_j> param arrays or a cov_tri assay
      terms <- colnames(out$coefficients) %||% coef_terms
      if (!is.null(terms) && length(terms) > 0) {
        K <- length(terms)
        tri_try <- matrix(NA_real_, nrow = K * (K + 1) / 2, ncol = nrow(out$coefficients))
        t_idx <- 1L
        for (i in seq_len(K)) {
          for (j in i:K) {
            A <- .gds_safe_assay(res, paste0("cov:", terms[i], ":", terms[j]))
            if (!is.null(A)) {
              tri_try[t_idx, ] <- if (is.array(A)) as.numeric(A[, 1, 1, drop = TRUE]) else as.numeric(A)
            }
            t_idx <- t_idx + 1L
          }
        }
        if (all(is.finite(tri_try[, 1]))) cov_tri <- tri_try
      }
      if (is.null(cov_tri)) {
        cov_assay <- .gds_safe_assay(res, "cov_tri")
        if (!is.null(cov_assay)) {
          cov_tri <- if (is.array(cov_assay)) {
            t(cov_assay[, 1, , drop = TRUE])
          } else t(as.matrix(cov_assay))
        }
      }
    }
    if (is.null(cov_tri)) {
      mf <- tryCatch(compute_manual_fit(), error = function(e) NULL)
      cov_tri <- mf$cov_tri
    }
    if (!is.null(cov_tri)) {
      out$cov <- list(type = "tri", tri = cov_tri, coef_names = colnames(out$coefficients))
    }
  }

  # Optional contrasts if reducer produced them
  con_est <- .gds_safe_assay(res, "con_est")
  con_se  <- .gds_safe_assay(res, "con_se")
  con_z   <- .gds_safe_assay(res, "con_z")
  if (!is.null(con_est) && !is.null(con_se)) {
    out$contrasts <- list(
      names    = rownames(con_est) %||% paste0("c", seq_len(nrow(con_est))),
      estimate = t(con_est),
      se       = t(con_se),
      z        = if (!is.null(con_z)) t(con_z) else t(con_est / pmax(con_se, .Machine$double.eps))
    )
  }

  # Reducer support for fit-time contrasts is backend-dependent. Compute the
  # declared contrast directly from the same subject-level beta/variance data
  # so the public fmri_meta(..., contrasts=) contract never silently degrades.
  if (!is.null(contrasts)) {
    if (!identical(combine, NULL)) {
      stop("Fit-time coefficient contrasts are unavailable for t-only combination models.",
           call. = FALSE)
    }
    if (identical(robust, "t")) {
      stop("GDS-backed fit-time contrasts do not support robust = 't'.",
           call. = FALSE)
    }

    gcomp_con <- require_gcomp()
    beta_con <- .gds_safe_assay(gcomp_con, "beta")
    var_con <- .gds_safe_assay(gcomp_con, "var")
    se_con <- .gds_safe_assay(gcomp_con, "se")
    if (is.null(beta_con) || (is.null(var_con) && is.null(se_con))) {
      stop("Fit-time contrasts require beta plus variance or standard-error assays.",
           call. = FALSE)
    }

    X_con <- design_matrix(dim(beta_con)[2])
    design_names <- colnames(X_con)
    if (is.null(design_names) && ncol(X_con) == 1L) {
      design_names <- "(Intercept)"
      colnames(X_con) <- design_names
    }
    if (is.matrix(contrasts)) {
      if (ncol(contrasts) != ncol(X_con)) {
        stop("contrasts columns must match predictors", call. = FALSE)
      }
      if (!is.null(colnames(contrasts)) &&
          !identical(colnames(contrasts), design_names)) {
        stop("contrast matrix columns must exactly match design-matrix columns",
             call. = FALSE)
      }
      Cmat <- t(contrasts)
      contrast_names <- rownames(contrasts) %||%
        paste0("contrast", seq_len(nrow(contrasts)))
    } else if (is.numeric(contrasts)) {
      if (!is.null(names(contrasts))) {
        unknown <- setdiff(names(contrasts), design_names)
        if (length(unknown)) {
          stop("Unknown contrast coefficient(s): ", paste(unknown, collapse = ", "),
               call. = FALSE)
        }
        weights_con <- setNames(numeric(ncol(X_con)), design_names)
        weights_con[names(contrasts)] <- contrasts
      } else if (length(contrasts) == ncol(X_con)) {
        weights_con <- contrasts
      } else {
        stop("Unnamed contrast vector must have length equal to number of predictors",
             call. = FALSE)
      }
      Cmat <- matrix(weights_con, ncol = 1L)
      contrast_names <- "contrast1"
    } else {
      stop("Unsupported 'contrasts' type; provide a matrix or numeric vector",
           call. = FALSE)
    }

    # Same subject-level data and model as the coefficient fit: reuse it.
    contrast_core <- compute_manual_fit()
    contrast_fit <- .meta_contrasts_from_tri(contrast_core$beta, contrast_core$cov_tri, Cmat)
    out$contrasts <- list(
      names = contrast_names,
      estimate = t(contrast_fit$estimate),
      se = t(contrast_fit$se),
      z = t(contrast_fit$z)
    )
  }

  class(out) <- "fmri_meta"
  out
}
