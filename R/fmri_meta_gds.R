# Implementation for fmri_meta.group_data_gds (wrapper lives in fmri_meta.R)
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

  method  <- match.arg(method)
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

  # Fallback: if var is absent but se present, compute meta directly
  gcomp <- try(fmrigds::compute(gd), silent = TRUE)
  if (inherits(gcomp, "try-error")) gcomp <- NULL
  an_try <- .gds_safe_assay_names(gcomp %||% gd)
  # Pre-alloc locals for direct path
  beta_local <- se_local <- z_local <- p_local <- NULL
  tau2_local <- I2_local <- Q_local <- df_local <- NULL

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
    z_combined <- vapply(seq_len(ncol(t_mat)), function(j) {
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
    beta_local <- matrix(z_combined, nrow = 1L)
    se_local <- matrix(1, nrow = 1L, ncol = length(z_combined))
    z_local <- beta_local
    p_local <- 2 * stats::pnorm(abs(z_local), lower.tail = FALSE)
    tau2_local <- I2_local <- Q_local <- df_local <-
      rep(NA_real_, length(z_combined))
    res <- NULL
  } else if (length(an_try) > 0 && ("beta" %in% an_try) &&
      (("se" %in% an_try) || ("var" %in% an_try))) {
    if (is.null(gcomp)) gcomp <- fmrigds::compute(gd)
    B <- fmrigds::assay(gcomp, "beta")
    SE <- .gds_safe_assay(gcomp, "se")
    VAR <- .gds_safe_assay(gcomp, "var")
    X <- .gds_safe_model_matrix(gd, formula, fallback_n = dim(B)[2])
    # Assume single contrast or loop across contrasts if present
    if (length(dim(B)) == 3L) {
      Kc <- dim(B)[3]
    } else {
      Kc <- 1L
    }
    # Fit per-contrast and then verify consistency; for Kc==1, just use directly
    Y <- t(first_contrast_matrix(B))
    V <- if (!is.null(VAR)) {
      t(first_contrast_matrix(VAR))
    } else {
      t(first_contrast_matrix(SE)^2)
    }
    if (identical(weights, "equal")) {
      V[] <- 1
    } else if (identical(weights, "custom")) {
      custom <- .fmri_ttest_validate_weights_custom(
        weights_custom, nrow(V), ncol(V)
      )
      V <- if (is.vector(custom)) {
        matrix(1 / custom, nrow = nrow(V), ncol = ncol(V))
      } else {
        1 / custom
      }
    }
    if (identical(return_cov, "tri") && method %in% c("pm","reml")) {
      res_fit <- fmrireg::fmri_meta_fit_cov(Y = Y, V = V, X = X, method = method, robust = robust,
                                            n_threads = getOption("fmrireg.num_threads", 0))
    } else {
      res_fit <- fmrireg::fmri_meta_fit(Y = Y, V = V, X = X, method = method, robust = robust,
                                        n_threads = getOption("fmrireg.num_threads", 0))
    }
    beta_local <- res_fit$beta
    se_local   <- res_fit$se
    z_local    <- res_fit$z
    p_local    <- 2 * stats::pnorm(abs(res_fit$z), lower.tail = FALSE)
    tau2_local <- res_fit$tau2
    I2_local   <- res_fit$I2_fe
    Q_local    <- res_fit$Q_fe
    df_local   <- res_fit$df
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
      # Fallback: manual compute using available assays
      gcomp <- fmrigds::compute(gd)
      B <- fmrigds::assay(gcomp, 'beta')
      Varr <- .gds_safe_assay(gcomp, 'var')
      SEa  <- .gds_safe_assay(gcomp, 'se')
      X <- .gds_safe_model_matrix(gd, formula, fallback_n = dim(B)[2])
      Y <- t(first_contrast_matrix(B))
      V <- if (!is.null(Varr)) {
        t(first_contrast_matrix(Varr))
      } else {
        t(first_contrast_matrix(SEa)^2)
      }
      if (identical(return_cov, "tri") && method %in% c("pm","reml")) {
        rf <- fmrireg::fmri_meta_fit_cov(Y = Y, V = V, X = X, method = method, robust = robust,
                                          n_threads = getOption("fmrireg.num_threads", 0))
      } else {
        rf <- fmrireg::fmri_meta_fit(Y = Y, V = V, X = X, method = method, robust = robust,
                                      n_threads = getOption("fmrireg.num_threads", 0))
      }
      beta_local <- rf$beta
      se_local   <- rf$se
      z_local    <- rf$z
      p_local    <- 2 * stats::pnorm(abs(rf$z), lower.tail = FALSE)
      tau2_local <- rf$tau2
      I2_local   <- rf$I2_fe
      Q_local    <- rf$Q_fe
      df_local   <- rf$df
      res <- NULL
    } else {
      res <- res_try
    }
  }

  # Extract core assays if present
  # Initialize diagnostics to avoid unbound variables in certain paths
  tau2 <- NULL; I2 <- NULL; Q <- NULL; Q_df <- NULL

  # Prefer regression-style param arrays when available
  X_cols <- tryCatch(colnames(fmrigds::model_matrix(data, formula)), error = function(e) NULL)
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
  manual_fit <- NULL
  compute_manual_fit <- function(require_cov = FALSE) {
    if (!require_cov && !is.null(manual_fit)) return(manual_fit)
    if (require_cov && !is.null(manual_fit) && !is.null(manual_fit$cov_tri)) return(manual_fit)
    gcomp_local <- gcomp
    if (is.null(gcomp_local)) {
      gcomp_local <- fmrigds::compute(gd)
      gcomp <<- gcomp_local
    }
    B <- fmrigds::assay(gcomp_local, "beta")
    Varr <- .gds_safe_assay(gcomp_local, "var")
    SEa  <- .gds_safe_assay(gcomp_local, "se")
    n_subj <- dim(B)[2]
    X <- .gds_safe_model_matrix(gd, formula, fallback_n = n_subj)
    Y <- t(first_contrast_matrix(B))
    V <- if (!is.null(Varr)) {
      t(first_contrast_matrix(Varr))
    } else {
      t(first_contrast_matrix(SEa)^2)
    }
    mf <- if (require_cov && method %in% c("pm","reml")) {
      fmrireg::fmri_meta_fit_cov(Y = Y, V = V, X = X, method = method, robust = robust,
                                 n_threads = getOption("fmrireg.num_threads", 0))
    } else {
      fmrireg::fmri_meta_fit(Y = Y, V = V, X = X, method = method, robust = robust,
                             n_threads = getOption("fmrireg.num_threads", 0))
    }
    manual_fit <<- mf
    mf
  }
  diag_missing <- is.null(combine) && (
    (method %in% c("dl", "pm", "reml") && is.null(tau2)) ||
      is.null(I2) || is.null(Q) || is.null(Q_df)
  )
  cov_manual <- NULL
  if (diag_missing || (identical(return_cov, "tri") && is.null(res))) {
    mf <- compute_manual_fit(require_cov = identical(return_cov, "tri"))
    if (method %in% c("dl","pm","reml") && is.null(tau2)) tau2 <- mf$tau2
    if (is.null(I2) && !is.null(mf$I2_fe)) I2 <- mf$I2_fe
    if (is.null(Q) && !is.null(mf$Q_fe)) Q <- mf$Q_fe
    if (is.null(Q_df) && !is.null(mf$df)) Q_df <- mf$df
    if (identical(return_cov, "tri") && !is.null(mf$cov_tri)) cov_manual <- mf$cov_tri
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
  sample_labels <- .fmri_ttest_sample_labels(data, NULL)

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
    model        = list(X = tryCatch(fmrigds::model_matrix(data, formula), error = function(e) NULL),
                        formula = formula),
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

  # Optional covariance triangle
  if (identical(return_cov, "tri")) {
    # Try to reconstruct from param-array cov:<term_i>:<term_j>
    terms <- colnames(out$coefficients) %||% coef_terms
    if (!is.null(terms) && length(terms) > 0) {
      K <- length(terms)
      tsize <- K*(K+1)/2
      # Store as [tsize x P] so downstream matches legacy
      cov_tri <- matrix(NA_real_, nrow = tsize, ncol = nrow(out$coefficients))
      t_idx <- 1L
      for (i in seq_len(K)) {
        for (j in i:K) {
          nm <- paste0("cov:", terms[i], ":", terms[j])
          A <- .gds_safe_assay(res, nm)
          if (!is.null(A)) {
            v <- if (is.array(A)) as.numeric(A[,1,1, drop = TRUE]) else as.numeric(A)
            cov_tri[t_idx, ] <- v
          }
          t_idx <- t_idx + 1L
        }
      }
      if (all(is.finite(cov_tri[1:min(3, nrow(cov_tri)), 1]))) {
        out$cov <- list(type = "tri", tri = cov_tri, coef_names = colnames(out$coefficients))
      } else {
        # Fallback: attempt legacy assay
        cov_assay <- .gds_safe_assay(res, "cov_tri")
        if (!is.null(cov_assay)) {
          cov_assay <- if (is.array(cov_assay)) {
            d <- dim(cov_assay)
            t(cov_assay[, 1, , drop = TRUE])
          } else t(as.matrix(cov_assay))
          out$cov <- list(type = "tri", tri = cov_assay, coef_names = colnames(out$coefficients))
        } else if (!is.null(cov_manual)) {
          out$cov <- list(type = "tri", tri = cov_manual, coef_names = colnames(out$coefficients))
        }
      }
    } else if (!is.null(cov_manual)) {
      out$cov <- list(type = "tri", tri = cov_manual, coef_names = colnames(out$coefficients))
    } else {
      # As a final fallback, compute covariance directly if still missing
      mf <- compute_manual_fit(require_cov = TRUE)
      if (!is.null(mf$cov_tri)) {
        out$cov <- list(type = "tri", tri = mf$cov_tri, coef_names = colnames(out$coefficients))
      }
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

    gcomp_con <- gcomp %||% fmrigds::compute(gd)
    beta_con <- .gds_safe_assay(gcomp_con, "beta")
    var_con <- .gds_safe_assay(gcomp_con, "var")
    se_con <- .gds_safe_assay(gcomp_con, "se")
    if (is.null(beta_con) || (is.null(var_con) && is.null(se_con))) {
      stop("Fit-time contrasts require beta plus variance or standard-error assays.",
           call. = FALSE)
    }
    first_contrast <- function(a) {
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
    Y_con <- t(first_contrast(beta_con))
    V_con <- if (!is.null(var_con)) {
      t(first_contrast(var_con))
    } else {
      t(first_contrast(se_con)^2)
    }
    if (identical(weights, "equal")) V_con[] <- 1

    X_con <- .gds_safe_model_matrix(gd, formula, fallback_n = nrow(Y_con))
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

    contrast_fit <- fmri_meta_fit_contrasts(
      Y = Y_con,
      V = V_con,
      X = X_con,
      Cmat = Cmat,
      method = method,
      robust = robust,
      n_threads = n_threads
    )
    out$contrasts <- list(
      names = contrast_names,
      estimate = t(contrast_fit$c_beta),
      se = t(contrast_fit$c_se),
      z = t(contrast_fit$c_z)
    )
  }

  class(out) <- "fmri_meta"
  out
}
