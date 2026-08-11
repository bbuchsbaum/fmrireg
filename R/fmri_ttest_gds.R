#' fmri_ttest method for gds-backed group_data
#'
#' Delegates to fmrigds reducers for meta or classical (OLS) engines and wraps
#' outputs into an fmri_ttest_fit-compatible object.
#'
#' @inheritParams fmri_ttest
#' @keywords internal
fmri_ttest.group_data_gds <- function(gd,
                           formula = ~ 1,
                           engine = c("auto", "meta", "classic", "welch"),
                           paired = FALSE,
                           mu0 = 0,
                           contrast = NULL,
                           mc = NULL,
                           alpha = 0.05,
                           sign = c("AminusB", "BminusA"),
                           mask = NULL,
                           voxelwise_cov = NULL,
                           center_voxelwise = TRUE,
                           voxel_name = "voxel_cov",
                           weights = c("ivw", "equal", "custom"),
                           weights_custom = NULL,
                           combine = NULL) {

  # Ensure our fmrigds registrations are available at first use
  try(.ensure_fmrigds_registered(), silent = TRUE)

  engine <- match.arg(engine)
  sign   <- match.arg(sign)
  weights <- match.arg(weights)
  meta_method <- getOption("fmrireg.meta.method", "pm")
  if (!is.null(mc)) mc <- match.arg(mc, c("bh", "by", "spatial_fdr"))

  if (paired) {
    stop("paired=TRUE is not supported for gds-backed data", call. = FALSE)
  }
  if (!isTRUE(all.equal(mu0, 0))) {
    stop("mu0 is not supported for gds-backed data", call. = FALSE)
  }
  if (!is.null(mask)) {
    stop("mask is not supported for gds-backed data", call. = FALSE)
  }
  if (!is.null(voxelwise_cov)) {
    stop("voxelwise_cov is not supported for gds-backed data", call. = FALSE)
  }

  if (weights == "custom") {
    if (is.null(weights_custom)) {
      stop("weights='custom' requires 'weights_custom'", call. = FALSE)
    }
  }

  # Map fmrireg weight names to fmrigds names
  weights_gds <- switch(weights,
                        ivw = "1/var",
                        equal = "equal",
                        custom = "custom",
                        weights)  # fallback

  # Engine auto-detect: prefer meta when beta+se/var are present
  if (engine == "auto") {
    an <- .gds_safe_assay_names(gd)
    if (("beta" %in% an) && ("se" %in% an || "var" %in% an)) {
      engine <- "meta"
    } else {
      engine <- "classic"
    }
  }

  covars <- tryCatch(as.data.frame(fmrigds::col_data(gd)), error = function(e) NULL)
  if (is.null(covars)) {
    n_subjects <- tryCatch(length(fmrigds::subjects(gd)), error = function(e) 0L)
    covars <- data.frame(.row = seq_len(n_subjects))
  }
  n_subjects_gds <- tryCatch(
    length(fmrigds::subjects(gd)), error = function(e) NA_integer_
  )
  X <- .gds_safe_model_matrix(
    gd, formula,
    fallback_n = if (is.finite(n_subjects_gds)) n_subjects_gds else NULL
  )
  coef_terms <- colnames(X) %||% character(0)
  group_info <- if (!is.null(X)) .fmri_ttest_group_term(X, covars) else NULL
  feature_group <- .fmri_ttest_feature_group(gd, NULL)
  sample_labels <- .fmri_ttest_sample_labels(gd, NULL)
  exact_contrast <- NULL

  if (engine == "welch") {
    has_group <- "group" %in% colnames(covars) &&
      is.factor(covars$group) && nlevels(covars$group) == 2L
    if (!has_group || is.null(X) || ncol(X) > 2L) {
      stop("engine='welch' requires a two-level factor named 'group' and no additional covariates.",
           call. = FALSE)
    }
    if (!is.null(contrast)) {
      stop("GDS-backed Welch fits currently test the group coefficient directly; omit 'contrast'.",
           call. = FALSE)
    }
    Y <- .fmri_ttest_materialize_effects(gd, "Welch engine")
    w <- welch_t_cpp(Y, as.integer(covars$group))
    t_value <- as.numeric(w$t)
    df_value <- as.numeric(w$df)
    p_value <- 2 * stats::pt(abs(t_value), df = df_value, lower.tail = FALSE)
    z_value <- stats::qnorm(pmax(1e-300, 1 - p_value / 2)) * sign(t_value)
    out <- list(
      beta = rbind(
        "(Intercept)" = (w$muA + w$muB) / 2,
        group = w$muA - w$muB
      ),
      se = matrix(
        NA_real_, nrow = 2L, ncol = ncol(Y),
        dimnames = list(c("(Intercept)", "group"), NULL)
      ),
      t = rbind("(Intercept)" = NA_real_, group = t_value),
      z = rbind("(Intercept)" = NA_real_, group = z_value),
      p = rbind("(Intercept)" = NA_real_, group = p_value),
      df = rbind("(Intercept)" = NA_real_, group = df_value),
      formula = formula,
      engine = engine,
      n_subjects = nrow(Y),
      n_features = ncol(Y),
      roi_names = sample_labels
    )
  } else if (engine == "classic") {
    if (is.null(X)) {
      stop("Could not construct the GDS subject-level design matrix.", call. = FALSE)
    }
    Y <- .fmri_ttest_materialize_effects(gd, "classic engine")
    ols <- ols_t_cpp(Y, X)
    rownames(ols$beta) <- coef_terms
    rownames(ols$se) <- coef_terms
    rownames(ols$t) <- coef_terms
    df_value <- rep(ols$df, ncol(Y))
    p_value <- matrix(
      NA_real_, nrow = nrow(ols$t), ncol = ncol(ols$t),
      dimnames = dimnames(ols$t)
    )
    for (i in seq_len(nrow(p_value))) {
      p_value[i, ] <- 2 * stats::pt(
        abs(ols$t[i, ]), df = df_value, lower.tail = FALSE
      )
    }
    z_value <- matrix(
      NA_real_, nrow = nrow(ols$t), ncol = ncol(ols$t),
      dimnames = dimnames(ols$t)
    )
    for (i in seq_len(nrow(z_value))) {
      z_value[i, ] <- stats::qnorm(pmax(1e-300, 1 - p_value[i, ] / 2)) *
        sign(ols$t[i, ])
    }
    out <- list(
      beta = ols$beta,
      se = ols$se,
      t = ols$t,
      z = z_value,
      p = p_value,
      df = matrix(
        ols$df, nrow = nrow(ols$t), ncol = ncol(ols$t),
        byrow = TRUE, dimnames = dimnames(ols$t)
      ),
      formula = formula,
      engine = engine,
      n_subjects = nrow(Y),
      n_features = ncol(Y),
      roi_names = sample_labels
    )
  } else if (engine == "meta") {
    raw_contrast_weights <- NULL
    if (!is.null(contrast)) {
      canonical_contrast_weights <- .fmri_ttest_resolve_contrast(
        contrast,
        coef_names = .fmri_ttest_canonical_coef_names(coef_terms, group_info),
        group_info = group_info
      )
      raw_contrast_weights <- .fmri_ttest_raw_contrast_weights(
        canonical_contrast_weights,
        coef_names = coef_terms,
        group_info = group_info,
        target_sign = sign,
        source_sign = "BminusA"
      )
    }

    meta_fit <- fmri_meta(
      gd,
      formula = formula,
      method = meta_method,
      weights = weights,
      weights_custom = weights_custom,
      combine = combine,
      contrasts = if (!is.null(raw_contrast_weights)) {
        matrix(
          raw_contrast_weights,
          nrow = 1L,
          dimnames = list("contrast1", coef_terms)
        )
      } else NULL,
      verbose = FALSE
    )

    out <- list(
      beta = t(meta_fit$coefficients),
      se = t(meta_fit$se),
      z = t(meta_fit$coefficients / meta_fit$se),
      p = t(2 * stats::pnorm(abs(meta_fit$coefficients / meta_fit$se), lower.tail = FALSE)),
      df = matrix(Inf, nrow = ncol(meta_fit$coefficients), ncol = nrow(meta_fit$coefficients)),
      formula = formula,
      engine = engine,
      n_subjects = tryCatch(length(fmrigds::subjects(gd)), error = function(e) NA_integer_),
      n_features = nrow(meta_fit$coefficients),
      roi_names = sample_labels
    )

    if (!is.null(meta_fit$contrasts)) {
      exact_contrast <- list(
        estimate = as.numeric(meta_fit$contrasts$estimate[, 1]),
        se = as.numeric(meta_fit$contrasts$se[, 1]),
        z = as.numeric(meta_fit$contrasts$z[, 1]),
        p = 2 * stats::pnorm(abs(meta_fit$contrasts$z[, 1]), lower.tail = FALSE),
        df = rep(Inf, nrow(meta_fit$contrasts$estimate))
      )
    }
  } else {
    pl <- fmrigds::as_plan(gd)
    if (engine == "meta") {
      pl <- fmrigds::reduce(pl, method = "meta:fe", formula = formula, weights = weights_gds)
    } else {
      pl <- fmrigds::reduce(pl, method = "ols:voxelwise", formula = formula)
    }
    res <- fmrigds::compute(pl)

    beta <- .gds_coerce_matrix(.gds_safe_assay(res, "beta"))
    se   <- .gds_coerce_matrix(.gds_safe_assay(res, "se"))
    tmat <- .gds_coerce_matrix(.gds_safe_assay(res, "t"))
    z    <- .gds_coerce_matrix(.gds_safe_assay(res, "z"))
    p    <- .gds_coerce_matrix(.gds_safe_assay(res, "p"))
    df_assay <- .gds_safe_assay(res, "df")
    anames <- .gds_safe_assay_names(res)

    if (length(coef_terms) > 0 && all(paste0("coef:", coef_terms) %in% anames)) {
      get_term_matrix <- function(prefix) {
        mats <- lapply(coef_terms, function(tn) .gds_coerce_matrix(.gds_safe_assay(res, paste0(prefix, tn))))
        if (length(mats) == 0 || any(vapply(mats, is.null, logical(1)))) return(NULL)
        do.call(cbind, mats)
      }
      beta_coef <- get_term_matrix("coef:")
      se_coef   <- get_term_matrix("se_coef:")
      t_coef    <- get_term_matrix("t_coef:")
      p_coef    <- get_term_matrix("p_coef:")
      if (is.null(beta) && !is.null(beta_coef)) beta <- beta_coef
      if (is.null(se)   && !is.null(se_coef))   se   <- se_coef
      if (is.null(tmat) && !is.null(t_coef))    tmat <- t_coef
      if (is.null(p)    && !is.null(p_coef))    p    <- p_coef
    }
    if (is.null(z) && !is.null(tmat)) z <- tmat

    df_mat <- NULL
    if (!is.null(df_assay)) {
      v <- if (is.array(df_assay)) {
        d <- dim(df_assay)
        if (length(d) == 3) as.numeric(df_assay[, 1, 1]) else as.numeric(df_assay)
      } else as.numeric(df_assay)
      K_est <- tryCatch(length(coef_terms), error = function(e) NA_integer_)
      if (!is.finite(K_est) || is.na(K_est) || K_est < 1) {
        K_est <- if (!is.null(beta)) nrow(t(beta)) else 1L
      }
      P_est <- tryCatch({
        if (!is.null(beta)) nrow(beta) else if (!is.null(z)) nrow(z) else if (!is.null(p)) nrow(p) else length(v)
      }, error = function(e) length(v))
      df_mat <- matrix(v, nrow = K_est, ncol = P_est, byrow = TRUE)
    }

    out <- list(
      beta = if (!is.null(beta)) t(beta) else NULL,
      se   = if (!is.null(se))   t(se)   else NULL,
      t    = if (!is.null(tmat)) t(tmat) else NULL,
      z    = if (!is.null(z))    t(z)    else NULL,
      p    = if (!is.null(p))    t(p)    else NULL,
      df   = df_mat,
      formula = formula,
      engine  = engine,
      n_subjects = tryCatch(length(fmrigds::subjects(gd)), error = function(e) NA_integer_),
      n_features = tryCatch(ncol(beta %||% z %||% p), error = function(e) NA_integer_),
      roi_names = sample_labels
    )
  }

  # Apply row/col names for readability and symmetry
  if (!is.null(out$beta) && is.matrix(out$beta) && length(coef_terms) == nrow(out$beta)) rownames(out$beta) <- coef_terms
  if (!is.null(out$se)   && is.matrix(out$se)   && length(coef_terms) == nrow(out$se))   rownames(out$se)   <- coef_terms
  if (!is.null(out$t)    && is.matrix(out$t)    && length(coef_terms) == nrow(out$t))    rownames(out$t)    <- coef_terms
  if (!is.null(out$z)    && is.matrix(out$z)    && length(coef_terms) == nrow(out$z))    rownames(out$z)    <- coef_terms
  if (!is.null(out$p)    && is.matrix(out$p)    && length(coef_terms) == nrow(out$p))    rownames(out$p)    <- coef_terms
  if (!is.null(out$df)   && is.matrix(out$df)   && length(coef_terms) == nrow(out$df))   rownames(out$df)   <- coef_terms
  if (!is.null(sample_labels)) {
    if (!is.null(out$beta) && is.matrix(out$beta) && length(sample_labels) == ncol(out$beta)) colnames(out$beta) <- sample_labels
    if (!is.null(out$se)   && is.matrix(out$se)   && length(sample_labels) == ncol(out$se))   colnames(out$se)   <- sample_labels
    if (!is.null(out$t)    && is.matrix(out$t)    && length(sample_labels) == ncol(out$t))    colnames(out$t)    <- sample_labels
    if (!is.null(out$z)    && is.matrix(out$z)    && length(sample_labels) == ncol(out$z))    colnames(out$z)    <- sample_labels
    if (!is.null(out$p)    && is.matrix(out$p)    && length(sample_labels) == ncol(out$p))    colnames(out$p)    <- sample_labels
    if (!is.null(out$df)   && is.matrix(out$df)   && length(sample_labels) == ncol(out$df))   colnames(out$df)   <- sample_labels
  }

  out <- .fmri_ttest_normalize_group_rows(out, group_info)

  if (!is.null(contrast) && is.null(exact_contrast)) {
    weights <- .fmri_ttest_resolve_contrast(
      contrast,
      coef_names = rownames(out$beta),
      group_info = group_info
    )
    exact_contrast <- .fmri_ttest_single_coef_contrast(out, weights)
    if (is.null(exact_contrast)) {
      stop("gds-backed contrasts currently support only single coefficients", call. = FALSE)
    }
  }

  if (!is.null(exact_contrast)) {
    out <- .fmri_ttest_store_contrast(out, exact_contrast)
  }
  source_sign <- if (identical(engine, "welch")) "AminusB" else "BminusA"
  out <- .fmri_ttest_apply_group_sign(
    out, group_info, target_sign = sign, source_sign = source_sign
  )
  out <- .fmri_ttest_apply_mc(out, mc, alpha, feature_group)
  out$mc <- mc
  out$alpha <- alpha
  out$sign <- sign
  out$method <- if (identical(engine, "meta")) {
    if (!is.null(combine)) paste0("combine:", combine) else
      meta_method %||% getOption("fmrireg.meta.method", "pm")
  } else {
    engine
  }
  out$combine <- if (identical(engine, "meta")) combine else NULL
  out$weights <- if (identical(engine, "meta")) weights else NULL
  if (!is.null(group_info)) out$group_levels <- group_info$levels

  class(out) <- c("fmri_ttest_fit", "list")
  out
}
