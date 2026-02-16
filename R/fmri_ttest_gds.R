#' fmri_ttest method for gds-backed group_data
#'
#' Delegates to fmrigds reducers for meta or classical (OLS) engines and wraps
#' outputs into an fmri_ttest_fit-compatible object.
#'
#' @inheritParams fmri_ttest
#' @keywords internal
fmri_ttest.group_data_gds <- function(gd,
                           formula = ~ 1,
                           engine = c("auto", "meta", "classic"),
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
  if (!is.null(mc)) mc <- match.arg(mc, c("bh","by"))

  if (weights == "custom") {
    if (is.null(weights_custom)) {
      stop("weights='custom' requires 'weights_custom'", call. = FALSE)
    }
    stop("weights='custom' is not yet supported for gds-backed data; use weights='ivw' or 'equal'.", call. = FALSE)
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

  pl <- fmrigds::as_plan(gd)
  if (engine == "meta") {
    pl <- fmrigds::reduce(pl, method = "meta:fe", formula = formula, weights = weights_gds)
  } else {
    pl <- fmrigds::reduce(pl, method = "ols:voxelwise", formula = formula)
  }
  if (!is.null(mc)) {
    pl <- fmrigds::posthoc(pl, paste0("fdr:", mc), options = list(alpha = alpha))
  }
  res <- fmrigds::compute(pl)

  beta <- .gds_coerce_matrix(.gds_safe_assay(res, "beta"))
  se   <- .gds_coerce_matrix(.gds_safe_assay(res, "se"))
  tmat <- .gds_coerce_matrix(.gds_safe_assay(res, "t"))
  z    <- .gds_coerce_matrix(.gds_safe_assay(res, "z"))
  p    <- .gds_coerce_matrix(.gds_safe_assay(res, "p"))
  q    <- .gds_coerce_matrix(.gds_safe_assay(res, "q"))
  df_assay <- .gds_safe_assay(res, "df")

  # Names and labels
  coef_terms <- tryCatch(colnames(fmrigds::model_matrix(gd, formula)), error = function(e) NULL)
  sample_labels <- attr(gd, "fmrireg_sample_labels")
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

  # Compose df matrix robustly if present
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
    q    = if (!is.null(q))    t(q)    else NULL,
    df   = df_mat,
    formula = formula,
    engine  = engine,
    n_subjects = tryCatch(length(fmrigds::subjects(gd)), error = function(e) NA_integer_),
    n_features = tryCatch(ncol(beta %||% z %||% p), error = function(e) NA_integer_),
    roi_names = sample_labels
  )

  # Apply row/col names for readability and symmetry
  if (!is.null(out$beta) && is.matrix(out$beta) && length(coef_terms) == nrow(out$beta)) rownames(out$beta) <- coef_terms
  if (!is.null(out$se)   && is.matrix(out$se)   && length(coef_terms) == nrow(out$se))   rownames(out$se)   <- coef_terms
  if (!is.null(out$t)    && is.matrix(out$t)    && length(coef_terms) == nrow(out$t))    rownames(out$t)    <- coef_terms
  if (!is.null(out$z)    && is.matrix(out$z)    && length(coef_terms) == nrow(out$z))    rownames(out$z)    <- coef_terms
  if (!is.null(out$p)    && is.matrix(out$p)    && length(coef_terms) == nrow(out$p))    rownames(out$p)    <- coef_terms
  if (!is.null(out$q)    && is.matrix(out$q)    && length(coef_terms) == nrow(out$q))    rownames(out$q)    <- coef_terms
  if (!is.null(out$df)   && is.matrix(out$df)   && length(coef_terms) == nrow(out$df))   rownames(out$df)   <- coef_terms
  if (!is.null(sample_labels)) {
    if (!is.null(out$beta) && is.matrix(out$beta) && length(sample_labels) == ncol(out$beta)) colnames(out$beta) <- sample_labels
    if (!is.null(out$se)   && is.matrix(out$se)   && length(sample_labels) == ncol(out$se))   colnames(out$se)   <- sample_labels
    if (!is.null(out$t)    && is.matrix(out$t)    && length(sample_labels) == ncol(out$t))    colnames(out$t)    <- sample_labels
    if (!is.null(out$z)    && is.matrix(out$z)    && length(sample_labels) == ncol(out$z))    colnames(out$z)    <- sample_labels
    if (!is.null(out$p)    && is.matrix(out$p)    && length(sample_labels) == ncol(out$p))    colnames(out$p)    <- sample_labels
    if (!is.null(out$q)    && is.matrix(out$q)    && length(sample_labels) == ncol(out$q))    colnames(out$q)    <- sample_labels
    if (!is.null(out$df)   && is.matrix(out$df)   && length(sample_labels) == ncol(out$df))   colnames(out$df)   <- sample_labels
  }

  # Apply sign convention for group coefficient if present
  if (sign == "BminusA" && !is.null(out$beta) && any(rownames(out$beta) == "group")) {
    gi <- which(rownames(out$beta) == "group")
    if (length(gi) == 1L) {
      out$beta[gi, ] <- -out$beta[gi, ]
      if (!is.null(out$t)) out$t[gi, ] <- -out$t[gi, ]
      if (!is.null(out$z)) out$z[gi, ] <- -out$z[gi, ]
    }
  }

  class(out) <- c("fmri_ttest_fit", "list")
  out
}
