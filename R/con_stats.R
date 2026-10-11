#' @keywords internal
#' @noRd
qr.lm <- getFromNamespace("qr.lm", "stats")

#' Fit Contrasts
#'
#' @description
#' Generic function for fitting contrasts to model objects.
#'
#' @param object A fitted model object
#' @param ... Additional arguments passed to methods
#'
#' @return Contrast results (format depends on method)
#' @export
fit_contrasts <- function(object, ...) UseMethod("fit_contrasts")

#' Extract basic components from a linear model fit
#'
#' @param lmfit A fitted linear model object.
#' @return A list containing coefficient matrix (`betamat`), residual standard
#'   deviations (`sigma`), and residual degrees of freedom (`dfres`).
#' @keywords internal
#' @noRd
.lm_basic_stats <- function(lmfit) {
  betamat <- as.matrix(coef(lmfit))
  rss     <- colSums(as.matrix(lmfit$residuals^2))
  dfres   <- lmfit$df.residual
  sigma   <- sqrt(rss / dfres)
  list(betamat = betamat, sigma = sigma, dfres = dfres)
}

#' t statistic with a scale-free degenerate case
#'
#' `est / se`, except that an exactly zero standard error (a perfect fit) gives
#' 0. There is deliberately no absolute tolerance on `se`: a threshold such as
#' `sqrt(.Machine$double.eps)` silently zeroes valid statistics whenever the
#' data are on a small scale.
#' @keywords internal
#' @noRd
.safe_tstat <- function(est, se) {
  out <- est / se
  zero <- !is.na(se) & se == 0
  out[zero] <- 0
  out
}

#' Effective-df inputs for a possibly rank-deficient design
#'
#' The residual df is `n - rank`, so `n = df + rank` and the parameter count
#' that `calculate_effective_df()` must subtract is the rank, not `ncol(X)`.
#' @keywords internal
#' @noRd
.effective_df_from_rank <- function(df, p, aliased, robust_weights, ar_order) {
  if (is.null(robust_weights) && !(ar_order > 0)) {
    return(df)
  }
  rank <- p - length(aliased)
  calculate_effective_df(df + rank, rank, robust_weights, ar_order, method = "simple")
}

.lm_cov_unscaled <- function(lmfit, varnames = NULL, warn = TRUE) {
  qrx <- qr.lm(lmfit)
  coefs <- as.matrix(coef(lmfit))
  p <- nrow(coefs)
  rank <- qrx$rank
  pivot <- qrx$pivot %||% seq_len(p)
  estimable <- if (rank > 0L) pivot[seq_len(rank)] else integer(0)
  aliased <- if (rank < p) pivot[seq.int(rank + 1L, p)] else integer(0)

  cov.unscaled <- matrix(0, nrow = p, ncol = p)
  if (rank > 0L) {
    R <- qr.R(qrx)[seq_len(rank), seq_len(rank), drop = FALSE]
    cov_est <- chol2inv(R)
    cov.unscaled[estimable, estimable] <- cov_est
  }

  if (is.null(varnames)) {
    varnames <- rownames(coefs)
  }
  if (!is.null(varnames) && length(varnames) == p) {
    dimnames(cov.unscaled) <- list(varnames, varnames)
  }

  rank_info <- list(
    rank = rank,
    is_full_rank = rank == p,
    estimable = estimable,
    aliased = aliased,
    tol = NA_real_
  )
  cov.unscaled <- .attach_rank_attrs(cov.unscaled, rank_info)

  if (warn && length(aliased) > 0L) {
    warning(.rank_deficiency_message(rank, p, aliased, varnames), call. = FALSE)
  }

  cov.unscaled
}

.aliased_indices <- function(...) {
  vals <- list(...)
  for (val in vals) {
    aliased <- attr(val, "aliased", exact = TRUE)
    if (!is.null(aliased) && length(aliased) > 0L) {
      return(as.integer(aliased))
    }
  }
  integer(0)
}

.contrast_uses_aliased <- function(L, aliased, tol = sqrt(.Machine$double.eps)) {
  if (length(aliased) == 0L) {
    return(FALSE)
  }
  L <- as.matrix(L)
  aliased <- aliased[aliased <= ncol(L)]
  if (length(aliased) == 0L) {
    return(FALSE)
  }
  any(abs(L[, aliased, drop = FALSE]) > tol, na.rm = TRUE)
}

.warn_nonestimable_contrast <- function(contrast_name = NULL, aliased, coef_names = NULL) {
  aliased_names <- if (!is.null(coef_names) && length(coef_names) >= max(aliased, 0L)) {
    coef_names[aliased]
  } else {
    paste0("#", aliased)
  }
  if (length(aliased_names) > 8L) {
    aliased_names <- c(aliased_names[seq_len(8L)], "...")
  }
  prefix <- if (!is.null(contrast_name) && nzchar(contrast_name)) {
    paste0("Contrast '", contrast_name, "'")
  } else {
    "Contrast"
  }
  warning(
    paste0(
      prefix,
      " references non-estimable coefficients in a rank-deficient design: ",
      paste(aliased_names, collapse = ", ")
    ),
    call. = FALSE
  )
}

#' @keywords internal
#' @noRd
fit_Ftests <- function(object) {
  w <- object$weights
  ssr <- if (is.null(w)) {
    apply(object$residuals, 2, function(vals) sum(vals^2))
  } else {
    apply(object$residuals, 2, function(vals) sum((vals^2 *w)))
  }
  
  mss <- if (is.null(w)) {
    apply(object$fitted.values, 2, function(vals) sum(vals^2))
  } else {
    apply(object$fitted.values, 2, function(vals) sum(vals^2 * w))
  }
  
  #if (ssr < 1e-10 * mss) 
  #  warning("ANOVA F-tests on an essentially perfect fit are unreliable")
  
  dfr <- df.residual(object)
  p <- object$rank
  
  p1 <- 1L:p
  #comp <- object$effects[p1]
  asgn <- object$assign[qr.lm(object)$pivot][p1]
  
  
  #nmeffects <- c("(Intercept)", attr(object$terms, "term.labels"))
  #tlabels <- nmeffects[1 + unique(asgn)]
  
  df <- c(lengths(split(asgn, asgn)), dfr)
  
  I.p <- diag(nrow(coefficients(object)))
  nterms <- length(unique(asgn))
  hmat <- lapply(1:nterms, function(i) {
    subs <- which(asgn == i)
    hyp.matrix <- I.p[subs, , drop=FALSE]
    hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop = FALSE]
  })
  
  ret <- lapply(seq_along(ssr), function(i) {
    comp <- object$effects[,i]
    ss <- c(unlist(lapply(split(comp^2, asgn), sum)), ssr[i])
    ms <- ss/df
    f <- ms/(ssr[i]/dfr)
    
    P <- pf(f, df, dfr, lower.tail = FALSE)
    list(F = f, P = P)
  })
  
  FMat <- do.call(rbind, lapply(ret, "[[", "F"))
  PMat <- do.call(rbind, lapply(ret, "[[", "P"))
  
  list(F=FMat, P=PMat)
  
}



#' Beta Statistics for Linear Model
#'
#' @description
#' This function calculates beta statistics for a linear model fit.
#'
#' @param lmfit Fitted linear model object.
#' @param varnames Vector of variable names.
#' @param se Logical flag indicating whether to calculate standard errors (default: TRUE).
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{\code{estimate}: Estimated beta coefficients.}
#'   \item{\code{stat}: t-statistics of the beta coefficients (if \code{se = TRUE}).}
#'   \item{\code{se}: Standard errors of the beta coefficients (if \code{se = TRUE}).}
#'   \item{\code{prob}: Probabilities associated with the t-statistics (if \code{se = TRUE}).}
#'   \item{\code{stat_type}: Type of the statistics calculated.}
#' }
#' @examples
#' data(mtcars)
#' lm_fit <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' beta_stats(lm_fit, c("Intercept", "wt", "qsec", "am"))
#' @noRd
#' @autoglobal
beta_stats <- function(lmfit, varnames, se = TRUE) {
  basics  <- .lm_basic_stats(lmfit)
  betamat <- t(basics$betamat)
  sigma   <- basics$sigma
  rdf     <- basics$dfres
  colnames(betamat) <- varnames

  if (se) {
    cov.unscaled <- .lm_cov_unscaled(lmfit, varnames = varnames)
    aliased <- attr(cov.unscaled, "aliased", exact = TRUE) %||% integer(0)
    coef_se      <- sqrt(diag(cov.unscaled))
    if (length(aliased) > 0L) {
      coef_se[aliased] <- NA_real_
      betamat[, aliased] <- NA_real_
    }
    vc           <- outer(sigma, coef_se)
    colnames(vc) <- varnames
    tstat        <- .safe_tstat(betamat, vc)
    prob         <- 2 * pt(-abs(tstat), rdf)

    ret <- dplyr::tibble(
      type        = "beta",
      name        = "parameter_estimates",
      stat_type   = "tstat",
      df.residual = rdf,
      conmat      = list(NULL),
      colind      = list(NULL),
      data        = list(dplyr::tibble(
        estimate = list(betamat),
        se       = list(vc),
        stat     = list(tstat),
        prob     = list(prob),
        sigma    = list(sigma)
      ))
    )
  } else {
    ret <- dplyr::tibble(
      type        = "beta",
      name        = "parameter_estimates",
      stat_type   = "effect",
      df.residual = rdf,
      conmat      = list(NULL),
      colind      = list(NULL),
      data        = list(tibble(
        estimate = list(betamat),
        se       = list(NULL),
        stat     = list(NULL),
        prob     = list(NULL),
        sigma    = list(sigma)
      ))
    )
  }

  ret
}

#' Normalise an F-contrast weight matrix to hypotheses-in-rows orientation
#'
#' fmrireg accepts F-contrast weights against a term with `length(colind)`
#' coefficients in either orientation:
#'   * hypothesis-oriented: `n_hypotheses x length(colind)` (each row is a linear
#'     hypothesis), or
#'   * coefficient/cell-oriented: `length(colind) x r` (rows index cells,
#'     columns index hypotheses), as produced by `oneway_contrast()`.
#' This returns the weights as `n_hypotheses x length(colind)`, so callers derive both the
#' full contrast matrix and its numerator rank from `nrow(result)`.
#' A square matrix (both dims `== length(colind)`) is treated as already
#' hypothesis-oriented, preserving historical precedence.
#' @keywords internal
#' @noRd
.orient_fcontrast <- function(conmat, colind) {
  conmat <- as.matrix(conmat)
  k <- length(colind)
  if (ncol(conmat) == k) {
    conmat
  } else if (nrow(conmat) == k) {
    t(conmat)
  } else {
    stop(sprintf(
      "F contrast weight matrix dimensions %d x %d do not match length(colind) %d",
      nrow(conmat), ncol(conmat), k
    ), call. = FALSE)
  }
}

#' Fit F-contrasts for Linear Model
#'
#' @description
#' This function calculates F-contrasts for a fitted linear model.
#'
#' @param lmfit Fitted linear model object.
#' @param conmat Contrast matrix.
#' @param colind Column indices corresponding to the variables in the contrast matrix.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{\code{estimate}: Estimated contrasts.}
#'   \item{\code{se}: Residual variance.}
#'   \item{\code{stat}: F-statistics for the contrasts.}
#'   \item{\code{prob}: Probabilities associated with the F-statistics.}
#'   \item{\code{stat_type}: Type of the statistics calculated.}
#' }
#' @noRd
#' @keywords internal

fit_Fcontrasts <- function(lmfit, conmat, colind, cov.unscaled = NULL) {
  basics  <- .lm_basic_stats(lmfit)
  betamat <- basics$betamat
  sigma2  <- basics$sigma^2
  rdf     <- basics$dfres

  # Callers fitting several contrasts on one model can compute this once.
  if (is.null(cov.unscaled)) {
    cov.unscaled <- .lm_cov_unscaled(lmfit)
  }
  aliased <- attr(cov.unscaled, "aliased", exact = TRUE) %||% integer(0)

  # Orient to hypotheses-in-rows so the numerator rank `r` and the full
  # contrast matrix `cmat` are derived from a single source. Sizing `cmat` by
  # nrow(conmat) while leaving r = nrow(conmat) is wrong for cell-oriented
  # weights (e.g. oneway_contrast: k x (k-1)), where r must be (k-1).
  Cw <- .orient_fcontrast(conmat, colind)
  r  <- nrow(Cw)
  cmat <- matrix(0, r, nrow(betamat))
  cmat[, colind] <- Cw

  if (.contrast_uses_aliased(cmat, aliased)) {
    .warn_nonestimable_contrast(NULL, aliased, rownames(betamat))
    estimate <- rep(NA_real_, ncol(betamat))
    Fstat <- rep(NA_real_, ncol(betamat))
    return(structure(list(
      conmat      = conmat,
      estimate    = estimate,
      se          = sigma2,
      df.residual = rdf,
      stat        = Fstat,
      prob        = rep(NA_real_, ncol(betamat)),
      stat_type   = "Fstat"),
      class = c("Fstat", "result_stat")
    ))
  }

  M  <- cmat %*% cov.unscaled %*% t(cmat)
  cm <- tryCatch(solve(M), error = function(e) {
    warning(paste(
      "Singular matrix in F-contrast computation (C(X'X)^-1C'). Details:",
      e$message
    ))
    matrix(NaN, nrow(M), ncol(M))
  })

  # Exclude zero-weight coefficients so NA-valued aliased betas (ruled out of
  # this contrast by .contrast_uses_aliased() above) do not poison the
  # estimate via 0 * NA, mirroring the t-contrast path.
  nz <- which(colSums(abs(cmat)) > 0)
  cmat_nz <- cmat[, nz, drop = FALSE]
  beta_nz <- betamat[nz, , drop = FALSE]
  # Per-voxel quadratic forms u_v' cm u_v for all voxels at once (U is r x V).
  U <- cmat_nz %*% beta_nz
  estimate <- colSums(U * (cm %*% U)) / r

  Fstat <- estimate / sigma2

  structure(list(
    conmat      = conmat,
    estimate    = estimate,
    se          = sigma2,
    df.residual = rdf,
    stat        = Fstat,
    prob        = pf(Fstat, r, rdf, lower.tail = FALSE),
    stat_type   = "Fstat"),
    class = c("Fstat", "result_stat")
  )
}

#' @keywords internal
#' @noRd
#' @autoglobal
estimate_contrast.contrast <- function(x, fit, colind, ...) {
  ret <- fit_contrasts(fit, x$weights, colind, se = TRUE, cov.unscaled = list(...)$cov.unscaled)
  
  ## assumes se computed
  tibble(
    type="contrast",
    name=x$name,
    stat_type=ret$stat_type,
    df.residual=ret$df.residual,
    conmat=list(x$weights),
    colind=list(colind),
    data=list(tibble(
      estimate=ret$estimate,
      se=ret$se,
      stat=ret$stat,
      prob=ret$prob,
      sigma=ret$sigma)))
    
}


#' @export
#' @autoglobal
#' @keywords internal
#' @noRd
estimate_contrast.Fcontrast <- function(x, fit, colind, ...) {
  ret <- fit_Fcontrasts(fit, x$weights, colind, cov.unscaled = list(...)$cov.unscaled)
  tibble(
    type="Fcontrast",
    name=x$name,
    stat_type=ret$stat_type,
    df.residual=ret$df.residual,
    conmat=list(x$weights),
    colind=list(colind),
    data=list(tibble(
      estimate=ret$estimate,
      se=ret$se,
      stat=ret$stat,
      prob=ret$prob)))
}

#' Fit Contrasts for Linear Model (Default Method)
#'
#' @description
#' This function calculates contrasts for a fitted linear model.
#'
#' @param object The fitted linear model object.
#' @param conmat The contrast matrix or contrast vector.
#' @param colind The subset column indices in the design associated with the contrast.
#' @param se Whether to compute standard errors, t-statistics, and p-values (default: TRUE).
#' @param cov.unscaled Optional unscaled coefficient covariance \eqn{(X'X)^{-1}}
#'   of \code{object}, to reuse across several contrasts on the same fit
#'   (default: computed from \code{object}).
#' @param ... Additional arguments (unused)
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{\code{conmat}: Contrast matrix.}
#'   \item{\code{sigma}: Residual standard error.}
#'   \item{\code{df.residual}: Degrees of freedom for residuals.}
#'   \item{\code{estimate}: Estimated contrasts.}
#'   \item{\code{se}: Standard errors of the contrasts (if \code{se = TRUE}).}
#'   \item{\code{stat}: t-statistics for the contrasts (if \code{se = TRUE}).}
#'   \item{\code{prob}: Probabilities associated with the t-statistics (if \code{se = TRUE}).}
#'   \item{\code{stat_type}: Type of the statistics calculated.}
#' }
#' @method fit_contrasts default
#' @export
fit_contrasts.default <- function(object, conmat, colind, se = TRUE, cov.unscaled = NULL, ...) {

  conmat <- as.matrix(conmat)

  basics  <- .lm_basic_stats(object)
  betamat <- basics$betamat
  sigma   <- basics$sigma
  rdf     <- basics$dfres

  ncoef <- nrow(betamat)
  cmat  <- matrix(0, ncoef, 1)
  cmat[colind, ] <- conmat

  # Aliased (non-estimable) coefficients are NA in betamat. For an estimable
  # contrast they carry zero weight, but 0 * NA = NA in R would poison the
  # estimate, so restrict the projection to the non-zero-weight coefficients.
  # Contrasts that genuinely load on aliased coefficients are flagged NA below
  # via .contrast_uses_aliased().
  nz <- which(cmat[, 1] != 0)
  ct <- drop(t(cmat[nz, , drop = FALSE]) %*% betamat[nz, , drop = FALSE])

  if (se) {
    if (is.null(cov.unscaled)) cov.unscaled <- .lm_cov_unscaled(object)
    aliased <- attr(cov.unscaled, "aliased", exact = TRUE) %||% integer(0)

    if (.contrast_uses_aliased(t(cmat), aliased)) {
      .warn_nonestimable_contrast(NULL, aliased, rownames(betamat))
      ct <- rep(NA_real_, length(sigma))
      se_vec <- rep(NA_real_, length(sigma))
      stat <- rep(NA_real_, length(sigma))
      prob <- rep(NA_real_, length(sigma))
    } else {
      var_est      <- as.numeric(t(cmat) %*% cov.unscaled %*% cmat)
      se_vec       <- sqrt(var_est) * sigma
      stat         <- ct / se_vec
      prob         <- 2 * pt(-abs(stat), rdf)
    }

    structure(list(
      conmat      = cmat,
      sigma       = sigma,
      df.residual = rdf,
      estimate    = ct,
      se          = se_vec,
      stat        = stat,
      prob        = prob,
      stat_type   = "tstat"),
      class = c("tstat", "result_stat")
    )
  } else {
    structure(list(
      conmat      = cmat,
      sigma       = sigma,
      df.residual = rdf,
      estimate    = ct,
      stat_type   = "effects"),
      class = c("effect", "result_stat")
    )
  }
}

#' @keywords internal
#' @noRd
#' @autoglobal
.fast_t_contrast  <- function(B, sigma2, XtXinv, l, df, 
                               robust_weights = NULL, ar_order = 0,
                               aliased = integer(0), contrast_name = NULL) {
  # B:      p x V matrix (betas)
  # sigma2: V-vector (residual variance)
  # XtXinv: p x p matrix (inverse crossproduct of design)
  # l:      1 x p vector/matrix (contrast weights)
  # df:     scalar (residual degrees of freedom)
  # robust_weights: Optional vector of robust weights
  # ar_order: Order of AR model (0 if none)
  
  # Ensure l is a row vector (matrix)
  if (!is.matrix(l)) {
    l <- matrix(l, nrow = 1)
  }

  if (.contrast_uses_aliased(l, aliased)) {
    .warn_nonestimable_contrast(contrast_name, aliased, rownames(B) %||% colnames(XtXinv))
    nvox <- ncol(B)
    return(list(
      estimate = rep(NA_real_, nvox),
      se       = rep(NA_real_, nvox),
      stat     = rep(NA_real_, nvox),
      prob     = rep(NA_real_, nvox),
      sigma    = sqrt(pmax(0, sigma2)),
      stat_type = "tstat"
    ))
  }
  
  est  <- drop(l %*% B)                # 1 × V  (BLAS GEMV)
  s2   <- as.numeric(l %*% XtXinv %*% t(l)) # scalar Var(est) / sigma2
  
  # Avoid division by zero or NaNs if s2 or sigma2 are zero/negative
  se   <- sqrt(pmax(0, s2 * sigma2))            # 1 × V
  
  df_effective <- .effective_df_from_rank(df, nrow(B), aliased, robust_weights, ar_order)

  tval <- .safe_tstat(est, se)
  pval <- 2 * pt(-abs(tval), df_effective)

  list(estimate = est,
       se       = se,
       stat     = tval,
       prob     = pval,
       sigma    = sqrt(pmax(0, sigma2)), # Include sigma for potential compatibility
       stat_type = "tstat")
}

#' @keywords internal
#' @noRd
#' @autoglobal
.fast_F_contrast <- function(B, sigma2, XtXinv, L, df,
                              robust_weights = NULL, ar_order = 0,
                              aliased = integer(0), contrast_name = NULL) {
  # B:      p x V matrix (betas)
  # sigma2: V-vector (residual variance)
  # XtXinv: p x p matrix
  # L:      r x p matrix (contrast weights)
  # df:     scalar (residual degrees of freedom)
  # robust_weights: Optional vector of robust weights
  # ar_order: Order of AR model (0 if none)
  
  if (!is.matrix(L)) {
      stop(".fast_F_contrast requires L to be a matrix.")
  }

  if (.contrast_uses_aliased(L, aliased)) {
    .warn_nonestimable_contrast(contrast_name, aliased, rownames(B) %||% colnames(XtXinv))
    nvox <- ncol(B)
    return(list(
      estimate = rep(NA_real_, nvox),
      se       = sigma2,
      stat     = rep(NA_real_, nvox),
      prob     = rep(NA_real_, nvox),
      stat_type = "Fstat"
    ))
  }
  
  r   <- nrow(L)
  U   <- L %*% B                    # r × V   (BLAS GEMM) : L*beta
  M   <- L %*% XtXinv %*% t(L)      # r × r   : L (X'X)^-1 L'
  
  # Use tryCatch for solve() in case M is singular
  Cinv <- tryCatch(solve(M), error = function(e) {
    warning(paste("Singular matrix in F-contrast computation (L(X'X)^-1 L'). Details:", e$message))
    # Return a matrix of NaNs or zeros? Let's use NaNs.
    matrix(NaN, nrow = r, ncol = r)
  })

  # Per-voxel quadratic form qf_v = u_v' Cinv u_v = diag(t(U) %*% Cinv %*% U).
  # With U = r x V and tmp = t(U) %*% Cinv = V x r, the v-th quadratic form is
  # the row sum of (tmp * t(U)), so this must be rowSums (length V), not
  # colSums (which would sum over voxels and return a length-r vector).
  if (any(is.nan(Cinv))) {
      qf <- rep(NaN, ncol(U))
  } else {
      tmp  <- t(U) %*% Cinv           # V x r
      qf  <- rowSums(tmp * t(U))      # V-vector, each element is u_v' Cinv u_v
  }
  
  estimate <- qf / r # Numerator mean square: (LB)' (L(X'X)^-1 L')^-1 (LB) / r
  se <- sigma2 # Denominator mean square (residual variance)
  
  df_effective <- .effective_df_from_rank(df, nrow(B), aliased, robust_weights, ar_order)

  # A zero residual variance (perfect fit) has no F; no absolute tolerance,
  # which would discard small-scale data.
  Fval <- ifelse(!(se > 0) | is.nan(qf), NaN, estimate / se)
                 
  pval <- pf(Fval, r, df_effective, lower.tail = FALSE)

  list(estimate = estimate, # Numerator MS
       se       = se,       # Denominator MS (sigma2)
       stat     = Fval,
       prob     = pval,
       stat_type = "Fstat")
}

#' @keywords internal
#' @noRd
#' @autoglobal
fit_lm_contrasts_fast <- function(B, sigma2, XtXinv, conlist, fconlist, df,
                                  robust_weights = NULL, ar_order = 0,
                                  aliased = NULL) {
  # Validate inputs
  if (!is.matrix(B)) stop("B must be a matrix (p x V)")
  if (!is.matrix(XtXinv)) stop("XtXinv must be a matrix (p x p)")
  if (nrow(B) != nrow(XtXinv) || nrow(XtXinv) != ncol(XtXinv)) {
    stop("Dimension mismatch: B must be p x V, XtXinv must be p x p")
  }
  aliased <- aliased %||% .aliased_indices(B, XtXinv)
  
  # Process t-contrasts
  tcontrast_results <- process_t_contrasts(
    B = B, sigma2 = sigma2, XtXinv = XtXinv, 
    conlist = conlist, df = df,
    robust_weights = robust_weights, ar_order = ar_order,
    aliased = aliased
  )
  
  # Process F-contrasts  
  fcontrast_results <- process_f_contrasts(
    B = B, sigma2 = sigma2, XtXinv = XtXinv,
    fconlist = fconlist, df = df,
    robust_weights = robust_weights, ar_order = ar_order,
    aliased = aliased
  )
  
  # Combine results
  all_results <- c(tcontrast_results, fcontrast_results)
  
  # Ensure proper naming
  names(all_results) <- vapply(all_results, function(x) x$name[1], character(1))
  
  return(all_results)
}

#' Process t-contrasts for fast linear model fitting
#' @keywords internal
#' @noRd
#' @autoglobal
process_t_contrasts <- function(B, sigma2, XtXinv, conlist, df, 
                                robust_weights = NULL, ar_order = 0,
                                aliased = integer(0)) {
  if (length(conlist) == 0) return(list())
  
  results <- vector("list", length(conlist))
  names(results) <- names(conlist)
  
  for (i in seq_along(conlist)) {
    con_name <- names(conlist)[i]
    contrast_weights <- conlist[[i]]
    
    tryCatch({
      # Extract colind - this is required for proper contrast computation
      colind <- extract_colind(contrast_weights, con_name, "t-contrast")
      
      # Create full contrast vector
      full_contrast <- create_full_contrast_vector(contrast_weights, colind, nrow(XtXinv))
      
      # Compute contrast statistics
      stats <- .fast_t_contrast(B, sigma2, XtXinv, full_contrast, df, 
                               robust_weights, ar_order, aliased, con_name)
      
      # Package results in expected format
      results[[i]] <- package_tcontrast_result(
        con_name, contrast_weights, colind, df, stats
      )
      
    }, error = function(e) {
      stop(sprintf(
        "Failed to compute t-contrast '%s': %s",
        con_name, e$message
      ), call. = FALSE)
    })
  }
  
  results
}

#' Process F-contrasts for fast linear model fitting  
#' @keywords internal
#' @noRd
#' @autoglobal
process_f_contrasts <- function(B, sigma2, XtXinv, fconlist, df,
                                robust_weights = NULL, ar_order = 0,
                                aliased = integer(0)) {
  if (length(fconlist) == 0) return(list())
  
  results <- vector("list", length(fconlist))
  names(results) <- names(fconlist)
  
  for (i in seq_along(fconlist)) {
    con_name <- names(fconlist)[i]
    contrast_matrix <- fconlist[[i]]
    
    tryCatch({
      # Extract colind - this is required for proper contrast computation
      colind <- extract_colind(contrast_matrix, con_name, "F-contrast")
      
      # Create full contrast matrix
      full_contrast <- create_full_contrast_matrix(contrast_matrix, colind, nrow(XtXinv))
      
      # Compute contrast statistics
      stats <- .fast_F_contrast(B, sigma2, XtXinv, full_contrast, df,
                               robust_weights, ar_order, aliased, con_name)
      
      # Package results in expected format
      results[[i]] <- package_fcontrast_result(
        con_name, contrast_matrix, colind, df, stats
      )
      
    }, error = function(e) {
      stop(sprintf(
        "Failed to compute F-contrast '%s': %s",
        con_name, e$message
      ), call. = FALSE)
    })
  }
  
  results
}

#' Extract column indices from contrast object
#' @keywords internal
#' @noRd
extract_colind <- function(contrast_obj, con_name, contrast_type) {
  colind <- attr(contrast_obj, "colind")
  
  if (is.null(colind)) {
    stop(sprintf(
      "Missing 'colind' attribute for %s '%s'. Contrast objects must have column indices attached.",
      contrast_type, con_name
    ))
  }
  
  if (!is.numeric(colind) || any(colind < 1)) {
    stop(sprintf(
      "Invalid 'colind' attribute for %s '%s': must be positive integers",
      contrast_type, con_name  
    ))
  }
  
  return(as.integer(colind))
}

#' Create full contrast vector with proper dimensions
#' @keywords internal  
#' @noRd
create_full_contrast_vector <- function(contrast_weights, colind, p) {
  # Ensure contrast_weights is a vector
  if (is.matrix(contrast_weights)) {
    if (nrow(contrast_weights) == 1) {
      contrast_weights <- as.vector(contrast_weights)
    } else if (ncol(contrast_weights) == 1) {
      contrast_weights <- as.vector(contrast_weights) 
    } else {
      stop("t-contrast weights must be a vector or single-row/column matrix")
    }
  }
  
  # Validate dimensions
  if (length(contrast_weights) != length(colind)) {
    stop(sprintf(
      "Dimension mismatch: contrast weights length (%d) != colind length (%d)",
      length(contrast_weights), length(colind)
    ))
  }
  
  if (max(colind) > p) {
    stop(sprintf(
      "Column index out of bounds: max colind (%d) > design matrix columns (%d)",
      max(colind), p
    ))
  }
  
  # Create full contrast vector
  full_contrast <- matrix(0, nrow = 1, ncol = p)
  full_contrast[1, colind] <- contrast_weights
  
  return(full_contrast)
}

#' Create full contrast matrix with proper dimensions
#' @keywords internal
#' @noRd  
create_full_contrast_matrix <- function(contrast_matrix, colind, p) {
  if (!is.matrix(contrast_matrix)) {
    stop("F-contrast weights must be a matrix")
  }

  # Accept either orientation (e.g. oneway_contrast supplies cell-oriented
  # k x (k-1) weights); normalise to hypotheses-in-rows so the fast F engine
  # sees an r x p matrix and derives r = nrow(L) correctly.
  Cw <- .orient_fcontrast(contrast_matrix, colind)

  if (max(colind) > p) {
    stop(sprintf(
      "Column index out of bounds: max colind (%d) > design matrix columns (%d)",
      max(colind), p
    ))
  }

  # Create full contrast matrix
  full_contrast <- matrix(0, nrow = nrow(Cw), ncol = p)
  full_contrast[, colind] <- Cw

  return(full_contrast)
}

#' Package t-contrast results in expected format
#' @keywords internal
#' @noRd
package_tcontrast_result <- function(con_name, original_weights, colind, df, stats) {
  dplyr::tibble(
    type = "contrast",
    name = con_name,
    stat_type = stats$stat_type,
    df.residual = df,
    conmat = list(original_weights),
    colind = list(colind),
    data = list(dplyr::tibble(
      estimate = stats$estimate,
      se = stats$se,
      stat = stats$stat,
      prob = stats$prob,
      sigma = stats$sigma
    ))
  )
}

#' Package F-contrast results in expected format  
#' @keywords internal
#' @noRd
package_fcontrast_result <- function(con_name, original_weights, colind, df, stats) {
  dplyr::tibble(
    type = "Fcontrast", 
    name = con_name,
    stat_type = stats$stat_type,
    df.residual = df,
    conmat = list(original_weights),
    colind = list(colind),
    data = list(dplyr::tibble(
      estimate = stats$estimate,
      se = stats$se,
      stat = stats$stat,
      prob = stats$prob
    ))
  )
}

#' @keywords internal
#' @noRd
#' @autoglobal
beta_stats_matrix <- function(Betas, XtXinv, sigma, dfres, varnames, 
                              robust_weights = NULL, ar_order = 0) {
  # Betas:    p x V matrix
  # XtXinv:   p x p matrix
  # sigma:    V-vector (residual std dev)
  # dfres:    Scalar residual df
  # varnames: p-vector of coefficient names
  # robust_weights: Optional vector of robust weights
  # ar_order: Order of AR model (0 if none)
  
  p <- nrow(Betas)
  V <- ncol(Betas)
  sigma2 <- sigma^2 # V-vector
  aliased <- .aliased_indices(Betas, XtXinv)
  
  df_effective <- .effective_df_from_rank(dfres, p, aliased, robust_weights, ar_order)
  
  # Compute SE for each beta, each voxel
  # se(beta_i) = sqrt( [XtXinv]_ii * sigma2 )
  diag_XtXinv <- diag(XtXinv)
  if (any(diag_XtXinv < 0)) {
      # Handle potential numerical issues if diagonal is negative
      warning("Negative diagonal elements found in XtXinv during beta SE calculation.")
      diag_XtXinv[diag_XtXinv < 0] <- NaN
  }
  if (length(aliased) > 0L) {
    diag_XtXinv[aliased] <- NA_real_
  }
  se_scaling <- sqrt(diag_XtXinv) # p-vector
  
  # Outer product: V-vector sigma with p-vector se_scaling -> V x p matrix
  vc <- outer(sigma, se_scaling) # SEs for all betas (rows=voxels, cols=betas)
  
  betamat <- t(Betas) # V x p matrix (voxels x betas)
  colnames(betamat) <- varnames
  colnames(vc) <- varnames
  if (length(aliased) > 0L) {
    betamat[, aliased] <- NA_real_
    vc[, aliased] <- NA_real_
  }
  
  # Compute t-stats: element-wise division
  # Avoid division by zero
  tstat <- .safe_tstat(betamat, vc)
  
  # Compute p-values using effective df
  prob <- 2 * pt(-abs(tstat), df_effective)
  
  # Package into the same tibble structure as beta_stats
  ret <- dplyr::tibble(
    type = "beta",
    name = "parameter_estimates",
    stat_type = "tstat",
    df.residual = dfres,
    conmat = list(NULL), # Keep consistent structure
    colind = list(NULL),
    data = list(dplyr::tibble(
      estimate = list(betamat), # V x p matrix
      se = list(vc),       # V x p matrix
      stat = list(tstat),    # V x p matrix
      prob = list(prob),     # V x p matrix
      sigma = list(sigma)    # V-vector
    ))
  )
  
  return(ret)
}

#' @keywords internal
#' @noRd
.enforce_contrast_scope <- function(colind, scope, contrast_name) {
  if (is.null(scope) || is.null(scope$allowed_colind) || is.null(colind)) {
    return(invisible(TRUE))
  }

  allowed <- as.integer(scope$allowed_colind)
  used <- as.integer(colind)
  bad <- setdiff(used, allowed)
  if (length(bad) == 0L) {
    return(invisible(TRUE))
  }

  reason <- as.character(scope$reason %||% "")
  msg <- paste0(
    "Contrast '", contrast_name, "' references unsupported coefficient indices for this fit: ",
    paste(bad, collapse = ", ")
  )
  if (nzchar(reason)) {
    msg <- paste0(msg, ". ", reason)
  }

  stop(msg, call. = FALSE)
}

#' Fit Contrasts for fMRI Linear Model Objects
#'
#' @description
#' S3 method for computing contrasts on fitted fmri_lm objects.
#'
#' @param object An fmri_lm object
#' @param contrasts A list of contrast specifications
#' @param ... Additional arguments (unused)
#'
#' @return A list of contrast results
#' @export
#' @method fit_contrasts fmri_lm
fit_contrasts.fmri_lm <- function(object, contrasts, ...) {
  if (is.null(contrasts) || length(contrasts) == 0) {
    return(list())
  }

  sigma <- object$result$sigma
  if (is.null(sigma)) {
    sigma <- object$result$betas$data[[1]]$sigma[[1]]
  }
  sigma <- as.numeric(sigma)

  beta_mat_raw <- tryCatch(
    suppressWarnings(as.matrix(coef(object, type = "betas"))),
    error = function(e) object$result$betas$data[[1]]$estimate[[1]]
  )

  event_idx_hint <- as.integer(object$result$event_indices %||% integer(0))
  expected_p <- length(event_idx_hint)

  if (expected_p > 0L) {
    if (ncol(beta_mat_raw) == expected_p) {
      beta_mat_vxp <- beta_mat_raw
    } else if (nrow(beta_mat_raw) == expected_p) {
      beta_mat_vxp <- t(beta_mat_raw)
    } else if (length(sigma) > 1L && ncol(beta_mat_raw) == length(sigma)) {
      beta_mat_vxp <- t(beta_mat_raw)
    } else {
      stop("Unable to align beta matrix orientation with event coefficient count")
    }
  } else if (length(sigma) > 1L) {
    if (nrow(beta_mat_raw) == length(sigma)) {
      beta_mat_vxp <- beta_mat_raw
    } else if (ncol(beta_mat_raw) == length(sigma)) {
      beta_mat_vxp <- t(beta_mat_raw)
    } else {
      stop("Unable to align beta matrix orientation with voxel count from sigma")
    }
  } else {
    beta_mat_vxp <- beta_mat_raw
  }

  betas <- t(as.matrix(beta_mat_vxp))  # p x V
  p <- nrow(betas)
  V <- ncol(betas)

  beta_names <- colnames(beta_mat_vxp)
  if (is.null(beta_names)) {
    dm_names <- tryCatch(colnames(design_matrix(object$model)), error = function(e) NULL)
    event_idx <- as.integer(object$result$event_indices %||% integer(0))
    if (!is.null(dm_names) && length(event_idx) == ncol(beta_mat_vxp) && length(event_idx) > 0L &&
        max(event_idx) <= length(dm_names)) {
      beta_names <- dm_names[event_idx]
    } else {
      beta_names <- paste0("beta_", seq_len(ncol(beta_mat_vxp)))
    }
  }

  if (length(sigma) == 1L) {
    sigma <- rep(sigma, V)
  }
  if (length(sigma) != V) {
    stop("Length of sigma does not match number of voxels in fmri_lm result")
  }
  sigma2 <- sigma^2

  df_residual <- object$result$rdf
  if (is.null(df_residual)) {
    df_residual <- object$result$betas$df.residual[1]
  }
  df_residual <- as.numeric(df_residual)[1]
  contrast_scope <- object$result$contrast_scope %||% NULL

  XtXinv <- object$result$cov.unscaled
  if (is.null(XtXinv)) {
    XtXinv <- object$result$XtXinv
  }
  if (is.null(XtXinv)) {
    fmrimod <- object$model
    tmats <- term_matrices(fmrimod)
    data_env <- list2env(tmats)
    data_env[[".y"]] <- rep(0, nrow(tmats[[1]]))
    form <- get_formula(fmrimod)
    X <- model.matrix(form, data_env)
    XtXinv <- .fast_preproject(X)$XtXinv
  }
  xtx_aliased <- attr(XtXinv, "aliased", exact = TRUE) %||% integer(0)
  XtXinv <- as.matrix(XtXinv)
  if (nrow(XtXinv) != p || ncol(XtXinv) != p) {
    event_idx <- as.integer(object$result$event_indices %||% integer(0))
    if (length(event_idx) == p && length(event_idx) > 0L &&
        max(event_idx) <= nrow(XtXinv) && max(event_idx) <= ncol(XtXinv)) {
      XtXinv <- XtXinv[event_idx, event_idx, drop = FALSE]
      xtx_aliased <- stats::na.omit(match(xtx_aliased, event_idx))
    } else {
      stop(
        sprintf(
          "Incompatible covariance dimensions for post-hoc contrasts: XtXinv is %dx%d but beta matrix has %d coefficients",
          nrow(XtXinv), ncol(XtXinv), p
        ),
        call. = FALSE
      )
    }
  }
  if (!is.null(beta_names) && length(beta_names) == p) {
    rownames(betas) <- beta_names
    dimnames(XtXinv) <- list(beta_names, beta_names)
  }
  attr(betas, "aliased") <- as.integer(xtx_aliased)
  attr(XtXinv, "aliased") <- as.integer(xtx_aliased)

  vm <- object$result$variance_model
  structured_covariance <- NULL
  df_inference <- rep(df_residual, V)
  if (inherits(vm, "fmri_lm_variance_model")) {
    if (length(vm$df_inference) == 1L) {
      df_inference <- rep(vm$df_inference, V)
    } else if (length(vm$df_inference) == V) {
      df_inference <- as.numeric(vm$df_inference)
    }
    if (identical(vm$covariance_scope, "voxel") &&
        is.list(vm$covariance) && length(vm$covariance) == V) {
      event_idx <- as.integer(object$result$event_indices %||% integer(0))
      structured_covariance <- lapply(vm$covariance, function(cv) {
        cv <- as.matrix(cv)
        if (nrow(cv) == p && ncol(cv) == p) return(cv)
        if (length(event_idx) == p && max(event_idx) <= nrow(cv)) {
          return(cv[event_idx, event_idx, drop = FALSE])
        }
        stop("Structured covariance dimensions do not match post-hoc coefficients.",
             call. = FALSE)
      })
    }
  }

  compute_t <- function(l, name) {
    if (is.null(structured_covariance)) {
      return(.fast_t_contrast(
        betas, sigma2, XtXinv, l, df_residual,
        aliased = xtx_aliased, contrast_name = name
      ))
    }
    estimate <- drop(l %*% betas)
    se <- vapply(structured_covariance, function(cv) {
      sqrt(pmax(drop(crossprod(l, cv %*% l)), 0))
    }, numeric(1))
    stat <- ifelse(se > 0, estimate / se, NA_real_)
    list(estimate = estimate, se = se, stat = stat,
         prob = 2 * stats::pt(-abs(stat), df = df_inference),
         stat_type = "tstat")
  }

  compute_F <- function(L, name) {
    if (is.null(structured_covariance)) {
      return(.fast_F_contrast(
        betas, sigma2, XtXinv, L, df_residual,
        aliased = xtx_aliased, contrast_name = name
      ))
    }
    r <- qr(L)$rank
    stat <- vapply(seq_len(V), function(v) {
      u <- L %*% betas[, v]
      M <- L %*% structured_covariance[[v]] %*% t(L)
      tryCatch(drop(crossprod(u, solve(M, u))) / r,
               error = function(e) NA_real_)
    }, numeric(1))
    list(estimate = stat, se = rep(1, V), stat = stat,
         prob = stats::pf(stat, r, df_inference, lower.tail = FALSE),
         stat_type = "Fstat")
  }

  package <- function(con_name, type, stats) {
    list(
      name = con_name,
      type = type,
      estimate = as.numeric(stats$estimate),
      se = as.numeric(stats$se),
      stat = as.numeric(stats$stat),
      prob = as.numeric(stats$prob),
      df.residual = df_residual,
      df.inference = as.numeric(df_inference),
      stat_type = stats$stat_type
    )
  }

  contrast_results <- lapply(seq_along(contrasts), function(i) {
    con_spec <- contrasts[[i]]
    con_name <- names(contrasts)[i]
    if (is.null(con_name) || !nzchar(con_name)) {
      con_name <- if (is.list(con_spec) && is.character(con_spec$name) &&
                      length(con_spec$name) == 1L && nzchar(con_spec$name)) {
        con_spec$name
      } else {
        paste0("contrast_", i)
      }
    }

    resolved <- tryCatch(
      .posthoc_contrast_matrix(con_spec, object, beta_names, p,
                               scope_check = function(ci) {
                                 .enforce_contrast_scope(ci, contrast_scope, con_name)
                               }),
      error = function(e) {
        stop(sprintf("Cannot compute post-hoc contrast '%s': %s",
                     con_name, conditionMessage(e)), call. = FALSE)
      }
    )
    .enforce_contrast_scope(resolved$colind, contrast_scope, con_name)

    if (identical(resolved$type, "F")) {
      package(con_name, "Fcontrast", compute_F(resolved$L, con_name))
    } else {
      package(con_name, "contrast", compute_t(drop(resolved$L), con_name))
    }
  })

  contrast_results <- contrast_results[!vapply(contrast_results, is.null, logical(1))]
  names(contrast_results) <- vapply(contrast_results, `[[`, character(1), "name")
  contrast_results
}

#' Map names to post-hoc coefficient positions
#'
#' Exact match against the event design column names first; otherwise a name
#' matches a column whose name ends in `_<name>` (the term-tag prefix that
#' fmridesign adds), provided that match is unique. Anything else is an error.
#' @keywords internal
#' @noRd
.posthoc_match_names <- function(nms, beta_names) {
  if (is.null(beta_names)) {
    stop("the fit has no coefficient names to match contrast weights against",
         call. = FALSE)
  }
  idx <- match(nms, beta_names)
  for (j in which(is.na(idx))) {
    hit <- which(endsWith(beta_names, paste0("_", nms[j])))
    if (length(hit) == 1L) idx[j] <- hit
  }
  if (anyNA(idx)) {
    stop(sprintf(
      "weight names do not match any event coefficient: %s (available: %s)",
      paste(nms[is.na(idx)], collapse = ", "),
      paste(beta_names, collapse = ", ")
    ), call. = FALSE)
  }
  if (anyDuplicated(idx)) {
    stop("weight names map to the same coefficient more than once", call. = FALSE)
  }
  idx
}

#' Resolve a contrast_spec against the fitted event terms
#'
#' Tries every (non-feature) event term with `contrast_weights()`; the spec
#' must resolve against exactly one of them, with every weight row naming one
#' of that term's design columns (the same mapping `contrast_weights()` uses
#' for an event model during the fit).
#' @return list(W = cells x r weights with design-column rownames, term, is_f)
#' @keywords internal
#' @noRd
.posthoc_resolve_spec <- function(con_spec, object) {
  emod <- object$model$event_model
  event_terms <- terms(emod)
  dm <- design_matrix(emod)
  col_indices <- attr(dm, "col_indices")
  ev_names <- colnames(dm)
  tnames <- names(event_terms)
  if (is.null(col_indices) || is.null(tnames)) {
    stop("the event design has no column-index metadata", call. = FALSE)
  }

  # Variables the specification refers to (e.g. `cond` in ~ cond == "A").
  spec_vars <- unique(unlist(lapply(con_spec, function(el) {
    if (inherits(el, "formula")) all.vars(el) else NULL
  })))

  hits <- list()
  for (k in seq_along(event_terms)) {
    term <- event_terms[[k]]
    if (inherits(term, "feature_term")) next
    term_vars <- unique(c(names(term$event_table), names(term$events), term$varname))
    if (length(spec_vars) > 0L && !all(spec_vars %in% term_vars)) next
    warns <- character(0)
    cw <- tryCatch(
      withCallingHandlers(
        contrast_weights(con_spec, term),
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) NULL
    )
    if (is.null(cw) || is.null(cw$weights)) next
    W <- as.matrix(cw$weights)
    rn <- rownames(W)
    term_cols <- ev_names[col_indices[[tnames[k]]]]
    if (is.null(rn) || length(term_cols) == 0L) next
    m <- match(rn, term_cols)
    if (anyNA(m)) {
      prefix <- paste0(tnames[k], "_")
      if (all(startsWith(term_cols, prefix))) {
        m <- match(rn, substring(term_cols, nchar(prefix) + 1L))
      }
    }
    if (anyNA(m)) next
    rownames(W) <- term_cols[m]
    hits[[length(hits) + 1L]] <- list(W = W, term = tnames[k],
                                      is_f = inherits(cw, "Fcontrast"),
                                      warnings = warns)
  }

  if (length(hits) == 0L) {
    stop("the contrast specification does not resolve against any event term ",
         "of the fitted model", call. = FALSE)
  }
  if (length(hits) > 1L) {
    stop(sprintf("the contrast specification is ambiguous: it resolves against terms %s",
                 paste(vapply(hits, `[[`, character(1), "term"), collapse = ", ")),
         call. = FALSE)
  }
  # Warnings from the term the spec resolved against are the user's to see
  # (as they would be for the same contrast declared in the model formula).
  for (msg in hits[[1L]]$warnings) warning(msg, call. = FALSE)
  hits[[1L]]
}

#' Build the r x p post-hoc contrast matrix for one specification
#'
#' Accepts a `contrast_spec` (any type that `contrast_weights()` understands),
#' a computed contrast object (with `weights`/`offset_weights`), a named
#' weight vector, or a weight matrix (hypotheses in rows, or cells in rows as
#' produced by `oneway_contrast()`). Weights are placed by design column name,
#' or by an explicit `colind` attribute; unnamed weights without `colind` are
#' accepted only when they cover every event coefficient. There is no
#' positional fallback.
#' @return list(type = "t" | "F", L = r x p matrix, colind)
#' @keywords internal
#' @noRd
.posthoc_contrast_matrix <- function(con_spec, object, beta_names, p,
                                     scope_check = NULL) {
  try_names <- function(nms) {
    if (is.null(nms)) return(NULL)
    tryCatch(.posthoc_match_names(nms, beta_names), error = function(e) NULL)
  }
  # W: cells in rows (length(colind) x r)
  build <- function(W, colind, is_f) {
    colind <- as.integer(colind)
    # Scope (e.g. event-only coefficients) is reported before generic bounds.
    if (!is.null(scope_check)) scope_check(colind)
    if (length(colind) != nrow(W)) {
      stop("contrast dimensions do not match the selected coefficients", call. = FALSE)
    }
    if (anyNA(colind) || any(colind < 1L) || any(colind > p)) {
      stop(sprintf("contrast column indices must lie in 1..%d", p), call. = FALSE)
    }
    is_f <- is_f || ncol(W) > 1L
    L <- matrix(0, nrow = ncol(W), ncol = p)
    L[, colind] <- t(W)
    list(type = if (is_f) "F" else "t", L = L, colind = colind)
  }

  if (inherits(con_spec, "contrast_spec")) {
    hit <- .posthoc_resolve_spec(con_spec, object)
    return(build(hit$W, .posthoc_match_names(rownames(hit$W), beta_names), hit$is_f))
  }

  is_f <- inherits(con_spec, "Fcontrast")
  weights <- con_spec
  if (is.list(con_spec) && !is.null(con_spec$weights)) {
    # A computed contrast object: offset_weights are named by event design
    # column; plain weights by the term's columns.
    weights <- con_spec$offset_weights %||% con_spec$weights
    if (is.null(con_spec$offset_weights)) {
      attr(weights, "colind") <- attr(con_spec$weights, "colind") %||%
        attr(con_spec, "colind") %||% con_spec$colind
    }
  }
  if (is.list(weights) || !is.numeric(weights)) {
    stop("unsupported contrast specification of class ",
         paste(class(con_spec), collapse = "/"), call. = FALSE)
  }
  colind <- attr(weights, "colind")

  if (!is.matrix(weights)) {
    w <- as.numeric(weights)
    if (is.null(colind)) {
      colind <- if (!is.null(names(weights))) {
        .posthoc_match_names(names(weights), beta_names)
      } else if (length(w) == p) {
        seq_len(p)
      } else {
        stop(sprintf(paste(
          "%d unnamed weights cannot be placed among %d event coefficients:",
          "name them by design column or attach a 'colind' attribute"
        ), length(w), p), call. = FALSE)
      }
    }
    return(build(matrix(w, ncol = 1L), colind, is_f))
  }

  W <- matrix(as.numeric(weights), nrow(weights), ncol(weights),
              dimnames = dimnames(weights))
  if (!is.null(colind)) {
    return(build(t(.orient_fcontrast(W, colind)), colind, is_f))
  }
  if (!is.null(idx <- try_names(colnames(W)))) {
    return(build(t(W), idx, is_f))           # hypotheses in rows
  }
  if (!is.null(idx <- try_names(rownames(W)))) {
    return(build(W, idx, is_f))              # cells in rows
  }
  if (ncol(W) == p) {
    return(build(t(W), seq_len(p), is_f))
  }
  if (nrow(W) == p) {
    return(build(W, seq_len(p), is_f))
  }
  stop(sprintf(paste(
    "cannot place a %d x %d weight matrix among %d event coefficients:",
    "give it design column names or a 'colind' attribute"
  ), nrow(W), ncol(W), p), call. = FALSE)
}
