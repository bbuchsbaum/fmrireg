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
    cov.unscaled <- chol2inv(qr.lm(lmfit)$qr)
    coef_se      <- sqrt(diag(cov.unscaled))
    vc           <- outer(sigma, coef_se)
    colnames(vc) <- varnames
    tstat        <- ifelse(abs(vc) < .Machine$double.eps^0.5, 0, betamat / vc)
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

fit_Fcontrasts <- function(lmfit, conmat, colind) {
  basics  <- .lm_basic_stats(lmfit)
  betamat <- basics$betamat
  sigma2  <- basics$sigma^2
  rdf     <- basics$dfres

  cov.unscaled <- chol2inv(qr.lm(lmfit)$qr)

  cmat <- matrix(0, nrow(conmat), nrow(betamat))
  if (ncol(conmat) == length(colind)) {
    cmat[, colind] <- conmat
  } else if (nrow(conmat) == length(colind)) {
    cmat[, colind] <- t(conmat)
  } else {
    stop(sprintf(
      "F contrast weight matrix dimensions %d x %d do not match length(colind) %d",
      nrow(conmat), ncol(conmat), length(colind)
    ))
  }

  r  <- nrow(conmat)
  M  <- cmat %*% cov.unscaled %*% t(cmat)
  cm <- tryCatch(solve(M), error = function(e) {
    warning(paste(
      "Singular matrix in F-contrast computation (C(X'X)^-1C'). Details:",
      e$message
    ))
    matrix(NaN, nrow(M), ncol(M))
  })

  estimate <- purrr::map_dbl(1:ncol(betamat), function(i) {
    cb <- cmat %*% betamat[, i]
    drop(t(cb) %*% cm %*% cb) / r
  })

  Fstat <- estimate / sigma2

  structure(list(
    conmat      = conmat,
    estimate    = estimate,
    se          = sigma2,
    df.residual = rdf,
    stat        = Fstat,
    prob        = 1 - pf(Fstat, r, rdf),
    stat_type   = "Fstat"),
    class = c("Fstat", "result_stat")
  )
}

#' @keywords internal
#' @noRd
#' @autoglobal
estimate_contrast.contrast <- function(x, fit, colind, ...) {
  ret <- fit_contrasts(fit, x$weights, colind, se=TRUE)
  
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
  ret <- fit_Fcontrasts(fit, x$weights, colind)
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
fit_contrasts.default <- function(object, conmat, colind, se = TRUE, ...) {

  conmat <- as.matrix(conmat)

  basics  <- .lm_basic_stats(object)
  betamat <- basics$betamat
  sigma   <- basics$sigma
  rdf     <- basics$dfres

  ncoef <- nrow(betamat)
  cmat  <- matrix(0, ncoef, 1)
  cmat[colind, ] <- conmat

  ct <- drop(t(cmat) %*% betamat)

  if (se) {
    cov.unscaled <- chol2inv(qr.lm(object)$qr)
    var_est      <- as.numeric(t(cmat) %*% cov.unscaled %*% cmat)
    se_vec       <- sqrt(var_est) * sigma

    structure(list(
      conmat      = cmat,
      sigma       = sigma,
      df.residual = rdf,
      estimate    = ct,
      se          = se_vec,
      stat        = ct / se_vec,
      prob        = 2 * pt(-abs(ct / se_vec), rdf),
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
                               robust_weights = NULL, ar_order = 0) {
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
  
  est  <- drop(l %*% B)                # 1 × V  (BLAS GEMV)
  s2   <- as.numeric(l %*% XtXinv %*% t(l)) # scalar Var(est) / sigma2
  
  # Avoid division by zero or NaNs if s2 or sigma2 are zero/negative
  se   <- sqrt(pmax(0, s2 * sigma2))            # 1 × V
  
  # Calculate effective df if needed
  if (!is.null(robust_weights) || ar_order > 0) {
    p <- nrow(B)
    n <- df + p  # Total observations
    df_effective <- calculate_effective_df(n, p, robust_weights, ar_order, method = "simple")
  } else {
    df_effective <- df
  }
  
  tval <- ifelse(se < .Machine$double.eps^0.5, 0, est / se)
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
                              robust_weights = NULL, ar_order = 0) {
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
  
  r   <- nrow(L)
  U   <- L %*% B                    # r × V   (BLAS GEMM) : L*beta
  M   <- L %*% XtXinv %*% t(L)      # r × r   : L (X'X)^-1 L'
  
  # Use tryCatch for solve() in case M is singular
  Cinv <- tryCatch(solve(M), error = function(e) {
    warning(paste("Singular matrix in F-contrast computation (L(X'X)^-1 L'). Details:", e$message))
    # Return a matrix of NaNs or zeros? Let's use NaNs.
    matrix(NaN, nrow = r, ncol = r)
  })

  # qf = diag(t(U) %*% Cinv %*% U)
  # Efficient computation: colSums((t(U) %*% Cinv) * t(U))
  # Check if Cinv contains NaNs
  if (any(is.nan(Cinv))) {
      qf <- rep(NaN, ncol(U))
  } else {
      tmp  <- t(U) %*% Cinv           # V x r
      qf  <- colSums(tmp * t(U))      # V-vector, each element is u_v' Cinv u_v
  }
  
  estimate <- qf / r # Numerator mean square: (LB)' (L(X'X)^-1 L')^-1 (LB) / r
  se <- sigma2 # Denominator mean square (residual variance)
  
  # Avoid division by zero/NaNs
  # Calculate effective df if needed
  if (!is.null(robust_weights) || ar_order > 0) {
    p <- nrow(B)
    n <- df + p  # Total observations
    df_effective <- calculate_effective_df(n, p, robust_weights, ar_order, method = "simple")
  } else {
    df_effective <- df
  }
  
  Fval <- ifelse(abs(se) < .Machine$double.eps^0.5 | is.nan(qf), 
                 NaN, 
                 estimate / se)
                 
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
                                  robust_weights = NULL, ar_order = 0) {
  # Validate inputs
  if (!is.matrix(B)) stop("B must be a matrix (p x V)")
  if (!is.matrix(XtXinv)) stop("XtXinv must be a matrix (p x p)")
  if (nrow(B) != nrow(XtXinv) || nrow(XtXinv) != ncol(XtXinv)) {
    stop("Dimension mismatch: B must be p x V, XtXinv must be p x p")
  }
  
  # Process t-contrasts
  tcontrast_results <- process_t_contrasts(
    B = B, sigma2 = sigma2, XtXinv = XtXinv, 
    conlist = conlist, df = df,
    robust_weights = robust_weights, ar_order = ar_order
  )
  
  # Process F-contrasts  
  fcontrast_results <- process_f_contrasts(
    B = B, sigma2 = sigma2, XtXinv = XtXinv,
    fconlist = fconlist, df = df,
    robust_weights = robust_weights, ar_order = ar_order
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
                                robust_weights = NULL, ar_order = 0) {
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
                               robust_weights, ar_order)
      
      # Package results in expected format
      results[[i]] <- package_tcontrast_result(
        con_name, contrast_weights, colind, df, stats
      )
      
    }, error = function(e) {
      warning(sprintf("Failed to compute t-contrast '%s': %s", con_name, e$message))
      results[[i]] <- NULL
    })
  }
  
  # Filter out failed contrasts
  Filter(Negate(is.null), results)
}

#' Process F-contrasts for fast linear model fitting  
#' @keywords internal
#' @noRd
#' @autoglobal
process_f_contrasts <- function(B, sigma2, XtXinv, fconlist, df,
                                robust_weights = NULL, ar_order = 0) {
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
                               robust_weights, ar_order)
      
      # Package results in expected format
      results[[i]] <- package_fcontrast_result(
        con_name, contrast_matrix, colind, df, stats
      )
      
    }, error = function(e) {
      warning(sprintf("Failed to compute F-contrast '%s': %s", con_name, e$message))
      results[[i]] <- NULL  
    })
  }
  
  # Filter out failed contrasts
  Filter(Negate(is.null), results)
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
  
  # Validate dimensions
  if (ncol(contrast_matrix) != length(colind)) {
    stop(sprintf(
      "Dimension mismatch: contrast matrix columns (%d) != colind length (%d)",
      ncol(contrast_matrix), length(colind)
    ))
  }
  
  if (max(colind) > p) {
    stop(sprintf(
      "Column index out of bounds: max colind (%d) > design matrix columns (%d)", 
      max(colind), p
    ))
  }
  
  # Create full contrast matrix
  full_contrast <- matrix(0, nrow = nrow(contrast_matrix), ncol = p)
  full_contrast[, colind] <- contrast_matrix
  
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
  
  # Calculate effective degrees of freedom if needed
  if (!is.null(robust_weights) || ar_order > 0) {
    n <- dfres + p  # Total observations
    df_effective <- calculate_effective_df(n, p, robust_weights, ar_order, method = "simple")
  } else {
    df_effective <- dfres
  }
  
  # Compute SE for each beta, each voxel
  # se(beta_i) = sqrt( [XtXinv]_ii * sigma2 )
  diag_XtXinv <- diag(XtXinv)
  if (any(diag_XtXinv < 0)) {
      # Handle potential numerical issues if diagonal is negative
      warning("Negative diagonal elements found in XtXinv during beta SE calculation.")
      diag_XtXinv[diag_XtXinv < 0] <- NaN
  }
  se_scaling <- sqrt(diag_XtXinv) # p-vector
  
  # Outer product: V-vector sigma with p-vector se_scaling -> V x p matrix
  vc <- outer(sigma, se_scaling) # SEs for all betas (rows=voxels, cols=betas)
  
  betamat <- t(Betas) # V x p matrix (voxels x betas)
  colnames(betamat) <- varnames
  colnames(vc) <- varnames
  
  # Compute t-stats: element-wise division
  # Avoid division by zero
  tstat <- ifelse(abs(vc) < .Machine$double.eps^0.5, 0, betamat / vc)
  
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
  # Get the fitted model components from the fmri_lm object
  betas <- t(object$result$betas$data[[1]]$estimate[[1]])  # p x V matrix (transpose to get correct orientation)
  sigma <- object$result$sigma                              # V-vector
  df_residual <- object$result$rdf                         # scalar
  
  # We need the design matrix to compute XtXinv
  # Get this from the model object
  fmrimod <- object$model
  tmats <- term_matrices(fmrimod)
  data_env <- list2env(tmats)
  data_env[[".y"]] <- rep(0, nrow(tmats[[1]]))
  form <- get_formula(fmrimod)
  X <- model.matrix(form, data_env)
  
  # Compute XtXinv
  XtX <- crossprod(X)
  XtXinv <- tryCatch(
    chol2inv(chol(XtX)),
    error = function(e) {
      # Use SVD-based pseudoinverse for singular matrices
      svd_result <- svd(XtX)
      d <- svd_result$d
      tol <- max(dim(XtX)) * .Machine$double.eps * max(d)
      pos <- d > tol
      if (sum(pos) == 0) {
        stop("Completely singular design matrix in contrast computation")
      }
      svd_result$v[, pos] %*% diag(1/d[pos], nrow = sum(pos)) %*% t(svd_result$u[, pos])
    }
  )
  
  # Process each contrast
  contrast_results <- lapply(names(contrasts), function(con_name) {
    con_spec <- contrasts[[con_name]]
    
    tryCatch({
      # Get contrast weights - this should be computed from the contrast specification
      # For now, assume it's a pair_contrast and extract weights
      if (inherits(con_spec, "pair_contrast_spec")) {
        # Get the term for this contrast (assume first term for simplicity)
        event_terms <- terms(object$model$event_model)
        if (length(event_terms) > 0) {
          term <- event_terms[[1]]
          con_weights_obj <- contrast_weights(con_spec, term)
          con_weights <- as.vector(con_weights_obj$weights)
          
          # Find matching columns in design matrix
          # For simplicity, use the condition column indices
          condition_cols <- grep("condition", colnames(X))
          if (length(condition_cols) >= length(con_weights)) {
            colind <- condition_cols[1:length(con_weights)]
          } else {
            # Fallback to first few columns  
            colind <- 1:length(con_weights)
          }
        } else {
          stop("No event terms available for contrast computation")
        }
      } else {
        stop("Unsupported contrast type: ", class(con_spec))
      }
      
      # Compute contrast using simple matrix multiplication
      if (length(colind) > 0 && length(con_weights) > 0 && length(colind) == length(con_weights)) {
        # Simple contrast computation: weights %*% betas for relevant columns
        contrast_estimates <- drop(con_weights %*% betas[colind, , drop = FALSE])
        
        # For now, return a simplified result
        list(
          name = con_name,
          estimate = contrast_estimates,
          se = rep(1, length(contrast_estimates)),  # Placeholder
          stat = contrast_estimates,  # Placeholder
          prob = rep(0.05, length(contrast_estimates))  # Placeholder
        )
      } else {
        warning(paste("Dimension mismatch for contrast", con_name, 
                     "- colind:", length(colind), "weights:", length(con_weights)))
        NULL
      }
    }, error = function(e) {
      warning(paste("Error processing contrast", con_name, ":", e$message))
      NULL
    })
  })
  
  # Filter out NULL results and set names
  contrast_results <- contrast_results[!sapply(contrast_results, is.null)]
  names(contrast_results) <- sapply(contrast_results, function(x) x$name)
  
  contrast_results
}

