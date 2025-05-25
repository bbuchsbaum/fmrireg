#' @keywords internal
#' @noRd
qr.lm <- getFromNamespace("qr.lm", "stats")

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
beta_stats <- function(lmfit, varnames, se=TRUE) {
  cfs <- coef(lmfit)
  
  betamat <- if (is.vector(cfs)) {
    as.matrix(coef(lmfit))
  } else {
    betamat <- cfs
  }
  
  if (se) {
    Qr <- qr.lm(lmfit)
    cov.unscaled <- chol2inv(Qr$qr)
    rss <- colSums(as.matrix(lmfit$residuals^2))
    
    rdf <- lmfit$df.residual
    resvar <- rss/rdf
    sigma <- sqrt(resvar)
    
    vc <- sapply(1:ncol(betamat), function(i) {
      vcv <- cov.unscaled * sigma[i]^2
      sqrt(diag(vcv))
    })
    
    
    #prob <- 2 * (1 - pt(abs(betamat/vc), lmfit$df.residual))
    #tstat <- betamat/vc
    
    vc <- t(vc)
    colnames(vc) <- varnames
    
    betamat <- t(betamat)
    colnames(betamat) <- varnames
    #vc <- t(vc)
    
    
    ## assumes se computed
    ret <- dplyr::tibble(
      type="beta",
      name="parameter_estimates",
      stat_type="tstat",
      df.residual=rdf,
      conmat=list(NULL),
      colind=list(NULL),
      data=list(dplyr::tibble(
        estimate=list(betamat),
        se=list(vc),
        stat=list(betamat/vc),
        prob=list(2 * (1 - pt(abs(betamat/vc), lmfit$df.residual))),
        sigma=list(sigma)))
    )
  } else {
    stop("unimplemented")
    betamat <- t(betamat)
    colnames(betamat) <- varnames
    ret <- dplyr::tibble(
      type="beta",
      name="parameter_estimates",
      stat_type="effect",
      df.residual=rdf,
      conmat=list(NULL),
      colind=list(NULL),
      data=list(tibble(
        estimate=list(betamat),
        se=list(NULL),
        stat=list(NULL),
        prob=list(NULL),
        sigma=list(NULL)))
    )
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
fit_Fcontrasts <- function(lmfit, conmat, colind) {
  Qr <- qr.lm(lmfit)

  cov.unscaled <- try(chol2inv(Qr$qr))
  
  cmat <- matrix(0, nrow(conmat), ncol(cov.unscaled))
  if (ncol(conmat) == length(colind)) {
    cmat[, colind] <- conmat
  } else if (nrow(conmat) == length(colind)) {
    cmat[, colind] <- t(conmat)
  } else {
    stop(
      sprintf(
        "F contrast weight matrix dimensions %d x %d do not match length(colind) %d",
        nrow(conmat), ncol(conmat), length(colind)
      )
    )
  }
  
  cfs <- coef(lmfit)
  
  betamat <- if (is.vector(cfs)) {
    as.matrix(coef(lmfit))
  } else {
    cfs
  }
  
  rss <- colSums(as.matrix(lmfit$residuals^2))
  rdf <- lmfit$df.residual
  resvar <- rss/rdf
  
  sigma <- sqrt(resvar)
  sigma2 <- sigma^2
  
  #msr <- summary.lm(reg)$sigma  # == SSE / (n-p)
  r <- nrow(conmat)
  
  
  #                        -1
  #     (Cb - d)' ( C (X'X)   C' ) (Cb - d) / r
  # F = ---------------------------------------
  #                 SSE / (n-p)
  #
  
  cm <- solve((cmat %*% cov.unscaled %*% t(cmat)))
  
  #Fstat <- map_dbl(1:ncol(betamat), function(i) {
  #  b <- betamat[,i]
  #  cb <- cmat %*% b
  #  Fstat <- t(cb) %*% cm %*% (cb) / r / sigma2[i]
  #})
  
  estimate <- purrr::map_dbl(1:ncol(betamat), function(i) {
    b <- betamat[,i]
    cb <- cmat %*% b
    t(cb) %*% cm %*% (cb) / r 
  })
  
  se <- sigma2
  Fstat <- estimate/se
  

  return(
    structure(list(
      conmat=conmat,
      estimate=estimate,
      se=sigma2,
      df.residual=rdf,
      stat=Fstat,
      prob=1-pf(Fstat,r,rdf),
      stat_type="Fstat")
      ,class=c("Fstat", "result_stat"))
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

#' Fit Contrasts for Linear Model
#'
#' @description
#' This function calculates contrasts for a fitted linear model.
#'
#' @param lmfit The fitted linear model object.
#' @param conmat The contrast matrix or contrast vector.
#' @param colind The subset column indices in the design associated with the contrast.
#' @param se Whether to compute standard errors, t-statistics, and p-values (default: TRUE).
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
#' @noRd
fit_contrasts <- function(lmfit, conmat, colind, se=TRUE) {
  
  if (!is.matrix(conmat) && is.numeric(conmat)) {
    conmat <- t(matrix(conmat, 1, length(conmat)))
  }
  
  cfs <- coef(lmfit)
  
  if (is.vector(cfs)) {
    betamat <- as.matrix(coef(lmfit))
  } else {
    betamat <- cfs
  }
  
  ncoef <- nrow(betamat)
  cmat <- matrix(0, ncoef, 1)
  cmat[colind,] <- conmat
  
  ct <- as.vector(t(cmat) %*% betamat)
  rss <- colSums(as.matrix(lmfit$residuals^2))
  rdf <- lmfit$df.residual
  resvar <- rss/rdf
  sigma <- sqrt(resvar)
  
  p1 <- 1:lmfit$rank
  
  if (se) {
    Qr <- qr.lm(lmfit)
    cov.unscaled <- try(chol2inv(Qr$qr))
  
    if (inherits(cov.unscaled, "try-error")) {
      stop("fit_contrasts: error computing contrast covariance")
    }
  
    vc <- purrr::map_dbl(1:ncol(betamat), function(i) {
      vcv <- cov.unscaled * sigma[i]^2
      sqrt(diag(t(cmat) %*% vcv %*% cmat))
    })
  
    return(
      structure(list(
        conmat=cmat,
        sigma=sigma,
        df.residual=lmfit$df.residual,
        estimate=ct,
        se=vc,
        stat=ct/vc,
        prob=2 * (1 - pt(abs(ct/vc), lmfit$df.residual)),
        stat_type="tstat"),
        class=c("tstat", "result_stat"))
    )
  } else {
    return(
      structure(list(
        conmat=cmat,
        sigma=sigma,
        df.residual=lmfit$df.residual,
        estimate=ct,
        stat_type="effects"),
        class=c("effect", "result_stat"))
    )
  }
}

#' @keywords internal
#' @noRd
#' @autoglobal
.fast_t_contrast  <- function(B, sigma2, XtXinv, l, df) {
  # B:      p x V matrix (betas)
  # sigma2: V-vector (residual variance)
  # XtXinv: p x p matrix (inverse crossproduct of design)
  # l:      1 x p vector/matrix (contrast weights)
  # df:     scalar (residual degrees of freedom)
  
  # Ensure l is a row vector (matrix)
  if (!is.matrix(l)) {
    l <- matrix(l, nrow = 1)
  }
  
  est  <- drop(l %*% B)                # 1 × V  (BLAS GEMV)
  s2   <- as.numeric(l %*% XtXinv %*% t(l)) # scalar Var(est) / sigma2
  
  # Avoid division by zero or NaNs if s2 or sigma2 are zero/negative
  se   <- sqrt(pmax(0, s2 * sigma2))            # 1 × V
  
  tval <- ifelse(se < .Machine$double.eps^0.5, 0, est / se)
  pval <- 2 * pt(-abs(tval), df)

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
.fast_F_contrast <- function(B, sigma2, XtXinv, L, df) {
  # B:      p x V matrix (betas)
  # sigma2: V-vector (residual variance)
  # XtXinv: p x p matrix
  # L:      r x p matrix (contrast weights)
  # df:     scalar (residual degrees of freedom)
  
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
  Fval <- ifelse(abs(se) < .Machine$double.eps^0.5 | is.nan(qf), 
                 NaN, 
                 estimate / se)
                 
  pval <- pf(Fval, r, df, lower.tail = FALSE)

  list(estimate = estimate, # Numerator MS
       se       = se,       # Denominator MS (sigma2)
       stat     = Fval,
       prob     = pval,
       stat_type = "Fstat")
}

#' @keywords internal
#' @noRd
#' @autoglobal
fit_lm_contrasts_fast <- function(B, sigma2, XtXinv, conlist, fconlist, df) {
  # B:        p x V matrix (betas)
  # sigma2:   V-vector (residual variance)
  # XtXinv:   p x p matrix
  # conlist:  Named list of simple contrast vectors/matrices (1xp or px1)
  # fconlist: Named list of F contrast matrices (rxp)
  # df:       Scalar residual df

  # ----- simple contrasts (t) -----
  simples <- purrr::imap(conlist, function(l, nm) {
    # Ensure l has colind attribute attached in fmri_lm_fit
    colind <- attr(l, "colind")
    if (is.null(colind)) {
        warning(paste("Missing 'colind' attribute for simple contrast:", nm))
        # Cannot compute contrast without knowing which columns it applies to.
        # Need to decide how to handle this - skip contrast? Error?
        # For now, return NULL, will be filtered out by bind_rows.
        # Return NULL for this specific contrast, imap will handle it.
        return(NULL)
    }
    
    # Create full contrast vector/matrix padded with zeros
    p <- nrow(XtXinv)
    if (is.matrix(l)) { # Should be 1xp or px1
        if (ncol(l) == 1) l <- t(l) # Ensure it's 1xp
        if (ncol(l) != length(colind)) stop(paste("Contrast matrix columns mismatch colind for contrast:", nm))
        full_l <- matrix(0, nrow = 1, ncol = p)
        full_l[, colind] <- l
    } else { # Vector
        if (length(l) != length(colind)) stop(paste("Contrast vector length mismatch colind for contrast:", nm))
        full_l <- matrix(0, nrow = 1, ncol = p)
        full_l[, colind] <- l
    }
    
    res <- .fast_t_contrast(B, sigma2, XtXinv, full_l, df)
    # Package into tibble matching estimate_contrast.contrast output structure
    dplyr::tibble(type      = "contrast",
           name      = nm,
           stat_type = res$stat_type,
           df.residual = df,
           conmat    = list(l), # Store original contrast weights
           colind    = list(colind),
           data      = list(dplyr::tibble(
             estimate = res$estimate, # V-vector -- REMOVE list() wrapper
             se       = res$se,       # V-vector -- REMOVE list() wrapper
             stat     = res$stat,     # V-vector -- REMOVE list() wrapper
             prob     = res$prob,     # V-vector -- REMOVE list() wrapper
             sigma    = res$sigma     # V-vector -- REMOVE list() wrapper
           )))
  })

  # ----- F contrasts -----
  Fcons <- purrr::imap(fconlist, function(L, nm) {
    colind <- attr(L, "colind")
    if (is.null(colind)) {
        warning(paste("Missing 'colind' attribute for F contrast:", nm))
        return(NULL) # Return NULL for this specific contrast
    }
        
    # Create full contrast matrix padded with zeros
    p <- nrow(XtXinv)
    if (!is.matrix(L)) stop("F contrast weights must be a matrix")
    if (ncol(L) != length(colind)) stop(paste("F contrast matrix columns mismatch colind for contrast:", nm))
        
    full_L <- matrix(0, nrow = nrow(L), ncol = p)
    full_L[, colind] <- L
    
    res <- .fast_F_contrast(B, sigma2, XtXinv, full_L, df)
    # Package into tibble matching estimate_contrast.Fcontrast output structure
    dplyr::tibble(type      = "Fcontrast",
           name      = nm,
           stat_type = res$stat_type,
           df.residual = df,
           conmat    = list(L), # Store original contrast weights
           colind    = list(colind),
           data      = list(dplyr::tibble(
             estimate = res$estimate, # V-vector (Num MS) -- REMOVE list() wrapper
             se       = res$se,       # V-vector (Den MS = sigma2) -- REMOVE list() wrapper
             stat     = res$stat,     # V-vector (F stat) -- REMOVE list() wrapper
             prob     = res$prob      # V-vector -- REMOVE list() wrapper
           )))
  })

  # Return a list containing a single tibble, similar to fit_lm_contrasts output?
  # Original fit_lm_contrasts returned list(contrasts=conres, bstats=bstats, fit=fit)
  # The callers (chunkwise_lm, runwise_lm) process this.
  # Let's return just the combined tibble for now, callers need adjustment.
  # Or return a list structure mimicking the original?
  # Mimic structure: list(contrasts = list(combined_tibble), bstats = ..., fit=NULL)
  # For now, return just the tibble of contrasts.
  # Callers will need modification.
  
  # Let's return the structure expected by runwise/chunkwise logic if possible.
  # The original `fit_lm_contrasts` returned a list where each element was a contrast result (tibble).
  # Let's try returning a list of tibbles, one per contrast.
  # Combine the lists of tibbles, filtering out any NULLs from failed contrasts
  all_cons_list <- c(Filter(Negate(is.null), simples), Filter(Negate(is.null), Fcons))
  
  # Ensure the list is named correctly
  names(all_cons_list) <- unlist(lapply(all_cons_list, function(x) x$name[1]))
  
  # Return the list of tibbles
  all_cons_list
}

#' @keywords internal
#' @noRd
#' @autoglobal
beta_stats_matrix <- function(Betas, XtXinv, sigma, dfres, varnames) {
  # Betas:    p x V matrix
  # XtXinv:   p x p matrix
  # sigma:    V-vector (residual std dev)
  # dfres:    Scalar residual df
  # varnames: p-vector of coefficient names
  
  p <- nrow(Betas)
  V <- ncol(Betas)
  sigma2 <- sigma^2 # V-vector
  
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
  
  # Compute p-values
  prob <- 2 * pt(-abs(tstat), dfres)
  
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

