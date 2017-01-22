
fit_Ftests <- function(object, model) {
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
  asgn <- object$assign[stats:::qr.lm(object)$pivot][p1]
  
  
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
    comp <- object$effects[p1,i]
    ss <- c(unlist(lapply(split(comp^2, asgn), sum)), ssr[i])
    ms <- ss/df
    f <- ms/(ssr[i]/dfr)
    
    P <- pf(f, df, dfr, lower.tail = FALSE)
    list(F=f, P=p)
  })
  
  FMat <- do.call(rbind, lapply(ret, "[[", "F"))
  PMat <- do.call(rbind, lapply(ret, "[[", "P"))
  
  list(F=FMat, P=PMat)
  
  
}


fit_contrasts <- function(lmfit, conmat) {
  if (is.vector(conmat)) {
    conmat <- matrix(conmat, 1, length(conmat))
  }
  
  Qr <- stats:::qr.lm(lmfit)
  
  p1 <- 1:lmfit$rank
  cov.unscaled <- chol2inv(Qr$qr[p1,p1,drop=FALSE])
  betamat <- lmfit$coefficients
  ct <- conmat %*% betamat
  
  rss <- colSums(lmfit$residuals^2)
  rdf <- lmfit$df.residual
  resvar <- rss/rdf
  sigma <- sqrt(resvar)
  
  vc <- sapply(1:ncol(betamat), function(i) {
    vcv <- cov.unscaled * sigma[i]^2
    sqrt(diag(conmat %*% vcv %*% t(conmat)))
  })
  
  prob <- 2 * (1 - pt(abs(ct/vc), lmfit$df.residual))
  tstat <- ct/vc
  ret <- as_tibble(cbind(as.vector(ct), as.vector(vc), as.vector(tstat), as.vector(prob)))
  names(ret) <- c("estimate", "se", "tstat", "prob")
  ret

  
}


#' @export
#' @import assertthat
contrast_set <- function(...) {
  ret <- list(...)
  assertthat::assert_that(all(sapply(ret, inherits, "contrast_spec")))
  class(ret) <- c("contrast_set", "list")
  ret
}

contrast_formula <- function(form, where=TRUE, split_by=NULL, id=NULL) {
  ret <- list(A=form,
              where=substitute(where),
              split_by=substitute(split_by),
              id=id)
  
  class(ret) <- c("contrast_formula_spec", "list")
  ret
  
}

#' @export
contrast <- function(A, B=NULL, where=TRUE, split_by=NULL, id=NULL) {
  ret <- list(A=substitute(A),
              B=substitute(B),
              where=substitute(where),
              split_by=substitute(split_by),
              id=id)

  class(ret) <- c("contrast_spec", "list")
  ret
}


#' @export
poly_contrast <- function(A, where=TRUE, degree=1, value_map=NULL, id=NULL) {
  ret <- list(
    A=substitute(A),
    B=NULL,
    where=substitute(where),
    degree=degree,
    value_map=value_map,
    id=id)
  
  class(ret) <- c("poly_contrast_spec", "contrast_spec", "list")
  ret
}

#' @export
contrast_weights.contrast_formula_spec <- function(x, term) {
  browser()
  
  term.cells <- cells(term)
  row.names(term.cells) <- longnames(term)
  
  keep <- eval(x$where, envir=term.cells, enclos=parent.frame())	
  reduced.term.cells <- subset(term.cells, keep)
  
  keepA <- eval(x$A, envir=term.cells, enclos=parent.frame())
  
  vals <- if (is.null(x$value_map)) {
    as.numeric(as.character(reduced.term.cells[keepA,]))
  } else {
    unlist(x$value_map[as.character(reduced.term.cells[keepA,])])
  }
  
  weights <- matrix(0, NROW(term.cells), x$degree)	
  pvals <- stats::poly(vals, degree=x$degree)
  row.names(weights) <- row.names(term.cells)
  colnames(weights) <- paste("poly", 1:x$degree, sep="")
  
  weights[keep, ] <- pvals
  
  ret <- list(
    weights=weights,
    contrast_spec=x)
  
  class(ret) <- c("poly_contrast", "contrast", "list")
  ret
  
}


#' @export
contrast_weights.poly_contrast_spec <- function(x, term) {
 
  term.cells <- cells(term)
  row.names(term.cells) <- longnames(term)
 
  keep <- eval(x$where, envir=term.cells, enclos=parent.frame())	
  reduced.term.cells <- subset(term.cells, keep)
  
  keepA <- eval(x$A, envir=term.cells, enclos=parent.frame())
  
  vals <- if (is.null(x$value_map)) {
    as.numeric(as.character(reduced.term.cells[keepA,]))
  } else {
    unlist(x$value_map[as.character(reduced.term.cells[keepA,])])
  }
  
  weights <- matrix(0, NROW(term.cells), x$degree)	
  pvals <- stats::poly(vals, degree=x$degree)
  row.names(weights) <- row.names(term.cells)
  colnames(weights) <- paste("poly", 1:x$degree, sep="")
  
  weights[keep, ] <- pvals
  
  ret <- list(
    weights=weights,
    contrast_spec=x)
  
  class(ret) <- c("poly_contrast", "contrast", "list")
  ret
  
}



# 
# .computeWeights <- function(cmod, pat) {
#   patmat <- patternMatch(cmod, pat)
#   if (sum(patmat) <= 0) {
#     stop(paste("failed to match any conditions with pattern: ", pat))
#   }
#   
#   cvec <- ifelse(patmat, 1, 0)
#   cvec/sum(cvec)	
# }


makeWeights <- function(keepA, keepB=NULL) {
  weights <- matrix(0, length(keepA), 1)
  numA <- sum(keepA)
  
  weights[keepA,1] <- rep(1/numA, numA)
  if (!is.null(keepB)) {
    numB <- sum(keepB)
    weights[keepB,1] <- -rep(1/numB, numB)
  }
  
  weights
}

#' @export
contrast_weights.contrast_spec <- function(x, term) {
  term.cells <- cells(term)
  row.names(term.cells) <- longnames(term)
  
  count <- attr(term.cells, "count")		
  term.cells <- subset(term.cells, count > 0)
  
  keep <- eval(x$where, envir=term.cells, enclos=parent.frame())	
  keepA <- eval(x$A, envir=term.cells, enclos=parent.frame())
  keepB <- if (is.null(x$B)) NULL else eval(x$B, envir=term.cells, enclos=parent.frame())
  split_fac <- eval(x$split_by, envir=term.cells, enclos=parent.frame())
  
  weights <- if (!is.null(split_fac)) {
    split.levs <- levels(x$split_fac)
    do.call(cbind,lapply(split.levs, function(lev) {					
      keepC <- split.fac == lev 
      if (is.null(keepB)) {
        makeWeights(keep & keepA & keepC)
      } else {
        makeWeights(keep & keepA & keepC, keep & keepB & keepC)
      }
    }))
  } else { 
    if (is.null(keepB)) {
      makeWeights(keep & keepA)
    } else {
      makeWeights(keep & keepA, keep & keepB) 
    }
  }
    
  row.names(weights) <- row.names(term.cells)
  
  ret <- list(
    weights=weights,
    contrast_spec=x)
  
  class(ret) <- c("contrast", "list")
  ret  
}

#' @importFrom gmodels estimable
#' @export
estcon.contrast <- function(x, fit, indices) {
  wts <- numeric(length(fit$assign))
  wts[indices] <- x$weights
  
  gmodels::estimable(fit, wts)
}

#' @export
print.contrast_set <- function(x) {
  for (con in x) {
    print(con)
    cat("\n")
  }
}

#' @export
print.contrast_spec <- function(x) {
  cat("contrast:", "\n")
  cat(" A: ", Reduce(paste, deparse(x$A)), "\n")
  if (!is.null(x$B))
    cat(" B: ", Reduce(paste, deparse(x$B)), "\n")
  if (x$where[1] != TRUE && length(x$where) > 1)
    cat(" where: ", x$where, "\n")

}

#' @export
print.contrast <- function(x) {
  print(x$contrast_spec)
  cat(" weights: ", x$weights)
}

#' @export
print.poly_contrast_spec <- function(x) {
  cat("poly contrast", "\n")
  cat(" A: ", Reduce(paste, deparse(x$A)), "\n")
  cat(" degree: ", x$degree, "\n")
  cat(" values: ", unlist(x$value_map), "\n")
}


