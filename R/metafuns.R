#' Meta-analysis using Stouffer's method
#'
#' Combines one-sided (upper-tail) p-values from several runs into a single
#' weighted Stouffer z-score, \eqn{Z = \sum_i w_i z_i / \sqrt{\sum_i w_i^2}}
#' with \eqn{z_i = \Phi^{-1}(1 - p_i)} and \eqn{w_i = 1/se_i}. Under the null
#' hypothesis \eqn{Z} is standard normal whatever the number of runs.
#'
#' @param pval A numeric matrix of upper-tail p-values (voxels x runs).
#' @param se A numeric matrix of standard errors used to form the weights
#'   \eqn{1/se}; a matrix of ones gives the unweighted Stouffer z.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{estimate}{The combined Stouffer z-scores.}
#'   \item{se}{The standard error of the combined z-score (one).}
#'   \item{stat}{Same as the estimate.}
#'   \item{prob}{The combined upper-tail p-values.}
#'   \item{stat_type}{"zfstat".}
#' }
#'
#' @keywords internal
#' @noRd
meta_stouffer <- function(pval, se) {
  pval <- as.matrix(pval)
  se <- matrix(as.numeric(se), nrow = nrow(pval), ncol = ncol(pval))
  # Work in the upper tail and keep the z-scores finite: an exact 0 or 1
  # p-value would map to +/-Inf.
  pval <- pmin(pmax(pval, .Machine$double.xmin), 1 - .Machine$double.eps / 2)
  z <- stats::qnorm(pval, lower.tail = FALSE)
  w <- 1 / se
  ok <- is.finite(z) & is.finite(w)
  w[!ok] <- 0
  z[!ok] <- 0
  den <- sqrt(rowSums(w^2))
  zcomb <- ifelse(den > 0, rowSums(w * z) / den, NA_real_)

  list(
    estimate = zcomb,
    se = rep(1, length(zcomb)),
    stat = zcomb,
    prob = stats::pnorm(zcomb, lower.tail = FALSE),
    stat_type = "zfstat"
  )
}


#' Fixed-effects meta-analysis
#'
#' This function performs a fixed-effects meta-analysis on input beta estimates and standard errors.
#' It supports two types of weighting schemes: inverse-variance weighting and equal weighting.
#' 
#' @param beta A numeric matrix of beta estimates from multiple studies or tests.
#' @param se A numeric matrix of standard errors corresponding to the beta estimates.
#' @param weighting A character string specifying the weighting scheme to use. 
#'   Options are "inv_var" (default) for inverse-variance weighting and "equal" for equal weighting.
#' 
#' @return A list containing the following elements:
#' \itemize{
#'   \item{estimate}{A numeric vector of the combined beta estimates.}
#'   \item{se}{A numeric vector of the pooled standard errors.}
#'   \item{stat}{A numeric vector of the combined z-scores.}
#'   \item{prob}{A numeric vector of the combined p-values.}
#'   \item{stat_type}{A character string indicating the type of statistic used, in this case "zstat".}
#' }
#' 
#' @keywords internal
#' @examples
#' beta <- matrix(c(0.5, 1, 0.75), nrow=1)
#' se <- matrix(c(0.2, 0.15, 0.25), nrow=1)
#' result <- meta_fixef(beta, se, weighting = "inv_var")
#' print(result)
#' @autoglobal
#' @noRd
meta_fixef <- function(ctab, weighting=c("inv_var", "equal")) {
  weighting <- match.arg(weighting)
  
  # Extract the list-column, then unlist each element before cbind
  se_list <- ctab$data %>% purrr::map(~ .$se)
  beta_list <- ctab$data %>% purrr::map(~ .$estimate)
  
  # Check if the first element is a list (indicative of fast path output)
  # If so, unlist each element. Assume consistency across the list.
  if (length(se_list) > 0 && is.list(se_list[[1]])) {
    se_list <- lapply(se_list, function(x) x[[1]])
  }
  if (length(beta_list) > 0 && is.list(beta_list[[1]])) {
    beta_list <- lapply(beta_list, function(x) x[[1]])
  }
  
  # Now se_list and beta_list should contain numeric vectors/matrices
  se <- do.call(cbind, se_list)
  beta <- do.call(cbind, beta_list)
  
  # Handle cases where cbind might produce a vector if only one contrast/run
  if (is.vector(se)) se <- matrix(se, ncol = 1)
  if (is.vector(beta)) beta <- matrix(beta, ncol = 1)
  
  # Ensure dimensions match if possible (can be tricky if runs have different numbers of voxels)
  # This assumes all runs/contrasts being pooled have the same number of rows (voxels/observations)
  if (nrow(se) != nrow(beta)) {
      stop("Mismatch in number of observations between estimates and standard errors in meta_fixef.")
  }

  ret <- do_fixef(se, beta, weighting)
  
  # Package results
  dplyr::tibble(type=ctab$type[1], name=ctab$name[1], stat_type="meta_zstat", 
                conmat=list(ctab$conmat[[1]]),
         colind=list(ctab$colind[[1]]), data=list(ret)) # ret from do_fixef is already a tibble
  
}


#' @keywords internal
#' @noRd
#' @autoglobal
do_fixef <- function(se, beta, weighting) {
  if (weighting == "inv_var") {
    inv_var <- 1/(se^2)
    # A run with se == 0 (an exact fit, or residuals below the precision of
    # the residual sum of squares) has infinite precision, and Inf/Inf would
    # turn the pooled estimate into NaN. Take the limit of inverse-variance
    # weighting instead: equal weight on the zero-variance runs, none on the
    # others.
    zero_se <- !is.na(se) & se == 0
    exact_rows <- rowSums(zero_se) > 0
    if (any(exact_rows)) {
      inv_var[exact_rows, ] <- ifelse(zero_se[exact_rows, , drop = FALSE], 1, 0)
    }
    wts <- inv_var/rowSums(inv_var)
    wbeta <- beta * wts
    wbeta <- rowSums(wbeta)
    pooledse <- sqrt(rowSums(wts*wts*(se^2)))
  } else {
    wbeta <- rowSums(beta)
    pooledse <- sqrt(rowSums((se^2)))
  }
  
  tibble(estimate=wbeta,
         se=pooledse,
         stat=wbeta/pooledse,
         prob=stats::pchisq((wbeta/pooledse)^2, 1, lower.tail = FALSE),
         stat_type="zstat")
}
  
  



#' Meta-analysis of F-contrasts
#'
#' This function performs a meta-analysis on a list of F-contrasts results obtained from different studies or tests.
#'
#' @param fres A list of F-contrast results, where each element of the list contains the results for a particular study or test.
#'
#' @return A list containing meta-analysis results for each contrast. Each element of the list includes:
#' \itemize{
#'   \item{estimate}{A numeric vector of the combined z-scores.}
#'   \item{se}{A numeric vector of the pooled standard errors.}
#'   \item{stat}{A numeric vector of the combined z-scores.}
#'   \item{prob}{A numeric vector of the combined p-values.}
#'   \item{stat_type}{A character string indicating the type of statistic used, in this case "zfstat".}
#' }
#' 
#' @keywords internal
#' @noRd
meta_Fcontrasts <- function(ftab) {
  pval <- do.call(cbind, ftab$data %>% purrr::map(~ .$prob))
  # For F-contrasts we pool evidence from p-values and avoid using sigma2-like
  # quantities as inverse-variance weights.
  se <- matrix(1, nrow = nrow(pval), ncol = ncol(pval))

  ret <- meta_stouffer(pval, se)
  dplyr::tibble(type=ftab$type[1], name=ftab$name[1], stat_type="meta_zfstat", 
         conmat=list(ftab$conmat[[1]]),
         colind=list(ftab$colind[[1]]), data=list(tibble::as_tibble(ret)))
}


#' Meta-analysis of contrasts
#'
#' This function performs a meta-analysis on a list of contrast results obtained from different studies or tests.
#'
#' @param cres A list of contrast results, where each element of the list contains the results for a particular study or test.
#' @param weighting A character string specifying the weighting method to use for the meta-analysis. Options are "inv_var" (inverse variance weighting) or "equal" (equal weighting). Default is "inv_var".
#'
#' @return A list containing meta-analysis results for each contrast, including:
#' \itemize{
#'   \item{estimate}{A matrix of the combined contrast estimates.}
#'   \item{se}{A matrix of the pooled standard errors.}
#'   \item{stat}{A matrix of the combined z-scores.}
#'   \item{prob}{A matrix of the combined p-values.}
#'   \item{stat_type}{A character string indicating the type of statistic used, in this case "zstat".}
#' }
#' 
#' @keywords internal
#' @noRd
#' @global name
meta_contrasts <- function(cres, weighting=c("inv_var", "equal")) {
  weighting <- match.arg(weighting)
  ctab <- unlist(cres, recursive = FALSE) %>% dplyr::bind_rows()

  # No contrasts were requested, so there is nothing to pool. Return an empty
  # contrast table rather than a placeholder row: a synthesised contrast with a
  # NULL estimate is indistinguishable from a real one downstream, and writers
  # and accessors already handle a zero-row table correctly.
  if (nrow(ctab) == 0) {
    return(empty_contrast_table())
  }

  gsplit <- ctab %>%
    dplyr::group_by(.data$name, .data$type) %>%
    dplyr::group_split()
  lapply(gsplit, function(tab) {
    type <- tab$type[1]
    if (type == "Fcontrast") {
      meta_Fcontrasts(tab)
    } else if (type == "contrast") {
      meta_fixef(tab)
    }
  }) %>% dplyr::bind_rows()
}

# Internal: convert t-statistics to z-scores with matching tail probability.
# Works on the log scale so that very large |t| (two-sided p below ~1e-16)
# gives a finite z instead of overflowing to Inf. Keeps the shape of `t`.
.t_to_z <- function(t, df) {
  lp <- stats::pt(-abs(t), df = df, log.p = TRUE)
  -stats::qnorm(lp, log.p = TRUE) * sign(t)
}

# Internal: build the S x P weight matrix for combine_t_statistics.
.combine_weight_matrix <- function(tmat, weights, weights_custom, se_mat, context) {
  S <- nrow(tmat)
  P <- ncol(tmat)
  from_custom <- function(msg) {
    if (is.vector(weights_custom) && length(weights_custom) == S) {
      matrix(as.numeric(weights_custom), nrow = S, ncol = P)
    } else if (is.matrix(weights_custom) && all(dim(weights_custom) == dim(tmat))) {
      weights_custom
    } else {
      warning(msg, call. = FALSE)
      NULL
    }
  }
  W <- NULL
  if (weights == "ivw") {
    if (!is.null(se_mat)) {
      W <- 1 / (se_mat^2)
    } else if (!is.null(weights_custom)) {
      W <- from_custom("weights_custom shape does not match; falling back to equal weights")
    } else {
      warning(sprintf("IVW requested for %s but no SE available; using equal weights", context),
              call. = FALSE)
    }
  } else if (weights == "custom") {
    if (is.null(weights_custom)) {
      warning("weights='custom' but 'weights_custom' not provided; using equal weights", call. = FALSE)
    } else {
      W <- from_custom("weights_custom must be length S or SxP; using equal weights")
    }
  }
  if (is.null(W)) W <- matrix(1, nrow = S, ncol = P)
  matrix(as.numeric(W), nrow = S, ncol = P)
}

# Internal: Combine t-statistics across subjects per feature
# Supports equal, inverse-variance (if SE provided), or custom weights.
# A subject whose t-statistic (or weight) is missing for a feature is dropped
# from that feature: the result equals the combination over the remaining
# subjects.
combine_t_statistics <- function(tmat, df, method = c("stouffer", "fisher", "lancaster"),
                                 weights = c("equal", "ivw", "custom"),
                                 weights_custom = NULL, se_mat = NULL) {
  method <- match.arg(method)
  weights <- match.arg(weights)
  tmat <- as.matrix(tmat)
  S <- nrow(tmat)
  # Ensure df length matches subjects
  if (length(df) == 1L) df <- rep(df, S)
  if (length(df) != S) {
    stop("Length of 'df' must be 1 or equal to number of subjects", call. = FALSE)
  }
  dfm <- matrix(as.numeric(df), nrow = S, ncol = ncol(tmat))
  # Log of the one-sided tail of |t| per subject/feature
  log_p_half <- stats::pt(-abs(tmat), df = dfm, log.p = TRUE)
  has_t <- is.finite(tmat) & is.finite(dfm) & dfm > 0

  sign_of_mean <- function() {
    tt <- tmat
    tt[!has_t] <- NA_real_
    sgn <- sign(colMeans(tt, na.rm = TRUE))
    ifelse(is.na(sgn) | sgn == 0, 1, sgn)
  }
  finish <- function(log_p_comb, m) {
    # two-sided combined p -> |z|; the log scale keeps tiny p finite
    zmag <- -stats::qnorm(log_p_comb - log(2), log.p = TRUE)
    z <- zmag * sign_of_mean()
    z[m == 0] <- NA_real_
    z
  }

  if (method == "stouffer") {
    W <- .combine_weight_matrix(tmat, weights, weights_custom, se_mat, "combine")
    z_i <- -stats::qnorm(log_p_half, log.p = TRUE) * sign(tmat)
    use <- has_t & is.finite(W)
    W[!use] <- 0
    z_i[!use] <- 0
    num <- colSums(W * z_i)
    den <- sqrt(colSums(W^2))
    return(ifelse(den > 0, num / den, NA_real_))
  }

  if (method == "fisher") {
    if (weights != "equal") {
      warning("Weighted Fisher is not supported; using equal weights", call. = FALSE)
    }
    log_p_two <- log_p_half + log(2)
    log_p_two[!has_t] <- 0
    m <- colSums(has_t)
    X2 <- -2 * colSums(log_p_two)
    log_p_comb <- stats::pchisq(X2, df = 2 * pmax(1L, m), lower.tail = FALSE, log.p = TRUE)
    return(finish(log_p_comb, m))
  }

  # lancaster (weighted Fisher via df weights)
  W <- .combine_weight_matrix(tmat, weights, weights_custom, se_mat, "lancaster")
  use <- has_t & is.finite(W) & W > 0
  W[!use] <- NA_real_
  # Normalise weights per feature by the mean of the weights in use so that
  # the total degrees of freedom equal 2 * (number of subjects used).
  mu <- colMeans(W, na.rm = TRUE)
  mu[!is.finite(mu) | mu <= 0] <- 1
  nu <- 2 * pmax(sweep(W, 2L, mu, "/"), 1e-8)
  nu[!use] <- 0
  log_p_two <- log_p_half + log(2)
  Xj <- matrix(0, nrow = S, ncol = ncol(tmat))
  Xj[use] <- stats::qchisq(log_p_two[use], df = nu[use], lower.tail = FALSE, log.p = TRUE)
  m <- colSums(use)
  log_p_comb <- stats::pchisq(colSums(Xj), df = pmax(colSums(nu), 1e-8),
                              lower.tail = FALSE, log.p = TRUE)
  finish(log_p_comb, m)
}

#' @keywords internal
#' @noRd
#' @importFrom dplyr tibble
meta_betas <- function(bstats, colind, weighting=c("inv_var", "equal")) {
  weighting <- match.arg(weighting)
  
  beta_ncols <- vapply(bstats, function(x) {
    ncol(x$data[[1]]$estimate[[1]])
  }, integer(1))
  se_ncols <- vapply(bstats, function(x) {
    ncol(x$data[[1]]$se[[1]])
  }, integer(1))
  max_col <- min(beta_ncols, se_ncols)
  if (max_col < 1L) {
    stop("No shared coefficient columns found in meta_betas.", call. = FALSE)
  }
  
  # Filter to columns available in every run. Some preprocessing paths can
  # prune run-specific nuisance columns unequally, so the first run's width is
  # not a safe bound for all runs.
  valid_colind <- colind[colind <= max_col]
  
  if (length(valid_colind) == 0) {
    warning("No valid column indices found in meta_betas. Using all columns common to every run.")
    valid_colind <- 1:max_col
  }
  
  res <- lapply(valid_colind, function(i) {
    beta <- do.call(cbind, lapply(bstats, function(x) x$data[[1]]$estimate[[1]][,i]))
    se <- do.call(cbind, lapply(bstats, function(x) x$data[[1]]$se[[1]][,i]))
    do_fixef(se,beta, weighting)
  })
  
  estimate = do.call(cbind, lapply(res, function(x) x$estimate))
  se=do.call(cbind, lapply(res, function(x) x$se))
  stat=do.call(cbind, lapply(res, function(x) x$stat))
  prob=do.call(cbind, lapply(res, function(x) x$prob))
  stat_type="zstat"
  
  ret <- dplyr::tibble(type="beta", name="parameter_estimates", stat_type="meta_zstat", conmat=list(NULL),
         colind=list(valid_colind), data=list(tibble(
           estimate=list(estimate),
           se=list(se),
           stat=list(stat),
           prob=list(prob))))

}

#' Pool runwise beta statistics by global design column
#'
#' Unlike \code{meta_betas()}, which assumes run-local column \code{j} is the
#' same coefficient in every run, this pools each global column over exactly
#' the runs whose design contains it. Event columns are shared by all runs and
#' are pooled by inverse-variance weighting; a run-specific baseline column
#' (drift, block intercept) appears in one run and keeps that run's estimate.
#'
#' @param bstats List of per-run beta-statistic tibbles.
#' @param colmaps List of integer vectors: \code{colmaps[[r]][j]} is the global
#'   column of run \code{r}'s local column \code{j}.
#' @param colind Global columns to report, in output order.
#' @return A beta-statistics tibble whose matrices have one column per
#'   \code{colind}, and \code{colind} recorded as those global indices.
#' @keywords internal
#' @noRd
meta_betas_mapped <- function(bstats, colmaps, colind, weighting = c("inv_var", "equal")) {
  weighting <- match.arg(weighting)
  stopifnot(length(bstats) == length(colmaps))
  n_vox <- nrow(as.matrix(bstats[[1]]$data[[1]]$estimate[[1]]))
  pieces <- lapply(colind, function(g) {
    runs <- which(vapply(colmaps, function(m) g %in% m, logical(1)))
    if (length(runs) == 0L) {
      na <- rep(NA_real_, n_vox)
      return(list(estimate = na, se = na, stat = na, prob = na))
    }
    beta <- do.call(cbind, lapply(runs, function(r) {
      as.matrix(bstats[[r]]$data[[1]]$estimate[[1]])[, match(g, colmaps[[r]])]
    }))
    se <- do.call(cbind, lapply(runs, function(r) {
      as.matrix(bstats[[r]]$data[[1]]$se[[1]])[, match(g, colmaps[[r]])]
    }))
    do_fixef(se, beta, weighting)
  })
  collect <- function(field) {
    do.call(cbind, lapply(pieces, function(x) as.numeric(x[[field]])))
  }

  dplyr::tibble(type = "beta", name = "parameter_estimates", stat_type = "meta_zstat",
                conmat = list(NULL), colind = list(colind),
                data = list(tibble(
                  estimate = list(collect("estimate")),
                  se = list(collect("se")),
                  stat = list(collect("stat")),
                  prob = list(collect("prob")))))
}
