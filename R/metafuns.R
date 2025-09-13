#' Meta-analysis using Stouffer's method
#'
#' This function performs a meta-analysis on input p-values and standard errors using Stouffer's method.
#' Stouffer's method combines z-scores by weighting them by their inverse variance.
#' 
#' @param pval A numeric vector of p-values from multiple studies or tests.
#' @param se A numeric vector of standard errors corresponding to the p-values.
#' 
#' @return A list containing the following elements:
#' \itemize{
#'   \item{estimate}{A numeric vector of the weighted z-scores obtained using Stouffer's method.}
#'   \item{se}{A numeric vector of the combined standard errors.}
#'   \item{stat}{Same as the estimate, the weighted z-scores.}
#'   \item{prob}{A numeric vector of the combined p-values.}
#'   \item{stat_type}{A character string indicating the type of statistic used, in this case "zfstat".}
#' }
#' 
#' @keywords internal
#' @noRd
meta_stouffer <- function(pval, se) {
  if (any(pval == 0)) {
    pval[pval == 0] <- .Machine$double.eps
  }
  inv_var <- 1/(se^2)
  wts <- inv_var/rowSums(inv_var)
  zscore <- qnorm(1-pval)
  wzscore <- zscore * wts
  wzscore <- rowSums(wzscore)
  
  return(
    list(
      estimate=wzscore,
      se=sqrt(rowSums(wts*wts*(se^2))),
      stat=wzscore,
      prob=1-pnorm(wzscore),
      stat_type="zfstat"
    )
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
         prob=1-pchisq((wbeta/pooledse)^2,1),
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
  se <- do.call(cbind, ftab$data %>% purrr::map(~ .$se))
  
  ret <- meta_stouffer(pval,se)
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
  ctab <- unlist(cres, recursive=FALSE) %>% dplyr::bind_rows()
  if (nrow(ctab) == 0) {
    return(dplyr::tibble(type="contrast", name="meta_contrast", stat_type="meta_zstat", conmat=list(NULL),
         colind=list(NULL), data=list(tibble(estimate=list(NULL), se=list(NULL), stat=list(NULL), prob=list(NULL)))))
  }

  gsplit <- ctab %>% dplyr::group_by(name,type) %>% dplyr::group_split()
  
  lapply(gsplit, function(tab) {
    type <- tab$type[1]
    if (type == "Fcontrast") {
      meta_Fcontrasts(tab)
    } else if (type == "contrast") {
      meta_fixef(tab)
    }
  }) %>% dplyr::bind_rows()
  
}

# Internal: Combine t-statistics across subjects per feature
# Supports equal, inverse-variance (if SE provided), or custom weights.
combine_t_statistics <- function(tmat, df, method = c("stouffer", "fisher", "lancaster"),
                                 weights = c("equal", "ivw", "custom"),
                                 weights_custom = NULL, se_mat = NULL) {
  method <- match.arg(method)
  weights <- match.arg(weights)
  # Ensure df length matches subjects
  if (length(df) == 1L) df <- rep(df, nrow(tmat))
  if (length(df) != nrow(tmat)) {
    stop("Length of 'df' must be 1 or equal to number of subjects", call. = FALSE)
  }
  # Two-tailed p-values per subject/voxel
  p_two <- matrix(NA_real_, nrow = nrow(tmat), ncol = ncol(tmat))
  for (i in seq_len(nrow(tmat))) {
    ti <- tmat[i, ]
    p_two[i, ] <- 2 * stats::pt(abs(ti), df = df[i], lower.tail = FALSE)
  }
  if (method == "stouffer") {
    z_i <- stats::qnorm(pmax(1e-300, 1 - p_two/2)) * sign(tmat)
    W <- matrix(1, nrow = nrow(tmat), ncol = ncol(tmat))
    if (weights == "ivw") {
      if (!is.null(se_mat)) {
        W <- 1/(se_mat^2)
      } else if (!is.null(weights_custom)) {
        if (is.vector(weights_custom) && length(weights_custom) == nrow(tmat)) {
          W <- matrix(weights_custom, nrow = nrow(tmat), ncol = ncol(tmat))
        } else if (is.matrix(weights_custom) && all(dim(weights_custom) == dim(tmat))) {
          W <- weights_custom
        } else {
          warning("weights_custom shape does not match; falling back to equal weights", call. = FALSE)
        }
      } else {
        warning("IVW requested for combine but no SE available; using equal weights", call. = FALSE)
      }
    } else if (weights == "custom") {
      if (is.null(weights_custom)) {
        warning("weights='custom' but 'weights_custom' not provided; using equal weights", call. = FALSE)
      } else if (is.vector(weights_custom) && length(weights_custom) == nrow(tmat)) {
        W <- matrix(weights_custom, nrow = nrow(tmat), ncol = ncol(tmat))
      } else if (is.matrix(weights_custom) && all(dim(weights_custom) == dim(tmat))) {
        W <- weights_custom
      } else {
        warning("weights_custom must be length S or SxP; using equal weights", call. = FALSE)
      }
    }
    # Ensure matrix shape
    W <- matrix(as.numeric(W), nrow = nrow(tmat), ncol = ncol(tmat))
    num <- colSums(W * z_i, na.rm = TRUE)
    den <- sqrt(colSums(W^2, na.rm = TRUE))
    z_comb <- ifelse(den > 0, num/den, NA_real_)
    return(z_comb)
  } else if (method == "fisher") {
    if (weights != "equal") {
      warning("Weighted Fisher is not supported; using equal weights", call. = FALSE)
    }
    m <- colSums(is.finite(p_two))
    X2 <- -2 * colSums(log(p_two), na.rm = TRUE)
    p_comb <- stats::pchisq(X2, df = 2 * pmax(1L, m), lower.tail = FALSE)
    zmag <- stats::qnorm(pmax(1e-300, 1 - p_comb/2))
    sgn <- sign(colMeans(tmat, na.rm = TRUE))
    z_comb <- zmag * ifelse(sgn == 0, 1, sgn)
    return(z_comb)
  } else { # lancaster (weighted Fisher via df weights)
    # Build weights per subject/voxel
    if (weights == "equal") {
      W <- matrix(1, nrow = nrow(tmat), ncol = ncol(tmat))
    } else if (weights == "ivw") {
      if (!is.null(se_mat)) {
        W <- 1/(se_mat^2)
      } else if (!is.null(weights_custom)) {
        if (is.vector(weights_custom) && length(weights_custom) == nrow(tmat)) {
          W <- matrix(weights_custom, nrow = nrow(tmat), ncol = ncol(tmat))
        } else if (is.matrix(weights_custom) && all(dim(weights_custom) == dim(tmat))) {
          W <- weights_custom
        } else {
          warning("weights_custom shape does not match; falling back to equal weights", call. = FALSE)
          W <- matrix(1, nrow = nrow(tmat), ncol = ncol(tmat))
        }
      } else {
        warning("IVW requested for lancaster but no SE available; using equal weights", call. = FALSE)
        W <- matrix(1, nrow = nrow(tmat), ncol = ncol(tmat))
      }
    } else { # custom
      if (is.null(weights_custom)) {
        warning("weights='custom' but 'weights_custom' not provided; using equal weights", call. = FALSE)
        W <- matrix(1, nrow = nrow(tmat), ncol = ncol(tmat))
      } else if (is.vector(weights_custom) && length(weights_custom) == nrow(tmat)) {
        W <- matrix(weights_custom, nrow = nrow(tmat), ncol = ncol(tmat))
      } else if (is.matrix(weights_custom) && all(dim(weights_custom) == dim(tmat))) {
        W <- weights_custom
      } else {
        warning("weights_custom must be length S or SxP; using equal weights", call. = FALSE)
        W <- matrix(1, nrow = nrow(tmat), ncol = ncol(tmat))
      }
    }
    # Lancaster: df_j = 2 * normalized weights so sum df ≈ 2m
    # Normalize weights per voxel by mean of finite weights
    Wn <- W
    for (j in seq_len(ncol(W))) {
      ww <- W[, j]
      ok <- is.finite(ww) & apply(!is.na(p_two[, j, drop = FALSE]), 1, all)
      mu <- mean(ww[ok], na.rm = TRUE)
      if (!is.finite(mu) || mu <= 0) mu <- 1
      Wn[, j] <- ww / mu
    }
    # Degrees per subject/voxel; coerce back to matrix to avoid dim drop
    nu <- 2 * pmax(1e-8, Wn)
    nu <- matrix(as.numeric(nu), nrow = nrow(Wn), ncol = ncol(Wn))
    # Combine via chi-square with df per subject
    # Compute chi-square quantiles per subject/voxel
    Xj <- matrix(NA_real_, nrow = nrow(tmat), ncol = ncol(tmat))
    for (i in seq_len(nrow(tmat))) {
      Xj[i, ] <- stats::qchisq(pmax(1e-300, 1 - p_two[i, ]), df = nu[i, ])
    }
    Xsum <- colSums(Xj, na.rm = TRUE)
    nu_sum <- colSums(nu, na.rm = TRUE)
    p_comb <- stats::pchisq(Xsum, df = nu_sum, lower.tail = FALSE)
    zmag <- stats::qnorm(pmax(1e-300, 1 - p_comb/2))
    sgn <- sign(colMeans(tmat, na.rm = TRUE))
    z_comb <- zmag * ifelse(sgn == 0, 1, sgn)
    return(z_comb)
  }
}


#' @keywords internal
#' @noRd
#' @importFrom dplyr tibble
meta_betas <- function(bstats, colind, weighting=c("inv_var", "equal")) {
  weighting <- match.arg(weighting)

  len <- length(colind)
  
  # Check the dimensions of the first beta matrix to understand the structure
  first_beta <- bstats[[1]]$data[[1]]$estimate[[1]]
  max_col <- ncol(first_beta)
  
  # Filter colind to only include valid column indices
  valid_colind <- colind[colind <= max_col]
  
  if (length(valid_colind) == 0) {
    warning("No valid column indices found in meta_betas. Using all available columns.")
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
