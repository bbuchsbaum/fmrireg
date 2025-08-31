#' fmrireg t-tests for Group Analysis
#'
#' Performs group-level t-tests on fMRI data with support for one-sample, two-sample,
#' paired, and ANCOVA designs. Provides both meta-analysis and classical t-test engines,
#' mirroring AFNI 3dttest++ functionality within the fmrireg framework.
#'
#' @param gd A group_data object or data frame with subject data
#' @param formula R formula for between-subjects design. Examples:
#'   \itemize{
#'     \item ~ 1: One-sample t-test
#'     \item ~ 1 + group: Two-sample t-test
#'     \item ~ 1 + group + age: ANCOVA with covariate
#'   }
#' @param engine Character string; analysis engine to use:
#'   \itemize{
#'     \item "auto": Automatically select based on data (default)
#'     \item "meta": Use inverse-variance meta-analysis (requires SE)
#'     \item "classic": Use OLS/Student t-test
#'     \item "welch": Use Welch t-test for unequal variances
#'   }
#' @param paired Logical; if TRUE, perform paired t-test on differences (default: FALSE)
#' @param mu0 Numeric scalar; constant to test against for one-sample tests (default: 0)
#' @param contrast Optional named numeric vector for linear combination of coefficients
#' @param mc Character string or NULL; multiple comparisons correction:
#'   \itemize{
#'     \item NULL: No correction (default)
#'     \item "bh": Benjamini-Hochberg FDR
#'     \item "spatial_fdr": Spatially-aware FDR
#'   }
#' @param alpha Numeric scalar; FDR level if mc is not NULL (default: 0.05)
#' @param sign Character string; sign convention for group differences:
#'   \itemize{
#'     \item "AminusB": A - B (default)
#'     \item "BminusA": B - A
#'   }
#' @param mask Optional mask object or path
#' @param voxelwise_cov Optional S x P matrix of voxelwise covariates
#' @param center_voxelwise Logical; center voxelwise covariate per feature (default: TRUE)
#' @param voxel_name Character string; name for voxelwise coefficient (default: "voxel_cov")
#'
#' @return An fmri_ttest_fit object containing:
#'   \itemize{
#'     \item beta: Matrix of coefficients
#'     \item se: Matrix of standard errors (if available)
#'     \item t: Matrix of t-statistics (classic engine)
#'     \item z: Matrix of z-scores
#'     \item p: Matrix of p-values
#'     \item df: Matrix of degrees of freedom
#'     \item q: Matrix of FDR-adjusted p-values (if mc is used)
#'   }
#'
#' @examples
#' \dontrun{
#' # One-sample t-test
#' fit <- fmri_ttest(gd, formula = ~ 1)
#' 
#' # Two-sample t-test
#' fit <- fmri_ttest(gd, formula = ~ 1 + group)
#' 
#' # Paired t-test
#' fit <- fmri_ttest(gd_diff, formula = ~ 1, paired = TRUE)
#' 
#' # ANCOVA with age covariate
#' fit <- fmri_ttest(gd, formula = ~ 1 + group + age)
#' 
#' # With spatial FDR correction
#' fit <- fmri_ttest(gd, formula = ~ 1 + group, mc = "spatial_fdr", alpha = 0.05)
#' }
#'
#' @export
fmri_ttest <- function(gd,
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
                      voxel_name = "voxel_cov") {
  
  engine <- match.arg(engine)
  sign <- match.arg(sign)
  
  if (!is.null(mc)) {
    mc <- match.arg(mc, c("bh", "spatial_fdr"))
  }
  
  # 1) Normalize input to matrices Y (S x P), optional V (S x P), covars
  if (inherits(gd, "group_data")) {
    # Extract first block (assuming single contrast)
    if (length(gd$blocks) != 1) {
      stop("fmri_ttest expects a single contrast block", call. = FALSE)
    }
    B <- gd$blocks[[1]]
    Y <- B$Y  # S x P effects (betas or derived)
    V <- if (!is.null(B$V)) B$V else NULL  # S x P sampling variances if available
    covars <- B$covars
    feature_group <- if (!is.null(B$feature)) B$feature$group else NULL
  } else {
    # Convert data frame to group_data
    if (!inherits(gd, "data.frame")) {
      stop("gd must be a group_data object or data frame", call. = FALSE)
    }
    gd <- group_data(gd, mask = mask)
    B <- gd$blocks[[1]]
    Y <- B$Y
    V <- if (!is.null(B$V)) B$V else NULL
    covars <- B$covars
    feature_group <- NULL
  }
  
  # 2) Subtract mu0 from Y for testing against constant
  Y <- Y - mu0
  
  # 3) Build design matrix X from formula
  if (is.null(covars)) {
    covars <- data.frame(intercept = rep(1, nrow(Y)))
  }
  
  X <- stats::model.matrix(stats::as.formula(formula), data = covars)
  storage.mode(X) <- "double"
  K <- ncol(X)
  S <- nrow(X)
  P <- ncol(Y)
  
  if (S != nrow(Y)) {
    stop("Design matrix rows must match data rows", call. = FALSE)
  }
  
  # 4) Choose engine
  use_meta <- if (engine == "auto") {
    !is.null(V)
  } else {
    engine == "meta"
  }
  
  res <- list()
  
  if (paired) {
    # Paired: assume Y contains within-subject differences
    message("paired=TRUE: performing one-sample test on differences vs mu0")
    
    # One-sample test (intercept only)
    Xp <- cbind("(Intercept)" = 1)
    
    if (!is.null(voxelwise_cov)) {
      ols <- fmri_ols_fit(Y, Xp, voxelwise = voxelwise_cov,
                          center_voxelwise = center_voxelwise,
                          voxel_name = voxel_name)
    } else {
      ols <- ols_t_cpp(Y, Xp)
      rownames(ols$beta) <- rownames(ols$se) <- rownames(ols$t) <- "(Intercept)"
    }
    
    T <- ols$t[1, , drop = FALSE]
    df <- rep(ols$df, P)
    Pval <- 2 * stats::pt(abs(as.numeric(T)), df = df, lower.tail = FALSE)
    Z <- stats::qnorm(pmax(1e-300, 1 - Pval/2)) * sign(as.numeric(T))
    
    res$beta <- ols$beta
    res$se <- ols$se
    res$t <- T
    res$df <- matrix(df, nrow = 1, ncol = P)
    res$z <- matrix(Z, nrow = 1, ncol = P)
    res$p <- matrix(Pval, nrow = 1, ncol = P)
    
    rownames(res$beta) <- rownames(res$se) <- rownames(res$t) <- 
      rownames(res$df) <- rownames(res$z) <- rownames(res$p) <- "(Intercept)"
    
  } else if (use_meta) {
    # Meta-regression using existing infrastructure
    method <- getOption("fmrireg.meta.method", "pm")
    robust <- getOption("fmrireg.meta.robust", "none")
    
    if (!is.null(voxelwise_cov)) {
      # Update meta_fit wrapper to handle voxelwise
      out <- fmri_meta_fit_extended(
        Y = Y, V = V, X = X,
        method = method, robust = robust,
        voxelwise = voxelwise_cov,
        center_voxelwise = center_voxelwise,
        voxel_name = voxel_name,
        n_threads = getOption("fmrireg.num_threads", 0)
      )
    } else {
      out <- fmri_meta_fit(
        Y = Y, V = V, X = X,
        method = method, robust = robust,
        n_threads = getOption("fmrireg.num_threads", 0)
      )
    }
    
    res$beta <- out$beta
    res$se <- out$se
    res$z <- out$z
    res$p <- 2 * stats::pnorm(abs(out$z), lower.tail = FALSE)
    res$df <- matrix(Inf, nrow = nrow(out$z), ncol = ncol(out$z))
    
    # Handle contrast if specified
    if (!is.null(contrast)) {
      w <- as.numeric(contrast[colnames(out$beta)])
      w[is.na(w)] <- 0
      num <- colSums(out$beta * w, na.rm = TRUE)
      den <- sqrt(colSums((out$se * w)^2, na.rm = TRUE))
      zc <- num / den
      res$z_contrast <- zc
      res$p_contrast <- 2 * stats::pnorm(abs(zc), lower.tail = FALSE)
    }
    
  } else {
    # Classic Student t or Welch
    has_group <- "group" %in% colnames(covars) && 
                 is.factor(covars$group) && 
                 nlevels(covars$group) == 2
    
    if (engine == "welch" && has_group && K <= 2) {
      # Welch t-test
      g <- as.integer(covars$group)
      w <- welch_t_cpp(Y, g)
      
      T <- as.numeric(w$t)
      df <- as.numeric(w$df)
      Pval <- 2 * stats::pt(abs(T), df = df, lower.tail = FALSE)
      Z <- stats::qnorm(pmax(1e-300, 1 - Pval/2)) * sign(T)
      
      res$beta <- rbind(
        "(Intercept)" = (w$muA + w$muB) / 2,
        group = w$muA - w$muB
      )
      res$se <- matrix(NA_real_, nrow = 2, ncol = P,
                       dimnames = list(c("(Intercept)", "group"), NULL))
      res$t <- rbind("(Intercept)" = NA_real_, group = T)
      res$df <- rbind("(Intercept)" = NA_real_, group = df)
      res$z <- rbind("(Intercept)" = NA_real_, group = Z)
      res$p <- rbind("(Intercept)" = NA_real_, group = Pval)
      
    } else {
      # General OLS t-test/ANCOVA
      if (!is.null(voxelwise_cov)) {
        ols <- fmri_ols_fit(Y, X, voxelwise = voxelwise_cov,
                           center_voxelwise = center_voxelwise,
                           voxel_name = voxel_name)
      } else {
        ols <- ols_t_cpp(Y, X)
        rownames(ols$beta) <- rownames(ols$se) <- rownames(ols$t) <- colnames(X)
      }
      
      # Convert t to z and p
      df <- rep(ols$df, P)
      Pval <- apply(ols$t, 1, function(tt) {
        2 * stats::pt(abs(tt), df = df, lower.tail = FALSE)
      })
      Pval <- t(Pval)
      rownames(Pval) <- rownames(ols$t)
      
      Z <- matrix(NA_real_, nrow = nrow(ols$t), ncol = P,
                  dimnames = dimnames(ols$t))
      for (i in seq_len(nrow(Z))) {
        Zi <- stats::qnorm(pmax(1e-300, 1 - Pval[i,]/2)) * sign(ols$t[i,])
        Z[i,] <- Zi
      }
      
      res$beta <- ols$beta
      res$se <- ols$se
      res$t <- ols$t
      res$df <- matrix(ols$df, nrow = nrow(ols$t), ncol = ncol(ols$t),
                       byrow = TRUE, dimnames = dimnames(ols$t))
      res$z <- Z
      res$p <- Pval
    }
  }
  
  # 5) Apply sign convention if needed
  if (sign == "BminusA" && "group" %in% rownames(res$beta)) {
    group_idx <- which(rownames(res$beta) == "group")
    res$beta[group_idx, ] <- -res$beta[group_idx, ]
    if (!is.null(res$t)) res$t[group_idx, ] <- -res$t[group_idx, ]
    if (!is.null(res$z)) res$z[group_idx, ] <- -res$z[group_idx, ]
  }
  
  # 6) Multiple comparisons correction
  if (!is.null(mc)) {
    if (!is.null(res$z_contrast)) {
      # Correct contrast
      if (mc == "bh") {
        res$q_contrast <- stats::p.adjust(res$p_contrast, method = "BH")
      } else if (mc == "spatial_fdr" && !is.null(feature_group)) {
        out <- spatial_fdr(z = as.numeric(res$z_contrast),
                          group = feature_group,
                          alpha = alpha)
        res$q_contrast <- out$q
        res$reject_contrast <- out$reject
      }
    } else {
      # Correct each coefficient
      corrs <- vector("list", nrow(res$z))
      for (i in seq_len(nrow(res$z))) {
        if (mc == "bh") {
          corrs[[i]] <- stats::p.adjust(res$p[i,], method = "BH")
        } else if (mc == "spatial_fdr" && !is.null(feature_group)) {
          out <- spatial_fdr(z = as.numeric(res$z[i,]),
                            group = feature_group,
                            alpha = alpha)
          corrs[[i]] <- out$q
        } else {
          corrs[[i]] <- res$p[i,]
        }
      }
      res$q <- do.call(rbind, corrs)
      rownames(res$q) <- rownames(res$z)
    }
  }
  
  # Store metadata
  res$call <- match.call()
  res$formula <- formula
  res$engine <- if (use_meta) "meta" else engine
  res$n_subjects <- S
  res$n_features <- P
  
  class(res) <- c("fmri_ttest_fit", "list")
  res
}

#' Print method for fmri_ttest_fit
#' 
#' @param x An fmri_ttest_fit object
#' @param ... Additional print arguments
#' @export
print.fmri_ttest_fit <- function(x, ...) {
  cat("fmri_ttest Results\n")
  cat("==================\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Engine:", x$engine, "\n")
  cat("Subjects:", x$n_subjects, "\n")
  cat("Features:", x$n_features, "\n")
  
  if (!is.null(x$beta)) {
    cat("\nCoefficients:\n")
    print(summary(x$beta))
  }
  
  invisible(x)
}

#' Summary method for fmri_ttest_fit
#' 
#' @param object An fmri_ttest_fit object
#' @param ... Additional summary arguments
#' @export
summary.fmri_ttest_fit <- function(object, ...) {
  print(object)
  
  if (!is.null(object$q)) {
    cat("\nMultiple comparisons correction applied\n")
    cat("Significant features (FDR < 0.05):\n")
    for (i in seq_len(nrow(object$q))) {
      n_sig <- sum(object$q[i,] < 0.05, na.rm = TRUE)
      cat("  ", rownames(object$q)[i], ":", n_sig, "\n")
    }
  }
  
  invisible(object)
}