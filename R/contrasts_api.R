#' Compute linear model contrast statistics (public API)
#'
#' A polished wrapper around the internal fast contrast engine. Accepts
#' name-based contrast specifications and returns a consistent tibble by default.
#'
#' Inputs mirror a standard GLM contrast computation: provide `B` (p x V),
#' `XtXinv` (p x p), residual degrees of freedom `df`, and either `sigma` or
#' `sigma2` (per-voxel noise). Contrasts can be specified as a single mixed
#' `contrasts` list or separately as `t_contrasts` (vectors) and `f_contrasts`
#' (matrices). When contrast vectors/matrices are named (for vectors) or have
#' column names (for matrices), names are matched to `columns` to determine the
#' appropriate design indices; otherwise lengths must match the number of
#' parameters `p`.
#'
#' @param B Numeric matrix (p x V) of coefficients.
#' @param XtXinv Numeric matrix (p x p), inverse cross-product of the design.
#' @param df Residual degrees of freedom.
#' @param sigma Optional numeric vector (length V) of residual std. dev.; ignored if `sigma2` provided.
#' @param sigma2 Optional numeric vector (length V) of residual variances.
#' @param contrasts Optional named list mixing t- and F-contrasts; vectors are t, matrices are F.
#' @param t_contrasts Optional named list of numeric vectors (t-contrasts).
#' @param f_contrasts Optional named list of numeric matrices (F-contrasts).
#' @param columns Optional character vector (length p) naming coefficients; used for name matching.
#' @param output Either "stacked" (default; tibble) or "list" (raw list of tibbles).
#' @param robust_weights Optional numeric vector of robust weights or NULL. Used only for df adjustments.
#' @param ar_order Integer AR order; used only for effective df adjustments.
#' @param drop_failed Logical; drop contrasts that fail validation (default TRUE).
#'
#' @return A tibble with rows per contrast (default) or a named list if `output = "list"`.
#' @export
compute_lm_contrasts <- function(
  B,
  XtXinv,
  df,
  sigma = NULL,
  sigma2 = NULL,
  contrasts = NULL,
  t_contrasts = NULL,
  f_contrasts = NULL,
  columns = NULL,
  output = c("stacked", "list"),
  robust_weights = NULL,
  ar_order = 0,
  drop_failed = TRUE
) {
  output <- match.arg(output, c("stacked", "list"))

  if (!is.matrix(B)) stop("B must be a numeric matrix (p x V)")
  if (!is.matrix(XtXinv)) stop("XtXinv must be a numeric matrix (p x p)")
  mode(B) <- "numeric"; mode(XtXinv) <- "numeric"
  p <- nrow(B); V <- ncol(B)
  if (nrow(XtXinv) != p || ncol(XtXinv) != p) stop("XtXinv must be p x p")

  # Resolve sigma2
  if (!is.null(sigma2)) {
    s2 <- as.numeric(sigma2)
  } else if (!is.null(sigma)) {
    s2 <- as.numeric(sigma)^2
  } else {
    stop("Provide either `sigma` or `sigma2`")
  }
  if (length(s2) == 1L) s2 <- rep_len(s2, V)
  if (length(s2) != V) stop("Length of sigma/sigma2 must equal ncol(B)")

  # Derive coefficient names for name-based matching
  if (is.null(columns)) {
    columns <- rownames(B)
    if (is.null(columns) || anyNA(columns)) {
      columns <- colnames(XtXinv)
    }
  }
  if (!is.null(columns)) {
    if (length(columns) != p) stop("`columns` must be length p (nrow(B))")
  }

  norm <- .normalize_contrasts_for_engine(
    p = p,
    contrasts = contrasts,
    t_contrasts = t_contrasts,
    f_contrasts = f_contrasts,
    columns = columns,
    drop_failed = drop_failed
  )

  res <- fit_lm_contrasts_fast(
    B = B,
    sigma2 = s2,
    XtXinv = XtXinv,
    conlist = norm$t_list,
    fconlist = norm$f_list,
    df = df,
    robust_weights = robust_weights,
    ar_order = ar_order
  )

  if (output == "list") {
    return(res)
  }
  if (length(res)) dplyr::bind_rows(res) else tibble::tibble()
}


#' Compute contrast statistics from sufficient statistics (public API)
#'
#' Convenience wrapper that accepts design/data sufficient statistics, computes
#' betas and residual variance, and delegates to `compute_lm_contrasts()`.
#'
#' @param XtX Numeric (p x p) cross-product of the design.
#' @param XtS Numeric (p x V) cross-product of design with data.
#' @param StS Numeric length-V vector of sum of squares per voxel.
#' @param df Residual degrees of freedom.
#' @inheritParams compute_lm_contrasts
#'
#' @return A tibble with rows per contrast (default) or a named list if `output = "list"`.
#' @export
compute_lm_contrasts_from_suffstats <- function(
  XtX,
  XtS,
  StS,
  df,
  sigma = NULL,
  sigma2 = NULL,
  contrasts = NULL,
  t_contrasts = NULL,
  f_contrasts = NULL,
  columns = NULL,
  output = c("stacked", "list"),
  robust_weights = NULL,
  ar_order = 0,
  drop_failed = TRUE
) {
  output <- match.arg(output, c("stacked", "list"))

  XtX <- as.matrix(XtX); XtS <- as.matrix(XtS); StS <- as.numeric(StS)
  p <- nrow(XtX); V <- ncol(XtS)
  stopifnot(ncol(XtX) == p, nrow(XtS) == p, length(StS) == V)

  # Robust inversion for XtX
  XtXinv <- tryCatch(chol2inv(chol(XtX)), error = function(e) solve(XtX))
  B <- XtXinv %*% XtS

  # Resolve sigma2 either from inputs or from suffstats
  if (is.null(sigma2) && is.null(sigma)) {
    SSE <- StS - colSums(B * XtS)
    SSE <- pmax(SSE, 0)  # guard against tiny negative due to numerics
    s2 <- pmax(SSE / df, .Machine$double.eps)
  } else if (!is.null(sigma2)) {
    s2 <- as.numeric(sigma2)
  } else {
    s2 <- as.numeric(sigma)^2
  }
  if (length(s2) == 1L) s2 <- rep_len(s2, V)
  if (length(s2) != V) stop("Length of sigma/sigma2 must equal ncol(XtS)")

  # Derive coefficient names for name-based matching (prefer columns if given)
  if (is.null(columns)) {
    columns <- rownames(XtS)
    if (is.null(columns) || anyNA(columns)) {
      columns <- colnames(XtX)
    }
  }
  if (!is.null(columns)) {
    if (length(columns) != p) stop("`columns` must be length p (nrow(XtX))")
  }

  norm <- .normalize_contrasts_for_engine(
    p = p,
    contrasts = contrasts,
    t_contrasts = t_contrasts,
    f_contrasts = f_contrasts,
    columns = columns,
    drop_failed = drop_failed
  )

  res <- fit_lm_contrasts_fast(
    B = B,
    sigma2 = s2,
    XtXinv = XtXinv,
    conlist = norm$t_list,
    fconlist = norm$f_list,
    df = df,
    robust_weights = robust_weights,
    ar_order = ar_order
  )

  if (output == "list") {
    return(res)
  }
  if (length(res)) dplyr::bind_rows(res) else tibble::tibble()
}


# Normalize contrast specifications into the internal lists expected by
# `fit_lm_contrasts_fast`, attaching a `colind` attribute as needed.
#
# @keywords internal
.normalize_contrasts_for_engine <- function(p, contrasts, t_contrasts, f_contrasts,
                                            columns, drop_failed = TRUE) {
  build_t <- function(w, nm) {
    if (is.matrix(w)) {
      if (nrow(w) == 1L) w <- as.vector(w)
      else if (ncol(w) == 1L) w <- as.vector(w)
      else stop(sprintf("t-contrast '%s' must be a vector or 1xP / Px1 matrix", nm))
    }
    w <- as.numeric(w)
    colind <- attr(w, "colind", exact = TRUE)
    if (is.null(colind)) {
      if (!is.null(names(w)) && !is.null(columns)) {
        idx <- match(names(w), columns)
        if (anyNA(idx)) {
          stop(sprintf("Names in t-contrast '%s' not found in `columns`", nm))
        }
        colind <- as.integer(idx)
        ord <- order(colind)
        colind <- colind[ord]; w <- w[ord]
      } else if (length(w) == p) {
        colind <- seq_len(p)
      } else {
        stop(sprintf("t-contrast '%s' must be named to match `columns` or length p", nm))
      }
    }
    attr(w, "colind") <- as.integer(colind)
    w
  }

  build_f <- function(M, nm) {
    if (!is.matrix(M)) {
      M <- matrix(as.numeric(M), nrow = 1L)
    }
    mode(M) <- "numeric"
    colind <- attr(M, "colind", exact = TRUE)
    if (is.null(colind)) {
      if (!is.null(colnames(M)) && !is.null(columns)) {
        idx <- match(colnames(M), columns)
        if (anyNA(idx)) {
          stop(sprintf("Column names in F-contrast '%s' not found in `columns`", nm))
        }
        colind <- as.integer(idx)
        ord <- order(colind)
        M <- M[, ord, drop = FALSE]
        colind <- colind[ord]
      } else if (ncol(M) == p) {
        colind <- seq_len(p)
      } else {
        stop(sprintf("F-contrast '%s' must have colnames matching `columns` or ncol == p", nm))
      }
    }
    attr(M, "colind") <- as.integer(colind)
    M
  }

  t_list <- list(); f_list <- list()

  if (!is.null(contrasts)) {
    if (!is.list(contrasts)) stop("`contrasts` must be a named list")
    if (is.null(names(contrasts)) || any(!nzchar(names(contrasts)))) {
      names(contrasts) <- sprintf("con%02d", seq_along(contrasts))
    }
    for (nm in names(contrasts)) {
      obj <- contrasts[[nm]]
      ok <- try({
        if (is.matrix(obj) && nrow(obj) > 1L) {
          f_list[[nm]] <- build_f(obj, nm)
        } else {
          t_list[[nm]] <- build_t(obj, nm)
        }
        TRUE
      }, silent = TRUE)
      if (inherits(ok, "try-error") && !drop_failed) stop(attr(ok, "condition")$message)
    }
  }

  if (!is.null(t_contrasts)) {
    if (is.null(names(t_contrasts)) || any(!nzchar(names(t_contrasts)))) {
      names(t_contrasts) <- sprintf("t%02d", seq_along(t_contrasts))
    }
    for (nm in names(t_contrasts)) {
      ok <- try({ t_list[[nm]] <- build_t(t_contrasts[[nm]], nm); TRUE }, silent = TRUE)
      if (inherits(ok, "try-error") && !drop_failed) stop(attr(ok, "condition")$message)
    }
  }

  if (!is.null(f_contrasts)) {
    if (is.null(names(f_contrasts)) || any(!nzchar(names(f_contrasts)))) {
      names(f_contrasts) <- sprintf("F%02d", seq_along(f_contrasts))
    }
    for (nm in names(f_contrasts)) {
      ok <- try({ f_list[[nm]] <- build_f(f_contrasts[[nm]], nm); TRUE }, silent = TRUE)
      if (inherits(ok, "try-error") && !drop_failed) stop(attr(ok, "condition")$message)
    }
  }

  list(t_list = Filter(Negate(is.null), t_list),
       f_list = Filter(Negate(is.null), f_list))
}

