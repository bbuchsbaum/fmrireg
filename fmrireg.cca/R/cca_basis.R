# Cache for PCA+alpha basis by TR/span/library signature
.cca_basis_cache <- new.env(parent = emptyenv())

.cca_cache_key <- function(TR, span, lib_sig) paste0("TR=", TR, ";span=", span, ";lib=", lib_sig)

.cca_build_library <- function(library_fun, pgrid) {
  # Try to build a basis library using fmrihrf::hrf_library
  if (!is.null(library_fun) && is.function(library_fun) && !is.null(pgrid)) {
    return(fmrihrf::hrf_library(library_fun, pgrid))
  }
  NULL
}

.cca_pca_from_library <- function(TR, span, library_fun = NULL, pgrid = NULL) {
  key <- .cca_cache_key(TR, span, digest::digest(list(formals = if (!is.null(library_fun)) tryCatch(utils::capture.output(print(library_fun)), error = function(e) "fun") else "spmg2", pgrid)))
  if (!is.null(.cca_basis_cache[[key]])) return(.cca_basis_cache[[key]])

  tt <- seq(0, span, by = TR)
  lib <- .cca_build_library(library_fun, pgrid)
  if (!is.null(lib)) {
    L <- fmrihrf::evaluate(lib, tt)
    if (is.null(dim(L)) || nrow(L) != length(tt)) {
      # Unexpected shape; fall back
      L <- NULL
    }
  } else {
    L <- NULL
  }

  if (is.null(L)) {
    # Fallback: use SPMG2 as a lightweight proxy library
    sp <- fmrihrf::evaluate(fmrihrf::HRF_SPMG2, tt)
    if (is.null(dim(sp)) || ncol(sp) < 2L) stop("SPMG2 evaluation did not return two columns")
    L <- sp
  }

  # Mean shape and first PC of residuals
  y1 <- rowMeans(L)
  Lc <- L - matrix(y1, nrow(L), ncol(L))
  sv <- svd(Lc)
  y2 <- sv$u[, 1]
  # Normalize and orient PC1 to be mostly positive
  y1 <- y1 / sqrt(sum(y1^2) + 1e-12)
  sgn <- if (sum(y2) < 0) -1 else 1
  y2 <- sgn * (y2 / sqrt(sum(y2^2) + 1e-12))

  out <- list(y1 = y1, y2 = y2)
  .cca_basis_cache[[key]] <- out
  out
}

#' Two-column temporal basis for CCA engine (alpha-mix)
#'
#' Constructs two temporal basis functions via an alpha-mix of the
#' SPMG2 canonical HRF and its temporal derivative, mirroring the
#' `y1 ± alpha*y2` construction used in the plan. This serves as a
#' robust stand-in for a PCA+alpha basis while remaining lightweight
#' and TR-agnostic.
#'
#' @param alpha Numeric scalar; if NULL, auto-alpha is used (default 0.3 heuristic).
#' @param span Numeric. Support span in seconds for evaluation (default 30).
#' @param TR Numeric. Sampling interval (seconds). If NULL, defaults to 1.
#' @param library_fun Optional function for fmrihrf::hrf_library. If NULL, falls back to SPMG2 proxy.
#' @param pgrid Optional parameter grid (list/data.frame) for fmrihrf::hrf_library.
#' @return An HRF basis object understood by fmrihrf.
cca_basis <- function(alpha = NULL, span = 30, TR = 1,
                      library_fun = NULL, pgrid = NULL) {
  stopifnot(is.numeric(span), length(span) == 1L, is.finite(span), span > 0)
  if (is.null(TR) || !is.finite(TR) || TR <= 0) TR <- 1

  pca <- .cca_pca_from_library(TR = TR, span = span, library_fun = library_fun, pgrid = pgrid)
  a <- if (is.null(alpha)) 0.3 else as.numeric(alpha)
  if (!is.finite(a) || a < 0) a <- 0.3

  y1t <- pca$y1 + a * pca$y2
  y2t <- pca$y1 - a * pca$y2

  # Build evaluators from the sampled grid using linear interpolation
  tt <- seq(0, span, by = TR)
  f1 <- stats::approxfun(tt, y1t, yleft = 0, yright = 0)
  f2 <- stats::approxfun(tt, y2t, yleft = 0, yright = 0)

  h1 <- fmrihrf::as_hrf(f1, name = sprintf("cca2_y1_a%.2f", a), nbasis = 1L, span = span)
  h2 <- fmrihrf::as_hrf(f2, name = sprintf("cca2_y2_a%.2f", a), nbasis = 1L, span = span)
  fmrihrf::bind_basis(h1, h2)
}
