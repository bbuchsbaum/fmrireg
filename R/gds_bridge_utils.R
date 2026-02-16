.gds_safe_assay <- function(x, name) {
  if (is.null(x)) return(NULL)
  tryCatch(fmrigds::assay(x, name), error = function(e) NULL)
}

.gds_safe_assay_names <- function(x) {
  if (is.null(x)) return(character(0))
  tryCatch(fmrigds::assay_names(x), error = function(e) character(0))
}

.gds_coerce_matrix <- function(a) {
  if (is.null(a)) return(NULL)
  if (is.matrix(a)) return(a)
  if (is.data.frame(a)) return(as.matrix(a))
  if (is.array(a)) {
    d <- dim(a)
    if (is.null(d)) {
      vals <- as.vector(a)
      return(matrix(vals, nrow = length(vals), ncol = 1))
    }
    if (length(d) == 1L) {
      vals <- as.vector(a)
      return(matrix(vals, nrow = d[1], ncol = 1))
    }
    vals <- as.vector(a)
    return(matrix(vals, nrow = d[1], ncol = prod(d[-1])))
  }
  if (is.vector(a)) {
    vals <- as.vector(a)
    return(matrix(vals, nrow = length(vals), ncol = 1))
  }
  a
}

.gds_terms_by_features <- function(a) {
  m <- .gds_coerce_matrix(a)
  if (is.null(m)) return(NULL)
  t(m)
}

.gds_features_by_terms <- function(a) {
  m <- .gds_coerce_matrix(a)
  if (is.null(m)) return(NULL)
  m
}

.gds_safe_model_matrix <- function(gd, formula, fallback_n = NULL) {
  mm <- tryCatch(fmrigds::model_matrix(gd, formula), error = function(e) NULL)
  if (!is.null(mm)) return(mm)
  if (is.null(fallback_n)) return(NULL)
  stats_terms <- tryCatch(stats::terms(formula), error = function(e) NULL)
  intercept <- !is.null(stats_terms) &&
    length(attr(stats_terms, "term.labels")) == 0 &&
    isTRUE(attr(stats_terms, "intercept") == 1)
  if (!intercept) return(NULL)
  matrix(1, nrow = fallback_n, ncol = 1, dimnames = list(NULL, "(Intercept)"))
}
