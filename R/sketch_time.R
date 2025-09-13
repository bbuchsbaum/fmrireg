#' Build a time sketch matrix S (m x T)
#'
#' @param Tlen Integer time length
#' @param ctrl List(method = "gaussian"|"countsketch", m, iters)
#' @return Dense or sparse sketch matrix S
#' @keywords internal
make_time_sketch <- function(Tlen, ctrl) {
  stopifnot(is.list(ctrl))
  method <- match.arg(ctrl$method %||% "gaussian", c("gaussian", "countsketch", "srht", "ihs"))
  m <- ctrl$m
  stopifnot(!is.null(m), m > 0L, m <= Tlen)
  if (method == "gaussian") {
    S <- matrix(stats::rnorm(m * Tlen), nrow = m, ncol = Tlen) / sqrt(m)
  } else if (method == "countsketch") {
    rows <- sample.int(m, Tlen, replace = TRUE)
    signs <- sample(c(-1, 1), Tlen, replace = TRUE)
    S <- Matrix::sparseMatrix(i = rows, j = seq_len(Tlen), x = signs, dims = c(m, Tlen))
  } else {
    # SRHT / IHS do not materialize S; callers should use srht_apply()/ihs_latent_solve()
    S <- NULL
  }
  S
}

# --- SRHT helpers ---

#' Build SRHT plan
#' @keywords internal
make_srht_plan <- function(Tlen, m) {
  stopifnot(m > 0L, m <= Tlen)
  list(T = Tlen, m = m,
       signs = sample(c(-1, 1), Tlen, replace = TRUE),
       perm  = sample.int(Tlen) - 1L,   # zero-based for C++
       rows  = sort(sample.int(Tlen, m, replace = FALSE)) - 1L,
       scale = sqrt(Tlen / m))
}

#' Apply SRHT to a T x k matrix using a plan
#' @keywords internal
srht_apply <- function(M, plan) {
  .Call('_fmrireg_cpp_srht_apply', PACKAGE = 'fmrireg', M,
        as.integer(plan$rows),
        as.numeric(plan$signs),
        as.integer(plan$perm),
        as.numeric(plan$scale))
}

#' Iterative Hessian Sketch (multi-RHS), returning M and Ginv
#' @keywords internal
ihs_latent_solve <- function(X, Z, m, iters = 3L) {
  .Call('_fmrireg_cpp_ihs_latent', PACKAGE = 'fmrireg', X, Z, as.integer(m), as.integer(iters))
}
