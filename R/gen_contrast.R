

#' Fast factorial contrast generators
#'
#' Returns a matrix **N_cells x N_contrasts** - *row r is the design cell in
#' row r of `des`* (for `expand.grid()` input, the first factor varies fastest,
#' which is also the column order of fmridesign crossed terms such as
#' `hrf(A, B)`),
#' columns are independent contrasts (difference-coded for the factors you ask
#' for, grand-mean for the rest).  Suitable for `tcrossprod(dm, C)` or
#' `lm.fit(design, y)` followed by `%*% coef` in the usual way.
#'
#' @param des      data.frame with one column per factor (must be `factor`)
#' @param factors  character vector: which factor(s) get **difference coding**.
#'                 - `generate_main_effect_contrast()` takes a **single**
#'                   factor name.<br>
#'                 - `generate_interaction_contrast()` takes >= 2 for an
#'                   interaction (or 1 to reproduce a main-effect matrix).
#'
#' @return numeric matrix **nrow = prod levels(f) , ncol = prod (Li - 1)** for the
#'         chosen factors.
#'
#' @examples
#' des <- expand.grid(Time = factor(1:4),
#'                    Cond = factor(c("face","scene")))
#'
#' # Main effect of Time (4-1 = 3 contrasts)
#' M <- generate_main_effect_contrast(des, "Time")
#'
#' # Full TimexCond interaction ( (4-1)*(2-1) = 3 contrasts )
#' I <- generate_interaction_contrast(des, c("Time","Cond"))
#' dim(I)   # 8 rows (cells) x 3 columns (contrasts)
#' @export
generate_interaction_contrast <- function(des, factors) {

  stopifnot(all(factors %in% names(des)))
  fac_names <- names(des)

  build_block <- function(f, diff_needed) {
    L <- nlevels(f)
    if (diff_needed)             # (L x (L-1)) difference coding
      -t(diff(diag(L)))          # rows = levels, cols = contrasts
    else
      matrix(1, nrow = L, ncol = 1)
  }

  if (!all(vapply(des, is.factor, logical(1)))) {
    stop("every column of 'des' must be a factor", call. = FALSE)
  }
  n_cells <- prod(vapply(des, nlevels, 1L))
  if (nrow(des) != n_cells || anyDuplicated(des)) {
    stop("'des' must list every combination of factor levels exactly once ",
         "(e.g. expand.grid() of the levels)", call. = FALSE)
  }

  blocks <- Map(build_block, des, fac_names %in% factors)

  # Row r of the result is the contrast weight of the cell in row r of `des`:
  # a row-wise Kronecker product of the per-factor blocks, indexed by each
  # factor's level in that row. (A plain Reduce(kronecker, blocks) produces
  # rows with the FIRST factor varying slowest, the opposite of expand.grid()
  # and of fmridesign's crossed-term columns.) Columns are ordered with the
  # first factor's contrast index varying fastest.
  C_matrix <- matrix(1, nrow = nrow(des), ncol = 1)
  for (k in seq_along(blocks)) {
    Rk <- blocks[[k]][as.integer(des[[k]]), , drop = FALSE]
    C_matrix <- do.call(cbind, lapply(seq_len(ncol(Rk)), function(j) C_matrix * Rk[, j]))
  }

  C_matrix
}

#' @param factor Single factor name for the main effect.
#' @rdname generate_interaction_contrast
#' @export
generate_main_effect_contrast <- function(des, factor) {
  if (length(factor) != 1L)
    stop("main-effect contrast expects exactly one factor name")
  generate_interaction_contrast(des, factor)
}