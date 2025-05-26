

#' Fast factorial contrast generators
#'
#' Returns a matrix **N_cells × N_contrasts** – *each row is a design cell*,
#' columns are independent contrasts (difference‑coded for the factors you ask
#' for, grand‑mean for the rest).  Suitable for `tcrossprod(dm, C)` or
#' `lm.fit(design, y)` followed by `%*% coef` in the usual way.
#'
#' @param des      data.frame with one column per factor (must be `factor`)
#' @param factors  character vector: which factor(s) get **difference coding**.
#'                 • `generate_main_effect_contrast()` takes a **single**
#'                   factor name.<br>
#'                 • `generate_interaction_contrast()` takes ≥ 2 for an
#'                   interaction (or 1 to reproduce a main‑effect matrix).
#'
#' @return numeric matrix **nrow = ∏ levels(f) , ncol = ∏ (Lᵢ − 1)** for the
#'         chosen factors.
#'
#' @examples
#' des <- expand.grid(Time = factor(1:4),
#'                    Cond = factor(c("face","scene")))
#'
#' # Main effect of Time (4‑1 = 3 contrasts)
#' M <- generate_main_effect_contrast(des, "Time")
#'
#' # Full Time×Cond interaction ( (4‑1)*(2‑1) = 3 contrasts )
#' I <- generate_interaction_contrast(des, c("Time","Cond"))
#' dim(I)   # 8 rows (cells) × 3 columns (contrasts)
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

  blocks   <- Map(build_block, des, fac_names %in% factors)
  C_matrix <- Reduce(kronecker, blocks)

  # Assert rows = design cells
  n_cells <- prod(vapply(des, nlevels, 1L))
  stopifnot(nrow(C_matrix) == n_cells)

  C_matrix
}

#' @param factor Single factor name for the main effect.
#' @rdname generate_interaction_contrast
#' @export
generate_main_effect_contrast <- function(des, factor) {
  if (length(factor) != 1L)
    stop("main‑effect contrast expects exactly one factor name")
  generate_interaction_contrast(des, factor)
}