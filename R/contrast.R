#' Translate legacy contrast regex patterns
#'
#' Convert older column-naming patterns to the current naming scheme.
#'
#' @param pattern Character string with the legacy regex.
#' @return Updated regex string.
#' @keywords internal
#' @export
#' @name translate_legacy_pattern
#' @rdname translate_legacy_pattern

translate_legacy_pattern <- function(pattern) {
  # 1. Replace Var[Level] -> Var.Level (Do this first)
  # Handles VarName[LevelName] -> VarName.LevelName
  pattern <- gsub("([A-Za-z0-9_\\.]+)\\[([^]]+)\\]", "\\1.\\2", pattern, perl = TRUE)

  # 2. Replace :basis[digits] -> _b<digits>
  # Handles :basis[3] -> _b3
  #pattern <- gsub(":basis\\[(\\d+)\\]", "_b\\1", pattern, perl = TRUE) # CORRECTED LINE
  pattern <- gsub(":basis\\[(\\d+)\\](\\$?)$", "_b\\1\\2", pattern, perl = TRUE)
  # 3. Replace standalone : -> _ (interaction separator)
  # Uses lookarounds to avoid replacing potential future :: syntax
  pattern <- gsub("(?<!:):(?!:)", "_", pattern, perl = TRUE)
  
  pattern
}

#' Get condition names for a term
#' Wraps conditions() with standard arguments for internal use.
#' @param term An event_term object.
#' @param expanded Logical, whether to return basis-expanded names.
#' @return Character vector of condition names.
#' @keywords internal
#' @noRd
.condnames <- function(term, expanded = TRUE) {
  # Assumes conditions() is robust and handles drop.empty=FALSE internally
  tryCatch(conditions(term, drop.empty = FALSE, expand_basis = expanded),
           error = function(e) {
             stop(paste("Error retrieving condition names for term '", 
                        term$varname %||% "<unknown>", # Use varname if available
                        "' (expanded=", expanded, "): ", e$message), call. = FALSE)
           })
}

#' Calculate contrast weights from logical masks
#' Returns a *named* numeric vector.
#' @param names Character vector of all condition names.
#' @param A_mask Logical vector (same length as names) indicating TRUE for conditions in group A.
#' @param B_mask Logical vector (same length as names) indicating TRUE for conditions in group B. Optional.
#' @param tol Tolerance for sum-to-zero check.
#' @return Named numeric vector of weights.
#' @keywords internal
#' @noRd
.mask_to_weights <- function(names, A_mask, B_mask = NULL, tol = 1e-8) {
  # Input validation
  stopifnot(is.character(names)) # names can be empty if upstream condnames is empty, but length checks handle this.
  stopifnot(is.logical(A_mask), length(A_mask) == length(names))
  
  nA <- sum(A_mask)
  nB <- 0 # Initialize nB for the case where B_mask is NULL

  if (!is.null(B_mask)) {
    stopifnot(is.logical(B_mask), length(B_mask) == length(names))
    if (any(A_mask & B_mask)) {
      stop(".mask_to_weights: Masks for group A and group B overlap.", call. = FALSE)
    }
    nB <- sum(B_mask)
  }

  # Check for completely empty selection: this should be an error.
  if (nA == 0 && (is.null(B_mask) || nB == 0)) {
    stop("Cannot calculate contrast weights: No conditions were selected by the provided mask(s). This usually indicates that the patterns or formulas used to define the contrast did not match any existing conditions. Please check your contrast specification.", call. = FALSE)
  }

  # Initialize weights vector
  w <- numeric(length(names))
  if (length(names) > 0) { # Only assign names if 'names' is not empty
      names(w) <- names
  }
  
  # Assign weights if masks are not empty
  if (nA > 0) {
    w[A_mask] <- 1 / nA
  }
  # Ensure B_mask is not NULL before checking nB for assignment
  if (!is.null(B_mask) && nB > 0) {
    w[B_mask] <- -1 / nB
  }
  
  # Warnings for partially defined A-vs-B contrasts (that won't sum to zero as expected)
  if (!is.null(B_mask)) { # These warnings only make sense for A-vs-B type contrasts
    if (nA == 0 && nB > 0) { # A empty, B not (and B_mask was provided)
      warning(".mask_to_weights: For A-vs-B contrast, Mask A is empty but Mask B is not. Weights will not sum to zero as expected for a balanced comparison.", call. = FALSE)
    } else if (nB == 0 && nA > 0) { # B empty, A not (and B_mask was provided)
      warning(".mask_to_weights: For A-vs-B contrast, Mask B is empty but Mask A is not. Weights will not sum to zero as expected for a balanced comparison.", call. = FALSE)
    } else if (nA > 0 && nB > 0) { # Both A and B are defined and non-empty for A-vs-B
        if (abs(sum(w)) > tol) {
           # This scenario (both non-empty, but sum != 0) should be rare with 1/nA and -1/nB
           warning(".mask_to_weights: Weights for A-vs-B contrast (both groups non-empty) do not sum to zero (Sum: ", sum(w), "). This is unexpected.", call. = FALSE)
        }
    }
  } 
  # If B_mask is NULL (single group contrast, e.g. from unit_contrast or pattern_A only column_contrast):
  # - If nA == 0, we would have errored out above.
  # - If nA > 0, weights are 1/nA for A_mask elements. Sum-to-zero is not expected here. No warning needed.
  
  w
}

#' Robust weight calculation using logical masks
#' Returns a *named* numeric vector (names are condition names)
#' @keywords internal
# makeWeights <- function(keepA, keepB=NULL) {
#   weights <- matrix(0, length(keepA), 1)
#   numA <- sum(keepA)
#   weights[keepA,1] <- rep(1/numA, numA)
#   
#   if (!is.null(keepB)) {
#     numB <- sum(keepB)
#     weights[keepB,1] <- -rep(1/numB, numB)
#   }
#   
#   weights
# }

#' Robust weight calculation using logical masks
#' Returns a *named* numeric vector (names are condition names)
#' @keywords internal
#' @noRd
.make_weights <- function(cond_names, mask_A, mask_B = NULL) {
  # Check inputs
  stopifnot(length(cond_names) == length(mask_A))
  stopifnot(is.logical(mask_A))
  if (!is.null(mask_B)) {
      stopifnot(length(cond_names) == length(mask_B))
      stopifnot(is.logical(mask_B))
      # Ensure no overlap between A and B masks
      if (any(mask_A & mask_B)) {
          stop(".make_weights: Masks A and B overlap.", call. = FALSE)
      }
  }
  
  w <- numeric(length(mask_A))
  names(w) <- cond_names
  
  numA <- sum(mask_A)
  if (numA > 0) {
      w[mask_A] <-  1 / numA
  }
  
  if (!is.null(mask_B)) {
    numB <- sum(mask_B)
    if (numB > 0) {
        w[mask_B] <- -1 / numB
    } 
  }
  
  # Handle cases where one mask is empty (weights should still sum to 0)
  if (numA == 0 && !is.null(mask_B) && numB > 0) {
      warning(".make_weights: Mask A is empty but Mask B is not. Weights will not sum to zero.", call. = FALSE)
  } else if (numB == 0 && !is.null(mask_B) && numA > 0) {
       warning(".make_weights: Mask B is empty but Mask A is not. Weights will not sum to zero.", call. = FALSE)
  } else if (numA == 0 && (is.null(mask_B) || numB == 0)) {
       warning(".make_weights: Both masks A and B are empty.", call. = FALSE)
  }
  
  w
}


#' Get Indices of Design Matrix Columns Matching a Pattern
#'
#' Retrieves the full list of potential design matrix column names for an
#' event_term (using `conditions(term, drop.empty = FALSE)`) and returns the
#' numeric indices of columns whose names match the provided regex pattern.
#' Handles translation of legacy patterns for backward compatibility.
#'
#' @param term An object of class `event_term`.
#' @param pattern A character string containing a regular expression (can use legacy syntax).
#' @param ... Additional arguments passed to `grep()`.
#' @return An integer vector of the indices of matching column names.
#'         Returns `integer(0)` if no names match or if the term has no conditions.
#' @keywords internal
#' @noRd
.col_index <- function(term, pattern, ...) {
  
  # --- Removed legacy pattern translation block ---
  
  # Ensure conditions function is available (Should be implicitly available if package loaded)
  # if (!exists("conditions") || !is.function(conditions)) {
  #     stop(".col_index requires the 'conditions' function.")
  # }
  
  # --- Removed duplicated basis expansion logic --- 
  # The caller should decide if expanded names are needed and pass them or call conditions appropriately.
  # For now, assume the core use case (column_contrast) needs expanded names.
  # TODO: Refactor later to use a condnames() helper or accept names directly.
  expand_basis <- TRUE # Default assumption for now, needs refinement
  hrfspec <- attr(term, "hrfspec")
  if (!is.null(hrfspec) && !is.null(hrfspec$hrf)) {
    hrf_fun <- hrfspec$hrf
    # Check nbasis generic exists and is applicable
    # We still need to know if *any* expansion is possible for the default
    if (exists("nbasis") && is.function(nbasis) && !inherits(try(nbasis(hrf_fun), silent=TRUE), "try-error") && nbasis(hrf_fun) > 1) { 
      # Expansion is possible, keep expand_basis = TRUE
    } else {
      expand_basis <- FALSE # No basis functions or only 1
    }
  } else {
     expand_basis <- FALSE # No HRF spec
  }

  # Get all potential column names (expanded or not)
  # Assume conditions() is robust and will stop() on genuine error
  all_colnames <- conditions(term, drop.empty = FALSE, expand_basis = expand_basis)

  # --- Removed try() wrapper and error check --- 
  # if (inherits(all_colnames, "try-error") || length(all_colnames) == 0) {
  #   # Stop if conditions() failed, return integer(0) if empty result (grep handles this)
  #   stop(paste("Could not retrieve column names via conditions() for term",
  #                 term$varname, "in .col_index.")) 
  # }

  # Perform grep using the pattern directly
  grep(pattern, all_colnames, value = FALSE, ...)
}

#' Contrast Specification
#'
#' @description
#' Define a linear contrast using a formula expression.
#'
#' @param form A formula describing the contrast.
#' @param name A character label for the contrast.
#' @param where An expression defining the subset over which the contrast is applied (default: NULL).
#'
#' @return A list containing the contrast specification.
#'
#' @examples
#' # A minus B contrast
#' contrast(~ A - B, name="A_B")
#'
#' @export
contrast <- function(form, name, where=NULL) {
  assert_that(rlang::is_formula(form))
  if (!is.null(where)) {
    assert_that(rlang::is_formula(where))
  }
  ret <- list(A=form,
              B=NULL,
              where=where,
              name=name)
  
  class(ret) <- c("contrast_formula_spec", "contrast_spec", "list")
  ret
  
}

#' Unit Contrast
#'
#' @description
#' Construct a contrast that sums to 1 and is used to define contrasts against the baseline.
#'
#' @param A A formula representing the contrast expression.
#' @param name A character string specifying the name of the contrast.
#' @param where An optional formula specifying the subset of conditions to apply the contrast to.
#'
#' @return A unit_contrast_spec object containing the contrast that sums to 1.
#'
#' @examples
#' con <- unit_contrast(~ Face, name="Main_face")
#'
#' @export
unit_contrast <- function(A, name, where=NULL) {
  assert_that(rlang::is_formula(A)) 
  
  if (!is.null(where)) {
    assert_that(rlang::is_formula(where))
  }
  
  structure(
    list(A=A,
         B=NULL,
         where=where,
         name=name),
    class=c("unit_contrast_spec", "contrast_spec", "list")
  )
  
}

#' One Against All Contrast
#'
#' @description
#' Construct contrasts comparing each factor level against the average of the other levels.
#'
#' @param levels A vector of factor levels to be compared.
#' @param facname A character string specifying the name of the factor containing the supplied levels.
#' @param where An optional formula specifying the subset over which the contrast is computed.
#'
#' @return A contrast_set object containing contrasts comparing each factor level against the average of the other levels.
#'
#' @examples
#' fac <- factor(rep(c("A", "B", "C"), 2))
#' con <- one_against_all_contrast(levels(fac), "fac")
#'
#' @export
one_against_all_contrast <- function(levels, facname, where=NULL) {
  if (!is.null(where)) {
    assert_that(rlang::is_formula(where))
  }
  
  ret <- lapply(1:length(levels), function(i) {
    lev1 <- levels[i]
    levother <- levels[-i]
    pair_contrast(as.formula(paste("~", facname, " == ", paste0('"', lev1, '"'))), 
                  as.formula(paste0("~", facname, "!= ", paste0('"', lev1, '"'))), 
                  where=where, name=paste0("con_", lev1, "_vs_", "other"))
  })
  
  do.call(contrast_set, ret)
  
}


#' Create a Set of Contrasts
#'
#' @description
#' Construct a list of contrast_spec objects.
#'
#' @param ... A variable-length list of contrast_spec objects.
#'
#' @return A list of contrast_spec objects with class "contrast_set".
#'
#' @examples
#' c1 <- contrast(~ A - B, name="A_B")
#' c2 <- contrast(~ B - C, name="B_C")
#' contrast_set(c1,c2)
#'
#' @export
#' @import assertthat
#' @importFrom purrr map_lgl
contrast_set <- function(...) {
  ret <- list(...)
  assertthat::assert_that(all(purrr::map_lgl(ret, inherits, "contrast_spec")))
  class(ret) <- c("contrast_set", "list")
  ret
}


#' Pairwise Contrasts
#'
#' @description
#' Construct pairwise contrasts for all combinations of factor levels.
#'
#' @param levels A vector of factor levels to be compared.
#' @param facname The name of the factor variable (column name in the design) these levels belong to.
#' @param where An optional formula specifying the subset over which the contrast is computed.
#' @param name_prefix A character string to prefix the generated contrast names (default: "con").
#'
#' @return A contrast_set object containing pairwise contrasts for all combinations of factor levels.
#'
#' @examples
#' # Assuming 'my_factor' is a column name
#' pairwise_contrasts(c("A", "B", "C"), facname = "my_factor")
#' pairwise_contrasts(c("A", "B", "C"), facname = "my_factor", name_prefix = "pair")
#'
#' @export
#' @importFrom utils combn
pairwise_contrasts <- function(levels, facname, where=NULL, name_prefix = "con") {
  assert_that(is.character(facname), length(facname) == 1, msg = "'facname' must be a single string.")
  if (!is.null(where)) {
    assert_that(rlang::is_formula(where))
  }
  
  if (length(levels) < 2) {
    stop("pairwise_contrasts requires at least two levels.")
  }
  
  cbns <- combn(length(levels), 2)
  ret <- lapply(1:ncol(cbns), function(i) {
    lev1 <- levels[cbns[1,i]]
    lev2 <- levels[cbns[2,i]]
    # Construct formulas using the factor name
    formula_A <- as.formula(paste("~", facname, "==", paste0('"', lev1, '"')))
    formula_B <- as.formula(paste("~", facname, "==", paste0('"', lev2, '"')))
    pair_contrast(formula_A, formula_B, 
                  where=where, name=paste0(name_prefix, "_", lev1, "_", lev2))
  })
  
  do.call(contrast_set, ret)
}




#' Pair Contrast
#'
#' @description
#' Construct a sum-to-zero contrast between two logical expressions. This function is
#' particularly useful for comparing specific conditions or combinations of conditions.
#'
#' @param A A formula representing the first logical expression in the contrast.
#' @param B A formula representing the second logical expression in the contrast.
#' @param name A character string specifying the name of the contrast (mandatory).
#' @param where An optional formula specifying the subset over which the contrast is computed.
#'
#' @return A pair_contrast_spec object containing:
#'   \item{A}{First logical expression}
#'   \item{B}{Second logical expression}
#'   \item{where}{Subsetting formula (if provided)}
#'   \item{name}{Contrast name}
#'
#' @details
#' The contrast is constructed as (A - B), where A and B are logical expressions that
#' evaluate to TRUE/FALSE for each observation. The resulting contrast weights sum to zero.
#'
#' @examples
#' # Compare faces vs scenes
#' pair_contrast(~ category == "face", ~ category == "scene", name = "face_vs_scene")
#'
#' # Compare with subsetting
#' pair_contrast(~ category == "face", ~ category == "scene",
#'              name = "face_vs_scene_block1",
#'              where = ~ block == 1)
#'
#' @seealso
#' \code{\link{pairwise_contrasts}} for all pairwise comparisons,
#' \code{\link{contrast_set}} for creating sets of contrasts
#'
#' @export
pair_contrast <- function(A, B, name, where = NULL) {
  assert_that(rlang::is_formula(A))
  assert_that(rlang::is_formula(B))
  

  if (!is.null(where)) {
    assert_that(rlang::is_formula(where))
  }
  
  ret <- list(A=A,
              B=B,
              where=where,
              name=name)
    
  class(ret) <- c("pair_contrast_spec", "contrast_spec", "list")
  ret
}


#' One-way Contrast
#'
#' @description
#' Create a one-way contrast specification
#'
#' @param A A formula specifying the contrast
#' @param name The name of the contrast
#' @param where Optional environment for evaluating the formula
#' @return A oneway_contrast_spec object that can be used to generate contrast weights
#' @examples
#' # Create a one-way contrast for a factor 'basis'
#' con <- oneway_contrast(~ basis, name = "Main_basis")
#'
#' # Create a one-way contrast with a specific environment
#' con <- oneway_contrast(~ basis, name = "Main_basis",
#'                       where = new.env())
#'
#' @seealso \code{\link{interaction_contrast}} for testing interactions,
#'          \code{\link{pair_contrast}} for pairwise comparisons
#' @export
oneway_contrast <- function(A, name, where = NULL) {
  assert_that(rlang::is_formula(A)) 
  
  if (!is.null(where)) {
    assert_that(rlang::is_formula(where))
  }
  
  structure(
    list(A=A,
         B=NULL,
         where=where,
         name=name),
    class=c("oneway_contrast_spec", "contrast_spec", "list")
  )
}

#' Interaction Contrast
#'
#' @description
#' Create an interaction contrast specification
#'
#' @param A A formula specifying the interaction contrast
#' @param name The name of the contrast
#' @param where Optional environment for evaluating the formula
#' @return An interaction_contrast_spec object containing the specification for
#'         generating interaction contrast weights
#' @examples
#' # Create an interaction contrast for factors A and B
#' con <- interaction_contrast(~ A * B, name = "A_by_B")
#'
#' # Create an interaction contrast with a specific environment
#' con <- interaction_contrast(~ A * B, name = "A_by_B",
#'                           where = new.env())
#'
#' @seealso \code{\link{oneway_contrast}} for main effects,
#'          \code{\link{pair_contrast}} for pairwise comparisons
#' @export
interaction_contrast <- function(A, name, where = NULL) {
  assert_that(rlang::is_formula(A)) 
  
  if (!is.null(where)) {
    assert_that(rlang::is_formula(where))
  }
  
  
  
  structure(
    list(A=A,
         B=NULL,
         where=where,
         name=name),
    class=c("interaction_contrast_spec", "contrast_spec", "list")
  )
}

#' Column Contrast Specification
#'
#' @description
#' Define a contrast by directly targeting design matrix columns using regex patterns.
#' This is useful for contrasts involving continuous variables or specific basis functions.
#'
#' @param pattern_A A character string containing a regex pattern to identify the
#'   columns for the positive (+) part of the contrast.
#' @param pattern_B Optional character string containing a regex pattern for the
#'   negative (-) part (for A-B type contrasts). If NULL, creates a contrast testing
#'   the average of columns matching `pattern_A` against baseline (0).
#' @param name A character string name for the contrast (mandatory).
#' @param where Currently unused for column_contrast, but kept for API consistency.
#'
#' @return A `column_contrast_spec` object containing the specification.
#'
#' @details
#' This contrast type operates by finding design matrix columns whose names match
#' the provided patterns (`pattern_A`, `pattern_B`). It calculates weights such that
#' the average effect of the 'A' columns is compared to the average effect of the
#' 'B' columns (or baseline if `pattern_B` is NULL). Weights are assigned as +1/nA
#' for 'A' columns and -1/nB for 'B' columns, ensuring the contrast sums to zero
#' if both A and B groups are present.
#'
#' Use standard R regex syntax for the patterns. Remember to escape special
#' characters (e.g., `\\[`, `\\.`, `\\*`).
#'
#' @examples
#' # Test the main effect of a continuous modulator 'RT'
#' # Assumes RT is a column name, e.g., from columns(Scale(RT))
#' cc1 <- column_contrast(pattern_A = "^z_RT$", name = "Main_RT")
#'
#' # Compare Condition.A vs Condition.B for the 'RT' modulator effect
#' # Assumes condition names like "Condition.A_z_RT", "Condition.B_z_RT"
#' cc2 <- column_contrast(pattern_A = "^Condition\\.A_z_RT$",
#'                        pattern_B = "^Condition\\.B_z_RT$",
#'                        name = "CondA_vs_CondB_for_RT")
#'
#' # Test a specific basis function (e.g., basis spline #3)
#' # Assumes column names like "TermName_Condition.Tag_b03"
#' cc3 <- column_contrast(pattern_A = "_b03$", name = "Basis_3_Effect")
#'
#' @export
column_contrast <- function(pattern_A, pattern_B = NULL, name, where = NULL) {
  assert_that(is.character(pattern_A), length(pattern_A) == 1)
  if (!is.null(pattern_B)) {
    assert_that(is.character(pattern_B), length(pattern_B) == 1)
  }
  assert_that(is.character(name), length(name) == 1)
  if (!is.null(where)) {
      warning("'where' argument is currently ignored for column_contrast.")
      # assert_that(rlang::is_formula(where)) # Keep structure if needed later
  }

  ret <- list(
    pattern_A = pattern_A,
    pattern_B = pattern_B,
    where = where, # Store it even if unused for now
    name = name
  )

  class(ret) <- c("column_contrast_spec", "contrast_spec", "list")
  ret
}

#' Polynomial Contrast
#'
#' @description
#' Create polynomial contrasts for testing trends across ordered factor levels. This is
#' particularly useful for analyzing factors with a natural ordering (e.g., time, dose).
#'
#' @param A A formula specifying the ordered factor.
#' @param name A character string identifying the contrast.
#' @param where An optional formula for subsetting the data.
#' @param degree An integer specifying the degree of the polynomial (default: 1).
#' @param value_map An optional list mapping factor levels to numeric values.
#'
#' @return A poly_contrast_spec object containing the specification for generating
#'   polynomial contrast weights.
#'
#' @details
#' The function creates orthogonal polynomial contrasts up to the specified degree.
#' These contrasts can test for linear, quadratic, cubic, and higher-order trends
#' in the data. The value_map parameter allows for non-uniform spacing between levels.
#'
#' @examples
#' # Linear trend across time points
#' pcon <- poly_contrast(~ time, name = "linear_time", degree = 1)
#'
#' # Cubic trend with custom spacing
#' pcon <- poly_contrast(~ dose, name = "dose_cubic",
#'                      degree = 3,
#'                      value_map = list("low" = 0, "med" = 2, "high" = 5))
#'
#' @seealso
#' \code{\link{oneway_contrast}} for categorical contrasts,
#' \code{\link{interaction_contrast}} for interaction effects
#'
#' @export
poly_contrast <- function(A, name, where = NULL, degree = 1, value_map = NULL) {
  assert_that(rlang::is_formula(A))
  
  if (!is.null(where)) {
    assert_that(rlang::is_formula(where))
  }

  ret <- list(
    A=A,
    B=NULL,
    where=where,
    degree=degree,
    value_map=value_map,
    name=name)
  
  class(ret) <- c("poly_contrast_spec", "contrast_spec", "list")
  ret
}

#' Unit Contrast Weights
#'
#' @description
#' Compute the contrast weights for a unit_contrast_spec object.
#'
#' @param x A unit_contrast_spec object.
#' @param term A term object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the term, name, weights, condition names, and contrast specification.
#'
#' @export
contrast_weights.unit_contrast_spec <- function(x, term,...) {
  # Get cells (categorical only)
  term_cells <- cells(term)
  if (nrow(term_cells) == 0) {
       warning(paste("Contrast '", x$name, "': Term '", term$varname, "' has no categorical cells."), call. = FALSE)
       weights_out <- matrix(numeric(0), nrow = 0, ncol = 1)
       colnames(weights_out) <- x$name
       cell_names_out <- character(0)
  } else {
      # Apply 'where' clause
      keep <- if (!is.null(x$where)) {
        tryCatch(rlang::eval_tidy(rlang::f_rhs(x$where), data = term_cells), error = function(e) {
            warning(paste("Contrast '", x$name, "': Error evaluating 'where' clause: ", e$message), call. = FALSE)
            rep(FALSE, nrow(term_cells))
        })
      } else {
        rep(TRUE, nrow(term_cells))
      }
      
      relevant_cells <- term_cells[keep, , drop = FALSE]
      
      if (nrow(relevant_cells) == 0) {
          warning(paste("Contrast '", x$name, "' resulted in no relevant cells after applying the 'where' clause."), call. = FALSE)
          weights_out <- matrix(numeric(0), nrow = 0, ncol = 1)
          colnames(weights_out) <- x$name
          cell_names_out <- character(0)
      } else {
          # For unit_contrast, the formula A usually defines the set of cells to be averaged.
          # All cells within relevant_cells that are targeted by the formula should get a weight.
          # If x$A is ~MyFactor, all levels of MyFactor in relevant_cells are included.
          # Thus, keepA_rel should be TRUE for all rows in relevant_cells. 
          # If x$A is a logical expression, eval_tidy would handle it, but for unit_contrast 
          # the simple case implies averaging all specified cells.
          
          # The original eval_tidy might not be appropriate if x$A is just a variable name.
          # keepA_rel <- tryCatch(rlang::eval_tidy(rlang::f_rhs(x$A), data = relevant_cells), error = function(e) {
          #     warning(paste("Contrast '", x$name, "': Error evaluating formula A: ", e$message), call. = FALSE)
          #     rep(FALSE, nrow(relevant_cells))
          # })
          
          # For unit_contrast, we average all cells identified by the scope. 
          # The formula A in unit_contrast (e.g. ~Face) implicitly means all levels of Face.
          # Therefore, all `relevant_cells` (after `where` clause) should be included in the averaging.
          keepA_rel <- rep(TRUE, nrow(relevant_cells))
          
          # Create unique names/identifiers for the relevant cells
          cell_names_rel <- apply(relevant_cells, 1, paste, collapse = "_")
          
          # Calculate weights - average of cells identified by formula A
          # Use the new .mask_to_weights helper
          weights_rel_named <- .mask_to_weights(names = cell_names_rel, A_mask = keepA_rel, B_mask = NULL) 
          # .make_weights already averages (divides by nA). So weights_rel_named is correct.
          
          # Ensure output is matrix
          weights_out <- matrix(weights_rel_named, ncol = 1)
          rownames(weights_out) <- names(weights_rel_named)
          colnames(weights_out) <- x$name
          cell_names_out <- names(weights_rel_named)
      }
  }

  # Return structure focused on cell-based weights
  ret <- list(
    term=term,
    name=x$name,
    weights=weights_out, # Weights relative to relevant cells
    condnames=cell_names_out, # Names of relevant cells
    contrast_spec=x
  )
  
  class(ret) <- c("unit_contrast", "cell_contrast", "contrast", "list")
  ret
}


#' @export
`-.contrast_spec` <- function(e1, e2, ...){
  assert_that(inherits(e2, "contrast_spec"))
  structure(list(
    name=paste0(e1$name, ":", e2$name),
    con1=e1,
    con2=e2),
    class=c("contrast_diff_spec", "contrast_spec", "list")
  )
}


#' One-way Contrast Weights
#'
#' @description
#' Compute the contrast weights for an oneway_contrast_spec object.
#'
#' @param x An oneway_contrast_spec object.
#' @param term A term object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the term, name, weights, condition names, and contrast specification.
#'
#' @export
contrast_weights.oneway_contrast_spec <- function(x, term,...) {
  # Get cells (categorical only)
  term_cells <- cells(term)
  if (nrow(term_cells) == 0) {
       warning(paste("Contrast '", x$name, "': Term '", term$varname, "' has no categorical cells."), call. = FALSE)
       weights_out <- matrix(numeric(0), nrow = 0, ncol = 0) # F-contrast, ncol unknown yet
       cell_names_out <- character(0)
  } else { 
      # Apply 'where' clause
      keep <- if (!is.null(x$where)) {
        tryCatch(rlang::eval_tidy(rlang::f_rhs(x$where), data = term_cells), error = function(e) {
            warning(paste("Contrast '", x$name, "': Error evaluating 'where' clause: ", e$message), call. = FALSE)
            rep(FALSE, nrow(term_cells))
        })
      } else {
        rep(TRUE, nrow(term_cells))
      }
      
      relevant_cells <- term_cells[keep, , drop = FALSE]
      
      if (nrow(relevant_cells) == 0) {
          warning(paste("Contrast '", x$name, "' resulted in no relevant cells after applying the 'where' clause."), call. = FALSE)
          weights_out <- matrix(numeric(0), nrow = 0, ncol = 0)
          cell_names_out <- character(0)
      } else {
          # Identify the factor for the main effect
          fac_name <- all.vars(rlang::f_rhs(x$A))
          if (length(fac_name) > 1) {
              warning(paste("Contrast '", x$name, "': one-way contrast has >1 factor specified (", 
                            paste(fac_name, collapse=", "), "), using first: ", fac_name[1]), call.=FALSE)
              fac_name <- fac_name[1]
          }
          if (!(fac_name %in% names(relevant_cells))){
               stop(paste("Contrast '", x$name, "': factor ", fac_name, " not found in relevant cells."), call.=FALSE)
          }
          
          # Generate contrast matrix relative to the levels of the factor in relevant_cells
          # Assuming generate_main_effect_contrast works on the cell structure
          # It should return a matrix where rows correspond to relevant_cells
          cmat <- tryCatch(generate_main_effect_contrast(relevant_cells, fac_name), error = function(e) {
              stop(paste("Contrast '", x$name, "': Error generating main effect contrast for factor ", fac_name, ": ", e$message), call.=FALSE)
          })
          
          # Ensure rownames match cell identifiers (if helper doesn't add them)
          # Create unique names/identifiers for the relevant cells 
          cell_names_rel <- apply(relevant_cells, 1, paste, collapse = "_")
          rownames(cmat) <- cell_names_rel
          colnames(cmat) <- paste(x$name, seq_len(ncol(cmat)), sep="_") # Name F-contrast columns
          
          weights_out <- cmat
          cell_names_out <- cell_names_rel
      }
  }

  # Return structure focused on cell-based weights
  ret <- list(
    term = term,
    name = x$name,
    weights = weights_out, # Weights matrix relative to relevant cells
    condnames = cell_names_out, # Names of relevant cells
    contrast_spec = x
  )
  
  # Classify as Fcontrast if multiple columns, otherwise simple contrast
  base_class <- if(ncol(weights_out) > 1) "Fcontrast" else "contrast"
  class(ret) <- c("oneway_contrast", base_class, "cell_contrast", "contrast", "list")
  ret
}

#' Interaction Contrast Weights
#'
#' @description
#' Compute the contrast weights for an interaction_contrast_spec object.
#'
#' @param x An interaction_contrast_spec object.
#' @param term A term object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the term, name, weights, condition names, and contrast specification.
#'
#' @export
contrast_weights.interaction_contrast_spec <- function(x, term,...) {
  # Get cells (categorical only)
  term_cells <- cells(term)
  if (nrow(term_cells) == 0) {
       warning(paste("Contrast '", x$name, "': Term '", term$varname, "' has no categorical cells for interaction."), call. = FALSE)
       weights_out <- matrix(numeric(0), nrow = 0, ncol = 0)
       cell_names_out <- character(0)
  } else {
      # Apply 'where' clause
      keep <- if (!is.null(x$where)) {
        tryCatch(rlang::eval_tidy(rlang::f_rhs(x$where), data = term_cells), error = function(e) {
            warning(paste("Contrast '", x$name, "': Error evaluating 'where' clause: ", e$message), call. = FALSE)
            rep(FALSE, nrow(term_cells))
        })
      } else {
        rep(TRUE, nrow(term_cells))
      }
      
      relevant_cells <- term_cells[keep, , drop = FALSE]
      
      if (nrow(relevant_cells) == 0) {
          warning(paste("Contrast '", x$name, "' resulted in no relevant cells after applying the 'where' clause."), call. = FALSE)
          weights_out <- matrix(numeric(0), nrow = 0, ncol = 0)
          cell_names_out <- character(0)
      } else {
          # Identify factors for interaction
          factors <- all.vars(rlang::f_rhs(x$A))
          if (length(factors) < 2) {
              stop(paste("Contrast '", x$name, "': Interaction contrast requires at least two factors."), call.=FALSE)
          }
          if (!all(factors %in% names(relevant_cells))) {
              missing_facs <- factors[!factors %in% names(relevant_cells)]
              stop(paste("Contrast '", x$name, "': Factor(s)", paste(missing_facs, collapse=", "), "not found in relevant cells."), call.=FALSE)
          }
          
          # Generate interaction contrast matrix relative to relevant cells
          # Assuming generate_interaction_contrast works on cell structure
          cmat <- tryCatch(generate_interaction_contrast(relevant_cells, factors), error = function(e) {
              stop(paste("Contrast '", x$name, "': Error generating interaction contrast: ", e$message), call.=FALSE)
          })
          
          # Ensure rownames match cell identifiers
          cell_names_rel <- apply(relevant_cells, 1, paste, collapse = "_")
          rownames(cmat) <- cell_names_rel
          colnames(cmat) <- paste(x$name, seq_len(ncol(cmat)), sep="_") # Name F-contrast columns
          
          weights_out <- cmat
          cell_names_out <- cell_names_rel
      }
  }
  
  # Return structure focused on cell-based weights
  ret <- list(
    term = term,
    name = x$name,
    weights = weights_out, # Weights matrix relative to relevant cells
    condnames = cell_names_out, # Names of relevant cells
    contrast_spec = x
  )
  
  # Classify as Fcontrast (interactions usually are)
  base_class <- if(is.null(ncol(weights_out)) || ncol(weights_out) > 1) "Fcontrast" else "contrast"
  class(ret) <- c("interaction_contrast", base_class, "cell_contrast", "contrast", "list")
  ret
}

#' Polynomial Contrast Weights
#'
#' @description
#' Compute the contrast weights for a poly_contrast_spec object.
#'
#' @param x A poly_contrast_spec object.
#' @param term A term object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the term, name, weights, condition names, and contrast specification.
#'
#' @export
contrast_weights.poly_contrast_spec <- function(x, term,...) {
  # Get cells (categorical only)
  term_cells <- cells(term)
  if (nrow(term_cells) == 0) {
       warning(paste("Contrast '", x$name, "': Term '", term$varname, "' has no categorical cells for poly contrast."), call. = FALSE)
       weights_out <- matrix(numeric(0), nrow = 0, ncol = x$degree)
       colnames(weights_out) <- paste(x$name, 1:x$degree, sep="_")
       cell_names_out <- character(0)
  } else {
      # Apply 'where' clause
      keep <- if (!is.null(x$where)) {
        tryCatch(rlang::eval_tidy(rlang::f_rhs(x$where), data = term_cells), error = function(e) {
            warning(paste("Contrast '", x$name, "': Error evaluating 'where' clause: ", e$message), call. = FALSE)
            rep(FALSE, nrow(term_cells))
        })
      } else {
        rep(TRUE, nrow(term_cells))
      }
      
      relevant_cells <- term_cells[keep, , drop = FALSE]

      if (nrow(relevant_cells) == 0) {
          warning(paste("Contrast '", x$name, "' resulted in no relevant cells after applying the 'where' clause."), call. = FALSE)
          weights_out <- matrix(numeric(0), nrow = 0, ncol = x$degree)
          colnames(weights_out) <- paste(x$name, 1:x$degree, sep="_")
          cell_names_out <- character(0)
      } else {
          # Evaluate the formula A on relevant cells to get the factor levels
          vals_fac <- tryCatch(rlang::eval_tidy(rlang::f_rhs(x$A), data = relevant_cells), error = function(e) {
              stop(paste("Contrast '", x$name, "': Error evaluating formula A: ", e$message), call.=FALSE)
          })
          
          # Convert factor levels to numeric values (using value_map or direct coercion)
          vals_num <- if (is.null(x$value_map)) {
              tryCatch(as.numeric(as.character(vals_fac)), warning = function(w){
                  stop(paste("Contrast '", x$name, "': Cannot coerce factor levels from formula A to numeric for poly contrast. Use value_map? Error: ", w$message), call.=FALSE)
              })
          } else {
              mapped_vals <- x$value_map[as.character(vals_fac)]
              if(anyNA(mapped_vals) || length(mapped_vals) != length(vals_fac)){
                  stop(paste("Contrast '", x$name, "': value_map does not cover all factor levels present in relevant cells."), call.=FALSE)
              }
              unlist(mapped_vals)
          }
          
          # Check for sufficient unique points for the polynomial degree
          if (length(unique(vals_num)) <= x$degree) {
               stop(paste("Contrast '", x$name, "': Polynomial degree (", x$degree, 
                           ") is too high for the number of unique points (", length(unique(vals_num)), ") in relevant cells."), call.=FALSE)
          }

          # Calculate orthogonal polynomial weights relative to relevant cells
          pvals_mat <- tryCatch(stats::poly(vals_num, degree = x$degree), error = function(e){
               stop(paste("Contrast '", x$name, "': Error calculating polynomial weights: ", e$message), call.=FALSE)
          })
          
          # Create unique names/identifiers for the relevant cells 
          cell_names_rel <- apply(relevant_cells, 1, paste, collapse = "_")
          rownames(pvals_mat) <- cell_names_rel
          colnames(pvals_mat) <- paste(x$name, 1:x$degree, sep="_")
          
          weights_out <- pvals_mat
          cell_names_out <- cell_names_rel
      }
  }
  
  # Return structure focused on cell-based weights
  ret <- list(
    term = term,
    name = x$name,
    weights = weights_out, # Weights matrix relative to relevant cells
    condnames = cell_names_out, # Names of relevant cells
    contrast_spec = x
  )
  
  # Classify as Fcontrast (poly usually has multiple columns)
  base_class <- if(is.null(ncol(weights_out)) || ncol(weights_out) > 1) "Fcontrast" else "contrast"
  class(ret) <- c("poly_contrast", base_class, "cell_contrast", "contrast", "list")
  ret
}

#' Pair Contrast Weights
#'
#' @description
#' Compute the contrast weights for a pair_contrast_spec object.
#'
#' @param x A pair_contrast_spec object.
#' @param term A term object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the term, name, weights, condition names, and contrast specification.
#'
#' @export
contrast_weights.pair_contrast_spec <- function(x, term,...) {
  # Get cells (only categorical combinations) for evaluating formulas A & B
  term_cells <- cells(term)
  
  # Determine if basis expansion is needed
  expand_basis <- FALSE
  nbasis <- 1L
  hrfspec <- attr(term, "hrfspec")
  if (!is.null(hrfspec) && !is.null(hrfspec$hrf)) {
    hrf_fun <- hrfspec$hrf
    # Check nbasis using the nbasis generic
    if (!inherits(try(nbasis(hrf_fun), silent=TRUE), "try-error") && nbasis(hrf_fun) > 1) {
       nbasis <- nbasis(hrf_fun)
       expand_basis <- TRUE
    }
  }
  
  # Get base condition names (pre-expansion) - these are the targets for weights
  base_condnames_all <- try(conditions(term, drop.empty=FALSE, expand_basis=FALSE), silent=TRUE)
  if (inherits(base_condnames_all, "try-error") || length(base_condnames_all) == 0) {
      warning(paste("Contrast '", x$name, "': Failed to get base condition names for term '", term$varname, "'. Skipping."), call. = FALSE)
      return(NULL) # Return NULL if cannot proceed
  }
  
  if (nrow(term_cells) == 0) {
      # If term_cells is empty, we can still proceed if base_condnames_all exist,
      # but the contrast formulas A/B likely won't match anything.
      warning(paste("Contrast '", x$name, "': Term '", term$varname, "' has no observed categorical cells. Contrast weights might be all zero."), call. = FALSE)
  }
  
  # Evaluate 'where' clause on term_cells
  keep <- if (!is.null(x$where)) {
      tryCatch(rlang::eval_tidy(rlang::f_rhs(x$where), data = term_cells), error = function(e) {
          warning(paste("Contrast '", x$name, "': Error evaluating 'where' clause: ", e$message), call. = FALSE)
          rep(FALSE, nrow(term_cells)) # Default to FALSE on error
      })
  } else {
      rep(TRUE, nrow(term_cells)) # Keep all if no 'where'
  }
  relevant_cells <- term_cells[keep, , drop = FALSE]
  
  if (nrow(relevant_cells) == 0 && nrow(term_cells) > 0) {
       warning(paste("Contrast '", x$name, "' resulted in no relevant cells after applying the 'where' clause."), call. = FALSE)
       # Proceed, but weights will likely be zero
  }
  
  # Evaluate A and B formulas on the relevant cells
  keepA_rel <- if(nrow(relevant_cells) > 0) tryCatch(rlang::eval_tidy(rlang::f_rhs(x$A), data = relevant_cells), error = function(e) {
      warning(paste("Contrast '", x$name, "': Error evaluating formula A: ", e$message), call. = FALSE)
      logical(nrow(relevant_cells)) # Return all FALSE on error
  }) else logical(0)
  
  keepB_rel <- if (is.null(x$B) || nrow(relevant_cells) == 0) NULL else {
      tryCatch(rlang::eval_tidy(rlang::f_rhs(x$B), data = relevant_cells), error = function(e) {
          warning(paste("Contrast '", x$name, "': Error evaluating formula B: ", e$message), call. = FALSE)
          logical(nrow(relevant_cells)) # Return all FALSE on error
      })
  }
  
  # --- FIX: Map formula evaluation results to condition names --- 
  mask_A_full <- logical(length(base_condnames_all))
  mask_B_full <- if (!is.null(keepB_rel)) logical(length(base_condnames_all)) else NULL
  
  if (nrow(relevant_cells) > 0) {
    # Get the variable names involved in the term (must match cells colnames)
    term_vars <- names(relevant_cells)
    
    for (i in 1:nrow(relevant_cells)) {
      cell_row <- relevant_cells[i, , drop = FALSE]
      
      # Construct the corresponding condition name(s) for this cell row
      # This needs to handle interactions correctly.
      # For a single factor 'condition' with level 'a', it should produce 'condition.a'
      # For factors 'condition', 'task' with levels 'a', 'X', it should produce 'condition.a_task.X'
      tokens <- vector("character", length(term_vars))
      for (j in seq_along(term_vars)) {
         var_name <- term_vars[j]
         level_val <- cell_row[[var_name]] # Get the factor level value
         # Assume level_token handles sanitization if needed
         tokens[j] <- level_token(var_name, level_val)
      }
      target_cond_name <- make_cond_tag(tokens)
      
      # Find the index in the full list
      idx <- match(target_cond_name, base_condnames_all)
      
      if (!is.na(idx)) {
        if (i <= length(keepA_rel) && keepA_rel[i]) {
          mask_A_full[idx] <- TRUE
        }
        if (!is.null(mask_B_full) && i <= length(keepB_rel) && keepB_rel[i]) {
          mask_B_full[idx] <- TRUE
        }
      }
    }
  }
  # --- END FIX --- 

  # Calculate base weights relative to the full set of base conditions
  base_weights_named <- .make_weights(base_condnames_all, mask_A_full, mask_B_full)
  
  # --- Basis Expansion --- 
  if (expand_basis) {
     # Get expanded condition names
     expanded_condnames <- try(conditions(term, drop.empty = FALSE, expand_basis = TRUE), silent = TRUE)
     if (inherits(expanded_condnames, "try-error") || length(expanded_condnames) == 0) {
         warning(paste("Contrast '", x$name, "': Failed to get expanded condition names for term '", term$varname, "'. Skipping expansion."), call. = FALSE)
         weights_out <- matrix(base_weights_named, ncol = 1)
         rownames(weights_out) <- names(base_weights_named)
         colnames(weights_out) <- x$name
         cell_names_out <- names(base_weights_named)
     } else if (length(expanded_condnames) != length(base_weights_named) * nbasis) {
          warning(paste("Contrast '", x$name, "': Mismatch between expanded names (", length(expanded_condnames), 
                        ") and expected from base weights * nbasis (", length(base_weights_named) * nbasis, "). Skipping expansion."), call. = FALSE)
          weights_out <- matrix(base_weights_named, ncol = 1)
          rownames(weights_out) <- names(base_weights_named)
          colnames(weights_out) <- x$name
          cell_names_out <- names(base_weights_named)
     } else {
         # --- REMOVED MANUAL BROADCAST ---
         # weights_out_vec <- rep(base_weights_named, each = nbasis)
         # --- Instead, align using names ---
         weights_out <- matrix(0, nrow = length(expanded_condnames), ncol = 1)
         rownames(weights_out) <- expanded_condnames
         colnames(weights_out) <- x$name
         # Find expanded names corresponding to non-zero base weights
         base_names_with_weight <- names(base_weights_named)[base_weights_named != 0]
         # Create regex patterns to match expanded names
         # patterns_to_match <- paste0("^", base_names_with_weight, "(\_b\\d+)?$") # Match name + optional _b suffix
         
         for(i in seq_along(base_names_with_weight)){
             base_name <- base_names_with_weight[i]
             weight_val <- base_weights_named[base_name]
             # Find all expanded conditions that start with this base name
             # Use fixed=TRUE for exact match of base name part
             pattern <- paste0("^", base_name, "(_b\\d+)?$")   # note: _b\\d+
             matching_expanded_indices <- grep(pattern, expanded_condnames, perl = TRUE)
            if(length(matching_expanded_indices) > 0){
                 weights_out[matching_expanded_indices, 1] <- weight_val 
             }
         }
         # Ensure weights still sum to zero (or close to it)
         if (abs(sum(weights_out)) > 1e-8) {
              warning(paste("Contrast ", x$name, ": Weights do not sum to zero after basis expansion."), call.=FALSE)
         }
         cell_names_out <- expanded_condnames
     }
  } else {
     # No expansion needed
     weights_out <- matrix(base_weights_named, ncol = 1)
     rownames(weights_out) <- names(base_weights_named)
     colnames(weights_out) <- x$name
     cell_names_out <- names(base_weights_named)
  }

  # Return structure focused on cell-based weights
  ret <- list(
    term = term,
    name = x$name,
    weights = weights_out, # Weights relative to (potentially expanded) relevant cells
    condnames = cell_names_out, # Names of the relevant cells corresponding to rows of weights
    contrast_spec = x
  )
  
  class(ret) <- c("cell_contrast", "contrast", "list") 
  ret  
}

#' Column Contrast Weights
#'
#' @description
#' Compute contrast weights for a `column_contrast_spec` object by targeting
#' design matrix columns based on regex patterns.
#'
#' @param x A `column_contrast_spec` object.
#' @param term An `event_term` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the contrast details:
#'   \item{term}{The original `event_term` object.}
#'   \item{name}{The name of the contrast.}
#'   \item{weights}{A numeric matrix where rows correspond to the full design
#'                  matrix columns (from `.condnames(term, expanded = TRUE)`)
#'                  and columns represent the contrast(s). Usually one column.}
#'   \item{condnames}{Character vector of all potential *expanded* condition names from `term`.}
#'   \item{contrast_spec}{The original `column_contrast_spec` object.}
#'
#' @export
#' @import assertthat
contrast_weights.column_contrast_spec <- function(x, term, ...) {

  # --- Use .condnames helper to get expanded names --- 
  all_colnames <- .condnames(term, expanded = TRUE)
  # Note: .condnames() already includes error handling for conditions()
  if (length(all_colnames) == 0) {
      # It's possible conditions() returns empty if term has no levels/columns
      warning(paste("Column contrast '", x$name, "': Term '", term$varname %||% "<unknown>",
                    "' resulted in zero condition names. Weights will be empty."), call. = FALSE)
      # Return structure with empty weights
      weights_out <- matrix(numeric(0), nrow = 0, ncol = 1)
      colnames(weights_out) <- x$name
      ret <- list(
          term = term,
          name = x$name,
          weights = weights_out,
          condnames = character(0),
          contrast_spec = x
      )
      class(ret) <- c("column_contrast", "contrast", "list")
      return(ret)
  }
  num_all_conds <- length(all_colnames)

  # --- Find indices using grep directly --- 
  # No need for .col_index intermediary anymore if we have names
  idx_A <- grep(x$pattern_A, all_colnames, value = FALSE)
  nA <- length(idx_A)
  if (nA == 0) {
    warning(paste("Column contrast '", x$name, "': pattern_A ('", x$pattern_A,
                  "') did not match any design matrix columns for term '", term$varname %||% "<unknown>", "'."),
            call. = FALSE)
  }

  idx_B <- integer(0)
  nB <- 0
  if (!is.null(x$pattern_B)) {
    idx_B <- grep(x$pattern_B, all_colnames, value = FALSE)
    nB <- length(idx_B)
    if (nB == 0) {
      warning(paste("Column contrast '", x$name, "': pattern_B ('", x$pattern_B,
                    "') did not match any design matrix columns for term '", term$varname %||% "<unknown>", "'."),
              call. = FALSE)
    }
  }

  # --- Retain overlap check --- 
  if (nA > 0 && nB > 0 && any(idx_A %in% idx_B)) {
    stop(paste("Column contrast '", x$name, "': pattern_A and pattern_B match overlapping columns.",
               " Indices A: ", paste(idx_A, collapse=", "),
               "; Indices B: ", paste(idx_B, collapse=", ")), 
         call. = FALSE)
  }

  # --- Calculate weights using .mask_to_weights --- 
  mask_A <- seq_along(all_colnames) %in% idx_A
  mask_B <- if (!is.null(x$pattern_B)) seq_along(all_colnames) %in% idx_B else NULL
  
  # .mask_to_weights handles 1/nA, -1/nB, checks, and warnings
  weights_vec <- .mask_to_weights(all_colnames, mask_A, mask_B)
  
  # Ensure output is a matrix
  weights_mat <- matrix(weights_vec, ncol = 1)
  rownames(weights_mat) <- all_colnames
  colnames(weights_mat) <- x$name

  # Return structure
  ret <- list(
    term = term,
    name = x$name,
    weights = weights_mat,
    condnames = all_colnames, # These are the expanded names used for weights
    contrast_spec = x
  )

  # Classify appropriately
  class(ret) <- c("column_contrast", "contrast", "list")
  ret
}

#' Contrast Formula Weights
#'
#' @description
#' Compute the contrast weights for a contrast_formula_spec object.
#'
#' @param x A contrast_formula_spec object.
#' @param term A term object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the term, name, weights, condition names, and contrast specification.
#'
#' @export
contrast_weights.contrast_formula_spec <- function(x, term,...) {

  term.cells <- cells(term)
  condnames <- shortnames(term)
  count <- attr(term.cells, "count")		
  term.cells <- subset(term.cells, count > 0)
  
  if (!is.null(x$where)) {
    keep <- rlang::eval_tidy(rlang::f_rhs(x$where), data=term.cells)
    assert_that(sum(keep) > 0)
    term.cells <- term.cells[keep,]
  } else {
    keep <- rep(TRUE, nrow(term.cells))
  }
  
  if (is_continuous(term)) {
    ## hack to handle contrasts with continuous terms
    facs <- !sapply(term$events, is_continuous)
    term.cells[,!facs] <- 1
  } 

  # Create a simple environment with shortnames as variables
  # Each shortname gets a vector with 1 in its position and 0 elsewhere
  weights <- matrix(0, NROW(term.cells), 1)
  eval_env <- new.env(parent = rlang::f_env(x$A))
  
  # Create indicator variables for each condition using shortnames
  for (i in seq_along(condnames)) {
    indicator <- rep(0, length(condnames))
    indicator[i] <- 1
    assign(condnames[i], indicator, envir = eval_env)
  }
  
  # Evaluate the contrast formula
  res <- tryCatch(rlang::eval_tidy(rlang::f_rhs(x$A), env=eval_env), error = function(e) {
       stop(paste("Contrast formula evaluation failed:", e$message, "\nAvailable names:", paste(condnames, collapse=", ")), call.=FALSE)
  })
  
  # Apply results to weights matrix
  weights[keep,1] <- as.vector(res)
  
  # Use longnames for rownames (the new format)
  row.names(weights) <- longnames(term)[keep]

  # Return structure
  ret <- list(
    term=term,
    name=x$name,
    weights=weights,
    condnames=longnames(term),
    contrast_spec=x)
  
  class(ret) <- c("contrast", "list")
  ret  
}

#' Contrast Difference Weights
#'
#' @description
#' Compute the contrast weights for a contrast_diff_spec object.
#'
#' @param x A contrast_diff_spec object.
#' @param term A term object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the term, name, weights, condition names, and contrast specification.
#'
#' @export
contrast_weights.contrast_diff_spec <- function(x, term,...) {
  wts1 <- contrast_weights(x$con1, term)
  wts2 <- contrast_weights(x$con2, term)

  ret <- structure(
    list(
      term=term,
      name=x$name,
      weights=wts1$weights - wts2$weights,
      condnames=longnames(term),
      contrast_spec=x),
    class=c("contrast_diff", "contrast")
  )

  ret  
}

#' Convert Contrast to GLT
#'
#' @description
#' Convert a contrast to an AFNI 'GLT' format.
#'
#' @param x The contrast to convert.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the GLT string, name, and contrast specification.
#'
#' @export
to_glt.contrast <- function(x,...) {
  if (is.matrix(x$weights) && ncol(x$weights) > 1) {
    glts <- lapply(1:ncol(x$weights), function(i) {
      paste0(signif(x$weights[,i],4), "*", x$condnames, collapse=" ")
    })
    
    ret <- list(glt_str=glts,
                name=paste0("GLT_", x$name, "_", 1:ncol(x$weights)),
                con=x)
    
    class(ret) <- "glt_contrast_list"
    ret
  } else {
    glt <- paste0(signif(x$weights,4), "*", x$condnames, collapse=" ")
    ret <- list(glt_str=glt,
       name=paste0("GLT_", x$name),
       con=x)
  
    class(ret) <- "glt_contrast"
    ret
  }
}

#' Write GLT to File
#'
#' @description
#' Write a generic function for writing GLT contrast to a file.
#'
#' @param x The GLT contrast object
#' @param ... Additional arguments passed to methods
#' @export
write_glt <- function(x, ...) {
  UseMethod("write_glt")
}

#' Write GLT to File
#'
#' @description
#' Write a GLT contrast to a file.
#'
#' @param x The GLT contrast to write.
#' @param fname The file name to write to.
#'
#' @export
write_glt.glt_contrast <- function(x, fname=NULL) {
  con <- if (is.null(fname)) {
    file(fname, "w")
  } else {
    file(paste0(x$name, ".txt"), "w")
  }
  
  writeLines(x$glt_str, con=con)
  close(con)
}

#' Estimated Contrast
#'
#' @description
#' Compute the estimated contrast for a given fit and indices.
#'
#' @param x The contrast to estimate.
#' @param fit The fit object.
#' @param indices The indices to use.
#' @param ... Additional arguments (currently unused).
#'
#' @return The estimated contrast.
#'
#' @noRd
#' @keywords internal
estcon.contrast <- function(x, fit, indices, ...) {
  wts <- numeric(length(fit$assign))
  wts[indices] <- x$weights
  
  gmodels::estimable(fit, wts)
}

#' Print Contrast Set
#'
#' @description
#' Print a contrast set.
#'
#' @param x The contrast set to print.
#' @param ... Additional arguments (currently unused).
#'
#' @export
print.contrast_set <- function(x, ...) {
  n_contrasts <- length(x)
  
  # Header
  cat("\n═══ Contrast Set ═══\n")
  
  # Summary
  cat("\n Overview:\n")
  cat("  • Number of contrasts:", n_contrasts, "\n")
  
  # Group contrasts by type
  types <- sapply(x, function(con) class(con)[1])
  type_table <- table(types)
  if (length(type_table) > 0) {
    cat("  • Types of contrasts:\n")
    for (type in names(type_table)) {
      cat("    -", type, ":", type_table[type], "\n")
    }
  }
  
  # List all contrasts
  cat("\n🔍 Individual Contrasts:\n")
  for (i in seq_along(x)) {
    cat("\n[", i, "] ", x[[i]]$name, " (", class(x[[i]])[1], ")\n", sep="")
    cat("    Formula: ")
    if (!is.null(x[[i]]$A)) cat(deparse(x[[i]]$A))
    if (!is.null(x[[i]]$B)) cat(" vs ", deparse(x[[i]]$B))
    cat("\n")
    if (!is.null(x[[i]]$where)) {
      cat("    Subset: ", deparse(x[[i]]$where), "\n")
    }
  }
  
  cat("\n")
  invisible(x)
}

#' Print Contrast Specification
#'
#' @description
#' Print a contrast specification.
#'
#' @param x The contrast specification to print.
#' @param ... Additional arguments (currently unused).
#'
#' @export
print.contrast_spec <- function(x,...) {
  cat("contrast:", x$name, "\n")
  cat(" A: ", Reduce(paste, deparse(x$A)), "\n")
  if (!is.null(x$B))
    cat(" B: ", Reduce(paste, deparse(x$B)), "\n")
  if (!is.null(x$where))
    cat(" where: ", Reduce(paste, deparse(x$where)), "\n")
  

}

#' Print Contrast
#'
#' @description
#' Print a contrast.
#'
#' @param x The contrast to print.
#' @param ... Additional arguments (currently unused).
#'
#' @export
print.contrast <- function(x,...) {
  print(x$contrast_spec)
  cat(" term: ", x$term$varname, "\n")
  cat(" weights: ", "\n")
  print(x$weights)
  cat(" conditions: ", x$condnames)
  
}

#' Print Polynomial Contrast Specification
#'
#' @description
#' Print a polynomial contrast specification.
#'
#' @param x The polynomial contrast specification to print.
#' @param ... Additional arguments (currently unused).
#'
#' @export
print.poly_contrast_spec <- function(x,...) {
  cat("poly contrast:", "\n")
  cat(" A: ", Reduce(paste, deparse(x$A)), "\n")
  cat(" degree: ", x$degree, "\n")
  if (!is.null(x$where)) {
    cat(" where: ", deparse(x$where), "\n")
  }
  
  if (!is.null(x$value_map)) {
    cat(" values: ", unlist(x$value_map), "\n")
  }
}

#' Print Contrast Difference Specification
#'
#' @description
#' Print a contrast difference specification.
#'
#' @param x The contrast difference specification to print.
#' @param ... Additional arguments (currently unused).
#'
#' @export
print.contrast_diff_spec <- function(x,...) {
  cat("contrast difference:", "\n")
  cat("  ", x$con1$name, "-", x$con2$name, "\n")
}

#' plot_contrasts
#'
#' @description
#' Generic function for plotting contrasts.
#'
#' @param x Object containing contrast information
#' @param ... Additional arguments passed to methods
#' @export
plot_contrasts <- function(x, ...) {
  UseMethod("plot_contrasts")
}

#' plot_contrasts.event_model
#'
#' @description
#' Produces a heatmap of all contrasts defined for an \code{event_model}.
#' Rows = each contrast (or column of an F-contrast), columns = each regressor in
#' the full design matrix, and the fill color = the contrast weight.
#'
#' @param x An \code{event_model} with (lazily) defined contrasts.
#' @param absolute_limits Logical; if \code{TRUE}, the color scale is fixed at (-1,1).
#'   If \code{FALSE}, the range is set to (min, max) of the weights.
#' @param rotate_x_text Logical; if \code{TRUE}, rotate x-axis labels for readability.
#' @param ... Further arguments passed to \code{geom_tile}, e.g. \code{color="grey80"}.
#'
#' @return A \code{ggplot2} object (a heatmap).
#' @import ggplot2
#' @export
plot_contrasts.event_model <- function(
    x,
    absolute_limits = FALSE,
    rotate_x_text   = TRUE,
    scale_mode      = c("auto", "diverging", "one_sided"),
    coord_fixed     = TRUE,
    ...
) {
  # 1) Extract all the design-matrix column names
  dm <- design_matrix(x)
  regressor_names <- colnames(dm)
  
  # 2) Gather contrast weights (the nested list by term, then by contrast)
  cws <- contrast_weights(x)
  
  # Flatten everything into one big matrix of contrast weights
  big_mat  <- NULL
  rownames <- character(0)
  
  add_contrast_row <- function(vec, row_name) {
    if (is.null(big_mat)) {
      big_mat <<- matrix(vec, nrow = 1)
      colnames(big_mat) <<- regressor_names
      rownames <<- row_name
    } else {
      big_mat <<- rbind(big_mat, vec)
      rownames <<- c(rownames, row_name)
    }
  }
  
  for (term_nm in names(cws)) {
    term_level <- cws[[term_nm]]
    if (is.null(term_level)) next
    
    for (contrast_nm in names(term_level)) {
      cw_obj <- term_level[[contrast_nm]]
      # By default, we store offset_weights in cw_obj$offset_weights
      W <- cw_obj$offset_weights
      if (is.null(W)) next  # skip if no offset_weights
      
      # (#designCols x #contrastCols)
      ncols <- ncol(W)
      for (k in seq_len(ncols)) {
        this_col <- W[, k]
        if (ncols > 1) {
          row_label <- paste0(contrast_nm, "_component", k)
        } else {
          row_label <- contrast_nm
        }
        add_contrast_row(this_col, row_label)
      }
    }
  }
  
  if (is.null(big_mat)) {
    stop("No contrasts found in this event_model.")
  }
  
  rownames(big_mat) <- rownames
  
  # 3) Convert big_mat to a long data frame
  df_long <- as.data.frame(big_mat, check.names = FALSE)
  df_long$ContrastName <- rownames(big_mat)
  df_long <- reshape2::melt(
    df_long,
    id.vars      = "ContrastName",
    variable.name = "Regressor",
    value.name    = "Weight"
  )
  
  # 4) Build the ggplot
  plt <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x    = ReorderFactor(Regressor),
      y    = ReorderFactor(ContrastName, reverse = TRUE),
      fill = Weight
    )
  ) +
    ggplot2::geom_tile(...)
  
  # Decide on color scale
  scale_mode <- match.arg(scale_mode)
  wmin <- min(df_long$Weight, na.rm = TRUE)
  wmax <- max(df_long$Weight, na.rm = TRUE)
  
  # If user set absolute_limits=TRUE, we might forcibly use [-1,1] or [0,1]
  # but let's allow scale_mode to override as well.
  if (scale_mode == "diverging") {
    # Diverging scale, centered on 0
    lim_low  <- if (absolute_limits) -1 else wmin
    lim_high <- if (absolute_limits)  1 else wmax
    midpt <- 0
    
    plt <- plt + ggplot2::scale_fill_gradient2(
      limits   = c(lim_low, lim_high),
      midpoint = midpt,
      low      = "blue",
      mid      = "white",
      high     = "red"
    )
    
  } else if (scale_mode == "one_sided") {
    # One-sided scale for 0..1 or 0..something
    lim_low  <- if (absolute_limits) 0 else wmin
    lim_high <- if (absolute_limits) 1 else wmax
    
    plt <- plt + ggplot2::scale_fill_gradient(
      limits = c(lim_low, lim_high),
      low    = "white",
      high   = "red"
    )
    
  } else {
    # scale_mode == "auto"
    # If we detect any negative weight, do diverging; otherwise do one-sided
    if (wmin < 0) {
      # diverging
      lim_low  <- if (absolute_limits) -1 else wmin
      lim_high <- if (absolute_limits)  1 else wmax
      plt <- plt + ggplot2::scale_fill_gradient2(
        limits   = c(lim_low, lim_high),
        midpoint = 0,
        low      = "blue",
        mid      = "white",
        high     = "red"
      )
    } else {
      # one-sided
      lim_low  <- if (absolute_limits) 0 else wmin
      lim_high <- if (absolute_limits) 1 else wmax
      plt <- plt + ggplot2::scale_fill_gradient(
        limits = c(lim_low, lim_high),
        low    = "white",
        high   = "red"
      )
    }
  }
  
  # 5) Theming
  plt <- plt +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      x    = "Regressor",
      y    = "Contrast",
      fill = "Weight"
    ) +
    ggplot2::theme(
      panel.grid  = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank()
    )
  
  if (rotate_x_text) {
    plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  
  # 6) Optionally fix the coordinate ratio to keep tiles square
  if (coord_fixed) {
    plt <- plt +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::coord_fixed()
  }
  
  plt
}

#' A small utility to preserve factor order in ggplot
#' 
#' Makes a factor from a character vector but preserves the order of appearance.
#' If `reverse=TRUE`, it reverses that order.
#' @keywords internal
#' @noRd
ReorderFactor <- function(x, reverse=FALSE) {
  levs <- unique(as.character(x))
  if (reverse) levs <- rev(levs)
  factor(x, levels=levs)
}

#' Contrast Weights for a Contrast Set
#'
#' @description
#' Compute the contrast weights for each contrast specification within a contrast_set object.
#'
#' @param x A contrast_set object (a list of contrast_spec objects).
#' @param term A term object against which weights should be computed.
#' @param ... Additional arguments passed to individual contrast_weights methods.
#'
#' @return A named list where each element is the result of calling contrast_weights 
#'         on the corresponding contrast_spec in the set. The list names are the 
#'         names of the individual contrasts.
#'
#' @export
#' @importFrom purrr map set_names
contrast_weights.contrast_set <- function(x, term, ...) {
  # Ensure x is a list (contrast_set inherits from list)
  if (!is.list(x)) {
    stop("Input 'x' must be a contrast_set (list).")
  }
  
  # Iterate through each contrast spec in the set
  results_list <- purrr::map(x, function(contrast_spec) {
    # Check if the element is actually a contrast_spec
    if (!inherits(contrast_spec, "contrast_spec")) {
      warning(paste("Element", contrast_spec$name, "is not a contrast_spec object, skipping."))
      return(NULL)
    }
    # Compute weights for the individual contrast spec
    contrast_weights(contrast_spec, term, ...)
  })
  
  # Filter out any NULL results (if any elements weren't contrast_spec)
  results_list <- results_list[!sapply(results_list, is.null)]
  
  # Set the names of the results list based on the names of the contrasts
  contrast_names <- purrr::map_chr(results_list, "name")
  results_list <- purrr::set_names(results_list, contrast_names)
  
  return(results_list)
}



