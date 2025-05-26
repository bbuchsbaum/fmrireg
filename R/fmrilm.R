#' Get the formula representation of an fMRI model
#'
#' This function extracts the formula from an \code{fmri_model} object.
#'
#' @param x An \code{fmri_model} object.
#' @return A formula representing the model.
#' @export
get_formula.fmri_model <- function(x) {
  assert_that(inherits(x, "fmri_model"), msg = "'x' must be an 'fmri_model' object")
  term_names <- names(terms(x))
  form <- paste(".y ~", paste(term_names, collapse = " + "), "-1")
  return(as.formula(form))
}

#' Extract Term Matrices from an fMRI Model
#'
#' This function extracts the term matrices from an \code{fmri_model}, which consists of event-related terms
#' and baseline-related terms. These matrices are used to build the design matrix in fMRI data analysis.
#'
#' @param x An \code{fmri_model} object containing the event and baseline models.
#' @param blocknum (Optional) A numeric vector specifying the block numbers to include. Defaults to all blocks.
#' @return A named list of term matrices, with event terms followed by baseline terms.
#'         Attributes \code{"event_term_indices"} and \code{"baseline_term_indices"} store the indices of event and baseline terms,
#'         \code{"blocknum"} stores the block numbers, and \code{"varnames"} stores the variable names.
#' @export
term_matrices.fmri_model <- function(x, blocknum = NULL) {
  assert_that(inherits(x, "fmri_model"), msg = "'x' must be an 'fmri_model' object")
  
  if (is.null(blocknum)) {
    blocknum <- sort(unique(x$event_model$blockids))
  }
  
  # Get the full convolved design matrix from the event model
  event_dm <- design_matrix(x$event_model, blockid = blocknum)
  
  # Get the baseline design matrix
  baseline_dm <- design_matrix(x$baseline_model, blockid = blocknum)
  
  # Extract individual term matrices using the col_indices attribute
  col_indices <- attr(x$event_model$design_matrix, "col_indices")
  if (is.null(col_indices)) {
    stop("Event model design matrix missing 'col_indices' attribute needed to extract individual term matrices.")
  }
  
  # Extract event term matrices from the full convolved design matrix
  eterms <- lapply(names(col_indices), function(term_name) {
    indices <- col_indices[[term_name]]
    as.matrix(event_dm[, indices, drop = FALSE])
  })
  names(eterms) <- names(col_indices)
  
  # Extract baseline term matrices (baseline terms are simpler, one per term)
  bterms <- lapply(baseline_terms(x), function(term) as.matrix(design_matrix(term, blockid = blocknum)))
  
  # Compute indices for event and baseline terms
  num_event_cols <- ncol(event_dm)
  num_baseline_cols <- ncol(baseline_dm)
  
  eterm_indices <- 1:num_event_cols
  bterm_indices <- (num_event_cols + 1):(num_event_cols + num_baseline_cols)
  
  # Combine term matrices
  term_matrices <- c(eterms, bterms)
  names(term_matrices) <- names(terms(x))
  
  # Collect variable names
  vnames <- c(colnames(event_dm), colnames(baseline_dm))
  
  # Set attributes
  attr(term_matrices, "event_term_indices") <- eterm_indices
  attr(term_matrices, "baseline_term_indices") <- bterm_indices
  attr(term_matrices, "blocknum") <- blocknum
  attr(term_matrices, "varnames") <- vnames
  
  return(term_matrices)
}



#' @keywords internal
#' @noRd
is.formula <- function(x) {
  inherits(x, "formula")
}

#' @keywords internal
#' @noRd
.fast_preproject <- function(X) {
  # Ensure X is a matrix
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  XtX   <- crossprod(X)        # p × p
  # Add small ridge for stability if needed, but try direct first
  # Rchol <- tryCatch(chol(XtX), error = function(e) chol(XtX + diag(ncol(XtX)) * 1e-10))
  Rchol <- chol(XtX)           # p × p  upper‑triangular
  Pinv  <- backsolve(Rchol, t(X), transpose = TRUE)  # (Rchol^-1)' Xᵀ = (Rchol'^-1) Xᵀ = (XtX)^-1 Xᵀ -> p x n
                                                    # backsolve solves R'y = x for y if transpose=TRUE
                                                    # We want to solve R z = t(X) for z, where R is upper.
                                                    # Let R = chol(XtX). We want z = R^-1 t(X)
                                                    # Pinv should be (XtX)^-1 X^T = (R'R)^-1 X^T = R^-1 R'^-1 X^T
                                                    # Let's use solve(R) %*% solve(t(R)) %*% t(X) ? No.
                                                    # Let's use chol2inv(Rchol) %*% t(X) ? Yes.
  
  # Revisit Pinv calculation based on user's note: Pinv = backsolve(Rchol, t(X)) # (R⁻¹) Xᵀ  →  p × n
  # This assumes R is upper triangular. backsolve solves R x = b.
  # We want (XtX)^-1 Xt = (R'R)^-1 Xt = R^-1 R'^-1 Xt
  # Let's test:
  # X <- matrix(rnorm(30*5), 30, 5); XtX <- crossprod(X); Rchol <- chol(XtX)
  # Pinv_chol2inv <- chol2inv(Rchol) %*% t(X)
  # Pinv_backsolve <- backsolve(Rchol, t(X)) # solves R z = t(X) -> z = R^-1 t(X) - This is NOT (XtX)^-1 Xt
  # Pinv_correct <- solve(XtX, t(X))
  # Let's stick to chol2inv for clarity and correctness
  
  XtXinv <- chol2inv(Rchol)
  Pinv   <- XtXinv %*% t(X) # p x n : (X'X)^-1 X'
  
  list(Pinv = Pinv,           # (X'X)^-1 X'
       XtXinv = XtXinv,       # (X'X)^-1
       Rchol = Rchol,         # For potential future use
       dfres = nrow(X) - ncol(X)) # rank issue not handled here, assumes full rank
}

#' @keywords internal
#' @noRd
.fast_lm_matrix <- function(Y, proj) {
  # Y : n × V  (double matrix, one chunk)
  # proj$Pinv : p × n : (X'X)^-1 X'
  
  if (!is.matrix(Y)) {
     Y <- as.matrix(Y)
  }
  
  # Ensure dimensions match: ncol(Pinv) == nrow(Y)
  if (ncol(proj$Pinv) != nrow(Y)) {
      stop(paste(".fast_lm_matrix: Dimension mismatch between Pinv (",
                 nrow(proj$Pinv),"x",ncol(proj$Pinv), ") and Y (",
                 nrow(Y),"x",ncol(Y),"). Number of time points (n) must match.", sep=""))
  }
  
  Betas <- proj$Pinv %*% Y              # p × V   GEMM 1: (X'X)^-1 X' Y
  
  # Resid = Y - X %*% Betas = Y - X %*% (X'X)^-1 X' Y
  # Avoid forming X %*% Betas if possible.
  # X is needed. How to get X efficiently?
  # We have proj$Pinv = (X'X)^-1 X'. Can we get X back? Not easily.
  # The user note suggests Resid <- Y - crossprod(t(proj$Pinv), Betas)
  # Let's check: crossprod(t(Pinv), Betas) = Pinv' %*% Betas = ( ((X'X)^-1 X')' ) %*% Betas
  # = (X (X'X)^-1) %*% Betas. This is NOT X %*% Betas.
  # We need X itself. X must be passed or reconstructed.
  # Let's require X to be passed to .fast_lm_matrix or reconstruct it inside the caller.
  # For now, assume X is available where .fast_lm_matrix is called.
  # Modify signature: .fast_lm_matrix <- function(X, Y, proj)
  # Or pass X into proj list: proj$X <- X? No, X can be large.
  
  # Let's rethink Resid calculation from user note:
  # Resid <- Y - crossprod(t(proj$Pinv), Betas)
  # This seems incorrect. The user's formula later is correct: Resid = Y - X %*% Betas
  # We absolutely need X. The calling functions (chunkwise_lm, runwise_lm) MUST provide X.
  
  # Let's modify the callers to pass X when using the fast path.
  # For now, let's assume X is passed.
  # Modify signature: .fast_lm_matrix <- function(X, Y, proj)
  
  # Reverting signature for now, will adapt callers.
  
  # How to get X in chunkwise_lm? `modmat` is X.
  # How to get X in runwise_lm? Need to reconstruct `X_run` inside the loop.
  
  # Let's proceed assuming X will be available in the calling context.
  # The user example doesn't pass X explicitly to .fast_lm_matrix. How is Resid calculated?
  # Resid <- Y - crossprod(t(proj$Pinv), Betas) # Let's re-examine this.
  # Is it possible X = t(solve(crossprod(proj$Pinv))) ? No.
  # What if Pinv was defined as t(X)? No.
  
  # OK, the user's code MUST be wrong. `Resid = Y - X %*% Betas` is the definition.
  # We will adapt the callers to provide X to this function.
  
  # Updated plan: Modify .fast_lm_matrix signature and callers.
  # .fast_lm_matrix <- function(X, Y, proj) { ... Resid <- Y - X %*% Betas ... }
  
  # Sticking to the user's provided code for now, assuming there's a BLAS trick I'm missing
  # regarding Resid <- Y - crossprod(t(proj$Pinv), Betas). This needs verification.
  # If Pinv = X (e.g. orthogonal design), then t(Pinv) = X', crossprod(t(Pinv), Betas) = X %*% Betas.
  # But Pinv = (X'X)^-1 X', so t(Pinv) = X (X'X)^-1.
  # crossprod(t(Pinv), Betas) = t(Pinv)' Betas = (X (X'X)^-1)' Betas = (X'X)^-1 X' Betas
  # This is Pinv %*% Betas... which is Betas itself? No. Pinv * Betas gives p x V.
  
  # Let's assume the user meant: X %*% Betas = X %*% (X'X)^-1 X' Y = H Y (H is hat matrix)
  # Maybe crossprod(t(proj$Pinv), Betas) is meant to be calculated differently?
  # If X is available, X %*% Betas is direct.
  
  # Let's implement based on needing X explicitly.
  # Modify signature:
  
  # .fast_lm_matrix <- function(X, Y, proj) {
  #   if (!is.matrix(Y)) { Y <- as.matrix(Y) }
  #   if (!is.matrix(X)) { X <- as.matrix(X) }
  #
  #   Betas <- proj$Pinv %*% Y              # p × V
  #   Fitted <- X %*% Betas                 # n x V
  #   Resid <- Y - Fitted                   # n × V
  #
  #   rss   <- colSums(Resid^2)
  #   sigma2 <- rss / proj$dfres # Use sigma2 for consistency
  #   sigma <- sqrt(sigma2)
  #
  #   list(betas = Betas,   # p x V
  #        fitted = Fitted, # n x V (optional, maybe not needed downstream?)
  #        residuals = Resid, # n x V (optional)
  #        rss   = rss,     # V
  #        sigma = sigma,   # V
  #        sigma2 = sigma2) # V
  # }
  # This requires passing X.
  
  # Let's try the user's version and see if it works in practice or if we need to debug later.
  # Keep user's version for now:
  
  Betas <- proj$Pinv %*% Y              # p × V   GEMM 1
  # This calculation of residuals is mathematically suspect unless Pinv has a special structure
  # or crossprod behaves differently than expected. Let's assume it requires X implicitly or is a placeholder.
  # We will likely need to replace this with Resid = Y - X %*% Betas in the calling function.
  # For now, keep the provided structure.
  # Fitted <- X %*% Betas # Requires X
  # Resid <- Y - Fitted # Requires Fitted
  
  # Calculate RSS without explicitly forming residuals using matrix algebra:
  # RSS = Y'Y - beta'X'Y = Y'Y - beta'X'X beta
  # We have Betas (p x V). We need X'Y. X'Y = X'X (X'X)^-1 X'Y = X'X Betas
  # rss_v = y_v'y_v - beta_v'(X'X)beta_v
  yTy <- colSums(Y^2) # V-vector
  XtX <- solve(proj$XtXinv) # Recover X'X from its inverse
  
  # Need to compute t(beta_v) %*% XtX %*% beta_v for each voxel v
  # This involves V matrix multiplications (1xp * pxp * px1) -> scalar
  # Can optimize: B' (X'X B) -> B is p x V. XtX is p x p
  # tmp = XtX %*% Betas (p x V)
  # diag(t(Betas) %*% tmp) -> V-vector
  XtX_Betas <- XtX %*% Betas # p x V
  beta_XtX_beta <- colSums(Betas * XtX_Betas) # Element-wise product and sum columns -> V-vector
  
  rss <- yTy - beta_XtX_beta # V-vector
  # Check for negative RSS due to numerical precision
  rss[rss < 0 & rss > -1e-10] <- 0
  if (any(rss < -1e-10)) {
      warning("Negative residual sum of squares computed in .fast_lm_matrix")
  }
  
  sigma2 <- rss / proj$dfres
  sigma <- sqrt(sigma2)
  
  # Return only what's needed downstream based on user notes
  list(betas = Betas,   # p x V
       rss   = rss,     # V
       sigma = sigma,   # V
       sigma2 = sigma2) # V
}

#' Create an fMRI Model
#'
#' This function creates an \code{fmri_model} by combining an event model and a baseline model.
#' If a baseline model is not provided, a default one is created based on the dataset.
#'
#' @param formula The model formula for experimental events.
#' @param block The model formula for block structure.
#' @param baseline_model (Optional) A \code{baseline_model} object. If \code{NULL}, a default baseline model is created.
#' @param dataset An \code{fmri_dataset} containing the event table and sampling frame.
#' @param drop_empty Logical. Whether to remove factor levels with zero size. Default is \code{TRUE}.
#' @param durations A vector of event durations. Default is \code{0}.
#' @return An \code{fmri_model} object.
#' @keywords internal
create_fmri_model <- function(formula, block, baseline_model = NULL, dataset, drop_empty = TRUE, durations = 0) {
  assert_that(is.formula(formula), msg = "'formula' must be a formula")
  assert_that(is.formula(block), msg = "'block' must be a formula")
  assert_that(inherits(dataset, "fmri_dataset"), msg = "'dataset' must be an 'fmri_dataset'")
  assert_that(is.numeric(durations), msg = "'durations' must be numeric")
  
  if (is.null(baseline_model)) {
    baseline_model <- baseline_model(
      basis = "bs",
      degree = max(ceiling(median(dataset$sampling_frame$blocklens) / 100), 3),
      sframe = dataset$sampling_frame
    )
  } else {
    assert_that(inherits(baseline_model, "baseline_model"),
                msg = "'baseline_model' must have class 'baseline_model'")
  }
  
  ev_model <- event_model(
    formula_or_list = formula,
    block = block,
    data = dataset$event_table,
    sampling_frame = dataset$sampling_frame,
    drop_empty = drop_empty,
    durations = durations
  )
  
  fmri_model(ev_model, baseline_model)
}



#' Fit a Linear Regression Model for fMRI Data Analysis
#'
#' This function fits a linear regression model for fMRI data analysis using the specified model formula,
#' block structure, and dataset. The model can be fit using either a runwise or chunkwise data splitting strategy,
#' and robust fitting can be enabled if desired.
#'
#' @param formula The model formula for experimental events.
#' @param block The model formula for block structure.
#' @param baseline_model (Optional) A \code{baseline_model} object. Default is \code{NULL}.
#' @param dataset An \code{fmri_dataset} object containing the time-series data.
#' @param durations A vector of event durations. Default is \code{0}.
#' @param drop_empty Logical. Whether to remove factor levels with zero size. Default is \code{TRUE}.
#' @param robust Logical. Whether to use robust fitting. Default is \code{FALSE}.
#' @param strategy The data splitting strategy, either \code{"runwise"} or \code{"chunkwise"}. Default is \code{"runwise"}.
#' @param nchunks Number of data chunks when strategy is \code{"chunkwise"}. Default is \code{10}.
#' @param use_fast_path Logical. If \code{TRUE}, use matrix-based computation for speed. Default is \code{FALSE}.
#' @param progress Logical. Whether to display a progress bar during model fitting. Default is \code{FALSE}.
#' @param ... Additional arguments.
#' @return A fitted linear regression model for fMRI data analysis.
#' @export
#' @seealso \code{\link{fmri_dataset}}, \code{\link{fmri_lm_fit}}
#' @examples
#' 
#' facedes <- subset(read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), 
#' header=TRUE), face_gen != "n/a")
#' facedes$face_gen <- droplevels(factor(facedes$face_gen))
#' sframe <- sampling_frame(rep(430/2,6), TR=2)
#' ev <- event_model(onset ~ hrf(face_gen, basis="gaussian"), data=facedes, 
#' block= ~ run, sampling_frame=sframe)
#' globonsets <- global_onsets(sframe, facedes$onset, blockids(ev))
#' reg1_signal <- regressor(globonsets[facedes$face_gen == "male"], hrf=HRF_GAUSSIAN)
#' reg2_signal <- regressor(globonsets[facedes$face_gen == "female"], hrf=HRF_GAUSSIAN)
#' time <- samples(sframe, global=TRUE)
#' y1 <- evaluate(reg1_signal, time)*1.5
#' y2 <- evaluate(reg2_signal, time)*3.0
#' y <- y1+y2
#' ys1 <- y + rnorm(length(y), sd=.02)
#' ys2 <- y + rnorm(length(y), sd=.02)
#' 
#' h <<- gen_hrf(hrf_bspline, N=7, span=25)
#' dset <- matrix_dataset(cbind(ys1,ys2), TR=2, run_length=sframe$blocklens, event_table=facedes)
#' flm <- fmri_lm(onset ~ hrf(face_gen, basis=gen_hrf(hrf_bspline, N=7, span=25)), block = ~ run, 
#' strategy="chunkwise", nchunks=1, dataset=dset)
#' 
fmri_lm <- function(formula, block, baseline_model = NULL, dataset, durations = 0, drop_empty = TRUE, robust = FALSE,
                    strategy = c("runwise", "chunkwise"), nchunks = 10, use_fast_path = FALSE, progress = FALSE, ...) {
  
  strategy <- match.arg(strategy)
  
  # Error checking
  assert_that(is.formula(formula), msg = "'formula' must be a formula")
  assert_that(is.formula(block), msg = "'block' must be a formula")
  assert_that(inherits(dataset, "fmri_dataset"), msg = "'dataset' must be an 'fmri_dataset'")
  assert_that(is.numeric(durations), msg = "'durations' must be numeric")
  assert_that(is.logical(drop_empty), msg = "'drop_empty' must be logical")
  assert_that(is.logical(robust), msg = "'robust' must be logical")
  assert_that(is.logical(use_fast_path), msg = "'use_fast_path' must be logical")
  if (strategy == "chunkwise") {
    assert_that(is.numeric(nchunks) && nchunks > 0, msg = "'nchunks' must be a positive number")
  }
  if (robust && use_fast_path) {
      warning("Robust fitting ('robust=TRUE') is not currently implemented with 'use_fast_path=TRUE'. Ignoring 'robust=TRUE'.")
      robust <- FALSE # Force robust to FALSE if fast path is used
  }
  
  model <- create_fmri_model(formula, block, baseline_model, dataset, durations = durations, drop_empty = drop_empty)
  # Pass use_fast_path down
  ret <- fmri_lm_fit(model, dataset, strategy, robust, nchunks,
                     use_fast_path = use_fast_path, progress = progress, ...)
  return(ret)
}


#' Fit an fMRI Linear Regression Model with a Specified Fitting Strategy
#'
#' This function fits an fMRI linear regression model using the specified \code{fmri_model} object, dataset,
#' and data splitting strategy (either \code{"runwise"} or \code{"chunkwise"}). It is primarily an internal function
#' used by the \code{fmri_lm} function.
#'
#' @param fmrimod An \code{fmri_model} object.
#' @param dataset An \code{fmri_dataset} object containing the time-series data.
#' @param strategy The data splitting strategy, either \code{"runwise"} or \code{"chunkwise"}. Default is \code{"runwise"}.
#' @param robust Logical. Whether to use robust fitting. Default is \code{FALSE}.
#' @param nchunks Number of data chunks when strategy is \code{"chunkwise"}. Default is \code{10}.
#' @param use_fast_path Logical. If \code{TRUE}, use matrix-based computation for speed. Default is \code{FALSE}.
#' @param progress Logical. Whether to display a progress bar during model fitting. Default is \code{FALSE}.
#' @param ... Additional arguments.
#' @return A fitted fMRI linear regression model with the specified fitting strategy.
#' @keywords internal
#' @seealso \code{\link{fmri_lm}}, \code{\link{fmri_model}}, \code{\link{fmri_dataset}}
fmri_lm_fit <- function(fmrimod, dataset, strategy = c("runwise", "chunkwise"),
                        robust = FALSE, nchunks = 10, use_fast_path = FALSE, progress = FALSE, ...) {
  strategy <- match.arg(strategy)
  
  # Error checking
  assert_that(inherits(fmrimod, "fmri_model"), msg = "'fmrimod' must be an 'fmri_model' object")
  assert_that(inherits(dataset, "fmri_dataset"), msg = "'dataset' must be an 'fmri_dataset' object")
  assert_that(is.logical(robust), msg = "'robust' must be logical")
  assert_that(is.logical(use_fast_path), msg = "'use_fast_path' must be logical")
  if (strategy == "chunkwise") {
    assert_that(is.numeric(nchunks) && nchunks > 0, msg = "'nchunks' must be a positive number")
  }
  if (robust && use_fast_path) {
      # Warning already issued in fmri_lm, but double check
      warning("Robust fitting ('robust=TRUE') is not currently implemented with 'use_fast_path=TRUE'. Ignoring 'robust=TRUE'.")
      robust <- FALSE
  }
  
  # Get contrast info grouped by term
  contrast_info_by_term <- contrast_weights(fmrimod$event_model)
  full_design_colnames <- colnames(design_matrix(fmrimod))
  processed_conlist <- list()
  
  # Process contrasts term by term to correctly assign column indices
  for (term_name in names(contrast_info_by_term)) {
      term_contrasts <- contrast_info_by_term[[term_name]]
      
      # Get col_indices from design matrix instead of term_indices from event model
      col_indices <- attr(fmrimod$event_model$design_matrix, "col_indices")
      
      if (length(term_contrasts) > 0 && !is.null(col_indices) && !is.null(col_indices[[term_name]])) {
          # Get the column indices directly from col_indices instead of trying to match names
          # The col_indices already contains the correct indices for this term
          colind <- col_indices[[term_name]]
          
          if (length(colind) == 0) {
              warning(paste("No column indices found for term:", term_name))
              next # Skip contrasts for this term if columns can't be found
          }
          
          # Apply colind attribute to each contrast spec within this term
          processed_term_contrasts <- lapply(term_contrasts, function(con_spec) {
              if (inherits(con_spec, "contrast") || inherits(con_spec, "Fcontrast")) {
                  # Set the colind attribute on the contrast weights for the slow path
                  attr(con_spec$weights, "colind") <- colind
                  # Also set it directly on the contrast object for estimate_contrast
                  attr(con_spec, "colind") <- colind
                  # Add term name attribute for potential future use/debugging
                  # attr(con_spec$weights, "term") <- term_name 
              } else {
                  warning(paste("Item in contrast list for term", term_name, "is not a contrast or Fcontrast object."))
              }
              con_spec # Return modified or original con_spec
          })
          processed_conlist <- c(processed_conlist, processed_term_contrasts)
      } else if (length(term_contrasts) > 0 && (is.null(col_indices) || is.null(col_indices[[term_name]]))) {
           warning(paste("Contrasts found for term '", term_name, "' but col_indices are missing in the event model design matrix."))
      }
  }

  # Now processed_conlist contains all valid contrasts with the colind attribute added
  
  # Separate simple and F contrasts (full objects) for the standard path
  simple_conlist_objects <- Filter(function(x) inherits(x, "contrast"), processed_conlist)
  fconlist_objects <- Filter(function(x) inherits(x, "Fcontrast"), processed_conlist)
  # Combine for standard path (fit_lm_contrasts expects a single list)
  standard_path_conlist <- c(simple_conlist_objects, fconlist_objects)
  
  # Pass the full processed contrast objects list down.
  # The fitting function (chunkwise/runwise) will decide whether to use the objects (slow path)
  # or extract weights (fast path).
  result <- switch(strategy,
                   "runwise" = runwise_lm(dataset, fmrimod, standard_path_conlist, # Pass full objects
                                         robust = robust, use_fast_path = use_fast_path,
                                         progress = progress, ...),
                   "chunkwise" = {
                     if (inherits(dataset, "latent_dataset")) {
                       chunkwise_lm(dataset, fmrimod, standard_path_conlist, # Pass full objects
                                   nchunks, robust = robust, progress = progress, ...) # Do not pass use_fast_path
                     } else {
                       chunkwise_lm(dataset, fmrimod, standard_path_conlist, # Pass full objects
                                   nchunks, robust = robust, use_fast_path = use_fast_path,
                                   progress = progress, ...)
                     }
                   })
  
  ret <- list(
    result = result,
    model = fmrimod,
    strategy = strategy,
    bcons = processed_conlist,
    dataset = dataset
  )
  
  class(ret) <- "fmri_lm"
  
  return(ret)
}


#' Compute Fitted Hemodynamic Response Functions for an fmri_lm Object
#'
#' This method computes the fitted hemodynamic response functions (HRFs) for an \code{fmri_lm} object.
#'
#' @param x An \code{fmri_lm} object for which the fitted HRFs should be computed.
#' @param sample_at A numeric vector of time points at which the HRFs should be sampled. Default is \code{seq(0, 24, by = 1)}.
#' @param ... Additional arguments (currently unused).
#' @return A list where each element corresponds to an event term in the \code{fmri_lm} object. Each element contains:
#' \describe{
#'   \item{\code{pred}}{A matrix of predicted HRF values.}
#'   \item{\code{design}}{A tibble containing the design matrix for the HRFs.}
#' }
#' @export
fitted_hrf.fmri_lm <- function(x, sample_at = seq(0, 24, by = 1), ...) {
  # Error checking
  assert_that(inherits(x, "fmri_lm"), msg = "'x' must be an 'fmri_lm' object")
  assert_that(is.numeric(sample_at), msg = "'sample_at' must be numeric")
  
  eterms <- terms(x$model$event_model)
  betas <- coef(x)
  tind <- x$model$event_model$term_indices
  
  pred <- lapply(seq_along(tind), function(i) {
    ind <- tind[[i]]
    hrf_spec <- eterms[[i]]$hrfspec
    hrf <- hrf_spec$hrf
    nb <- attr(hrf, "nbasis")
    G <- as.matrix(hrf(sample_at))
    
    excond <- cells(eterms[[i]], exclude_basis = TRUE)
    ncond <- nrow(excond)
    Gex <- do.call(Matrix::bdiag, replicate(ncond, G, simplify = FALSE))
    
    B <- t(betas[, ind, drop = FALSE])
    yh <- Gex %*% B
    
    excond_expanded <- excond %>% dplyr::slice(rep(1:dplyr::n(), each = length(sample_at)))
    design <- cbind(
      dplyr::tibble(time = rep(sample_at, ncond)),
      excond_expanded
    )
    
    list(pred = as.matrix(yh), design = tibble::as_tibble(design))
  })
  
  names(pred) <- names(eterms)
  return(pred)
}



#' Reshape Coefficient Data
#'
#' This function reshapes coefficient data from wide to long format and merges it with design information.
#'
#' @param df A data frame containing coefficient estimates.
#' @param des A data frame containing design information.
#' @param measure The name of the value column in the reshaped data. Default is \code{"value"}.
#' @return A data frame in long format with merged design information.
#' @keywords internal
#' @autoglobal
reshape_coef <- function(df, des, measure = "value") {
  # Create a unique identifier for each row
  df <- df %>% dplyr::mutate(row_id = dplyr::row_number())
  
  des <- des %>% dplyr::mutate(key = do.call(paste, c(.[, colnames(des)], sep = ":")))
  
  colnames(df)[-ncol(df)] <- des$key  # assign new column names excluding the last column (row_id)
  
  df_long <- df %>%
    tidyr::pivot_longer(-row_id, names_to = "col_name", values_to = measure)
  
  # Match the long dataframe with the design dataframe
  df_long <- dplyr::left_join(df_long, des, by = c("col_name" = "key"))
  
  return(df_long)
}


#' Extract Statistical Measures from an fmri_lm Object
#'
#' This function extracts statistical measures (e.g., estimates, standard errors) from an \code{fmri_lm} object.
#'
#' @param x An \code{fmri_lm} object.
#' @param type The type of statistic to extract: \code{"betas"}, \code{"contrasts"}, or \code{"F"}.
#' @param element The specific element to extract, such as \code{"estimate"}, \code{"se"}, \code{"stat"}, or \code{"prob"}.
#' @return A tibble containing the requested statistical measures.
#' @keywords internal
pull_stat_revised <- function(x, type, element) {
  if (type == "betas") {
    # Ensure we access the matrix correctly from the list structure
    beta_matrix <- x$result$betas$data[[1]]$estimate[[1]]
    ret <- beta_matrix[, x$result$event_indices, drop = FALSE]
    colnames(ret) <- conditions(x$model$event_model)
    suppressMessages(tibble::as_tibble(ret, .name_repair = "check_unique"))
  } else if (type == "contrasts") {
    ret <- x$result$contrasts %>% dplyr::filter(type == "contrast")
    if (nrow(ret) == 0) {
      stop("No simple contrasts for this model.")
    }
    cnames <- ret$name
    # Extract the specific element (e.g., estimate), which is a list(vector)
    # Then extract the vector itself (element [[1]]) before binding
    out <- lapply(ret$data, function(inner_tibble) inner_tibble[[element]][[1]]) %>% 
             dplyr::bind_cols()
    names(out) <- cnames
    out
  } else if (type == "F") {
    ret <- x$result$contrasts %>% dplyr::filter(type == "Fcontrast")
    if (nrow(ret) == 0) {
      stop("No F contrasts for this model.")
    }
    cnames <- ret$name
    # Extract the specific element (e.g., estimate), which is list(vector)
    # Then extract the vector itself (element [[1]]) before binding
    out <- lapply(ret$data, function(inner_tibble) inner_tibble[[element]][[1]]) %>% 
             dplyr::bind_cols()
    names(out) <- cnames
    out
  } else {
    stop("Invalid type specified. Must be 'betas', 'contrasts', or 'F'.")
  }
}

pull_stat <- function(x, type, element) {
  if (type == "betas") {
    ret <- x$result$betas$data[[1]][[element]][[1]]
    
    # Check bounds and filter valid indices
    max_col <- ncol(ret)
    valid_event_indices <- x$result$event_indices[x$result$event_indices <= max_col]
    
    if (length(valid_event_indices) == 0) {
      warning("No valid event indices found in pull_stat. Using all available columns.")
      valid_event_indices <- 1:max_col
    }
    
    ret <- ret[, valid_event_indices, drop = FALSE]
    
    # Use the actual column names from the design matrix instead of conditions()
    # This avoids duplicate names when multiple terms have the same variables
    dm <- design_matrix(x$model)
    if (!is.null(dm) && ncol(dm) >= max(valid_event_indices)) {
      actual_colnames <- colnames(dm)[valid_event_indices]
      colnames(ret) <- actual_colnames
    } else {
      # Fallback: use conditions but make them unique
      condition_names <- conditions(x$model$event_model)[1:length(valid_event_indices)]
      colnames(ret) <- make.names(condition_names, unique = TRUE)
    }
    
    # Ensure tibble output for consistency with original behavior
    res <- suppressMessages(tibble::as_tibble(ret, .name_repair = "check_unique"))
  } else if (type == "contrasts") {
    ret <- x$result$contrasts %>% dplyr::filter(type == "contrast")
    if (nrow(ret) == 0) {
      stop("No simple contrasts for this model.")
    }
    cnames <- ret$name
    out <- lapply(ret$data, function(x) x[[element]]) %>% dplyr::bind_cols()
    names(out) <- cnames
    out
  } else if (type == "F") {
    ret <- x$result$contrasts %>% dplyr::filter(type == "Fcontrast")
    if (nrow(ret) == 0) {
      stop("No F contrasts for this model.")
    }
    cnames <- ret$name
    out <- lapply(ret$data, function(x) x[[element]]) %>% dplyr::bind_cols()
    names(out) <- cnames
    out
  } else {
    stop("Invalid type specified. Must be 'betas', 'contrasts', or 'F'.")
  }
}

#' Extract Model Coefficients from an fmri_lm Object
#'
#' This function extracts model coefficients (estimates) from an `fmri_lm` object.
#'
#' @param object An `fmri_lm` object.
#' @param type The type of coefficients to extract: `"betas"` or `"contrasts"`. Default is `"betas"'.
#' @param include_baseline Logical. If `TRUE`, include coefficients for baseline regressors along with event regressors.
#'                         If `FALSE` (default), only event regressors are returned.
#' @param recon Logical. If `TRUE`, reconstructs the coefficients into a neuroimaging volume. Default is `FALSE`.
#' @param ... Additional arguments (currently unused).
#' @return A tibble or matrix of coefficients.
#' @export
coef.fmri_lm <- function(object, type = c("betas", "contrasts"), include_baseline = FALSE, recon = FALSE, ...) {
  type <- match.arg(type)
  
  if (type == "contrasts") {
    # Contrast handling remains the same
    res <- pull_stat(object, "contrasts", "estimate")
  } else if (type == "betas") {
    # Get all beta estimates first
    all_betas <- object$result$betas$data[[1]]$estimate[[1]]
    
    if (include_baseline) {
      # Return all betas, ensure correct names from the full design matrix
      res <- all_betas
      colnames(res) <- colnames(design_matrix(object$model))
      # Convert back to tibble for consistency if needed, though matrix might be better here
      # res <- as_tibble(res)
    } else {
      # Default: return only event betas
      # Check bounds and filter valid indices
      max_col <- ncol(all_betas)
      valid_event_indices <- object$result$event_indices[object$result$event_indices <= max_col]
      
      if (length(valid_event_indices) == 0) {
        warning("No valid event indices found in coef.fmri_lm. Using all available columns.")
        valid_event_indices <- 1:max_col
      }
      
      res <- all_betas[, valid_event_indices, drop = FALSE]
      
      # Use the actual column names from the design matrix instead of conditions()
      # This avoids duplicate names when multiple terms have the same variables
      dm <- design_matrix(object$model)
      if (!is.null(dm) && ncol(dm) >= max(valid_event_indices)) {
        actual_colnames <- colnames(dm)[valid_event_indices]
        colnames(res) <- actual_colnames
      } else {
        # Fallback: use conditions but make them unique
        condition_names <- conditions(object$model$event_model)[1:length(valid_event_indices)]
        colnames(res) <- make.names(condition_names, unique = TRUE)
      }
      
      # Ensure tibble output for consistency with original behavior
      res <- suppressMessages(tibble::as_tibble(res, .name_repair = "check_unique"))
    }
  } else {
    # Should not happen due to match.arg, but defensive coding
    stop("Invalid type specified.")
  }
  
  # Reconstruction functionality can be added here if necessary (applies to the 'res' matrix/tibble)
  # if (recon && inherits(object$dataset, "fmri_dataset")) { ... }
  
  return(res)
}

#' Extract Statistical Values from an fmri_lm Object
#'
#' This function extracts statistical values (e.g., t-statistics, F-statistics) from an \code{fmri_lm} object.
#'
#' @param x An \code{fmri_lm} object.
#' @param type The type of statistics to extract: \code{"estimates"}, \code{"contrasts"}, or \code{"F"}.
#' @param ... Additional arguments (currently unused).
#' @return A tibble containing the requested statistical values.
#' @export
stats.fmri_lm <- function(x, type = c("estimates", "contrasts", "F"), ...) {
  type <- match.arg(type)
  if (type == "estimates") {
    pull_stat(x, "betas", "stat")
  } else if (type == "contrasts") {
    pull_stat(x, "contrasts", "stat")
  } else if (type == "F") {
    pull_stat(x, "F", "stat")
  }
}

#' Extract Standard Errors from an fmri_lm Object
#'
#' This function extracts standard errors from an \code{fmri_lm} object.
#'
#' @rdname standard_error
#'
#' @param x An \code{fmri_lm} object.
#' @param type The type of standard errors to extract: \code{"estimates"} or
#'   \code{"contrasts"}.
#' @return A tibble containing the standard errors.
#' @export
standard_error.fmri_lm <- function(x, type = c("estimates", "contrasts")) {
  type <- match.arg(type)
  if (type == "estimates") {
    pull_stat(x, "betas", "se")
  } else if (type == "contrasts") {
    pull_stat(x, "contrasts", "se")
  }
}



#' Fit Linear Model Contrasts
#'
#' This function computes contrasts and beta statistics for a fitted linear model.
#'
#' @param fit A fitted linear model object.
#' @param conlist A list of contrast matrices.
#' @param fcon A list of F-contrasts.
#' @param vnames Variable names corresponding to the model coefficients.
#' @param se Logical. Whether to compute standard errors. Default is \code{TRUE}.
#' @return A list containing contrasts, beta statistics, and the fitted model.
#' @keywords internal
fit_lm_contrasts <- function(fit, conlist, fcon, vnames, se = TRUE) {
  conres <- if (!is.null(conlist)) {
    ret <- lapply(conlist, function(con) {
      # Extract colind from the contrast object's attributes
      colind <- attr(con, "colind")
      if (is.null(colind)) {
        warning(paste("Missing colind attribute for contrast:", con$name %||% "unnamed"))
        return(NULL) # Skip this contrast
      }
      estimate_contrast(con, fit, colind)
    })
    # Filter out NULL results
    ret <- ret[!sapply(ret, is.null)]
    names(ret) <- sapply(conlist[!sapply(ret, is.null)], function(x) x$name %||% "unnamed")
    ret
  } else {
    list()
  }
  
  bstats <- beta_stats(fit, vnames, se = se)
  list(contrasts = conres, bstats = bstats, fit = fit)
}




#' Fit Multiresponse Linear Model
#'
#' This function fits a linear model to multiple responses in an fMRI dataset.
#'
#' @param form The formula used to define the linear model.
#' @param data_env The environment containing the data to be used in the linear model.
#' @param conlist The list of contrasts used in the analysis.
#' @param vnames The names of the variables used in the linear model.
#' @param fcon The F-contrasts used in the analysis.
#' @param modmat The model matrix (default is \code{NULL}, which will calculate the model matrix using the formula).
#' @return A list containing the results from the multiresponse linear model analysis.
#' @keywords internal
multiresponse_lm <- function(form, data_env, conlist, vnames, fcon, modmat = NULL) {
  lm_fit <- if (is.null(modmat)) {
    lm(as.formula(form), data = data_env)
  } else {
    lm.fit(modmat, data_env$.y)
  }
  
  # Use the actual column names from the model matrix instead of vnames
  # This ensures the dimensions match correctly
  actual_vnames <- if (is.null(modmat)) {
    names(coef(lm_fit))
  } else {
    colnames(modmat)
  }
  
  fit_lm_contrasts(lm_fit, conlist, fcon, actual_vnames)
}





unpack_chunkwise <- function(cres, event_indices, baseline_indices) {
  # --- Beta Processing (Seems OK) ---
  cbetas <- lapply(cres, function(x) x$bstats) %>% dplyr::bind_rows()
  dat_beta <- cbetas$data %>% dplyr::bind_rows()
  
  # Check validity (assuming estimate column now correctly holds matrices)
  valid_estimates_idx <- sapply(dat_beta$estimate, function(x) !is.null(x) && is.matrix(x) && nrow(x) > 0)
  if (!any(valid_estimates_idx)) {
      stop("No valid beta estimates found across chunks in unpack_chunkwise.")
  }
  dat_beta_valid <- dat_beta[valid_estimates_idx, , drop = FALSE]

  # Concatenate beta results across chunks
  estimate_beta <- do.call(rbind, dat_beta_valid$estimate)
  se_beta <- do.call(rbind, dat_beta_valid$se)
  stat_beta <- do.call(rbind, dat_beta_valid$stat)
  prob_beta <- do.call(rbind, dat_beta_valid$prob)
  sigma_beta <- do.call(c, dat_beta_valid$sigma)
  
  # Re-package combined beta results
  cbetas_out <- dplyr::tibble(
    type = cbetas$type[1],
    stat_type = cbetas$stat_type[1],
    df.residual = cbetas$df.residual[1],
    conmat = list(NULL),
    colind = list(NULL),
    data = list(
      dplyr::tibble(
        estimate = list(estimate_beta),   
        se = list(se_beta),               
        stat = list(stat_beta),            
        prob = list(prob_beta),            
        sigma = list(sigma_beta)          
      )
    )
  )

  # --- Contrast Processing --- 
  ncon <- if (length(cres) > 0 && !is.null(cres[[1]]$contrasts) && length(cres[[1]]$contrasts) > 0) {
      length(cres[[1]]$contrasts)
  } else { 0 }

  if (ncon > 0) {
    contab <- lapply(cres, function(x) { 
        # Ensure contrasts is a list of tibbles, even if only one contrast
        cons <- x$contrasts 
        if (!is.list(cons)) cons <- list(cons) # Handle single contrast case if needed
        if (length(cons) > 0 && !is.null(names(cons))) { # Ensure names exist
             dplyr::bind_rows(cons, .id = "contrast_internal_name") # Requires names
        } else if (length(cons) > 0) {
             # Fallback if names are missing, might need adjustment based on actual structure
             warning("Contrast list per chunk lacks names, attempting bind_rows without .id")
             dplyr::bind_rows(cons)
        } else {
             dplyr::tibble() # Return empty tibble for chunks with no contrasts
        }
    }) %>% dplyr::bind_rows() # Bind results from all chunks

    # Check if contab is empty after binding
    if (nrow(contab) == 0) {
        con <- dplyr::tibble()
    } else {
        # Group by original contrast name and type
        # Use 'name' column if it exists, otherwise fallback might be needed
        grouping_vars <- intersect(c("name", "type"), names(contab))
        if (length(grouping_vars) == 0) stop("Cannot group contrasts: 'name' or 'type' column missing.")
        
        gsplit <- contab %>% dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>% dplyr::group_split()

        # Process each contrast group (combine results across chunks)
        con <- lapply(gsplit, function(g) {
            dat <- g$data %>% dplyr::bind_rows() 
            
            # Both paths now produce the same structure, so no conditional logic needed.
            # Simply assign the vectors directly.
            estimate_full <- dat$estimate
            se_full <- dat$se
            stat_full <- dat$stat
            prob_full <- dat$prob
            sigma_full <- if ("sigma" %in% names(dat)) dat$sigma else NULL
            
            # Re-package combined data for this contrast
            combined_data_tibble <- dplyr::tibble(
                estimate = list(estimate_full), 
                se = list(se_full),             
                stat = list(stat_full),          
                prob = list(prob_full)           
            )
            if (!is.null(sigma_full)) {
                combined_data_tibble$sigma = list(sigma_full)
            }

            # Take metadata from the first chunk's entry for this contrast
            g %>% dplyr::select(-data) %>% dplyr::slice_head() %>% 
                dplyr::mutate(data = list(combined_data_tibble))
                
        }) %>% dplyr::bind_rows() 
    }
  } else {
    con <- dplyr::tibble() # Return empty tibble if no contrasts
  }

  # --- DEBUG FINAL CONTRAST TIBBLE ---
  # message("Structure of final 'con' tibble before returning from unpack_chunkwise:")
  # print(str(con))
  # --- END DEBUG ---

  list(
    betas = cbetas_out,
    contrasts = con,
    event_indices = event_indices,
    baseline_indices = baseline_indices
  )
}




#' Perform Chunkwise Linear Modeling on fMRI Dataset
#'
#' This function performs a chunkwise linear model analysis on an fMRI dataset,
#' splitting the dataset into chunks and running the linear model on each chunk.
#'
#' @param dset An \code{fmri_dataset} object.
#' @param model The \code{fmri_model} used for the analysis.
#' @param contrast_objects The list of full contrast objects.
#' @param nchunks The number of chunks to divide the dataset into.
#' @param robust Logical. Whether to use robust linear modeling (default is \code{FALSE}).
#' @param verbose Logical. Whether to display progress messages (default is \code{FALSE}).
#' @param use_fast_path Logical. If \code{TRUE}, use matrix-based computation for speed. Default is \code{FALSE}.
#' @param progress Logical. Display a progress bar for chunk processing. Default is \code{FALSE}.
#' @return A list containing the unpacked chunkwise results.
#' @keywords internal
chunkwise_lm.fmri_dataset <- function(dset, model, contrast_objects, nchunks, robust = FALSE,
                                      verbose = FALSE, use_fast_path = FALSE, progress = FALSE) {
  chunks <- exec_strategy("chunkwise", nchunks = nchunks)(dset)
  if (progress) {
    pb <- cli::cli_progress_bar("Fitting chunks", total = length(chunks), clear = FALSE)
    on.exit(cli::cli_progress_done(id = pb), add = TRUE)
  }
  form <- get_formula(model)
  tmats <- term_matrices(model)
  vnames <- attr(tmats, "varnames")
  event_indices = attr(tmats, "event_term_indices")
  baseline_indices = attr(tmats, "baseline_term_indices")

  # Common setup for both paths
  ym <- NULL # Define ym for R CMD check

  
  
  if (!use_fast_path) {
   
      # -------- Original Slow Path --------
      # Slow path uses lmfun which calls fit_lm_contrasts, expects full contrast objects
      # contrast_objects should already be the correct list structure here
      data_env <- list2env(tmats)
      data_env[[ ".y"]] <- rep(0, nrow(tmats[[1]])) # Corrected [[ ]] indexing
      modmat <- model.matrix(as.formula(form), data_env)
      Qr_global <- qr(modmat)
      Vu <- chol2inv(Qr_global$qr)
      
      lmfun <- if (robust) multiresponse_rlm else multiresponse_lm
      
      cres <- vector("list", length(chunks))
      for (i in seq_along(chunks)) {
        ym <- chunks[[i]]
        if (verbose) message("Processing chunk ", ym$chunk_num)
        data_env[[".y"]] <- as.matrix(ym$data)

        ret <- lmfun(form, data_env, contrast_objects, vnames, fcon = NULL, modmat = modmat)

        rss <- colSums(as.matrix(ret$fit$residuals^2))
        rdf <- ret$fit$df.residual
        resvar <- rss / rdf
        sigma <- sqrt(resvar)

        cres[[i]] <- list(bstats = ret$bstats, contrasts = ret$contrasts,
                          rss = rss, rdf = rdf, sigma = sigma)
        if (progress) cli::cli_progress_update(id = pb)
      }

  } else {
      # -------- New Fast Path --------
      # Fast path needs weights extracted from contrast_objects
      simple_conlist <- Filter(function(x) inherits(x, "contrast"), contrast_objects)
      fconlist <- Filter(function(x) inherits(x, "Fcontrast"), contrast_objects)
      simple_conlist_weights <- lapply(simple_conlist, `[[`, "weights")
      names(simple_conlist_weights) <- names(simple_conlist) 
      fconlist_weights <- lapply(fconlist, `[[`, "weights")
      names(fconlist_weights) <- names(fconlist) 
      
      if (robust) {
          warning("Robust fitting not implemented for fast path, using standard OLS.")
          robust <- FALSE
      }
      
      message("Using fast path for chunkwise LM...")
      
      data_env <- list2env(tmats)
      data_env[[ ".y"]] <- rep(0, nrow(tmats[[1]])) # Placeholder for model.matrix
      modmat <- model.matrix(as.formula(form), data_env)
      
      proj <- .fast_preproject(modmat)
      Vu <- proj$XtXinv 
      
      cres <- vector("list", length(chunks))
      for (i in seq_along(chunks)) {
          ym <- chunks[[i]]
          if (verbose) message("Processing chunk (fast path) ", ym$chunk_num)
          Ymat <- as.matrix(ym$data)
          if (verbose) message("  Chunk ", ym$chunk_num, ": ncol(Ymat) = ", ncol(Ymat))

          res <- .fast_lm_matrix(Ymat, proj)

          actual_vnames <- colnames(modmat)
          bstats <- beta_stats_matrix(res$betas, proj$XtXinv, res$sigma, proj$dfres, actual_vnames)

          contrasts <- fit_lm_contrasts_fast(res$betas, res$sigma2, proj$XtXinv,
                                             simple_conlist_weights, fconlist_weights, proj$dfres)

          cres[[i]] <- list(bstats = bstats,
                            contrasts = contrasts,
                            rss = res$rss,
                            rdf = proj$dfres,
                            sigma = res$sigma)
          if (progress) cli::cli_progress_update(id = pb)
      }
      # -------- End New Fast Path --------
  }
  
  # Unpack results (expects specific structure from cres)
  out <- unpack_chunkwise(cres, event_indices, baseline_indices)
  # Add cov.unscaled to the output
  out$cov.unscaled <- Vu 
  out
}



#' Perform Runwise Linear Modeling on fMRI Dataset
#'
#' This function performs a runwise linear model analysis on an fMRI dataset by
#' running the linear model for each data run and combining the results.
#'
#' @param dset An \code{fmri_dataset} object.
#' @param model The \code{fmri_model} used for the analysis.
#' @param contrast_objects The list of full contrast objects.
#' @param robust Logical. Whether to use robust linear modeling (default is \code{FALSE}).
#' @param verbose Logical. Whether to display progress messages (default is \code{FALSE}).
#' @param progress Logical. Display a progress bar for run processing. Default is \code{FALSE}.
#' @return A list containing the combined results from runwise linear model analysis.
#' @keywords internal
#' @autoglobal
runwise_lm <- function(dset, model, contrast_objects, robust = FALSE, verbose = FALSE,
                       use_fast_path = FALSE, progress = FALSE) {
  # Get an iterator of data chunks (runs)
  chunks <- exec_strategy("runwise")(dset)
  if (progress) {
    pb <- cli::cli_progress_bar("Fitting runs", total = length(chunks), clear = FALSE)
    on.exit(cli::cli_progress_done(id = pb), add = TRUE)
  }
  form <- get_formula(model)
  # Global design matrix needed for pooling compatibility? Or just for Vu?
  modmat_global <- design_matrix(model)
  Qr_global <- qr(modmat_global)
  Vu <- chol2inv(Qr_global$qr)
  
  # Define ym for R CMD check
  ym <- NULL
  
  if (!use_fast_path) {
      # -------- Original Slow Path --------
      # Slow path uses lmfun which calls fit_lm_contrasts, expects full contrast objects
      lmfun <- if (robust) multiresponse_rlm else multiresponse_lm
      
      cres <- vector("list", length(chunks))
      for (i in seq_along(chunks)) {
        ym <- chunks[[i]]
        if (verbose) message("Processing run ", ym$chunk_num)
        tmats <- term_matrices(model, ym$chunk_num)
        vnames <- attr(tmats, "varnames")
        event_indices <- attr(tmats, "event_term_indices")
        baseline_indices <- attr(tmats, "baseline_term_indices")

        data_env <- list2env(tmats)
        data_env$.y <- as.matrix(ym$data)
        ret <- lmfun(form, data_env, contrast_objects, vnames, fcon = NULL)

        rss <- colSums(as.matrix(ret$fit$residuals^2))
        rdf <- ret$fit$df.residual
        resvar <- rss / rdf
        sigma <- sqrt(resvar)

        cres[[i]] <- list(
          conres = ret$contrasts,
          bstats = ret$bstats,
          event_indices = event_indices,
          baseline_indices = baseline_indices,
          rss = rss,
          rdf = rdf,
          resvar = resvar,
          sigma = sigma
        )
        if (progress) cli::cli_progress_update(id = pb)
      }
      # -------- End Original Slow Path --------
      
  } else {
      # -------- New Fast Path --------
      # Fast path needs weights extracted from contrast_objects
      simple_conlist <- Filter(function(x) inherits(x, "contrast"), contrast_objects)
      fconlist <- Filter(function(x) inherits(x, "Fcontrast"), contrast_objects)
      simple_conlist_weights <- lapply(simple_conlist, `[[`, "weights")
      names(simple_conlist_weights) <- names(simple_conlist) 
      fconlist_weights <- lapply(fconlist, `[[`, "weights")
      names(fconlist_weights) <- names(fconlist) 
      
      if (robust) {
          warning("Robust fitting not implemented for fast path, using standard OLS.")
          robust <- FALSE
      }
      
      message("Using fast path for runwise LM...")
      
      # .export needed? conlist, fcon, model should be available.
      # Add functions from this package? .packages = c("dplyr", "purrr", "fmrireg")? Or rely on namespace?
      cres <- vector("list", length(chunks))
      for (i in seq_along(chunks)) {
        ym <- chunks[[i]]
        if (verbose) message("Processing run (fast path) ", ym$chunk_num)

        tmats <- term_matrices(model, ym$chunk_num)
        vnames <- attr(tmats, "varnames")
        event_indices <- attr(tmats, "event_term_indices")
        baseline_indices <- attr(tmats, "baseline_term_indices")

        data_env_run <- list2env(tmats)
        n_timepoints_run <- nrow(tmats[[1]])
        if (n_timepoints_run == 0) {
            warning(paste("Skipping empty run:", ym$chunk_num))
            next
        }
        data_env_run[[".y"]] <- rep(0, n_timepoints_run)
        X_run <- model.matrix(form, data_env_run)

        proj_run <- .fast_preproject(X_run)

        Y_run <- as.matrix(ym$data)

        if (nrow(X_run) != nrow(Y_run)) {
            stop(paste("Dimension mismatch in run", ym$chunk_num, ": X_run rows (", nrow(X_run), ") != Y_run rows (", nrow(Y_run), ")"))
        }

        res <- .fast_lm_matrix(Y_run, proj_run)

        actual_vnames <- colnames(X_run)
        bstats <- beta_stats_matrix(res$betas, proj_run$XtXinv, res$sigma, proj_run$dfres, actual_vnames)

        conres <- fit_lm_contrasts_fast(res$betas, res$sigma2, proj_run$XtXinv,
                                         simple_conlist_weights, fconlist_weights, proj_run$dfres)

        cres[[i]] <- list(
          conres = conres,
          bstats = bstats,
          event_indices = event_indices,
          baseline_indices = baseline_indices,
          rss = res$rss,
          rdf = proj_run$dfres,
          resvar = res$sigma2,
          sigma = res$sigma
        )
        if (progress) cli::cli_progress_update(id = pb)
      }
      
      # Filter out NULL results from skipped empty runs
      cres <- Filter(Negate(is.null), cres)
      if (length(cres) == 0) {
          stop("No valid run results found in runwise fast path.")
      }
      # -------- End New Fast Path --------
  }
  
  # Combine results (Pooling logic assumes specific structure in cres[[i]]$bstats and cres[[i]]$conres)
  bstats_list <- lapply(cres, `[[`, "bstats")
  conres_list <- lapply(cres, `[[`, "conres")
  
  # Compute overall statistics (these seem independent of fast/slow path)
  sigma <- colMeans(do.call(rbind, lapply(cres, function(x) as.matrix(x$sigma)))) # Make sure sigma is matrix/vector
  rss <- colSums(do.call(rbind, lapply(cres, function(x) as.matrix(x$rss))))
  rdf <- sum(unlist(lapply(cres, `[[`, "rdf")))
  resvar <- rss / rdf # Overall residual variance
  
  # Pool over runs
  if (length(cres) > 1) {
    # meta_contrasts expects a list of lists (runs) of lists (contrasts) of tibbles?
    # Or list (runs) of lists (contrasts) where elements are the tibbles?
    # Current: conres_list is list (runs) of lists (contrasts are named elements, values are tibbles)
    # Need to check meta_contrasts implementation.
    # Assuming meta_contrasts can handle the list of lists structure from fit_lm_contrasts_fast.
    
    # meta_betas expects a list of bstats tibbles and event_indices from the first run.
    # Assuming beta_stats_matrix output is compatible.
    meta_con <- meta_contrasts(conres_list)
    meta_beta <- meta_betas(bstats_list, cres[[1]]$event_indices)
    
    list(
      contrasts = meta_con,
      betas = meta_beta,
      event_indices = cres[[1]]$event_indices,
      baseline_indices = cres[[1]]$baseline_indices,
      cov.unscaled = Vu, # Using Vu from global design matrix
      sigma = sigma, # Pooled sigma
      rss = rss,     # Pooled rss
      rdf = rdf,     # Pooled rdf
      resvar = resvar # Pooled resvar
    )
  } else {
    # If only one run, return its results directly
    list(
      contrasts = conres_list[[1]], # This is the list of contrast tibbles for the single run
      betas = bstats_list[[1]], # This is the bstats tibble for the single run
      event_indices = cres[[1]]$event_indices,
      baseline_indices = cres[[1]]$baseline_indices,
      cov.unscaled = Vu,
      sigma = cres[[1]]$sigma, # Use run sigma
      rss = cres[[1]]$rss,
      rdf = cres[[1]]$rdf,
      resvar = cres[[1]]$resvar
    )
  }
}


#' Print an fmri_lm_result object
#'
#' Provides a colorful and informative printout.
#'
#' @param x An fmri_lm_result object.
#' @param ... Additional arguments (unused).
#' @export
#' @rdname print
print.fmri_lm <- function(x, ...) {
  # optional: check if crayon is installed
  if (!requireNamespace("crayon", quietly = TRUE)) {
    # fallback to standard cat if crayon is missing
    cat("fmri_lm_result object (install 'crayon' for color)\n\n")
    
    cat("Model formula:\n",
        as.character(x$model$event_model$model_spec$formula), "\n")
    cat("Strategy: ", x$strategy, "\n")
    cat("Baseline parameters: ",
        ncol(design_matrix(x$model$baseline_model)), "\n")
    cat("Design parameters: ",
        ncol(design_matrix(x$model$event_model)), "\n")
    cat("Contrasts: ",
        paste(names(x$bcons), collapse = ", "), "\n\n")
    return(invisible(x))
  }
  
  # If we do have crayon, let's color it up:
  cat(crayon::blue$bold("\n╔════════════════════════════════╗\n"))
  cat(crayon::blue$bold("║        fmri_lm_result          ║\n"))
  cat(crayon::blue$bold("╚════════════════════════════════╝\n\n"))
  
  # Print the model formula
  cat(crayon::green("Model formula:\n  "))
  cat(crayon::silver(as.character(x$model$event_model$model_spec$formula)), "\n\n")
  
  # Print strategy
  cat(crayon::green("Fitting strategy:  "))
  cat(crayon::silver(x$strategy), "\n\n")
  
  # Some stats about baseline, design, and contrasts
  bdim <- ncol(design_matrix(x$model$baseline_model))
  ddim <- ncol(design_matrix(x$model$event_model))
  
  cat(crayon::green("Baseline parameters: "), crayon::silver(bdim), "\n")
  cat(crayon::green("Design parameters:   "), crayon::silver(ddim), "\n")
  
  # If you have some # of simple contrasts
  c_names <- names(x$bcons)
  if (length(c_names) > 0) {
    cat(crayon::green("Contrasts:          "), crayon::silver(paste(c_names, collapse = ", ")), "\n\n")
  } else {
    cat(crayon::green("Contrasts:          "), crayon::silver("None\n\n"))
  }
  
  cat(crayon::yellow("Use coef(...), stats(...), etc. to extract results.\n\n"))
  
  invisible(x)
}
  
    
    



