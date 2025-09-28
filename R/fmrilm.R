#' @method get_formula fmri_model
#' @rdname get_formula
#' @export
get_formula.fmri_model <- function(x,...) {
  assert_that(inherits(x, "fmri_model"), msg = "'x' must be an 'fmri_model' object")
  term_names <- names(terms(x))
  form <- paste(".y ~", paste(term_names, collapse = " + "), "-1")
  return(as.formula(form))
}



#' @keywords internal
#' @noRd
is.formula <- function(x) {
  inherits(x, "formula")
}

#' @keywords internal
#' @noRd
.fast_preproject <- function(X) {
  # Ensure X is a numeric matrix
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  # Validate no missing or non-finite values (tests expect an error on NA input)
  if (anyNA(X) || any(!is.finite(X))) {
    stop(".fast_preproject: design matrix contains NA/Inf values", call. = FALSE)
  }
  # Ensure the matrix is numeric
  if (!is.numeric(X)) {
    orig_colnames <- colnames(X)
    X <- matrix(as.numeric(X), nrow = nrow(X), ncol = ncol(X))
    # Preserve column names if they exist
    if (!is.null(orig_colnames)) {
      colnames(X) <- orig_colnames
    }
  }

  qr_decomp <- qr(X, LAPACK = TRUE)
  rank <- qr_decomp$rank
  p <- ncol(X)
  n <- nrow(X)

  if (rank == p) {
    XtX <- crossprod(X)
    Rchol <- tryCatch(chol(XtX),
                      error = function(e) chol(XtX + diag(ncol(XtX)) * 1e-8))
    XtXinv <- chol2inv(Rchol)
    Pinv <- XtXinv %*% t(X)
  } else {
    svd_result <- svd(X)
    d <- svd_result$d
    tol <- max(dim(X)) * .Machine$double.eps * max(d)
    pos <- d > tol
    U <- svd_result$u[, pos, drop = FALSE]
    V <- svd_result$v[, pos, drop = FALSE]
    D_inv <- diag(1/d[pos], nrow = sum(pos))
    Pinv <- V %*% D_inv %*% t(U)
    XtXinv <- V %*% D_inv^2 %*% t(V)
  }

  list(
    qr = qr_decomp,
    Pinv = Pinv,
    XtXinv = XtXinv,
    dfres = n - rank,
    rank = rank,
    is_full_rank = (rank == p)
  )
}

#' @keywords internal
#' @noRd
.fast_lm_matrix <- function(X, Y, proj, return_fitted = FALSE) {
  ctx <- glm_context(X = X, Y = Y, proj = proj)
  res <- solve_glm_core(ctx, return_fitted = return_fitted)
  res$sigma <- sqrt(res$sigma2)
  res
}

#' Fast row-wise robust regression for a single run
#'
#' Wrapper around robust_iterative_fitter for backward compatibility.
#' This function implements an IRLS algorithm using Huber or Tukey bisquare 
#' weights on the time-point residuals.
#'
#' @param X Design matrix (time points \eqn{\times} predictors)
#' @param Y Data matrix (time points \eqn{\times} voxels)
#' @param proj Preprojection list from \code{.fast_preproject(X)}
#' @param psi Psi function for weighting, either \code{"huber"} or
#'   \code{"bisquare"}
#' @param k_huber Tuning constant for Huber weights
#' @param c_tukey Tuning constant for Tukey bisquare weights
#' @param max_it Maximum number of IRLS iterations
#' @param sigma_fixed Optional fixed sigma value (for global scale estimation)
#' @keywords internal
#' @noRd
fast_rlm_run <- function(X, Y, proj,
                         psi = c("huber", "bisquare"),
                         k_huber = 1.345,
                         c_tukey = 4.685,
                         max_it = 2L,
                         sigma_fixed = NULL) {

  psi <- match.arg(psi)

  # Validate inputs: robust path should not accept NA/Inf
  if (anyNA(X) || any(!is.finite(X)) || anyNA(Y) || any(!is.finite(Y))) {
    stop("fast_rlm_run: X/Y contain NA/Inf values", call. = FALSE)
  }


  
  # Create initial GLM context
  if (missing(proj) || is.null(proj)) {
    proj <- .fast_preproject(X)
  }
  
  initial_ctx <- glm_context(X = X, Y = Y, proj = proj)
  
  # Create robust options
  cfg_robust_options <- list(
    type = psi,
    k_huber = k_huber,
    c_tukey = c_tukey,
    max_iter = max_it,
    scale_scope = if (is.null(sigma_fixed)) "run" else "global"
  )
  
  # Call robust_iterative_fitter
  result <- robust_iterative_fitter(
    initial_glm_ctx = initial_ctx,
    cfg_robust_options = cfg_robust_options,
    X_orig_for_resid = X,
    sigma_fixed = sigma_fixed
  )
  
  # Calculate standard errors
  se_beta <- sqrt(diag(result$XtWXi_final)) * result$sigma_robust_scale_final
  
  # Return in the expected format
  list(
    betas = result$betas_robust,
    se = se_beta,
    sigma = result$sigma_robust_scale_final,
    dfres = result$dfres,
    XtXinv = result$XtWXi_final,
    weights = result$robust_weights_final
  )
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
      degree = max(ceiling(median(fmrihrf::blocklens(dataset$sampling_frame)) / 100), 3),
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
  
  fmri_model(ev_model, baseline_model, dataset)
}



#' Fit a Linear Regression Model for fMRI Data Analysis
#'
#' `fmri_lm` is a generic for fitting fMRI regression models. The
#' default interface accepts a model formula and dataset. An
#' alternative method can be used with a preconstructed
#' \code{fmri_model} object that already contains the design and data.
#'
#' @param formula A model formula describing the event structure or an
#'   \code{fmri_model} object.
#' @param ... Additional arguments passed to the chosen method.
#' @return An object of class \code{fmri_lm}.
#' @export
fmri_lm <- function(formula, ...) {
  UseMethod("fmri_lm")
}

#' @rdname fmri_lm
#' @export
#' @param block The model formula for block structure.
#' @param baseline_model (Optional) A \code{baseline_model} object. Default is \code{NULL}.
#' @param dataset An \code{fmri_dataset} object containing the time-series data.
#' @param durations A vector of event durations. Default is \code{0}.
#' @param drop_empty Logical. Whether to remove factor levels with zero size. Default is \code{TRUE}.
#' @param robust Logical or character. Either \code{FALSE} (no robust fitting), 
#'   \code{TRUE} (use Huber), or one of \code{"huber"} or \code{"bisquare"}. Default is \code{FALSE}.
#' @param robust_options List of robust fitting options. See Details.
#' @param ar_options List of autoregressive modeling options. See Details.
#' @param strategy The data splitting strategy, either \code{"runwise"} or \code{"chunkwise"}. Default is \code{"runwise"}.
#' @param nchunks Number of data chunks when strategy is \code{"chunkwise"}. Default is \code{10}.
#' @param use_fast_path Logical. If \code{TRUE}, use matrix-based computation for speed. Default is \code{FALSE}.
#' @param progress Logical. Whether to display a progress bar during model fitting. Default is \code{FALSE}.
#' @param ar_voxelwise Logical. Estimate AR parameters voxel-wise (overrides \code{ar_options$voxelwise}).
#' @param parallel_voxels Logical. Parallelize across voxels where supported.
#' @param cor_struct Character. Shorthand for \code{ar_options$struct} (e.g., "ar1", "ar2", "arp").
#' @param cor_iter Integer. Shorthand for \code{ar_options$iter_gls}.
#' @param cor_global Logical. Shorthand for \code{ar_options$global}.
#' @param ar1_exact_first Logical. Shorthand for \code{ar_options$exact_first}.
#' @param ar_p Integer. Shorthand for \code{ar_options$p}.
#' @param robust_psi Character. Shorthand for \code{robust_options$type} (e.g., "huber", "bisquare").
#' @param robust_max_iter Integer. Shorthand for \code{robust_options$max_iter}.
#' @param robust_scale_scope Character. Shorthand for \code{robust_options$scale_scope} ("run" or "global").
#' @return A fitted linear regression model for fMRI data analysis.
#' 
#' @details
#' \code{robust_options} may contain:
#' \itemize{
#'   \item \code{type}: Character or logical. Type of robust fitting (\code{FALSE}, \code{"huber"}, \code{"bisquare"})
#'   \item \code{k_huber}: Numeric. Tuning constant for Huber's psi (default: 1.345)
#'   \item \code{c_tukey}: Numeric. Tuning constant for Tukey's bisquare psi (default: 4.685)
#'   \item \code{max_iter}: Integer. Maximum IRLS iterations (default: 2)
#'   \item \code{scale_scope}: Character. Scope for scale estimation (\code{"run"} or \code{"global"})
#'   \item \code{reestimate_phi}: Logical. Whether to re-estimate AR parameters after robust fitting
#' }
#' 
#' \code{ar_options} may contain:
#' \itemize{
#'   \item \code{struct}: Character. Correlation structure (\code{"iid"}, \code{"ar1"}, \code{"ar2"}, \code{"arp"})
#'   \item \code{p}: Integer. AR order when \code{struct = "arp"}
#'   \item \code{iter_gls}: Integer. Number of GLS iterations (default: 1)
#'   \item \code{global}: Logical. Use global AR coefficients (default: FALSE)
#'   \item \code{voxelwise}: Logical. Estimate AR parameters voxel-wise (default: FALSE)
#'   \item \code{exact_first}: Logical. Apply exact AR(1) scaling to first sample (default: FALSE)
#' }
#' 
#' @export
#' @seealso \code{\link{fmri_dataset}}, \code{\link{fmri_lm_fit}}, \code{\link{fmri_lm_control}}
#' @examples
#' 
#' facedes <- subset(read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), 
#' header=TRUE), face_gen != "n/a")
#' facedes$face_gen <- droplevels(factor(facedes$face_gen))
#' sframe <- sampling_frame(rep(430/2,6), TR=2)
#' ev <- event_model(onset ~ hrf(face_gen, basis="gaussian"), data=facedes, 
#' block= ~ run, sampling_frame=sframe)
#' globonsets <- fmrihrf::global_onsets(sframe, facedes$onset, facedes$run)
#' reg1_signal <- regressor(globonsets[facedes$face_gen == "male"], hrf=fmrihrf::HRF_GAUSSIAN)
#' reg2_signal <- regressor(globonsets[facedes$face_gen == "female"], hrf=fmrihrf::HRF_GAUSSIAN)
#' time <- samples(sframe, global=TRUE)
#' y1 <- fmrihrf::evaluate(reg1_signal, time)*1.5
#' y2 <- fmrihrf::evaluate(reg2_signal, time)*3.0
#' y <- y1+y2
#' ys1 <- y + rnorm(length(y), sd=.02)
#' ys2 <- y + rnorm(length(y), sd=.02)
#' 
#' h <<- gen_hrf(fmrihrf::hrf_bspline, N=7, span=25)
#' dset <- matrix_dataset(cbind(ys1,ys2), TR=2, 
#'                        run_length=fmrihrf::blocklens(sframe), 
#'                        event_table=facedes)
#' flm <- fmri_lm(onset ~ hrf(face_gen, 
#'                            basis=gen_hrf(fmrihrf::hrf_bspline, N=7, span=25)), 
#'                block = ~ run, 
#'                strategy="chunkwise", nchunks=1, dataset=dset)
#' 
fmri_lm.formula <- function(formula, block, baseline_model = NULL, dataset, durations = 0, drop_empty = TRUE,
                         robust = FALSE, robust_options = NULL, ar_options = NULL,
                         strategy = c("runwise", "chunkwise"), nchunks = 10, use_fast_path = FALSE, progress = FALSE,
                         ar_voxelwise = FALSE,
                         parallel_voxels = FALSE,
                    # Individual AR parameters for backward compatibility
                    cor_struct = NULL, cor_iter = NULL, cor_global = NULL, 
                    ar1_exact_first = NULL, ar_p = NULL,
                    # Individual robust parameters for backward compatibility
                    robust_psi = NULL, robust_max_iter = NULL, robust_scale_scope = NULL,
                    ...) {

  # --- low-rank / sketch engine fast-path (opt-in via ...$engine) ---
  dots <- list(...)
  engine  <- dots$engine
  lowrank <- dots$lowrank
  engine_ar_options <- dots$ar_options
  engine_robust_options <- dots$robust_options
  engine_cfg <- dots$cfg
  # Remove recognised keys so downstream unexpected-arg check passes
  dots$engine <- NULL
  dots$lowrank <- NULL
  dots$ar_options <- NULL
  dots$robust_options <- NULL
  dots$cfg <- NULL

  strategy <- match.arg(strategy)
  
  # Check for any unexpected arguments in ...
  if (length(dots) > 0) {
    stop("Unexpected arguments: ", paste(names(dots), collapse = ", "))
  }
  
  # Error checking
  assert_that(is.formula(formula), msg = "'formula' must be a formula")
  assert_that(is.formula(block), msg = "'block' must be a formula")
  assert_that(inherits(dataset, "fmri_dataset"), msg = "'dataset' must be an 'fmri_dataset'")
  assert_that(is.numeric(durations), msg = "'durations' must be numeric")
  assert_that(is.logical(drop_empty), msg = "'drop_empty' must be logical")
  assert_that(is.logical(robust) || robust %in% c("huber", "bisquare"), 
              msg = "'robust' must be logical or one of 'huber', 'bisquare'")
  assert_that(is.logical(use_fast_path), msg = "'use_fast_path' must be logical")
  assert_that(is.logical(ar_voxelwise), msg = "'ar_voxelwise' must be logical")
  if (strategy == "chunkwise") {
    assert_that(is.numeric(nchunks) && nchunks > 0, msg = "'nchunks' must be a positive number")
  }
  
  # Convert robust parameter to type for config
  if (is.logical(robust)) {
    robust_type <- if (robust) "huber" else FALSE
  } else {
    robust_type <- robust
  }
  
  # Merge engine-supplied option overrides into the explicit arguments
  if (!is.null(engine_robust_options)) {
    if (is.null(robust_options)) {
      robust_options <- engine_robust_options
    } else {
      robust_options <- utils::modifyList(robust_options, engine_robust_options)
    }
  }

  if (!is.null(engine_ar_options)) {
    if (is.null(ar_options)) {
      ar_options <- engine_ar_options
    } else {
      ar_options <- utils::modifyList(ar_options, engine_ar_options)
    }
  }

  if (is.null(robust_options)) {
    robust_options <- list()
  }
  if (!is.null(engine_robust_options)) {
    robust_options <- utils::modifyList(robust_options, engine_robust_options)
  }
  if (!is.null(robust_type) && !("type" %in% names(robust_options))) {
    robust_options$type <- robust_type
  }
  
  # Merge individual robust parameters for backward compatibility
  if (!is.null(robust_psi) && !("type" %in% names(robust_options))) {
    robust_options$type <- robust_psi
  }
  if (!is.null(robust_max_iter) && !("max_iter" %in% names(robust_options))) {
    robust_options$max_iter <- robust_max_iter
  }
  if (!is.null(robust_scale_scope) && !("scale_scope" %in% names(robust_options))) {
    robust_options$scale_scope <- robust_scale_scope
  }
  
  # Allow low-rank shorthand `order` to map into struct/p
  if (!is.null(ar_options$order) && is.null(ar_options$struct)) {
    ar_order_tmp <- as.integer(ar_options$order[1])
    if (!is.finite(ar_order_tmp) || ar_order_tmp <= 0L) {
      ar_options$struct <- "iid"
    } else if (ar_order_tmp == 1L) {
      ar_options$struct <- "ar1"
    } else if (ar_order_tmp == 2L) {
      ar_options$struct <- "ar2"
    } else {
      ar_options$struct <- "arp"
      ar_options$p <- ar_options$p %||% ar_order_tmp
    }
  }
  ar_options$order <- NULL
  
  if (is.null(ar_options)) {
    ar_options <- list()
  }
  if (!is.null(engine_ar_options)) {
    ar_options <- utils::modifyList(ar_options, engine_ar_options)
  }

  # Merge individual AR parameters for backward compatibility
  if (!is.null(cor_struct) && !("struct" %in% names(ar_options))) {
    ar_options$struct <- cor_struct
  }
  if (!is.null(cor_iter) && !("iter_gls" %in% names(ar_options))) {
    ar_options$iter_gls <- cor_iter
  }
  if (!is.null(cor_global) && !("global" %in% names(ar_options))) {
    ar_options$global <- cor_global
  }
  if (!is.null(ar1_exact_first) && !("exact_first" %in% names(ar_options))) {
    ar_options$exact_first <- ar1_exact_first
  }
  if (!is.null(ar_p) && !("p" %in% names(ar_options))) {
    ar_options$p <- ar_p
  }
  if (!("voxelwise" %in% names(ar_options))) {
    ar_options$voxelwise <- ar_voxelwise
  }
  
  # Create configuration object
  cfg <- if (!is.null(engine_cfg) && inherits(engine_cfg, "fmri_lm_config")) {
    engine_cfg
  } else {
    fmri_lm_control(robust_options = robust_options, ar_options = ar_options)
  }

  # If both were supplied, prefer the merged configuration but keep engine_cfg metadata
  if (!is.null(engine_cfg) && inherits(engine_cfg, "fmri_lm_config")) {
    cfg$robust <- engine_cfg$robust
    cfg$ar <- utils::modifyList(cfg$ar, engine_cfg$ar)
  }
  
  model <- create_fmri_model(formula, block, baseline_model, dataset, durations = durations, drop_empty = drop_empty)

  if (!is.null(engine)) {
    res <- fmri_lm_lowrank_dispatch(
      formula_or_model = model,
      dataset = dataset,
      engine = engine,
      lowrank = lowrank,
      cfg = cfg
    )
    if (!is.null(res)) {
      attr(res, "config") <- cfg
      attr(res, "strategy") <- "sketch"
      return(res)
    }
  }
  
  # Pass configuration object down
  # Note: We don't pass ... here because all parameters have been processed
  # and included in the cfg object
  ret <- fmri_lm_fit(model, dataset, strategy, cfg, nchunks,
                     use_fast_path = use_fast_path, progress = progress,
                     parallel_voxels = parallel_voxels)
  return(ret)
}

#' @rdname fmri_lm
#' @export
fmri_lm.fmri_model <- function(formula, dataset = NULL,
                               robust = FALSE, robust_options = NULL,
                               ar_options = NULL,
                               strategy = c("runwise", "chunkwise"), nchunks = 10,
                               use_fast_path = FALSE, progress = FALSE,
                               ar_voxelwise = FALSE, parallel_voxels = FALSE,
                               cor_struct = NULL, cor_iter = NULL,
                               cor_global = NULL, ar1_exact_first = NULL,
                               ar_p = NULL,
                               robust_psi = NULL, robust_max_iter = NULL,
                               robust_scale_scope = NULL,
                               ...) {
  # --- low-rank / sketch engine fast-path (opt-in via ...$engine) ---
  dots <- list(...)
  engine  <- dots$engine
  lowrank <- dots$lowrank
  engine_ar_options <- dots$ar_options
  engine_robust_options <- dots$robust_options
  engine_cfg <- dots$cfg
  dots$engine <- NULL
  dots$lowrank <- NULL
  dots$ar_options <- NULL
  dots$robust_options <- NULL
  dots$cfg <- NULL
  strategy <- match.arg(strategy)
  assert_that(inherits(formula, "fmri_model"))

  dataset <- dataset %||% formula$dataset %||% attr(formula, "dataset")
  if (is.null(dataset)) {
    stop("No dataset found in 'formula' and none supplied.")
  }
  assert_that(inherits(dataset, "fmri_dataset"))

  if (is.logical(robust)) {
    robust_type <- if (robust) "huber" else FALSE
  } else {
    robust_type <- robust
  }

  if (is.null(robust_options)) {
    robust_options <- list()
  }
  if (!is.null(robust_type) && !("type" %in% names(robust_options))) {
    robust_options$type <- robust_type
  }
  if (!is.null(robust_psi) && !("type" %in% names(robust_options))) {
    robust_options$type <- robust_psi
  }
  if (!is.null(robust_max_iter) && !("max_iter" %in% names(robust_options))) {
    robust_options$max_iter <- robust_max_iter
  }
  if (!is.null(robust_scale_scope) && !("scale_scope" %in% names(robust_options))) {
    robust_options$scale_scope <- robust_scale_scope
  }

  if (is.null(ar_options)) {
    ar_options <- list()
  }
  if (!is.null(cor_struct) && !("struct" %in% names(ar_options))) {
    ar_options$struct <- cor_struct
  }
  if (!is.null(cor_iter) && !("iter_gls" %in% names(ar_options))) {
    ar_options$iter_gls <- cor_iter
  }
  if (!is.null(cor_global) && !("global" %in% names(ar_options))) {
    ar_options$global <- cor_global
  }
  if (!is.null(ar1_exact_first) && !("exact_first" %in% names(ar_options))) {
    ar_options$exact_first <- ar1_exact_first
  }
  if (!is.null(ar_p) && !("p" %in% names(ar_options))) {
    ar_options$p <- ar_p
  }
  if (!("voxelwise" %in% names(ar_options))) {
    ar_options$voxelwise <- ar_voxelwise
  }

  if (!is.null(ar_options$order) && is.null(ar_options$struct)) {
    ar_order_tmp <- as.integer(ar_options$order[1])
    if (!is.finite(ar_order_tmp) || ar_order_tmp <= 0L) {
      ar_options$struct <- "iid"
    } else if (ar_order_tmp == 1L) {
      ar_options$struct <- "ar1"
    } else if (ar_order_tmp == 2L) {
      ar_options$struct <- "ar2"
    } else {
      ar_options$struct <- "arp"
      ar_options$p <- ar_options$p %||% ar_order_tmp
    }
  }
  ar_options$order <- NULL

  cfg <- if (!is.null(engine_cfg) && inherits(engine_cfg, "fmri_lm_config")) {
    engine_cfg
  } else {
    fmri_lm_control(robust_options = robust_options, ar_options = ar_options)
  }

  if (!is.null(engine_cfg) && inherits(engine_cfg, "fmri_lm_config")) {
    cfg$robust <- engine_cfg$robust
    cfg$ar <- utils::modifyList(cfg$ar, engine_cfg$ar)
  }

  if (!is.null(engine)) {
    res <- fmri_lm_lowrank_dispatch(
      formula_or_model = formula,
      dataset = dataset %||% formula$dataset %||% attr(formula, "dataset"),
      engine = engine,
      lowrank = lowrank,
      cfg = cfg
    )
    if (!is.null(res)) {
      attr(res, "config") <- cfg
      attr(res, "strategy") <- "sketch"
      return(res)
    }
  }

  ret <- fmri_lm_fit(formula, dataset, strategy, cfg, nchunks,
                     use_fast_path = use_fast_path, progress = progress,
                     parallel_voxels = parallel_voxels)
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
#' @param cfg An \code{fmri_lm_config} object containing all fitting options. See \code{\link{fmri_lm_control}}.
#' @param nchunks Number of data chunks when strategy is \code{"chunkwise"}. Default is \code{10}.
#' @param use_fast_path Logical. If \code{TRUE}, use matrix-based computation for speed. Default is \code{FALSE}.
#' @param progress Logical. Whether to display a progress bar during model fitting. Default is \code{FALSE}.
#' @param parallel_voxels Logical. If TRUE, voxelwise AR processing within runs
#'   is parallelised using `future.apply`. Default is \code{FALSE}.
#' @param ... Additional arguments.
#' @return A fitted fMRI linear regression model with the specified fitting strategy.
#' @keywords internal
#' @seealso \code{\link{fmri_lm}}, \code{\link{fmri_model}}, \code{\link{fmri_dataset}}
fmri_lm_fit <- function(fmrimod, dataset, strategy = c("runwise", "chunkwise"),
                        cfg, nchunks = 10, use_fast_path = FALSE, progress = FALSE,
                        parallel_voxels = FALSE, ...) {
  strategy <- match.arg(strategy)
  
  # Validate config
  assert_that(inherits(cfg, "fmri_lm_config"), msg = "'cfg' must be an 'fmri_lm_config' object")
  
  # Error checking
  assert_that(inherits(fmrimod, "fmri_model"), msg = "'fmrimod' must be an 'fmri_model' object")
  assert_that(inherits(dataset, "fmri_dataset"), msg = "'dataset' must be an 'fmri_dataset' object")
  assert_that(is.logical(use_fast_path), msg = "'use_fast_path' must be logical")
  if (strategy == "chunkwise") {
    assert_that(is.numeric(nchunks) && nchunks > 0, msg = "'nchunks' must be a positive number")
  }
  
  # Get contrast info (now flattened from contrast_weights.event_model)
  contrast_info_flattened <- contrast_weights(fmrimod$event_model)
  full_design_colnames <- colnames(design_matrix(fmrimod))
  processed_conlist <- list()
  
  # Get column indices from design matrix for each term
  col_indices <- attr(fmrimod$event_model$design_matrix, "col_indices")
  
  if (length(contrast_info_flattened) > 0 && !is.null(col_indices)) {
    # Process each contrast in the flattened structure
    for (contrast_name in names(contrast_info_flattened)) {
      con_spec <- contrast_info_flattened[[contrast_name]]
      
      # Extract term name from flattened name (e.g., "term#contrast" or "term.contrast")
      term_name <- trimws(contrast_name)
      if (grepl("#", term_name, fixed = TRUE)) {
        term_name <- sub("#.*$", "", term_name)
      } else {
        term_name <- sub("\\..*$", "", term_name)
      }
      
      if (!is.null(col_indices[[term_name]])) {
        # Get the column indices for this term
        colind <- col_indices[[term_name]]
        
        if (length(colind) == 0) {
          warning(paste("No column indices found for term:", term_name))
          next # Skip this contrast if columns can't be found
        }
        
        # Apply colind attribute to the contrast spec
        if (inherits(con_spec, "contrast") || inherits(con_spec, "Fcontrast")) {
          # Set the colind attribute on the contrast weights for the slow path
          attr(con_spec$weights, "colind") <- colind
          # Also set it directly on the contrast object for estimate_contrast
          attr(con_spec, "colind") <- colind
          # Add to processed list
          processed_conlist[[contrast_name]] <- con_spec
        } else {
          warning(paste("Item '", contrast_name, "' is not a contrast or Fcontrast object."))
        }
      } else {
        warning(paste("Contrast '", contrast_name, "' refers to term '", term_name, "' but col_indices are missing."))
      }
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
  phi_global <- NULL
  sigma_global <- NULL
  if (cfg$ar$global && cfg$ar$struct != "iid") {
    ar_order <- switch(cfg$ar$struct,
                       ar1 = 1L,
                       ar2 = 2L,
                       arp = cfg$ar$p)

    chunk_iter <- exec_strategy("runwise")(dataset)
    run_chunks <- collect_chunks(chunk_iter)
    
    form <- get_formula(fmrimod)
    resid_vec <- numeric(0)
    for (rch in run_chunks) {
      tmats_run <- term_matrices(fmrimod, rch$chunk_num)
      data_env_run <- list2env(tmats_run)
      n_time_run <- nrow(tmats_run[[1]])
      data_env_run[[".y"]] <- rep(0, n_time_run)
      X_run <- model.matrix(form, data_env_run)
      proj_run <- .fast_preproject(X_run)
      Y_run <- as.matrix(rch$data)
      ols <- .fast_lm_matrix(X_run, Y_run, proj_run, return_fitted = TRUE)
      resid_vec <- c(resid_vec, rowMeans(Y_run - ols$fitted))
    }
    phi_global <- estimate_ar_parameters(resid_vec, ar_order)
    cfg$ar$iter_gls <- 1L
  }

  if (cfg$robust$type != FALSE && cfg$robust$scale_scope == "global") {
    chunk_iter <- exec_strategy("runwise")(dataset)
    run_chunks <- collect_chunks(chunk_iter)
    form <- get_formula(fmrimod)
    row_med_all <- numeric(0)
    for (rch in run_chunks) {
      tmats_run <- term_matrices(fmrimod, rch$chunk_num)
      data_env_run <- list2env(tmats_run)
      n_time_run <- nrow(tmats_run[[1]])
      data_env_run[[".y"]] <- rep(0, n_time_run)
      X_run <- model.matrix(form, data_env_run)
      proj_run <- .fast_preproject(X_run)
      Y_run <- as.matrix(rch$data)
      ols <- .fast_lm_matrix(X_run, Y_run, proj_run, return_fitted = TRUE)
      row_med_all <- c(row_med_all,
                       matrixStats::rowMedians(abs(Y_run - ols$fitted)))
    }
    sigma_global <- 1.4826 * median(row_med_all)
    if (sigma_global <= .Machine$double.eps) sigma_global <- .Machine$double.eps
  }

  result <- switch(strategy,
                   "runwise" = runwise_lm(dataset, fmrimod, standard_path_conlist, # Pass full objects
                                         cfg = cfg, use_fast_path = use_fast_path,
                                         progress = progress,
                                         phi_fixed = phi_global,
                                         sigma_fixed = sigma_global,
                                         parallel_voxels = parallel_voxels
                                         ),
                   "chunkwise" = {
                    if (inherits(dataset, "latent_dataset")) {
                      chunkwise_lm(dataset, fmrimod, standard_path_conlist, # Pass full objects
                                   nchunks, cfg, verbose = FALSE, use_fast_path = FALSE,
                                   progress = progress,
                                   phi_fixed = phi_global,
                                   sigma_fixed = sigma_global
                                   ) # Do not pass use_fast_path
                    } else {
                      chunkwise_lm(dataset, fmrimod, standard_path_conlist, # Pass full objects
                                   nchunks, cfg = cfg, use_fast_path = use_fast_path,
                                   progress = progress,
                                   phi_fixed = phi_global,
                                   sigma_fixed = sigma_global
                                   )
                    }
                  })
  
  ret <- list(
    result = result,
    model = fmrimod,
    strategy = strategy,
    bcons = processed_conlist,
    dataset = dataset,
    ar_coef = result$ar_coef
  )
  
  class(ret) <- "fmri_lm"
  
  # Attach config as attribute for testing
  attr(ret, "config") <- cfg
  attr(ret, "strategy") <- strategy
  
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
  
  # Get term indices from the fmri_model's term_matrices
  tmats <- term_matrices(x$model)
  event_indices <- attr(tmats, "event_term_indices")
  
  # If no event terms, return empty list
  if (is.null(eterms) || length(eterms) == 0) {
    return(list())
  }
  
  # For each event term, compute its indices based on the term's structure
  # The event_indices should correspond to the columns for event terms
  pred <- lapply(seq_along(eterms), function(i) {
    eterm <- eterms[[i]]
    
    # Get the HRF specification (stored as an attribute in fmridesign)
    hrf_spec <- attr(eterm, "hrfspec")
    if (is.null(hrf_spec) && !is.null(eterm$hrfspec)) {
      # Backward-compatibility: some versions may store as a list element
      hrf_spec <- eterm$hrfspec
    }
    
    # Fallback HRF if spec is missing
    hrf <- if (!is.null(hrf_spec) && !is.null(hrf_spec$hrf)) hrf_spec$hrf else fmrihrf::HRF_SPMG1
    # Derive nbasis using fmrihrf helper when possible
    nb <- tryCatch({ fmrihrf::nbasis(hrf) }, error = function(e) NULL)
    
    # If nbasis is NULL, assume it's 1 (single basis function)
    if (is.null(nb)) {
      nb <- 1
    }
    
    # Get the conditions (cells) for this term
    excond <- cells(eterm, exclude_basis = TRUE)
    ncond <- nrow(excond)
    
    # Get the column indices for this term directly from the event_model's col_indices
    col_indices <- attr(x$model$event_model$design_matrix, "col_indices")
    term_name <- names(eterms)[i]
    ind <- col_indices[[term_name]]
    
    # Create the HRF basis matrix at sample points
    G <- as.matrix(hrf(sample_at))
    
    # Create block diagonal matrix for all conditions
    Gex <- do.call(Matrix::bdiag, replicate(ncond, G, simplify = FALSE))
    
    # Get the relevant betas for this term
    B <- t(betas[, ind, drop = FALSE])
    
    # Compute predicted HRF values
    yh <- Gex %*% B
    
    # Create expanded design info
    excond_expanded <- excond %>% dplyr::slice(rep(1:dplyr::n(), each = length(sample_at)))
    design <- cbind(
      dplyr::tibble(time = rep(sample_at, ncond)),
      excond_expanded
    )
    
    list(pred = as.matrix(yh), design = tibble::as_tibble(design))
  })
  
  # Set names from event terms
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

#' @method coef fmri_lm
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
      
      # Return as matrix - transpose to get conditions x voxels
      res <- t(res)
    }
  } else {
    # Should not happen due to match.arg, but defensive coding
    stop("Invalid type specified.")
  }
  
  # Reconstruction functionality can be added here if necessary (applies to the 'res' matrix/tibble)
  # if (recon && inherits(object$dataset, "fmri_dataset")) { ... }
  
  return(res)
}

#' @method stats fmri_lm
#' @rdname stats
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

#' @method standard_error fmri_lm
#' @rdname standard_error
#' @export
standard_error.fmri_lm <- function(x, type = c("estimates", "contrasts"),...) {
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
#' @param cfg An \code{fmri_lm_config} object containing all fitting options.
#' @param verbose Logical. Whether to display progress messages (default is \code{FALSE}).
#' @param use_fast_path Logical. If \code{TRUE}, use matrix-based computation for speed. Default is \code{FALSE}.
#' @param progress Logical. Display a progress bar for chunk processing. Default is \code{FALSE}.
#' @param phi_fixed Optional fixed AR parameters.
#' @param sigma_fixed Optional fixed robust scale estimate.
#' @return A list containing the unpacked chunkwise results.
#' @keywords internal
chunkwise_lm.fmri_dataset_old <- function(dset, model, contrast_objects, nchunks, cfg,
                                      verbose = FALSE, use_fast_path = FALSE, progress = FALSE,
                                      phi_fixed = NULL,
                                      sigma_fixed = NULL) {
  # Validate config
  assert_that(inherits(cfg, "fmri_lm_config"), msg = "'cfg' must be an 'fmri_lm_config' object")
  chunk_iter <- exec_strategy("chunkwise", nchunks = nchunks)(dset)
  chunks <- collect_chunks(chunk_iter)
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
      proj_global <- .fast_preproject(modmat)
      Vu <- proj_global$XtXinv
      
      lmfun <- if (cfg$robust$type != FALSE) multiresponse_rlm else multiresponse_lm
      
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
                          rss = rss, rdf = rdf, sigma = sigma,
                          ar_coef = NULL)
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
      
      if (cfg$robust$type != FALSE) {
          # Combined AR + Robust path for chunkwise
          ar_modeling <- cfg$ar$struct != "iid"
          ar_order <- switch(cfg$ar$struct,
                             ar1 = 1L,
                             ar2 = 2L,
                             arp = cfg$ar$p,
                             iid = 0L)
          
          if (ar_modeling) {
              message("Using fast path with AR + robust weighting...")
          } else {
              message("Using fast path with robust weighting...")
          }

          chunk_iter <- exec_strategy("runwise")(dset)
          run_chunks <- collect_chunks(chunk_iter)
          run_row_inds <- lapply(run_chunks, `[[`, "row_ind")

          run_info <- vector("list", length(run_chunks))

          # Pass-0: Per-run precomputation
          for (ri in seq_along(run_chunks)) {
              rch <- run_chunks[[ri]]
              tmats_run <- term_matrices(model, rch$chunk_num)
              data_env_run <- list2env(tmats_run)
              n_time_run <- nrow(tmats_run[[1]])
              data_env_run[[".y"]] <- rep(0, n_time_run)
              X_run_orig <- model.matrix(as.formula(form), data_env_run)
              proj_run_orig <- .fast_preproject(X_run_orig)
              Y_run_full <- as.matrix(rch$data)

              if (ar_modeling) {
                  # Step 1: Estimate AR parameters from initial OLS
                  glm_ctx_run_orig <- glm_context(X = X_run_orig, Y = Y_run_full, proj = proj_run_orig)
                  
                  phi_hat_run <- NULL
                  if (is.null(phi_fixed)) {
                      initial_fit <- solve_glm_core(glm_ctx_run_orig, return_fitted = TRUE)
                      resid_ols <- Y_run_full - initial_fit$fitted
                      phi_hat_run <- estimate_ar_parameters(rowMeans(resid_ols), ar_order)
                  } else {
                      phi_hat_run <- phi_fixed
                  }
                  
                  # Step 2: AR whitening
                  tmp <- ar_whiten_transform(X_run_orig, Y_run_full, phi_hat_run, cfg$ar$exact_first)
                  X_run_w <- tmp$X
                  Y_run_full_w <- tmp$Y
                  proj_run_w <- .fast_preproject(X_run_w)
                  
                  # Step 3: Create whitened GLM context
                  glm_ctx_run_whitened <- glm_context(X = X_run_w, Y = Y_run_full_w, proj = proj_run_w)
                  
                  # Determine sigma_fixed based on scope
                  sigma_fixed_for_run <- if (cfg$robust$scale_scope == "global" && !is.null(sigma_fixed)) {
                      sigma_fixed
                  } else {
                      NULL
                  }
                  
                  # Step 4: Robust fitting on whitened data
                  robust_fit_details_run <- robust_iterative_fitter(
                      initial_glm_ctx = glm_ctx_run_whitened,
                      cfg_robust_options = cfg$robust,
                      X_orig_for_resid = X_run_w,
                      sigma_fixed = sigma_fixed_for_run
                  )
                  
                  # Step 5: Optional re-estimation of AR parameters
                  if (!is.null(cfg$robust$reestimate_phi) && cfg$robust$reestimate_phi && is.null(phi_fixed)) {
                      # Calculate robust residuals on whitened data
                      resid_robust_w <- Y_run_full_w - X_run_w %*% robust_fit_details_run$betas_robust
                      
                      # Re-estimate AR parameters from robust residuals
                      phi_hat_run_updated <- estimate_ar_parameters(rowMeans(resid_robust_w), ar_order)
                      
                      # Update phi_hat_run for subsequent processing
                      phi_hat_run <- phi_hat_run_updated
                  }
                  
                  # Store results
                  run_info[[ri]] <- list(
                      phi_hat = phi_hat_run,
                      weights = robust_fit_details_run$robust_weights_final,
                      sqrtw = sqrt(robust_fit_details_run$robust_weights_final),
                      sigma = robust_fit_details_run$sigma_robust_scale_final,
                      X_orig = X_run_orig,
                      Y_orig = Y_run_full,
                      row_indices = run_row_inds[[ri]]
                  )
              } else {
                  # Robust-only path (no AR)
                  glm_ctx_run_orig <- glm_context(X = X_run_orig, Y = Y_run_full, proj = proj_run_orig)
                  
                  # Determine sigma_fixed based on scope
                  sigma_fixed_for_run <- if (cfg$robust$scale_scope == "global" && !is.null(sigma_fixed)) {
                      sigma_fixed
                  } else {
                      NULL
                  }
                  
                  # Call robust_iterative_fitter
                  robust_fit <- robust_iterative_fitter(
                      initial_glm_ctx = glm_ctx_run_orig,
                      cfg_robust_options = cfg$robust,
                      X_orig_for_resid = X_run_orig,
                      sigma_fixed = sigma_fixed_for_run
                  )

                  sqrtw <- sqrt(robust_fit$robust_weights_final)
                  
                  run_info[[ri]] <- list(
                      weights = robust_fit$robust_weights_final,
                      sqrtw = sqrtw,
                      sigma = robust_fit$sigma_robust_scale_final,
                      X_orig = X_run_orig,
                      Y_orig = Y_run_full,
                      row_indices = run_row_inds[[ri]]
                  )
              }
          }

          # Build global transformed matrices
          if (ar_modeling) {
              # For AR+Robust: First whiten, then weight
              X_global_list <- vector("list", length(run_chunks))
              
              for (ri in seq_along(run_chunks)) {
                  # First apply AR whitening
                  X_orig <- run_info[[ri]]$X_orig
                  phi_hat <- run_info[[ri]]$phi_hat
                  
                  # Whiten X (with dummy Y since we only need X transformation)
                  dummyY <- matrix(0, nrow(X_orig), 0)
                  X_whitened <- ar_whiten_transform(X_orig, dummyY, phi_hat, cfg$ar$exact_first)$X
                  
                  # Then apply robust weights
                  X_whitened_weighted <- X_whitened * run_info[[ri]]$sqrtw
                  X_global_list[[ri]] <- X_whitened_weighted
              }
              
              X_global_final_w <- do.call(rbind, X_global_list)
          } else {
              # For Robust-only: Just apply weights
              X_weighted_list <- vector("list", length(run_chunks))
              for (ri in seq_along(run_chunks)) {
                  X_weighted_list[[ri]] <- run_info[[ri]]$X_orig * run_info[[ri]]$sqrtw
              }
              X_global_final_w <- do.call(rbind, X_weighted_list)
          }
          
          proj_global_final_w <- .fast_preproject(X_global_final_w)
          
          Vu <- .fast_preproject(design_matrix(model))$XtXinv
          nvox <- ncol(dset$datamat)

          cres <- vector("list", length(run_chunks))

          chunk_iter <- exec_strategy("chunkwise", nchunks = nchunks)(dset)
          chunks <- collect_chunks(chunk_iter)
          if (progress) {
            pb <- cli::cli_progress_bar("Fitting chunks", total = length(chunks), clear = FALSE)
            on.exit(cli::cli_progress_done(id = pb), add = TRUE)
          }

          for (i in seq_along(chunks)) {
              ym <- chunks[[i]]
              if (verbose) message("Processing chunk (fast robust) ", ym$chunk_num)
              Ymat <- as.matrix(ym$data)
              
              # Transform Y data (whiten then weight)
              Y_transformed_list <- vector("list", length(run_chunks))
              
              if (ar_modeling) {
                  # For AR+Robust: First whiten, then weight
                  for (ri in seq_along(run_chunks)) {
                      rows <- run_info[[ri]]$row_indices
                      Y_chunk_segment <- Ymat[rows, , drop=FALSE]
                      phi_hat <- run_info[[ri]]$phi_hat
                      
                      # Whiten Y
                      dummyX <- matrix(0, length(rows), 0)
                      Y_whitened <- ar_whiten_transform(dummyX, Y_chunk_segment, phi_hat, cfg$ar$exact_first)$Y
                      
                      # Then apply robust weights
                      Y_whitened_weighted <- sweep(Y_whitened, 1, run_info[[ri]]$sqrtw, `*`)
                      Y_transformed_list[[ri]] <- Y_whitened_weighted
                  }
              } else {
                  # For Robust-only: Just apply weights
                  for (ri in seq_along(run_chunks)) {
                      rows <- run_info[[ri]]$row_indices
                      Y_weighted <- sweep(Ymat[rows, , drop=FALSE], 1, run_info[[ri]]$sqrtw, `*`)
                      Y_transformed_list[[ri]] <- Y_weighted
                  }
              }
              
              Y_chunk_final_w <- do.call(rbind, Y_transformed_list)
              
              # Create GLM context with pre-transformed matrices
              glm_ctx_final <- glm_context(
                  X = X_global_final_w, 
                  Y = Y_chunk_final_w, 
                  proj = proj_global_final_w
              )
              
              # Solve using transformed context
              res <- solve_glm_core(glm_ctx_final)
              
              # Calculate statistics
              actual_vnames <- colnames(X_global_final_w)
              
              # For sigma, we need to map run-specific sigmas to voxels
              # This is a simplification - in practice, we might need per-run statistics
              sigma_vec <- sqrt(res$sigma2)
              
              bstats <- beta_stats_matrix(res$betas, 
                                         proj_global_final_w$XtXinv, 
                                         sigma_vec,
                                         proj_global_final_w$dfres, 
                                         actual_vnames)
              
              contrasts <- fit_lm_contrasts_fast(res$betas, 
                                               res$sigma2, 
                                               proj_global_final_w$XtXinv,
                                               simple_conlist_weights, 
                                               fconlist_weights, 
                                               proj_global_final_w$dfres)
              
              cres[[i]] <- list(bstats = bstats,
                               contrasts = contrasts,
                               rss = res$rss,
                               rdf = proj_global_robustly_weighted$dfres,
                               sigma = sigma_vec,
                               ar_coef = NULL)
                               
              if (progress) cli::cli_progress_update(id = pb)
          }

          # Unpack results (expects specific structure from cres)
          out <- unpack_chunkwise(cres, event_indices, baseline_indices)
          # Add cov.unscaled to the output
          out$cov.unscaled <- Vu
          
          return(out)
      }

      message("Using fast path for chunkwise LM...")

      data_env <- list2env(tmats)
      data_env[[ ".y"]] <- rep(0, nrow(tmats[[1]])) # Placeholder for model.matrix
      modmat_orig <- model.matrix(as.formula(form), data_env)

      ar_modeling <- cfg$ar$struct != "iid"
      ar_order <- switch(cfg$ar$struct,
                         ar1 = 1L,
                         ar2 = 2L,
                         arp = cfg$ar$p,
                         iid = 0L)

      chunk_iter <- exec_strategy("runwise")(dset)
      run_chunks <- collect_chunks(chunk_iter)
      run_row_inds <- lapply(run_chunks, `[[`, "row_ind")

      if (ar_modeling && cfg$ar$iter_gls > 1L) {
          warning("iter_gls > 1 not supported in chunkwise fast path; using 1")
          cfg$ar$iter_gls <- 1L
      }

      phi_hat_list <- vector("list", length(run_chunks))

      if (ar_modeling) {
          X_w_list <- vector("list", length(run_chunks))
          for (ri in seq_along(run_chunks)) {
              rch <- run_chunks[[ri]]

              tmats_run <- term_matrices(model, rch$chunk_num)
              data_env_run <- list2env(tmats_run)
              n_time_run <- nrow(tmats_run[[1]])
              data_env_run[[".y"]] <- rep(0, n_time_run)
              X_run <- model.matrix(as.formula(form), data_env_run)
              proj_run <- .fast_preproject(X_run)
              Y_run <- as.matrix(rch$data)

              if (is.null(phi_fixed)) {
                  # Create GLM context for initial OLS
                  glm_ctx_run <- glm_context(X = X_run, Y = Y_run, proj = proj_run)
                  ols <- solve_glm_core(glm_ctx_run, return_fitted = TRUE)
                  resid_ols <- Y_run - ols$fitted
                  phi_hat_run <- estimate_ar_parameters(rowMeans(resid_ols), ar_order)
              } else {
                  phi_hat_run <- phi_fixed
              }
              dummyY <- matrix(0, nrow(X_run), 0)
              X_run <- ar_whiten_transform(X_run, dummyY, phi_hat_run, cfg$ar$exact_first)$X

              phi_hat_list[[ri]] <- phi_hat_run
              X_w_list[[ri]] <- X_run
          }
          modmat <- do.call(rbind, X_w_list)
      } else {
          modmat <- modmat_orig
      }
      proj <- .fast_preproject(modmat)
      Vu <- proj$XtXinv

      cres <- vector("list", length(chunks))
      for (i in seq_along(chunks)) {
          ym <- chunks[[i]]
          if (verbose) message("Processing chunk (fast path) ", ym$chunk_num)
          Ymat <- as.matrix(ym$data)

          if (ar_modeling) {
              for (ri in seq_along(run_chunks)) {
                  rows <- run_row_inds[[ri]]
                  phi <- phi_hat_list[[ri]]
                  if (!is.null(phi)) {
                      dummyX <- matrix(0, length(rows), 0)
                      subY <- Ymat[rows, , drop = FALSE]
                      subY <- ar_whiten_transform(dummyX, subY, phi, cfg$ar$exact_first)$Y
                      Ymat[rows, ] <- subY
                  }
              }
          }

          if (verbose) message("  Chunk ", ym$chunk_num, ": ncol(Ymat) = ", ncol(Ymat))

          # Create GLM context for chunk processing
          glm_ctx_chunk <- glm_context(X = modmat, Y = Ymat, proj = proj)
          res <- solve_glm_core(glm_ctx_chunk)

          actual_vnames <- colnames(modmat)
          sigma_vec <- sqrt(res$sigma2)
          bstats <- beta_stats_matrix(res$betas, proj$XtXinv, sigma_vec, proj$dfres, actual_vnames,
                                       robust_weights = NULL, ar_order = ar_order)

          contrasts <- fit_lm_contrasts_fast(res$betas, res$sigma2, proj$XtXinv,
                                             simple_conlist_weights, fconlist_weights, proj$dfres,
                                             robust_weights = NULL, ar_order = ar_order)

          cres[[i]] <- list(bstats = bstats,
                            contrasts = contrasts,
                            rss = res$rss,
                            rdf = proj$dfres,
                            sigma = sigma_vec,
                            ar_coef = NULL)
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
#' @param cfg An \code{fmri_lm_config} object containing all fitting options.
#' @param verbose Logical. Whether to display progress messages (default is \code{FALSE}).
#' @param use_fast_path Logical. Whether to use fast path computation (default is \code{FALSE}).
#' @param progress Logical. Display a progress bar for run processing. Default is \code{FALSE}.
#' @param phi_fixed Optional fixed AR parameters.
#' @param sigma_fixed Optional fixed robust scale estimate.
#' @param parallel_voxels Logical. If TRUE, process voxels in parallel using
#'   `future.apply`. Default is \code{FALSE}.
#' @return A list containing the combined results from runwise linear model analysis.
#' @keywords internal
#' @autoglobal
runwise_lm <- function(dset, model, contrast_objects, cfg, verbose = FALSE,
                       use_fast_path = FALSE, progress = FALSE,
                       phi_fixed = NULL,
                       sigma_fixed = NULL,
                       parallel_voxels = FALSE) {
  # Validate config
  assert_that(inherits(cfg, "fmri_lm_config"), msg = "'cfg' must be an 'fmri_lm_config' object")
  # Get an iterator of data chunks (runs)
  chunk_iter <- exec_strategy("runwise")(dset)
  chunks <- collect_chunks(chunk_iter)
  if (progress) {
    pb <- cli::cli_progress_bar("Fitting runs", total = length(chunks), clear = FALSE)
    on.exit(cli::cli_progress_done(id = pb), add = TRUE)
  }
  form <- get_formula(model)
  # Global design matrix for pooling across runs
  modmat_global <- as.matrix(design_matrix(model))
  Qr_global <- qr(modmat_global)
  proj_global <- .fast_preproject(modmat_global)
  if (Qr_global$rank < ncol(modmat_global)) {
    bad_cols <- colnames(modmat_global)[
      Qr_global$pivot[(Qr_global$rank + 1):ncol(modmat_global)]
    ]
    warning(
      sprintf(
        "Global design matrix is rank deficient; problematic columns: %s",
        paste(bad_cols, collapse = ", ")
      )
    )
    Vu <- MASS::ginv(crossprod(modmat_global))
    colnames(Vu) <- rownames(Vu) <- colnames(modmat_global)
  } else {
    Vu <- proj_global$XtXinv
    colnames(Vu) <- rownames(Vu) <- colnames(modmat_global)
  }
  
  # Define ym for R CMD check
  ym <- NULL
  
  if (!use_fast_path) {
      # -------- Slow Path using Integrated Solver --------
      # Check if integrated solver is available
      has_integrated_solver <- exists("solve_integrated_glm", 
                                     where = asNamespace("fmrireg"), 
                                     mode = "function") &&
                              exists("process_run_integrated", 
                                     where = asNamespace("fmrireg"),
                                     mode = "function")
      
      if (has_integrated_solver && !cfg$ar$voxelwise) {
          # Use integrated solver for non-voxelwise cases
          if (verbose) message("Using integrated solver for slow path...")
          
          cres <- vector("list", length(chunks))
          for (i in seq_along(chunks)) {
              ym <- chunks[[i]]
              if (verbose) message("Processing run ", ym$chunk_num)
              tmats <- term_matrices(model, ym$chunk_num)
              vnames <- attr(tmats, "varnames")
              event_indices <- attr(tmats, "event_term_indices")
              baseline_indices <- attr(tmats, "baseline_term_indices")
              
              data_env <- list2env(tmats)
              data_env$.y <- rep(0, nrow(tmats[[1]]))
              X_run <- model.matrix(form, data_env)
              Y_run <- as.matrix(ym$data)
              
              # Extract contrast weights
              simple_conlist <- Filter(function(x) inherits(x, "contrast"), contrast_objects)
              fconlist <- Filter(function(x) inherits(x, "Fcontrast"), contrast_objects)
              simple_conlist_weights <- lapply(simple_conlist, `[[`, "weights")
              names(simple_conlist_weights) <- names(simple_conlist)
              fconlist_weights <- lapply(fconlist, `[[`, "weights")
              names(fconlist_weights) <- names(fconlist)
              
              # Use integrated solver
              run_result <- process_run_integrated(
                  X_run = X_run, 
                  Y_run = Y_run, 
                  cfg = cfg,
                  phi_fixed = phi_fixed,
                  sigma_fixed = sigma_fixed,
                  conlist_weights = simple_conlist_weights,
                  fconlist_weights = fconlist_weights,
                  vnames = vnames
              )
              
              # Store results in expected format
              cres[[i]] <- list(
                  conres = run_result$contrasts,
                  bstats = beta_stats_matrix(run_result$betas, run_result$XtXinv, 
                                            run_result$sigma, run_result$dfres, 
                                            colnames(X_run),
                                            robust_weights = run_result$robust_weights, 
                                            ar_order = if (!is.null(run_result$phi_hat)) length(run_result$phi_hat) else 0),
                  event_indices = event_indices,
                  baseline_indices = baseline_indices,
                  rss = run_result$rss,
                  rdf = run_result$dfres,
                  resvar = mean(run_result$sigma^2),
                  sigma = run_result$sigma,
                  ar_coef = run_result$phi_hat
              )
              if (progress) cli::cli_progress_update(id = pb)
          }
      } else if (cfg$ar$voxelwise && cfg$ar$struct != "iid") {
          # Voxelwise AR path using modular solver
          ar_order <- switch(cfg$ar$struct,
                             ar1 = 1L,
                             ar2 = 2L,
                             arp = cfg$ar$p,
                             iid = 0L)

          simple_conlist <- Filter(function(x) inherits(x, "contrast"), contrast_objects)
          fconlist <- Filter(function(x) inherits(x, "Fcontrast"), contrast_objects)
          simple_conlist_weights <- lapply(simple_conlist, `[[`, "weights")
          names(simple_conlist_weights) <- names(simple_conlist)
          fconlist_weights <- lapply(fconlist, `[[`, "weights")
          names(fconlist_weights) <- names(fconlist)

          if (verbose) message("Using voxelwise AR(", ar_order, ") modeling...")

          cres <- vector("list", length(chunks))
          for (i in seq_along(chunks)) {
              ym <- chunks[[i]]
              if (verbose) message("Processing run ", ym$chunk_num, " with voxelwise AR")
              tmats <- term_matrices(model, ym$chunk_num)
              vnames <- attr(tmats, "varnames")
              event_indices <- attr(tmats, "event_term_indices")
              baseline_indices <- attr(tmats, "baseline_term_indices")

              data_env <- list2env(tmats)
              Y_run <- as.matrix(ym$data)
              data_env$.y <- Y_run
              X_run <- model.matrix(form, data_env)
              proj_run <- .fast_preproject(X_run)  # Compute once outside the loop

              n_voxels <- ncol(Y_run)
              p <- ncol(X_run)
              dfres <- nrow(X_run) - p

              betas_voxelwise <- matrix(NA_real_, p, n_voxels)
              sigma_voxelwise <- numeric(n_voxels)
              rss_voxelwise <- numeric(n_voxels)
              XtXinv_voxel_list <- vector("list", n_voxels)
              robust_w_list <- if (cfg$robust$type != FALSE) vector("list", n_voxels) else NULL

              process_voxel <- function(v) {
                  y_voxel <- matrix(Y_run[, v], ncol = 1)
                  ctx_vox <- glm_context(X_run, y_voxel, proj = proj_run)
                  ols <- solve_glm_core(ctx_vox, return_fitted = TRUE)
                  resid_ols <- y_voxel - ols$fitted

                  phi_voxel <- if (is.null(phi_fixed)) {
                      estimate_ar_parameters(drop(resid_ols), ar_order)
                  } else {
                      phi_fixed
                  }

                  tmp <- ar_whiten_transform(X_run, y_voxel, phi_voxel, cfg$ar$exact_first)
                  proj_w <- .fast_preproject(tmp$X)
                  ctx_vox_w <- glm_context(X = tmp$X, Y = tmp$Y, proj = proj_w, phi_hat = phi_voxel)

              if (cfg$robust$type != FALSE) {
                      sigma_fixed_for_vox <- if (cfg$robust$scale_scope == "global" && !is.null(sigma_fixed)) sigma_fixed else NULL
                      rfit <- robust_iterative_fitter(ctx_vox_w, cfg$robust, tmp$X, sigma_fixed_for_vox)
                      list(beta = rfit$betas_robust,
                           XtXinv = rfit$XtWXi_final,
                           rss = sum((tmp$Y - tmp$X %*% rfit$betas_robust)^2),
                           sigma = rfit$sigma_robust_scale_final,
                           rw = rfit$robust_weights_final,
                           phi = phi_voxel)
              } else {
                      fit_v <- solve_glm_core(ctx_vox_w)
                      list(beta = fit_v$betas,
                           XtXinv = proj_w$XtXinv,
                           rss = fit_v$rss,
                           sigma = sqrt(fit_v$sigma2),
                           rw = NULL,
                           phi = phi_voxel)
              }
          }

              voxel_indices <- seq_len(n_voxels)
              voxel_results <- if (parallel_voxels) {
                  with_package("future.apply")
                  if (progress) {
                      with_package("progressr")
                      progressr::with_progress({
                          p <- progressr::progressor(along = voxel_indices)
                          future.apply::future_lapply(voxel_indices, function(v) { p(); process_voxel(v) })
                      })
                  } else {
                      future.apply::future_lapply(voxel_indices, process_voxel)
                  }
              } else {
                  lapply(voxel_indices, process_voxel)
              }

              for (v in voxel_indices) {
                  res_v <- voxel_results[[v]]
                  betas_voxelwise[, v] <- res_v$beta
                  XtXinv_voxel_list[[v]] <- res_v$XtXinv
                  rss_voxelwise[v] <- res_v$rss
                  sigma_voxelwise[v] <- res_v$sigma
                  if (!is.null(robust_w_list)) robust_w_list[[v]] <- res_v$rw
              }

              bstats <- beta_stats_matrix_voxelwise(betas_voxelwise, XtXinv_voxel_list,
                                                    sigma_voxelwise, dfres, vnames,
                                                    robust_w_list, ar_order)

              conres <- fit_lm_contrasts_voxelwise(betas_voxelwise, sigma_voxelwise^2,
                                                    XtXinv_voxel_list, simple_conlist_weights,
                                                    fconlist_weights, dfres,
                                                    robust_w_list, ar_order)

              phi_voxel_list <- lapply(voxel_results, `[[`, "phi")

              cres[[i]] <- list(
                  conres = conres,
                  bstats = bstats,
                  event_indices = event_indices,
                  baseline_indices = baseline_indices,
                  rss = rss_voxelwise,
                  rdf = dfres,
                  resvar = sigma_voxelwise^2,
                  sigma = sigma_voxelwise,
                  ar_coef = phi_voxel_list
              )
              if (progress) cli::cli_progress_update(id = pb)
          }
      } else {
          # Standard slow path (no voxelwise AR)
          lmfun <- if (cfg$robust$type != FALSE) multiresponse_rlm else multiresponse_lm
          
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
                  sigma = sigma,
                  ar_coef = NULL
              )
              if (progress) cli::cli_progress_update(id = pb)
          }
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
      
      message("Using fast path for runwise LM...")

      ar_modeling <- cfg$ar$struct != "iid"
      ar_order <- switch(cfg$ar$struct,
                         ar1 = 1L,
                         ar2 = 2L,
                         arp = cfg$ar$p,
                         iid = 0L)
      
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

        # Check if we should use the integrated solver
        # These are internal functions, so check in the package namespace
        has_integrated_solver <- exists("solve_integrated_glm", 
                                       where = asNamespace("fmrireg"), 
                                       mode = "function") &&
                                exists("process_run_integrated", 
                                       where = asNamespace("fmrireg"),
                                       mode = "function")
        
        # Always use integrated solver if available (for consistency)
        use_integrated_solver <- has_integrated_solver
        
        if (getOption("fmrireg.debug", FALSE)) {
          message("Checking integrated solver: robust=", cfg$robust$type, 
                  ", ar_modeling=", ar_modeling,
                  ", has_integrated_solver=", has_integrated_solver,
                  ", use_integrated_solver=", use_integrated_solver)
        }
        
        if (use_integrated_solver) {
            # Use the integrated solver for AR, Robust, or AR+Robust
            run_result <- process_run_integrated(
                X_run = X_run, 
                Y_run = Y_run, 
                cfg = cfg,
                phi_fixed = phi_fixed,
                sigma_fixed = sigma_fixed,
                conlist_weights = simple_conlist_weights,
                fconlist_weights = fconlist_weights,
                vnames = vnames
            )
            
            # Extract results in format expected by runwise_lm
            betas_final <- run_result$betas
            sigma_vec <- run_result$sigma
            rss_final <- run_result$rss
            dfres_final <- run_result$dfres
            XtXinv_final <- run_result$XtXinv
            actual_vnames <- colnames(X_run)
            
            # Calculate beta stats
            bstats <- beta_stats_matrix(betas_final, XtXinv_final, sigma_vec, 
                                        dfres_final, actual_vnames,
                                        robust_weights = run_result$robust_weights, 
                                        ar_order = ar_order)
            
            # Contrasts were already computed by process_run_integrated
            conres <- run_result$contrasts
            
            # Store results
            cres[[i]] <- list(
              conres = conres,
              bstats = bstats,
              event_indices = event_indices,
              baseline_indices = baseline_indices,
              rss = rss_final,
              rdf = dfres_final,
              resvar = rss_final / dfres_final,
              sigma = sigma_vec,
              ar_coef = run_result$phi_hat
            )
            
            if (progress) cli::cli_progress_update(id = pb)
        } else {
            # Non-robust path - create GLM context for the run
            glm_ctx_run_orig <- glm_context(X = X_run, Y = Y_run, proj = proj_run)
            
            if (ar_modeling) {
                # Estimate phi_hat_run if needed
                phi_hat_run <- NULL
                if (is.null(phi_fixed)) {
                    initial_fit <- solve_glm_core(glm_ctx_run_orig, return_fitted = TRUE)
                    resid_ols <- Y_run - initial_fit$fitted
                    phi_hat_run <- estimate_ar_parameters(rowMeans(resid_ols), ar_order)
                } else {
                    phi_hat_run <- phi_fixed
                }
                
                # Iterative GLS
                for (iter in seq_len(cfg$ar$iter_gls)) {
                    # Whiten the data
                    tmp <- ar_whiten_transform(X_run, Y_run, phi_hat_run, cfg$ar$exact_first)
                    X_w <- tmp$X
                    Y_w <- tmp$Y
                    proj_w <- .fast_preproject(X_w)
                    
                    # Create whitened GLM context
                    glm_ctx_run_whitened <- glm_context(X = X_w, Y = Y_w, proj = proj_w, phi_hat = phi_hat_run)
                    
                    # Solve using whitened context
                    gls <- solve_glm_core(glm_ctx_run_whitened)
                    
                    # Update phi_hat if needed for next iteration
                    if (is.null(phi_fixed) && iter < cfg$ar$iter_gls) {
                        resid_gls <- Y_w - X_w %*% gls$betas
                        phi_hat_run <- estimate_ar_parameters(rowMeans(resid_gls), ar_order)
                    }
                }
                
                # Use whitened matrices for final results
                X_iter <- X_w
                Y_iter <- Y_w
                proj_iter <- proj_w
            } else {
                # IID case - use original context
                gls <- solve_glm_core(glm_ctx_run_orig)
                X_iter <- X_run
                Y_iter <- Y_run
                proj_iter <- proj_run
            }

            actual_vnames <- colnames(X_iter)
            sigma_vec <- sqrt(gls$sigma2)
            bstats <- beta_stats_matrix(gls$betas, proj_iter$XtXinv, sigma_vec, proj_iter$dfres, actual_vnames,
                                        robust_weights = NULL, ar_order = ar_order)

            conres <- fit_lm_contrasts_fast(gls$betas, gls$sigma2, proj_iter$XtXinv,
                                             simple_conlist_weights, fconlist_weights, proj_iter$dfres,
                                             robust_weights = NULL, ar_order = ar_order)

            cres[[i]] <- list(
              conres = conres,
              bstats = bstats,
              event_indices = event_indices,
              baseline_indices = baseline_indices,
              rss = gls$rss,
              rdf = proj_iter$dfres,
              resvar = gls$sigma2,
              sigma = sigma_vec,
              ar_coef = if (ar_modeling) list(phi_hat_run) else NULL
            )
        }
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
  rdf_vals <- unlist(lapply(cres, `[[`, "rdf"))
  sigma_mat <- do.call(rbind, lapply(seq_along(cres), function(i) {
    as.matrix(cres[[i]]$sigma^2) * rdf_vals[i]
  }))
  sigma <- sqrt(colSums(sigma_mat) / sum(rdf_vals))
  rss <- colSums(do.call(rbind, lapply(cres, function(x) as.matrix(x$rss))))
  rdf <- sum(rdf_vals)
  resvar <- rss / rdf # Overall residual variance
  
  # Pool over runs
  ar_coef_list <- lapply(cres, function(x) x$ar_coef)
  ar_coef_list <- Filter(function(x) !is.null(x) && length(unlist(x)) > 0, ar_coef_list)

  if (length(cres) > 1) {
    # meta_contrasts expects a list of lists (runs) of lists (contrasts) of tibbles?
    # Or list (runs) of lists (contrasts) where elements are the tibbles?
    # Current: conres_list is list (runs) of lists (contrasts are named elements, values are tibbles)
    # Need to check meta_contrasts implementation.
    # Assuming meta_contrasts can handle the list of lists structure from fit_lm_contrasts_fast.
    
    # meta_betas expects a list of bstats tibbles and event_indices from the first run.
    # Assuming beta_stats_matrix output is compatible.
    meta_con <- meta_contrasts(conres_list)
    # Include all beta indices (event + baseline)
    all_indices <- c(cres[[1]]$event_indices, cres[[1]]$baseline_indices)
    meta_beta <- meta_betas(bstats_list, all_indices)
    
    list(
      contrasts = meta_con,
      betas = meta_beta,
      event_indices = cres[[1]]$event_indices,
      baseline_indices = cres[[1]]$baseline_indices,
      cov.unscaled = Vu, # Using Vu from global design matrix
      sigma = sigma, # Pooled sigma
      rss = rss,     # Pooled rss
      rdf = rdf,     # Pooled rdf
      resvar = resvar, # Pooled resvar
      ar_coef = if (length(ar_coef_list) > 0) ar_coef_list else NULL
    )
  } else {
    # Single run - need to combine contrasts into single tibble format
    # The conres_list[[1]] is a list of contrast tibbles, but we need a single tibble
    single_contrasts <- conres_list[[1]]
    
    if (length(single_contrasts) > 0) {
      # Combine the list of contrast tibbles into a single tibble
      combined_contrasts <- dplyr::bind_rows(single_contrasts)
    } else {
      # Empty contrasts
      combined_contrasts <- tibble::tibble()
    }
    
    list(
      contrasts = combined_contrasts,  # Single tibble with all contrasts
      betas = bstats_list[[1]], # This is the bstats tibble for the single run
      event_indices = cres[[1]]$event_indices,
      baseline_indices = cres[[1]]$baseline_indices,
      cov.unscaled = Vu,
      sigma = cres[[1]]$sigma, # Use run sigma
      rss = cres[[1]]$rss,
      rdf = cres[[1]]$rdf,
      resvar = cres[[1]]$resvar,
      ar_coef = if (length(ar_coef_list) > 0) ar_coef_list[[1]] else NULL
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
  cat(crayon::blue$bold("\n==================================\n"))
  cat(crayon::blue$bold("        fmri_lm_result          \n"))
  cat(crayon::blue$bold("==================================\n\n"))
  
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
  
    
    
