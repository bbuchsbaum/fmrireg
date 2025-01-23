


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
  
  # Extract design matrices for event and baseline terms
  eterms <- lapply(event_terms(x), function(term) as.matrix(design_matrix(term, blocknum)))
  bterms <- lapply(baseline_terms(x), function(term) as.matrix(design_matrix(term, blocknum)))
  
  # Compute indices for event and baseline terms
  num_event_cols <- sum(sapply(eterms, ncol))
  num_baseline_cols <- sum(sapply(bterms, ncol))
  
  eterm_indices <- 1:num_event_cols
  bterm_indices <- (num_event_cols + 1):(num_event_cols + num_baseline_cols)
  
  # Combine term matrices
  term_matrices <- c(eterms, bterms)
  names(term_matrices) <- names(terms(x))
  
  # Collect variable names
  vnames <- unlist(lapply(term_matrices, colnames))
  
  # Set attributes
  attr(term_matrices, "event_term_indices") <- eterm_indices
  attr(term_matrices, "baseline_term_indices") <- bterm_indices
  attr(term_matrices, "blocknum") <- blocknum
  attr(term_matrices, "varnames") <- vnames
  
  return(term_matrices)
}


#' #' Extract term matrices for an fMRI model
#' #'
#' #' This function extracts the term matrices for an fMRI model, which consists of event-related terms
#' #' and baseline-related terms. The term matrices are used for building the design matrix in fMRI data analysis.
#' #'
#' #' @param x An object of class "fmri_model" containing the event and baseline models.
#' #' @param blocknum (Optional) A numeric vector specifying the block numbers to be included in the term matrices.
#' #'                 By default, all unique block numbers in the event model are included.
#' #' @return A named list of term matrices, with event terms followed by baseline terms.
#' #'         Attributes "event_term_indices" and "baseline_term_indices" store the indices of event and baseline terms,
#' #'         "blocknum" stores the block numbers, and "varnames" stores the variable names.
#' #' @export
#' #' @seealso fmri_model
#' term_matrices.fmri_model <- function(x, blocknum=NULL) {
#'   eterms <- lapply(event_terms(x), 
#'                    function(x) as.matrix(design_matrix(x, blocknum)))
#'   
#'   bterms <- lapply(baseline_terms(x), 
#'                    function(x) as.matrix(design_matrix(x, blocknum)))
#'   
#'   if (is.null(blocknum)) {
#'     blocknum <- sort(unique(x$event_model$blockids))
#'   }
#'   
#'   start <- 1
#'   eterm_indices <- 1:sum(map_int(eterms, ncol))
#'   start <- length(eterm_indices) +1
#'   bterm_indices <- start:(start+sum(map_int(bterms, ncol)))
#'   #browser()
#'   term_matrices <- c(eterms, bterms)
#'   names(term_matrices) <- names(terms(x))
#'   
#'   vnames <- unlist(lapply(term_matrices, colnames))
#'   
#'   attr(term_matrices, "event_term_indices") <- eterm_indices
#'   attr(term_matrices, "baseline_term_indices") <- bterm_indices
#'   attr(term_matrices, "blocknum") <- blocknum 
#'   attr(term_matrices, "varnames") <- vnames
#'   term_matrices
#' }


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
    formula = formula,
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
                    strategy = c("runwise", "chunkwise"), nchunks = 10, ...) {
  
  strategy <- match.arg(strategy)
  
  # Error checking
  assert_that(is.formula(formula), msg = "'formula' must be a formula")
  assert_that(is.formula(block), msg = "'block' must be a formula")
  assert_that(inherits(dataset, "fmri_dataset"), msg = "'dataset' must be an 'fmri_dataset'")
  assert_that(is.numeric(durations), msg = "'durations' must be numeric")
  assert_that(is.logical(drop_empty), msg = "'drop_empty' must be logical")
  assert_that(is.logical(robust), msg = "'robust' must be logical")
  if (strategy == "chunkwise") {
    assert_that(is.numeric(nchunks) && nchunks > 0, msg = "'nchunks' must be a positive number")
  }
  
  model <- create_fmri_model(formula, block, baseline_model, dataset, durations = durations, drop_empty = drop_empty)
  ret <- fmri_lm_fit(model, dataset, strategy, robust, nchunks, ...)
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
#' @param ... Additional arguments.
#' @return A fitted fMRI linear regression model with the specified fitting strategy.
#' @keywords internal
#' @seealso \code{\link{fmri_lm}}, \code{\link{fmri_model}}, \code{\link{fmri_dataset}}
fmri_lm_fit <- function(fmrimod, dataset, strategy = c("runwise", "chunkwise"), 
                        robust = FALSE, nchunks = 10, ...) {
  strategy <- match.arg(strategy)
  
  # Error checking
  assert_that(inherits(fmrimod, "fmri_model"), msg = "'fmrimod' must be an 'fmri_model' object")
  assert_that(inherits(dataset, "fmri_dataset"), msg = "'dataset' must be an 'fmri_dataset' object")
  assert_that(is.logical(robust), msg = "'robust' must be logical")
  if (strategy == "chunkwise") {
    assert_that(is.numeric(nchunks) && nchunks > 0, msg = "'nchunks' must be a positive number")
  }
  
  conlist <- unlist(contrast_weights(fmrimod$event_model), recursive = FALSE)
  fcons <- Fcontrasts(fmrimod$event_model)
  
  result <- switch(strategy,
                   "runwise" = runwise_lm(dataset, fmrimod, conlist, fcons, robust = robust, ...),
                   "chunkwise" = chunkwise_lm(dataset, fmrimod, conlist, fcons, nchunks, robust = robust, ...))
  
  ret <- list(
    result = result,
    model = fmrimod,
    strategy = strategy,
    bcons = conlist,
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
    
    list(pred = as.matrix(yh), design = as_tibble(design))
  })
  
  names(pred) <- names(eterms)
  return(pred)
}

# fitted_errors <- function(SE, B, nb) {
#   se_list <- lapply(seq_len(ncol(SE)), function(j) {
#     cov_matrix <- B[, j, drop = FALSE] %*% t(B[, j, drop = FALSE]) * SE[, j]^2
#     se <- sqrt(pmax(diag(G %*% cov_matrix %*% t(G)), 0)) # Dim: n_timepoints
#     se
#   })
# 
# }
  


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
pull_stat <- function(x, type, element) {
  if (type == "betas") {
    ret <- x$result$betas$data[[1]][[element]][[1]]
    ret <- ret[, x$result$event_indices, drop = FALSE]
    colnames(ret) <- conditions(x$model$event_model)
    suppressMessages(as_tibble(ret, .name_repair = "check_unique"))
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
#' This function extracts model coefficients (estimates) from an \code{fmri_lm} object.
#'
#' @param object An \code{fmri_lm} object.
#' @param type The type of coefficients to extract: \code{"betas"} or \code{"contrasts"}. Default is \code{"betas"}.
#' @param recon Logical. If \code{TRUE}, reconstructs the coefficients into a neuroimaging volume. Default is \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#' @return A tibble or matrix of coefficients.
#' @export
coef.fmri_lm <- function(object, type = c("betas", "contrasts"), recon = FALSE, ...) {
  type <- match.arg(type)
  res <- if (type == "betas") {
    pull_stat(object, "betas", "estimate")
  } else if (type == "contrasts") {
    pull_stat(object, "contrasts", "estimate")
  }
  
  # Reconstruction functionality can be added here if necessary
  # if (recon && inherits(object$dataset, "fmri_dataset")) {
  #   m <- get_mask(object$dataset)
  #   sp <- space(m)
  #   SparseNeuroVec(as.matrix(res), neuroim2::add_dim(sp, ncol(res)), mask = m)
  # } else {
  #   res
  # }
  
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
#' @param x An \code{fmri_lm} object.
#' @param type The type of standard errors to extract: \code{"estimates"} or \code{"contrasts"}.
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


#' Print an fmri_lm Object
#'
#' This function prints a summary of an \code{fmri_lm} object.
#'
#' @param x An \code{fmri_lm} object.
#' @param ... Additional arguments (currently unused).
#' @export
print.fmri_lm <- function(x, ...) {
  cat("fmri_lm model:\n", as.character(x$model$event_model$model_spec$formula), "\n")
  cat("  Baseline parameters:", ncol(design_matrix(x$model$baseline_model)), "\n")
  cat("  Design parameters:", ncol(design_matrix(x$model$event_model)), "\n")
  cat("  Contrasts:", paste(names(x$bcons), collapse = ", "), "\n")
}

# summary.fmri_lm <- function(x, type=c("coef", "contrasts", "Fcontrasts")) {
#   type <- match.arg(type)
#   if (type == "coef") {
#     betas=x$result$betas
#     list(
#       estimate=betas$estimate(),
#       se=betas$se(),
#       stat=betas$stat(),
#       prob=betas$prob())
#   } else if (type == "contrasts") {
#     x$result$contrasts
#   } else if (type == "Fcontrasts") {
#     x$result$Fcontrasts
#   } else {
#     stop()
#   }
# }



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
      estimate_contrast(con, fit, attr(con, "term_indices"))
    })
    names(ret) <- names(conlist)
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
  
  fit_lm_contrasts(lm_fit, conlist, fcon, vnames)
}



#' Unpack Chunkwise Results
#'
#' This function processes and unpacks the results of chunkwise analysis in fMRI data.
#'
#' @param cres The results of the chunkwise analysis.
#' @param event_indices The indices of the event-related effects.
#' @param baseline_indices The indices of the baseline-related effects.
#' @return A list containing the betas, contrasts, and statistical information from the chunkwise results.
#' @keywords internal
unpack_chunkwise <- function(cres, event_indices, baseline_indices) {
  cbetas <- lapply(cres, function(x) x$bstats) %>% dplyr::bind_rows()
  dat <- cbetas$data %>% dplyr::bind_rows()
  estimate <- dplyr::as_tibble(do.call(rbind, dat$estimate))
  se <- dplyr::as_tibble(do.call(rbind, dat$se))
  stat <- dplyr::as_tibble(do.call(rbind, dat$stat))
  prob <- dplyr::as_tibble(do.call(rbind, dat$prob))
  sigma <- dplyr::tibble(sigma = unlist(dat$sigma))
  cbetas <- dplyr::tibble(
    type = cbetas$type[1],
    stat_type = cbetas$stat_type[1],
    df.residual = cbetas$df.residual[1],
    conmat = list(NULL),
    colind = list(NULL),
    data = list(
      dplyr::tibble(
        estimate = list(estimate),
        se = list(se),
        stat = list(stat),
        prob = list(prob),
        sigma = list(sigma)
      )
    )
  )
  
  ncon <- length(cres[[1]]$contrasts)
  
  if (ncon > 0) {
    contab <- lapply(cres, function(x) x$contrasts %>% dplyr::bind_rows()) %>% dplyr::bind_rows()
    gsplit <- contab %>% dplyr::group_by(name, type) %>% dplyr::group_split()
    con <- lapply(gsplit, function(g) {
      dat <- g$data %>% dplyr::bind_rows()
      g %>% dplyr::select(-data) %>% dplyr::slice_head() %>% dplyr::mutate(data = list(dat))
    }) %>% dplyr::bind_rows()
  } else {
    con <- dplyr::tibble()
  }
  
  list(
    betas = cbetas,
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
#' @param conlist The list of contrasts used in the analysis.
#' @param fcon The F-contrasts used in the analysis.
#' @param nchunks The number of chunks to divide the dataset into.
#' @param robust Logical. Whether to use robust linear modeling (default is \code{FALSE}).
#' @param verbose Logical. Whether to display progress messages (default is \code{FALSE}).
#' @return A list containing the unpacked chunkwise results.
#' @keywords internal
chunkwise_lm.fmri_dataset <- function(dset, model, conlist, fcon, nchunks, robust = FALSE, verbose = FALSE) {
  chunks <- exec_strategy("chunkwise", nchunks = nchunks)(dset)
  form <- get_formula(model)
  tmats <- term_matrices(model)
  data_env <- list2env(tmats)
  data_env[[".y"]] <- rep(0, nrow(tmats[[1]]))
  modmat <- model.matrix(as.formula(form), data_env)
  Qr <- qr(modmat)
  Vu <- chol2inv(Qr$qr)
  
  lmfun <- if (robust) multiresponse_rlm else multiresponse_lm
  
  ym <- NULL
  cres <- foreach(ym = chunks, .combine = 'c', .verbose = verbose) %dopar% {
    if (verbose) message("Processing chunk ", ym$chunk_num)
    data_env[[".y"]] <- as.matrix(ym$data)
    ret <- lmfun(form, data_env, conlist, attr(tmats, "varnames"), fcon, modmat = modmat)
    ret$rss <- colSums(as.matrix(ret$fit$residuals^2))
    ret$rdf <- ret$fit$df.residual
    ret$resvar <- ret$rss / ret$rdf
    ret$sigma <- sqrt(ret$resvar)
    ret
  }
  
  event_indices = attr(tmats, "event_term_indices")
  baseline_indices = attr(tmats, "baseline_term_indices")
  
  out <- unpack_chunkwise(cres, event_indices, baseline_indices)
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
#' @param conlist The list of contrasts used in the analysis.
#' @param fcon The F-contrasts used in the analysis.
#' @param robust Logical. Whether to use robust linear modeling (default is \code{FALSE}).
#' @param verbose Logical. Whether to display progress messages (default is \code{FALSE}).
#' @return A list containing the combined results from runwise linear model analysis.
#' @keywords internal
runwise_lm <- function(dset, model, conlist, fcon, robust = FALSE, verbose = FALSE) {
  lmfun <- if (robust) multiresponse_rlm else multiresponse_lm
  
  # Get an iterator of data chunks
  chunks <- exec_strategy("runwise")(dset)
  
  form <- get_formula(model)
  modmat <- design_matrix(model)
  Qr <- qr(modmat)
  Vu <- chol2inv(Qr$qr)
  
  # Iterate over each data chunk
  cres <- foreach(ym = chunks, .combine = 'c', .verbose = verbose) %dopar% {
    if (verbose) message("Processing run ", ym$chunk_num)
    tmats <- term_matrices(model, ym$chunk_num)
    
    data_env <- list2env(tmats)
    data_env[[".y"]] <- as.matrix(ym$data)
    ret <- lmfun(form, data_env, conlist, attr(tmats, "varnames"), fcon)
    
    rss <- colSums(as.matrix(ret$fit$residuals^2))
    rdf <- ret$fit$df.residual
    resvar <- rss / rdf
    sigma <- sqrt(resvar)
    
    list(
      conres = ret$contrasts,
      bstats = ret$bstats,
      event_indices = attr(tmats, "event_term_indices"),
      baseline_indices = attr(tmats, "baseline_term_indices"),
      rss = rss,
      rdf = rdf,
      resvar = resvar,
      sigma = sigma
    )
  }
  
  # Combine results
  bstats <- lapply(cres, `[[`, "bstats")
  conres <- lapply(cres, `[[`, "conres")
  
  # Compute overall statistics
  sigma <- colMeans(do.call(rbind, lapply(cres, `[[`, "sigma")))
  rss <- colSums(do.call(rbind, lapply(cres, `[[`, "rss")))
  rdf <- colSums(do.call(rbind, lapply(cres, `[[`, "rdf")))
  
  # Pool over runs
  if (length(cres) > 1) {
    meta_con <- meta_contrasts(conres)
    meta_beta <- meta_betas(bstats, cres[[1]]$event_indices)
    list(
      contrasts = meta_con,
      betas = meta_beta,
      event_indices = cres[[1]]$event_indices,
      baseline_indices = cres[[1]]$baseline_indices,
      cov.unscaled = Vu,
      sigma = sigma,
      rss = rss,
      rdf = rdf,
      resvar = sigma^2
    )
  } else {
    list(
      contrasts = conres[[1]],
      betas = bstats[[1]],
      event_indices = cres[[1]]$event_indices,
      baseline_indices = cres[[1]]$baseline_indices,
      cov.unscaled = Vu,
      sigma = sigma,
      rss = rss,
      rdf = rdf,
      resvar = sigma^2
    )
  }
}
  
    
    



