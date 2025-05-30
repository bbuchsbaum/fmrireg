#' Multiresponse bootstrap linear model
#'
#' @description
#' Performs block bootstrap resampling for multiresponse linear models, particularly
#' useful for fMRI time series data where temporal dependencies exist. The function
#' implements a block bootstrap approach to maintain the temporal correlation structure
#' within the data.
#'
#' @details
#' The function performs the following steps:
#' 1. Fits the original linear model
#' 2. Implements block bootstrap resampling of residuals
#' 3. Reconstructs response variables using fitted values and resampled residuals
#' 4. Computes contrasts for each bootstrap sample
#'
#' The block bootstrap approach helps preserve temporal dependencies in the data by
#' resampling blocks of consecutive observations rather than individual observations.
#'
#' @param form Formula for the linear model. Required if modmat is NULL.
#' @param data_env Environment containing the data for the linear model.
#' @param conlist List of contrasts to be computed for each bootstrap sample.
#' @param vnames Vector of variable names.
#' @param fcon Contrasts for fixed effects.
#' @param modmat Optional pre-computed model matrix. If provided, form is ignored.
#' @param block_size Size of the blocks for the bootstrap (default: 30).
#'        Should be large enough to capture temporal dependencies but small enough
#'        to allow sufficient randomization.
#' @param boot_rows Logical flag indicating whether to bootstrap rows (default: FALSE).
#' @param nboot Number of bootstrap iterations (default: 100).
#' @param event_indices Indices of events for computing beta covariances.
#'
#' @return A list containing:
#' \itemize{
#'   \item original: The fitted original model with contrasts
#'   \item con_cov: Covariance matrices for contrasts (if contrasts provided)
#'   \item beta_cov: Covariance matrices for beta estimates
#'   \item nboot: Number of bootstrap iterations performed
#'   \item bootstrap: Logical indicating this is a bootstrap result
#' }
#'
#' @examples
#' \donttest{
#' # Simple example with synthetic data
#' X <- model.matrix(~ x1 + x2, data = data.frame(x1 = rnorm(100), x2 = rnorm(100)))
#' y <- matrix(rnorm(100 * 3), 100, 3)  # 3 response variables
#' result <- multiresponse_bootstrap_lm(modmat = X, data_env = list(.y = y),
#'                                     nboot = 100, block_size = 20)
#' }
#' @keywords internal
#' @importFrom foreach foreach %dopar%
multiresponse_bootstrap_lm <- function(form, data_env, 
                                     conlist, 
                                     vnames, 
                                     fcon, modmat=NULL, 
                                     block_size=30,
                                     boot_rows=FALSE,
                                     nboot=100,
                                     event_indices) {
  
  # Input validation
  if (is.null(modmat) && missing(form)) {
    stop("Either 'form' or 'modmat' must be provided")
  }
  
  if (!is.numeric(block_size) || block_size < 1) {
    stop("block_size must be a positive integer")
  }
  
  if (!is.numeric(nboot) || nboot < 1) {
    stop("nboot must be a positive integer")
  }
  
  if (!is.null(modmat)) {
    if (!inherits(modmat, "matrix")) {
      stop("modmat must be a matrix")
    }
    if (is.null(data_env$.y)) {
      stop("data_env$.y must contain response variables when using modmat")
    }
    if (nrow(modmat) != nrow(data_env$.y)) {
      stop("Number of rows in modmat must match number of rows in response variables")
    }
  }
  
  # Fit original model
  lm.orig <- if (is.null(modmat)) {
    lm(as.formula(form), data=data_env)
  } else {
    lm.fit(modmat, data_env$.y)
  }
  
  # Prepare for bootstrap
  yhat <- fitted(lm.orig)
  # Ensure yhat is a matrix
  if (is.vector(yhat)) {
    yhat <- matrix(yhat, ncol = 1)
  }
  rows <- 1:nrow(yhat)
  nblocks <- as.integer(length(rows)/block_size)
  
  # Create blocks for resampling
  blocks <- split(rows, rep(1:nblocks, each=length(rows)/nblocks, length.out=length(rows)))
  maxind <- max(blocks[[length(blocks)]])
  if (maxind < nrow(yhat)) {
    last_block <- (maxind+1):nrow(yhat)
    blocks <- c(blocks, list(last_block))
  }
  
  # Perform bootstrap
  boot_res <- vector("list", nboot)
  for (b in 1:nboot) {
    sam_blocks <- sample(1:length(blocks), replace=TRUE)
    samples <- unlist(blocks[sam_blocks])
    
    # Adjust sample size if necessary
    if (length(samples) > nrow(yhat)) {
      samples <- samples[1:nrow(yhat)]
    } else if (length(samples) < nrow(yhat)) {
      delta <- nrow(yhat) - length(samples) 
      samples <- c(samples, sample(rows, delta, replace=TRUE))
    }
    
    # Generate bootstrap sample
    if (is.matrix(lm.orig$residuals)) {
      rstar <- lm.orig$residuals[samples, , drop = FALSE]
    } else {
      rstar <- matrix(lm.orig$residuals[samples], ncol = 1)
    }
    ynew <- yhat + rstar
    
    # Fit bootstrap model and compute contrasts
    lm.boot <- lm.fit(modmat, ynew)
    boot_res[[b]] <- fit_lm_contrasts(lm.boot, conlist, fcon, vnames, se=FALSE) 
  }
  
  # Compute covariance matrices
  con_cov <- if (length(conlist) > 0) {
    lapply(names(boot_res[[1]]$contrasts), function(bname) {
      estimates <- lapply(boot_res, function(br) {
        # Access the nested structure: data[[1]]$estimate
        if (!is.null(br$contrasts[[bname]]) && nrow(br$contrasts[[bname]]) > 0) {
          br$contrasts[[bname]]$data[[1]]$estimate
        } else {
          NULL
        }
      })
      # Filter out NULL values and combine
      estimates <- estimates[!sapply(estimates, is.null)]
      if (length(estimates) > 0) {
        bm <- do.call(rbind, estimates)
        cov(bm)
      } else {
        NULL
      }
    })
  } else {
    NULL
  }
  
  beta_cov <- lapply(event_indices, function(i) {
    estimates <- lapply(boot_res, function(br) {
      # Access the nested structure: data[[1]]$estimate[[1]]
      # beta_stats returns estimates wrapped in an additional list
      if (!is.null(br$bstats) && nrow(br$bstats) > 0) {
        estimate_matrix <- br$bstats$data[[1]]$estimate[[1]]
        # Extract only the row corresponding to event_indices[i] if it exists
        if (is.matrix(estimate_matrix) && nrow(estimate_matrix) >= i) {
          estimate_matrix[i, , drop = FALSE]
        } else {
          NULL
        }
      } else {
        NULL
      }
    })
    # Filter out NULL values and combine
    estimates <- estimates[!sapply(estimates, is.null)]
    if (length(estimates) > 0) {
      bm <- do.call(rbind, estimates)
      cov(bm)
    } else {
      NULL
    }
  })
  
  # Return results
  orig <- fit_lm_contrasts(lm.orig, conlist, fcon, vnames) %>%
    purrr::list_modify(
      con_cov = con_cov,
      beta_cov = beta_cov,
      nboot = nboot,
      bootstrap = TRUE
    )

  orig
}
