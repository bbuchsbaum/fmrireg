# Spatial FDR: Structure-Adaptive Weighted BH for Multiple Comparisons
# Fast, spatially-aware FDR control for fMRI group analyses

#' Spatially-Aware Multiple Comparisons Correction
#'
#' Performs spatially-aware FDR control using structure-adaptive weighted 
#' Benjamini-Hochberg (SABHA-style) procedure. This method leverages spatial 
#' structure in the data to increase power while controlling the false discovery rate.
#'
#' @param z Numeric vector of Z-values (one per feature). Provide either z or p, not both.
#' @param p Numeric vector of p-values (two-sided). Provide either z or p, not both.
#' @param group Integer, character or factor vector of group IDs for each feature
#'   (e.g., parcel IDs, block IDs). Must have same length as z or p. Groups are
#'   numbered 1..G by sorted unique value, or by factor level (unused levels
#'   included) when \code{group} is a factor. This numbering indexes
#'   \code{neighbors}, \code{pi0_raw}, \code{pi0_smooth} and \code{m_g}.
#'   Features with \code{NA} group are excluded from \eqn{\pi_0} estimation and
#'   receive weight 1 before normalisation.
#' @param alpha Numeric scalar; FDR level to control (default: 0.05)
#' @param tau Numeric scalar; Storey threshold in (0,1) for \eqn{\pi_0} estimation (default: 0.5). 
#'   Higher values are more conservative.
#' @param lambda Numeric scalar; smoothing strength for groupwise \eqn{\pi_0} across neighbors (default: 1.0).
#'   Set to 0 for no smoothing, higher values for more smoothing.
#' @param neighbors Optional list of length G (number of groups) where element g
#'   is an integer vector of 1-based neighbor IDs of group g, both in the group
#'   numbering described under \code{group}. Used for spatial smoothing of \eqn{\pi_0}.
#' @param min_pi0 Numeric scalar; lower bound for \eqn{\pi_0} to stabilize weights (default: 0.05).
#'   Prevents infinite weights.
#' @param empirical_null Logical; if TRUE, estimate null distribution parameters
#'   (\eqn{\mu_0}, \eqn{\sigma_0}) from the z-values (default: TRUE). See Details.
#' @param verbose Logical; print progress messages (default: FALSE)
#'
#' @return Object of class "spatial_fdr_result" containing:
#'   \item{reject}{Logical vector indicating rejected hypotheses (discoveries)}
#'   \item{q}{Numeric vector of FDR-adjusted p-values (q-values)}
#'   \item{p}{Numeric vector of two-sided p-values used for testing}
#'   \item{weights}{Numeric vector of normalized weights used in weighted BH}
#'   \item{pi0_raw}{Numeric vector of raw \eqn{\pi_0} estimates per group}
#'   \item{pi0_smooth}{Numeric vector of smoothed \eqn{\pi_0} estimates per group}
#'   \item{threshold}{Numeric scalar; BH threshold used for rejection}
#'   \item{k}{Integer scalar; number of rejections}
#'   \item{mu0}{Numeric scalar; estimated null mean (if empirical_null = TRUE)}
#'   \item{sigma0}{Numeric scalar; estimated null SD (if empirical_null = TRUE)}
#'   \item{groups}{Integer vector of group IDs (1..G) in the numbering described
#'     under \code{group}}
#'   \item{G}{Integer scalar; number of groups}
#'   \item{m_g}{Integer vector; number of features with a finite p-value per group}
#'   \item{alpha}{Numeric scalar; FDR level used}
#'   \item{coef_name}{Character scalar; name of coefficient (for S3 method)}
#'
#' @details
#' This function implements a spatially-aware multiple testing procedure that:
#' 1. Estimates the proportion of null hypotheses (pi0) within spatial groups
#' 2. Optionally smooths these estimates across neighboring groups
#' 3. Uses the pi0 estimates to weight the Benjamini-Hochberg procedure
#' 4. Provides more power in regions with true signal while maintaining FDR control
#'
#' Weights are \eqn{1/\pi_0} of each feature's group, normalised to sum to the
#' number of finite p-values, and the weighted BH procedure of Genovese,
#' Roeder & Wasserman (2006) is applied: reject \eqn{p_i \le w_i \alpha k / n}.
#'
#' \strong{Empirical null.} When \code{z} is supplied with
#' \code{empirical_null = TRUE}, \eqn{\mu_0} and \eqn{\sigma_0} are estimated by
#' Efron's (2004) central matching: a quadratic is fitted by Poisson regression
#' to the log histogram counts of z within \eqn{\pm 1.5} robust SDs of the
#' median. Because a minority of non-null features only changes the intercept
#' of the log density near the centre and not its curvature, this estimate is
#' consistent under \eqn{N(0,1)} and stays close to it when, for example, 10\%
#' of features are true signals. With fewer than 500 finite z-values (or if the
#' fit fails) the median and normal-consistent MAD of all finite z-values are
#' used instead; with fewer than 10 the theoretical null is used.
#'
#' The computation is done by \code{fmrigds::spatial_fdr_sabha()} (also
#' available in fmrigds pipelines as the post-hoc \code{"fdr:sabha"}, and as
#' \code{"fmrireg:fdr:spatial"}); this function keeps fmrireg's argument names,
#' the \code{fmri_meta} method and the \code{spatial_fdr_result} class.
#'
#' The method is particularly effective for:
#' - Voxelwise analyses with spatial clustering of signal
#' - Parcel-based analyses with anatomical or functional grouping
#' - Any scenario where hypotheses can be grouped spatially or functionally
#'
#' @references
#' \itemize{
#'   \item Benjamini & Hochberg (1995). Controlling the false discovery rate.
#'   \item Storey (2002). A direct approach to false discovery rates.
#'   \item Hu et al. (2010). False discovery rate control with groups (SABHA).
#'   \item Genovese, Roeder & Wasserman (2006). False discovery control with
#'     p-value weighting.
#'   \item Efron (2004). Large-scale simultaneous hypothesis testing: the choice
#'     of a null hypothesis.
#' }
#'
#' @examples
#' # Simple synthetic example
#' set.seed(123)
#' n <- 1000
#' z_scores <- c(rnorm(800), rnorm(200, mean = 2))  # 200 true signals
#' group_ids <- rep(1:10, each = 100)  # 10 groups of 100 features each
#' 
#' # Basic usage without spatial smoothing
#' result <- spatial_fdr(z = z_scores, group = group_ids, alpha = 0.05)
#' summary(result)
#' 
#' # Create simple neighbor structure (each group neighbors with adjacent groups)
#' neighbors <- lapply(1:10, function(i) {
#'   c(if(i > 1) i-1, if(i < 10) i+1)
#' })
#' 
#' # With spatial smoothing
#' result_smooth <- spatial_fdr(z = z_scores, group = group_ids, 
#'                             neighbors = neighbors, lambda = 1.0)
#' print(result_smooth)
#' 
#' @seealso \code{\link{create_3d_blocks}}
#'
#' @export
spatial_fdr <- function(z = NULL, 
                       p = NULL, 
                       group = NULL,
                       alpha = 0.05, 
                       tau = 0.5,
                       lambda = 1.0, 
                       neighbors = NULL,
                       min_pi0 = 0.05,
                       empirical_null = TRUE,
                       verbose = FALSE) {
  
  # Check if first argument is an fmri_meta object
  if (!is.null(z) && inherits(z, "fmri_meta")) {
    # Handle fmri_meta objects specially
    # Arguments shift: z is object, p becomes coef, group stays group
    object <- z
    coef <- if (!is.null(p)) p else 1
    
    # Call the fmri_meta handler
    return(spatial_fdr_fmri_meta(object, coef = coef, group = group, 
                                 alpha = alpha, tau = tau, lambda = lambda,
                                 neighbors = neighbors, min_pi0 = min_pi0,
                                 empirical_null = empirical_null, 
                                 verbose = verbose))
  }
  
  # The procedure lives in fmrigds (spatial_fdr_sabha(), also the post-hoc
  # "fdr:sabha"); it validates the inputs with the same messages.
  if (verbose && !is.null(group)) {
    cat("Spatial FDR: Processing", length(z %||% p), "features in",
        if (is.factor(group)) nlevels(group) else length(unique(group[!is.na(group)])),
        "groups\n")
  }
  res <- fmrigds::spatial_fdr_sabha(
    z = if (!is.null(z)) as.numeric(z),
    p = if (!is.null(p)) as.numeric(p),
    group = group, alpha = alpha, tau = tau, lambda = lambda,
    neighbors = neighbors, min_pi0 = min_pi0, empirical_null = empirical_null
  )
  if (verbose) {
    if (res$mu0 != 0 || res$sigma0 != 1) {
      cat("  Empirical null: mu0 = ", round(res$mu0, 3),
          ", sigma0 = ", round(res$sigma0, 3), "\n", sep = "")
    }
    cat("  Results: ", sum(res$reject, na.rm = TRUE), " discoveries at FDR = ",
        alpha, "\n", sep = "")
  }
  fields <- c("reject", "q", "p", "weights", "pi0_raw", "pi0_smooth",
              "threshold", "k", "mu0", "sigma0", "groups", "G", "m_g",
              "alpha", "tau", "lambda")
  structure(res[fields], class = "spatial_fdr_result")
}

#' Print Spatial FDR Results
#'
#' @param x A spatial_fdr_result object
#' @param ... Additional print arguments
#' @return Invisibly returns the input object x
#' @export
print.spatial_fdr_result <- function(x, ...) {
  cat("Spatial FDR Results\n")
  cat("===================\n")
  cat("Features:", length(x$p), "\n")
  cat("Groups:", x$G, "\n")
  cat("FDR level:", x$alpha, "\n")
  cat("Discoveries:", sum(x$reject, na.rm = TRUE), 
      "(", round(100 * mean(x$reject, na.rm = TRUE), 2), "%)\n")
  cat("Threshold:", x$threshold, "\n")
  
  if (x$mu0 != 0 || x$sigma0 != 1) {
    cat("Empirical null: mu =", round(x$mu0, 3), ", sigma =", round(x$sigma0, 3), "\n")
  }
  
  cat("\nPi0 summary:\n")
  cat("  Range:", round(min(x$pi0_smooth), 3), "-", round(max(x$pi0_smooth), 3), "\n")
  cat("  Mean:", round(mean(x$pi0_smooth), 3), "\n")
  cat("  Median:", round(median(x$pi0_smooth), 3), "\n")
  
  invisible(x)
}

#' Summary of Spatial FDR Results
#'
#' @param object A spatial_fdr_result object
#' @param ... Additional summary arguments
#' @return Invisibly returns the input object
#' @export
summary.spatial_fdr_result <- function(object, ...) {
  print(object)
  
  # Group-level summary
  cat("\nGroup-level discoveries:\n")
  
  disc_by_group <- numeric(object$G)
  total_by_group <- object$m_g
  
  for (i in seq_along(object$reject)) {
    if (!is.na(object$groups[i]) && object$reject[i]) {
      g <- object$groups[i]
      if (g >= 1 && g <= object$G) {
        disc_by_group[g] <- disc_by_group[g] + 1
      }
    }
  }
  
  prop_disc <- disc_by_group / pmax(total_by_group, 1)
  
  cat("  Groups with discoveries:", sum(disc_by_group > 0), "/", object$G, "\n")
  cat("  Max proportion in group:", round(max(prop_disc), 3), "\n")
  cat("  Mean proportion in groups with signal:", 
      round(mean(prop_disc[disc_by_group > 0]), 3), "\n")
  
  invisible(object)
}

#' Create 3D Blocks for Voxelwise Analysis
#'
#' Creates spatial blocks from a 3D mask for use with spatial_fdr.
#' This is useful for voxelwise analyses where you want to group
#' nearby voxels together for more powerful multiple comparisons correction.
#'
#' @param mask Numeric 3D array or NeuroVol object defining the brain mask.
#'   Non-zero values indicate voxels to include.
#' @param block_size Integer vector of length 3 specifying block dimensions
#'   in voxels (default: c(10, 10, 10))
#' @param connectivity Integer scalar; type of connectivity for neighbors: 
#'   6 (face connectivity) or 26 (face, edge, and corner connectivity). Default: 26
#'
#' @return List with components:
#'   \item{group_id}{Integer vector of group IDs for each voxel in the mask}
#'   \item{neighbors}{List of length n_groups where element i contains 
#'     integer vector of 1-based neighbor IDs for group i}
#'   \item{n_groups}{Integer scalar; total number of groups created}
#'   \item{block_size}{Block size used}
#'
#' @examples
#' mask <- array(c(1L, 1L, 0L,
#'                 1L, 1L, 0L,
#'                 0L, 0L, 0L), dim = c(3, 3, 1))
#' blocks <- create_3d_blocks(mask, block_size = c(2, 2, 1))
#' blocks$n_groups
#'
#' @export
create_3d_blocks <- function(mask, block_size = c(10, 10, 10), connectivity = 26) {
  
  # Convert mask to array if needed (also accept file path)
  if (is.character(mask)) {
    if (!file.exists(mask)) {
      stop("mask path does not exist: ", mask, call. = FALSE)
    }
    vol <- try(neuroim2::read_vol(mask), silent = TRUE)
    if (inherits(vol, "try-error")) {
      stop("Failed to read mask NIfTI from path: ", mask, call. = FALSE)
    }
    mask_array <- as.array(vol)
  } else if (inherits(mask, "NeuroVol")) {
    mask_array <- as.array(mask)
  } else if (is.array(mask) && length(dim(mask)) == 3) {
    mask_array <- mask
  } else {
    stop("mask must be a 3D array, NeuroVol, or file path", call. = FALSE)
  }
  
  dims <- dim(mask_array)
  
  # Create block grid
  n_blocks <- ceiling(dims / block_size)
  
  # Assign voxels to blocks
  voxel_indices <- which(mask_array > 0)
  n_voxels <- length(voxel_indices)
  
  # Convert linear indices to subscripts
  coords <- arrayInd(voxel_indices, dims)
  
  # Compute block IDs
  block_coords <- ceiling(coords / matrix(block_size, nrow = n_voxels, ncol = 3, byrow = TRUE))
  
  # Convert block coordinates to linear block ID
  group_id <- (block_coords[, 3] - 1) * n_blocks[1] * n_blocks[2] +
             (block_coords[, 2] - 1) * n_blocks[1] +
              block_coords[, 1]
  
  # Compress to consecutive IDs
  group_id <- as.integer(factor(group_id))
  n_groups <- max(group_id)
  
  # Build neighbor list
  if (connectivity %in% c(6, 26)) {
    neighbors <- build_block_neighbors(n_blocks, connectivity)
    
    # Map to compressed IDs
    old_to_new <- integer(prod(n_blocks))
    for (i in seq_along(group_id)) {
      block_linear <- (block_coords[i, 3] - 1) * n_blocks[1] * n_blocks[2] +
                     (block_coords[i, 2] - 1) * n_blocks[1] +
                      block_coords[i, 1]
      old_to_new[block_linear] <- group_id[i]
    }
    
    # Remap neighbors
    neighbors_compressed <- vector("list", n_groups)
    for (g in seq_len(n_groups)) {
      old_id <- which(old_to_new == g)[1]
      if (!is.na(old_id) && old_id <= length(neighbors)) {
        old_neighbors <- neighbors[[old_id]]
        new_neighbors <- unique(old_to_new[old_neighbors])
        new_neighbors <- new_neighbors[new_neighbors > 0 & new_neighbors != g]
        neighbors_compressed[[g]] <- new_neighbors
      }
    }
    
    neighbors <- neighbors_compressed
  } else {
    neighbors <- NULL
  }
  
  list(
    group_id = group_id,
    neighbors = neighbors,
    n_groups = n_groups,
    block_size = block_size,
    n_voxels = n_voxels
  )
}

#' Build Block Neighbor Adjacency
#'
#' @param n_blocks Grid dimensions
#' @param connectivity 6 or 26
#' @return List of neighbor indices
#' @keywords internal
build_block_neighbors <- function(n_blocks, connectivity = 26) {
  n_total <- prod(n_blocks)
  neighbors <- vector("list", n_total)
  
  # Define offsets for connectivity
  if (connectivity == 6) {
    # Face neighbors only
    offsets <- rbind(
      c(-1, 0, 0), c(1, 0, 0),
      c(0, -1, 0), c(0, 1, 0),
      c(0, 0, -1), c(0, 0, 1)
    )
  } else {
    # All 26 neighbors
    offsets <- as.matrix(expand.grid(
      x = c(-1, 0, 1),
      y = c(-1, 0, 1),
      z = c(-1, 0, 1)
    ))
    offsets <- offsets[rowSums(abs(offsets)) > 0, ]
  }
  
  # Build adjacency
  for (k in seq_len(n_blocks[3])) {
    for (j in seq_len(n_blocks[2])) {
      for (i in seq_len(n_blocks[1])) {
        idx <- (k - 1) * n_blocks[1] * n_blocks[2] + 
               (j - 1) * n_blocks[1] + i
        
        nb <- integer(0)
        
        for (o in seq_len(nrow(offsets))) {
          ni <- i + offsets[o, 1]
          nj <- j + offsets[o, 2]
          nk <- k + offsets[o, 3]
          
          if (ni >= 1 && ni <= n_blocks[1] &&
              nj >= 1 && nj <= n_blocks[2] &&
              nk >= 1 && nk <= n_blocks[3]) {
            
            nb_idx <- (nk - 1) * n_blocks[1] * n_blocks[2] + 
                     (nj - 1) * n_blocks[1] + ni
            nb <- c(nb, nb_idx)
          }
        }
        
        neighbors[[idx]] <- nb
      }
    }
  }
  
  neighbors
}
