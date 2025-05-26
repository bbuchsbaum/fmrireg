###############################################################################
#         CODE FOR BUILDING A HRF LIBRARY, FINDING BEST HRF, and CLUSTERING
###############################################################################

#' Setup an HRF Library
#'
#' This function builds or loads a pre-defined "library" of candidate HRFs, then
#' performs a PCA/truncation to yield a basis set. For instance, if no \code{hrflib}
#' is provided, it uses a range of parameters (h1, h2, h3, h4, f1, f2) to generate
#' half-cosine-based HRFs via \code{hrf_half_cosine}, calls \code{gen_hrf_library},
#' then does a PCA. 
#'
#' @param hrflib An optional function that, when given a vector of times, returns
#'   a matrix of candidate HRFs as columns. If \code{NULL}, a built-in set of
#'   half-cosine parameter combinations is used.
#' @param rsam A numeric vector of time points at which HRFs are sampled (default: seq(0,24,1))
#' @param ncomp Number of principal components to keep (default 5)
#'
#' @return A list with elements:
#'   \item{L}{A matrix of HRFs (columns) from the full library}
#'   \item{Lk}{A matrix of dimension \code{nrow(L) x ncomp} (the first ncomp PCs)}
#'   \item{basis_set}{A list of HRF functions (one for each PC) created by \code{gen_empirical_hrf}}
#'
#' @keywords internal
setup_hrf_library <- function(hrflib = NULL, rsam = seq(0, 24, by = 1), ncomp = 5) {
  if (is.null(hrflib)) {
    # If no library provided, create a half-cosine param grid
    params <- expand.grid(
      h1 = seq(1, 3, by = 0.33),
      h2 = seq(3, 7, by = 0.33),
      h3 = seq(6, 8, by = 1),
      h4 = seq(5, 9, by = 1),
      f1 = seq(0, 0.2, by = 0.1),
      f2 = seq(0, 0.2, by = 0.1)
    )
    hrflib <- gen_hrf_library(hrf_half_cosine, params)  # from your code base
  }
  
  # Evaluate library at 'rsam' time points
  L <- hrflib(rsam)   # L is a matrix: rows= length(rsam), cols = # of candidate HRFs
  
  # PCA the library of HRFs
  pca_L <- multivarious::pca(L, preproc = multivarious::pass())
  Lk <- pca_L$u[, 1:ncomp, drop = FALSE]
  
  # Build a small list of actual R functions from the PCA basis
  basis_set <- lapply(seq_len(ncol(Lk)), function(i) {
    gen_empirical_hrf(rsam, Lk[, i])  # from your code base
  })
  
  list(L = L, Lk = Lk, basis_set = basis_set)
}


#' Generic S3 method for `find_best_hrf`
#'
#' This function is a generic for finding best HRFs across a dataset.
#' Usually, you implement \code{find_best_hrf.matrix_dataset} or
#' \code{find_best_hrf.fmri_mem_dataset} depending on your data type.
#'
#' @param dataset The dataset
#' @param onset_var The name of the onset variable
#' @param hrflib Possibly a function returning a matrix of candidate HRFs
#' @param rsam Time points at which HRFs are sampled
#' @param ncomp Number of basis components
#' @param nsplits Number of splits
#' @param block Possibly a formula specifying block structure
#' @param basemod Possibly a baseline model
#' @param cluster_series If TRUE, cluster the assigned HRFs
#' @param cluster_quantile The quantile used in the clustering
#' @param dist_method The distance measure
#' @param spatial_smoothing If TRUE, smooth HRF assignments in space
#' @param lambda Smoothing param
#' @param n_iterations Number of smoothing iterations
#' @param k_neighbors How many neighbors to consider
#' @param nboot Number of bootstrap iterations
#' @param ... Other args
#' 
#' @export
find_best_hrf <- function(dataset, onset_var, hrflib = NULL, rsam = seq(0, 24, by = 1),
                          ncomp = 5, nsplits = 3, block = NULL, basemod = NULL, 
                          cluster_series = TRUE, 
                          cluster_quantile = .5,
                          dist_method = "euclidean",
                          spatial_smoothing = FALSE,
                          lambda = 0.1,
                          n_iterations = 5,
                          k_neighbors = 6,
                          nboot = 1,
                          ...) {
  UseMethod("find_best_hrf")
}


#' Compute Best HRFs
#'
#' This is an internal helper that tries to match each column of X_proj
#' (representing data-based PCs) to columns of L_proj (the library's PCs),
#' and optionally performs a hierarchical clustering to group them.
#'
#' @param X_proj Typically the projection of the data
#' @param L_proj The projection of the library
#' @param L The full library matrix (rows= time, cols= candidates)
#' @param cluster_series If TRUE, we cluster with cutree
#' @param cluster_quantile The distance threshold quantile
#' @param dist_method The distance measure (by default "euclidean")
#'
#' @keywords internal
compute_best_hrfs <- function(X_proj, L_proj, L,
                              cluster_series = TRUE,
                              cluster_quantile = .5,
                              dist_method = "euclidean") 
{
  # Distances between each library PC and each data PC
  basis_euc <- proxy::dist(L_proj, X_proj, method = dist_method)
  basis_euc_mat <- as.matrix(basis_euc)
  
  # For each data col, find the library col with min distance
  best_match_indices <- apply(basis_euc_mat, 2, which.min)
  
  if (!cluster_series) {
    unique_indices <- unique(best_match_indices)
    return(list(
      best_hrf_indices = best_match_indices,
      clusters = NULL,
      cluster_representatives = unique_indices,
      reassigned_hrfs = best_match_indices,
      L_best_hrfs = L[, best_match_indices, drop = FALSE],
      L_unique = L[, unique_indices, drop = FALSE],
      L_proj = L_proj,
      X_proj = X_proj
    ))
  }
  
  # If cluster_series=TRUE, we do hierarchical clustering of the unique best HRFs
  unique_best_hrfs <- unique(best_match_indices)
  L_best_hrfs <- L[, unique_best_hrfs, drop=FALSE]
  
  # cluster
  hrf_dist <- dist(t(L_best_hrfs))
  hclust_result <- hclust(hrf_dist, method="average")
  # cutree at the threshold from 'cluster_quantile'
  threshold <- quantile(hrf_dist, cluster_quantile)
  clusters <- cutree(hclust_result, h=threshold)
  
  # define a "medoid" function
  cluster_medoid <- function(cluster_indices) {
    sub_dist <- as.matrix(hrf_dist)[cluster_indices, cluster_indices, drop=FALSE]
    # pick the row with minimal sum-of-distances
    medoid_rel_idx <- which.min(rowSums(sub_dist))
    cluster_indices[medoid_rel_idx]
  }
  
  unique_clusters <- unique(clusters)
  cluster_representatives <- sapply(unique_clusters, function(cl) {
    cluster_indices <- which(clusters == cl)
    unique_best_hrfs[cluster_medoid(cluster_indices)]
  })
  
  # Create a mapping from original indices to cluster representatives
  index_to_cluster <- clusters[match(best_match_indices, unique_best_hrfs)]
  reassigned_hrfs <- cluster_representatives[index_to_cluster]
  
  # Get the actual HRF shapes for the unique representatives
  L_unique <- L[, cluster_representatives, drop=FALSE]
  
  list(
    best_hrf_indices = unique_best_hrfs,
    clusters = clusters,
    cluster_representatives = cluster_representatives,
    reassigned_hrfs = reassigned_hrfs,
    L_best_hrfs = L_best_hrfs,
    L_unique = L_unique,
    L_proj = L_proj,
    X_proj = X_proj
  )
}


#' Smooth HRF Assignments
#'
#' This function tries to smooth out a voxel-wise label map of HRF assignments
#' by iterating over each voxel, checking its neighbors, picking the "most
#' common" neighbor label if that reduces cost, etc.
#'
#' @param reassigned_hrfs An integer vector of length n_voxels
#' @param L_best_hrfs A matrix of reference HRFs for each label
#' @param voxel_coords A matrix of voxel coordinates (n_voxels x 3), or similar
#' @param lambda Not actually used, but left here as a param
#' @param n_iterations Number of smoothing passes
#' @param k Number of neighbors
#'
#' @return A vector of the same length as `reassigned_hrfs`, updated (smoothed)
#' @keywords internal
smooth_hrf_assignments <- function(reassigned_hrfs, L_best_hrfs, voxel_coords,
                                   lambda = 0.1, n_iterations = 5, k = 6) 
{
  n_voxels <- length(reassigned_hrfs)
  if (n_voxels == 0) return(reassigned_hrfs)
  
  # build neighbor graph
  dists <- as.matrix(dist(voxel_coords))
  
  # for each voxel, the top k neighbors
  neighbor_indices <- apply(dists, 1, function(row) {
    # remove self
    self_idx <- which.min(row)
    neigh <- setdiff(seq_along(row), self_idx)
    # sort ascending
    neigh <- neigh[order(row[neigh])]
    # keep top k
    head(neigh, k)
  })
  
  get_hrf_shape <- function(idx) L_best_hrfs[, idx, drop=FALSE]
  
  # do repeated label smoothing
  for (iter in seq_len(n_iterations)) {
    new_assignments <- reassigned_hrfs
    
    for (v in seq_len(n_voxels)) {
      v_hrf_idx <- reassigned_hrfs[v]
      neigh_idx <- neighbor_indices[, v]
      neigh_hrfs <- reassigned_hrfs[neigh_idx]
      
      # pick the most frequent label in neighbors
      hrf_freq <- sort(table(neigh_hrfs), decreasing = TRUE)
      candidate_hrf <- as.numeric(names(hrf_freq)[1]) # the mode
      
      if (candidate_hrf != v_hrf_idx) {
        v_shape <- get_hrf_shape(v_hrf_idx)
        neigh_shapes <- L_best_hrfs[, neigh_hrfs, drop=FALSE]
        current_cost <- mean(apply(neigh_shapes, 2, function(x) sqrt(sum((v_shape - x)^2))))
        
        cand_shape <- get_hrf_shape(candidate_hrf)
        candidate_cost <- mean(apply(neigh_shapes, 2, function(x) sqrt(sum((cand_shape - x)^2))))
        
        if (candidate_cost < current_cost) {
          new_assignments[v] <- candidate_hrf
        }
      }
    }
    
    reassigned_hrfs <- new_assignments
  }
  
  reassigned_hrfs
}


###############################################################################
#         SPECIFIC METHODS FOR matrix_dataset AND fmri_mem_dataset
###############################################################################

#' find_best_hrf for a matrix_dataset
#'
#' This method tries to find the best-fitting HRFs for each voxel by:
#'  1) building/using an HRF library with PCA
#'  2) creating splits in the event_table
#'  3) convolving a "dummy" design
#'  4) projecting data + library
#'  5) picking best HRFs, optionally clustering + smoothing
#'
#' @inheritParams find_best_hrf
#' @return A list with details about the best HRFs
#' @export
find_best_hrf.matrix_dataset <- function(dataset, onset_var, hrflib = NULL, rsam = seq(0, 24, by = 1),
                                         ncomp = 5, nsplits = 3, block = NULL, basemod = NULL, 
                                         cluster_series = TRUE, 
                                         cluster_quantile = .5,
                                         dist_method = "euclidean",
                                         spatial_smoothing = FALSE,
                                         lambda = 0.1,
                                         n_iterations = 5,
                                         k_neighbors = 6,
                                         nboot = 1,
                                         ...) 
{
  # Helper function to do 1 "pipeline" iteration
  run_pipeline <- function(ds) {
    # Step 1: Setup HRF library
    hrf_setup <- setup_hrf_library(hrflib, rsam, ncomp)
    hrf_list  <- do.call(gen_hrf_set, hrf_setup$basis_set)  # merges them into a single function
    
    # Step 2: baseline model
    if (is.null(basemod)) {
      bmod <- baseline_model("constant", sframe = ds$sampling_frame)
    } else {
      bmod <- basemod
    }
    
    # Step 3: Create an artificial "splitvar"
    split_ids <- rep(seq_len(nsplits), length.out = nrow(ds$event_table))
    fac_var <- factor(paste0("split_", split_ids))
    ds$event_table <- ds$event_table %>% mutate(.splitvar = fac_var)
    
    # Step 4: event_model
    form <- as.formula(paste0(onset_var, " ~ hrf(.splitvar, basis=hrf_list)"))
    environment(form) <- environment()
    emod <- event_model(
      x = form,
      data = ds$event_table,
      block = block,
      sampling_frame = ds$sampling_frame
    )

    
    # Step 5: data & residual
    X <- get_data_matrix(ds)
    X_base <- as.matrix(design_matrix(bmod))
    X_resid <- resid(lsfit(X_base, X, intercept = FALSE))
    
    # PCA on X_resid
    pca_X <- multivarious::pca(X_resid, preproc = multivarious::center())
    ncomp_X <- min(ncomp, ncol(pca_X$u))   # clamp to the data's rank
    Mk <- pca_X$u[, seq_len(ncomp_X), drop=FALSE]
    
    # Also clamp the library side if ncomp_X < ncomp
    if (ncomp_X < ncomp) {
      hrf_setup$Lk <- hrf_setup$Lk[, seq_len(ncomp_X), drop=FALSE]
    }
    
    # design & projection
    Rk <- design_matrix(emod)
    C <- t(Rk) %*% Mk
    svd_C <- svd(C)
    
    X_proj <- t(X_resid) %*% Mk %*% svd_C$v
    L_proj <- t(hrf_setup$L) %*% hrf_setup$Lk %*% svd_C$v
    
    # Step 6: find best HRFs
    res <- compute_best_hrfs(X_proj, L_proj, hrf_setup$L, 
                             cluster_series = cluster_series,
                             cluster_quantile = cluster_quantile,
                             dist_method = dist_method)
    
    list(res = res, L_best_hrfs  = res$L_best_hrfs, dataset = ds)
  }

  
  # No bootstrap => 1 pass
  if (nboot == 1) {
    out <- run_pipeline(dataset)
    final_assignments <- if (!is.null(out$res$reassigned_hrfs)) {
      out$res$reassigned_hrfs
    } else {
      out$res$best_hrf_indices
    }
    # optional smoothing
    if (spatial_smoothing && !is.null(dataset$voxel_coords)) {
      final_assignments <- smooth_hrf_assignments(
        final_assignments,
        out$res$L_best_hrfs,
        dataset$voxel_coords,
        lambda = lambda,
        n_iterations = n_iterations,
        k = k_neighbors
      )
      out$res$reassigned_hrfs_smoothed <- final_assignments
    }
    # Wrap into a best_hrf_result object
    return(
      new_best_hrf_result(
        best_hrf_indices = out$res$best_hrf_indices,
        clusters = out$res$clusters,
        cluster_representatives = out$res$cluster_representatives,
        reassigned_hrfs = out$res$reassigned_hrfs,
        L_best_hrfs = out$res$L_best_hrfs,
        L_unique = out$res$L_unique,
        L_proj = out$res$L_proj,
        X_proj = out$res$X_proj,
        reassigned_hrfs_smoothed = out$res$reassigned_hrfs_smoothed
      )
    )
  }
  
  # Else bootstrap
  n_voxels <- nrow(get_data_matrix(dataset))
  bootstrap_assignments <- matrix(NA, nrow=n_voxels, ncol=nboot)
  last_run <- NULL
  
  for (b in seq_len(nboot)) {
    boot_events <- dataset$event_table[sample(nrow(dataset$event_table), replace=TRUE), ]
    ds_b <- dataset
    ds_b$event_table <- boot_events
    
    run_b <- run_pipeline(ds_b)
    final_assignments_b <- if (!is.null(run_b$res$reassigned_hrfs)) {
      run_b$res$reassigned_hrfs
    } else {
      run_b$res$best_hrf_indices
    }
    bootstrap_assignments[, b] <- final_assignments_b
    last_run <- run_b
  }
  
  # majority vote
  final_assignments <- apply(bootstrap_assignments, 1, function(x) {
    tt <- sort(table(x), decreasing=TRUE)
    as.numeric(names(tt)[1])
  })
  
  if (spatial_smoothing && !is.null(dataset$voxel_coords)) {
    final_assignments <- smooth_hrf_assignments(
      final_assignments,
      last_run$res$L_best_hrfs,
      dataset$voxel_coords,
      lambda=lambda,
      n_iterations=n_iterations,
      k=k_neighbors
    )
  }
  
  final_res <- last_run$res
  final_res$reassigned_hrfs <- final_assignments
  if (spatial_smoothing) {
    final_res$reassigned_hrfs_smoothed <- final_assignments
  }
  
  # Return S3
  new_best_hrf_result(
    best_hrf_indices = final_res$best_hrf_indices,
    clusters = final_res$clusters,
    cluster_representatives = final_res$cluster_representatives,
    reassigned_hrfs = final_res$reassigned_hrfs,
    L_best_hrfs = final_res$L_best_hrfs,
    L_unique = final_res$L_unique,
    L_proj = final_res$L_proj,
    X_proj = final_res$X_proj,
    reassigned_hrfs_smoothed = final_res$reassigned_hrfs_smoothed
  )
}

#' find_best_hrf for fmri_mem_dataset
#'
#' A simple pass-through that calls \code{find_best_hrf.matrix_dataset}.
#' 
#' @inheritParams find_best_hrf
#' @export
find_best_hrf.fmri_mem_dataset <- function(dataset, onset_var, hrflib = NULL, rsam = seq(0, 24, by = 1),
                                           ncomp = 5, nsplits = 3, block = NULL, basemod = NULL, 
                                           cluster_series = TRUE, 
                                           cluster_quantile = .5,
                                           dist_method = "euclidean",
                                           spatial_smoothing = FALSE,
                                           lambda = 0.1,
                                           n_iterations = 5,
                                           k_neighbors = 6,
                                           nboot = 1,
                                           ...) {
  find_best_hrf.matrix_dataset(
    dataset         = dataset,
    onset_var       = onset_var,
    hrflib          = hrflib,
    rsam            = rsam,
    ncomp           = ncomp,
    nsplits         = nsplits,
    block           = block,
    basemod         = basemod,
    cluster_series  = cluster_series,
    cluster_quantile= cluster_quantile,
    dist_method     = dist_method,
    spatial_smoothing = spatial_smoothing,
    lambda          = lambda,
    n_iterations    = n_iterations,
    k_neighbors     = k_neighbors,
    nboot           = nboot,
    ...
  )
}


#' Create a new "best_hrf_result" object
#'
#' This constructor wraps up the result from `find_best_hrf` 
#' into a proper S3 class with a pretty print method.
#'
#' @param best_hrf_indices numeric vector of unique best HRFs
#' @param clusters integer vector of clusters assigned to each unique HRF
#' @param cluster_representatives integer vector of cluster medoids
#' @param reassigned_hrfs integer vector of HRF IDs assigned to each voxel
#' @param L_best_hrfs matrix of all best HRFs found
#' @param L_unique matrix of unique HRFs (one per cluster/representative)
#' @param L_proj matrix of the library projected into data space
#' @param X_proj matrix of the data projection
#' @param reassigned_hrfs_smoothed optionally the smoothed HRF assignments
#'
#' @return An object of class \code{"best_hrf_result"}.
#' @keywords internal
new_best_hrf_result <- function(
    best_hrf_indices = NULL,
    clusters = NULL,
    cluster_representatives = NULL,
    reassigned_hrfs = NULL,
    L_best_hrfs = NULL,
    L_unique = NULL,
    L_proj = NULL,
    X_proj = NULL,
    reassigned_hrfs_smoothed = NULL
) {
  structure(
    list(
      best_hrf_indices = best_hrf_indices,
      clusters = clusters,
      cluster_representatives = cluster_representatives,
      reassigned_hrfs = reassigned_hrfs,
      L_best_hrfs = L_best_hrfs,
      L_unique = L_unique,
      L_proj = L_proj,
      X_proj = X_proj,
      reassigned_hrfs_smoothed = reassigned_hrfs_smoothed
    ),
    class = "best_hrf_result"
  )
}


#' Print method for "best_hrf_result"
#'
#' Uses crayon for colorful output
#'
#' @param x A \code{best_hrf_result} object
#' @param ... Not used
#' @export
print.best_hrf_result <- function(x, ...) {
  cat(crayon::blue("\n═══ Best HRF Result ═══\n\n"))
  
  # Basic info
  if (!is.null(x$reassigned_hrfs)) {
    cat(crayon::green(" Number of voxels assigned:"), length(x$reassigned_hrfs), "\n")
  }
  
  if (!is.null(x$clusters)) {
    cat(crayon::green(" Number of clusters:"), length(unique(x$clusters)), "\n")
  }
  
  if (!is.null(x$best_hrf_indices)) {
    cat(crayon::green(" Number of unique best HRFs:"), length(x$best_hrf_indices), "\n")
  }
  
  if (!is.null(x$L_unique)) {
    cat(crayon::green(" Number of representative HRFs:"), ncol(x$L_unique), "\n")
  }
  
  # If there's a smoothed assignment
  if (!is.null(x$reassigned_hrfs_smoothed)) {
    cat(crayon::yellow(" Smoothing was applied. See `$reassigned_hrfs_smoothed`.\n"))
  }
  
  cat(crayon::silver("\n(Use `$` to inspect fields: best_hrf_indices, L_best_hrfs, L_unique, etc.)\n"))
}
