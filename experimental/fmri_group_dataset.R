#' @import R6
NULL


#' R6 Class for Managing Multi-Subject fMRI Datasets
#' 
#' @description
#' The FMRIGroupDataset class provides a container and interface for managing multiple
#' fMRI datasets, typically representing different subjects in a study. It supports
#' operations like subsetting, applying functions across subjects, and accessing
#' individual subject data.
#' 
#' The class is designed to work with the BIDS (Brain Imaging Data Structure) format
#' and integrates with the bidser package for BIDS data handling.
#'
#' @details
#' # Key Features
#' * Manages multiple fMRI datasets with consistent interface
#' * Supports subject-wise operations and filtering
#' * Integrates with BIDS format and bidser package
#' * Provides method chaining for complex operations
#' 
#' # Class Structure
#' The class maintains:
#' * A list of individual fMRI datasets
#' * Subject IDs corresponding to each dataset
#' * Optional task and run information
#' * Methods for data access and manipulation
#'
#' @examples
#' \dontrun{
#' # Create individual datasets
#' ds1 <- fmri_dataset("sub-01_task-rest_bold.nii.gz", mask="mask.nii.gz", TR=2)
#' ds2 <- fmri_dataset("sub-02_task-rest_bold.nii.gz", mask="mask.nii.gz", TR=2)
#' 
#' # Create group dataset
#' group_ds <- FMRIGroupDataset$new(
#'   datasets = list(ds1, ds2),
#'   subject_ids = c("sub-01", "sub-02"),
#'   tasks = c("rest")
#' )
#' 
#' # Get specific subject's data
#' sub1_data <- group_ds$get_subject("sub-01")
#' 
#' # Apply function to all subjects
#' results <- group_ds$apply(function(ds, subj) {
#'   # Process individual dataset
#'   list(
#'     subject = subj,
#'     mean_signal = colMeans(get_data(ds))
#'   )
#' })
#' 
#' # Subset to specific subjects
#' subset_ds <- group_ds$subset(c("sub-01"))
#' }
#'
#' @seealso 
#' * [fmri_dataset()] for individual dataset creation
#' * [bidser::bids_project()] for BIDS project handling
#' 
#' @importFrom R6 R6Class
#' @export
FMRIGroupDataset <- R6Class("FMRIGroupDataset",
  public = list(
    #' @field datasets List of individual fMRI datasets
    #' Each element should be an object inheriting from class "fmri_dataset"
    datasets = NULL,
    
    #' @field subject_ids Vector of subject IDs
    #' Character vector of BIDS-compliant subject identifiers (e.g., "sub-01")
    subject_ids = NULL,
    
    #' @field tasks Vector of task names
    #' Optional character vector of task names present in the datasets
    tasks = NULL,
    
    #' @field runs Vector of run IDs
    #' Optional character vector of run identifiers
    runs = NULL,
    
    #' @description
    #' Create a new FMRIGroupDataset object
    #' 
    #' @param datasets List of individual fmri_dataset objects
    #' @param subject_ids Vector of subject IDs corresponding to datasets
    #' @param tasks Optional vector of task names
    #' @param runs Optional vector of run IDs
    #' 
    #' @return A new FMRIGroupDataset object
    #' 
    #' @examples
    #' \dontrun{
    #' # Basic creation
    #' group_ds <- FMRIGroupDataset$new(
    #'   datasets = list(ds1, ds2),
    #'   subject_ids = c("sub-01", "sub-02")
    #' )
    #' 
    #' # With task and run information
    #' group_ds <- FMRIGroupDataset$new(
    #'   datasets = list(ds1, ds2),
    #'   subject_ids = c("sub-01", "sub-02"),
    #'   tasks = c("rest", "task"),
    #'   runs = c("run-01", "run-02")
    #' )
    #' }
    initialize = function(datasets, subject_ids, tasks = NULL, runs = NULL) {
      if (length(datasets) != length(subject_ids)) {
        stop("Number of datasets must match number of subject IDs")
      }
      
      # Validate all datasets are fmri_dataset objects
      valid <- sapply(datasets, function(d) inherits(d, "fmri_dataset"))
      if (!all(valid)) {
        stop("All datasets must be fmri_dataset objects")
      }
      
      self$datasets <- datasets
      self$subject_ids <- subject_ids
      self$tasks <- tasks
      self$runs <- runs
      
      invisible(self)
    },
    
    #' @description
    #' Retrieve dataset for a specific subject
    #' 
    #' @param subject_id Character string specifying the subject ID
    #' @return The fmri_dataset object for the specified subject
    #' 
    #' @examples
    #' \dontrun{
    #' # Get single subject's dataset
    #' sub1_ds <- group_ds$get_subject("sub-01")
    #' 
    #' # Use in pipeline
    #' sub1_ds %>% 
    #'   get_data() %>%
    #'   process_data()
    #' }
    get_subject = function(subject_id) {
      idx <- match(subject_id, self$subject_ids)
      if (is.na(idx)) {
        stop("Subject ", subject_id, " not found in dataset")
      }
      self$datasets[[idx]]
    },
    
    #' @description
    #' Get all unique tasks across datasets
    #' 
    #' @return Character vector of task names
    #' 
    #' @details
    #' Returns explicitly set tasks if available, otherwise
    #' discovers unique tasks across all datasets
    get_tasks = function() {
      if (!is.null(self$tasks)) {
        return(self$tasks)
      }
      unique(unlist(lapply(self$datasets, function(d) d$tasks)))
    },
    
    #' @description
    #' Get all unique runs across datasets
    #' 
    #' @return Character vector of run IDs
    #' 
    #' @details
    #' Returns explicitly set runs if available, otherwise
    #' discovers unique runs across all datasets
    get_runs = function() {
      if (!is.null(self$runs)) {
        return(self$runs)
      }
      unique(unlist(lapply(self$datasets, function(d) d$runs)))
    },
    
    #' @description
    #' Apply a function to each dataset
    #' 
    #' @param fun Function to apply to each dataset
    #' @param ... Additional arguments passed to fun
    #' @return List of results from applying fun to each dataset
    #' 
    #' @details
    #' The function `fun` should accept at least two arguments:
    #' * dataset: The individual fmri_dataset object
    #' * subject_id: The ID of the current subject
    #' 
    #' @examples
    #' \dontrun{
    #' # Simple mean calculation
    #' means <- group_ds$apply(function(ds, subj) {
    #'   colMeans(get_data(ds))
    #' })
    #' 
    #' # More complex analysis
    #' results <- group_ds$apply(function(ds, subj, threshold) {
    #'   list(
    #'     subject = subj,
    #'     data = analyze_subject(ds, threshold = threshold)
    #'   )
    #' }, threshold = 0.05)
    #' }
    apply = function(fun, ...) {
      mapply(fun, self$datasets, self$subject_ids, ..., SIMPLIFY = FALSE)
    },
    
    #' @description
    #' Create new dataset containing only specified subjects
    #' 
    #' @param subject_ids Character vector of subject IDs to keep
    #' @return A new FMRIGroupDataset containing only the specified subjects
    #' 
    #' @examples
    #' \dontrun{
    #' # Basic subsetting
    #' sub_ds <- group_ds$subset(c("sub-01", "sub-02"))
    #' 
    #' # Chain operations
    #' results <- group_ds$
    #'   subset(c("sub-01", "sub-02"))$
    #'   apply(analyze_subject)
    #' }
    subset = function(subject_ids) {
      idx <- match(subject_ids, self$subject_ids)
      if (any(is.na(idx))) {
        stop("Some subject IDs not found: ", 
             paste(subject_ids[is.na(idx)], collapse=", "))
      }
      
      FMRIGroupDataset$new(
        datasets = self$datasets[idx],
        subject_ids = self$subject_ids[idx],
        tasks = self$tasks,
        runs = self$runs
      )
    },
    
    #' @description
    #' Print summary of group dataset
    #' 
    #' @return Invisibly returns self (for method chaining)
    print = function() {
      cat("FMRIGroupDataset\n")
      cat("================\n")
      cat("Subjects:", length(self$subject_ids), "\n")
      cat("Subject IDs:", paste(self$subject_ids, collapse=", "), "\n")
      
      tasks <- self$get_tasks()
      if (!is.null(tasks)) {
        cat("Tasks:", paste(tasks, collapse=", "), "\n")
      }
      
      runs <- self$get_runs()
      if (!is.null(runs)) {
        cat("Runs:", paste(runs, collapse=", "), "\n")
      }
      
      if (length(self$datasets) > 0) {
        cat("\nExample dataset (", self$subject_ids[1], "):\n", sep="")
        print(self$datasets[[1]])
      }
      
      invisible(self)
    }
  )
)
