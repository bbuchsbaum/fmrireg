#' Generate an AFNI Linear Model from an fmri_config object
#'
#' This function takes an fmri_config object and generates an AFNI linear model
#' by reading the necessary files, creating an fmri_dataset, event_model,
#' baseline_model, and fmri_model, and fitting the model using afni_lm.
#'
#' @param x An fmri_config object containing the necessary configuration for the analysis
#' @param ... Additional arguments passed to the function
#' @return An afni_lm object representing the fitted linear model
#' @importFrom utils read.table
#' @keywords internal
#' @noRd
gen_afni_lm.fmri_config <- function(x, ...) {
  
  nuisance_list <- if (!is.null(x$baseline_model$nuisance_files)) {
    lapply(x$baseline_model$nuisance_files, read.table, header=TRUE)
  }
  
  dset <- fmri_dataset(scans=x$scans, 
                       mask=x$mask, 
                       TR=x$TR, 
                       run_length=x$run_length, 
                       event_table=x$design, 
                       base_path=x$base_path,
                       censor=if (is.null(x$censor_file)) NULL else scan(paste0(x$base_path, "/", x$censor_file)))

  emodel <- event_model(x$event_model, data=x$design, block=as.formula(paste("~", x$block_column)),
                        sampling_frame=dset$sampling_frame)
  
  bmodel <- baseline_model(basis=x$baseline_model$basis, 
                           degree=x$baseline_model$degree, 
                           sframe=dset$sampling_frame, 
                           nuisance_list=nuisance_list)
  
  fmodel <- fmri_model(emodel, bmodel, dset)
  alm <- afni_lm(fmodel, dset, censor=dset$censor)
  
}

  

#' Set up an fMRI linear model for AFNI's 3dDeconvolve
#'
#' This function prepares an fMRI linear model for AFNI's 3dDeconvolve tool.
#' It takes an fmri_model object, an fmri_dataset object, and various options
#' to control the fitting process.
#'
#' @param fmri_mod An fmri_model object containing the event and baseline models
#' @param dataset An fmri_dataset object containing the scan data and other necessary information
#' @param working_dir The working directory (default is the current directory)
#' @param polort The number of polynomial baseline regressors (default is to suppress 'polort')
#' @param jobs The number of jobs to use with '3dDeconvolve' (default is 1)
#' @param censor A list of censoring vectors, one per run, or a single vector equal to the total number of scans (default is NULL)
#' @param options A list of options to be sent to 3dDeconvolve (default is an empty list)
#'
#' @return An afni_lm_spec object containing the fitted model, dataset, working directory, options, and command
#'
#' @examples 
#' etab <- data.frame(onset=c(1,30,15,25), fac=factor(c("A", "B", "A", "B")), 
#' run=c(1,1,2,2))
#' dset <- fmri_dataset(scans=c("s1.nii", "s2.nii"), mask="mask.nii", TR=1, 
#' run_length=c(50,50), event_table=etab)
#'
#' emodel <- event_model(onset ~ hrf(fac), block = ~ run, data=etab, 
#' sampling_frame=dset$sampling_frame)
#' bmodel <- baseline_model("bs", degree=4, sframe=dset$sampling_frame)
#' fmod <- fmri_model(emodel, bmodel)
#' alm <- afni_lm(fmod, dset, jobs=2, options=list(tout=TRUE, errts="residuals.nii.gz"))
#' @export
afni_lm <- function(fmri_mod, dataset, working_dir=".", polort=-1, jobs=1, 
                    censor=NULL, options=list()) {
  
  if (!is.null(censor)) {
    censor <- unlist(censor)
    assert_that(length(censor) == sum(fmrihrf::blocklens(dataset$sampling_frame)))
  }
  
  defopts <- list(noFDR=FALSE, 
                       fout=FALSE, 
                       rout=FALSE, 
                       tout=TRUE, 
                       float=TRUE, 
                       nocond=FALSE,
                       x1D_stop=FALSE,
                       cbucket="coefout", 
                       bucket="statout", 
                       jobs=jobs, 
                       polort=polort, 
                       censor=censor,
                       iresp=FALSE, 
                       noFDR=TRUE,
                       nofullf_atall=TRUE,
                       TR_times=1)
  
  for (optname in names(options)) {
    defopts[[optname]] <- options[[optname]]
  }
  
  cmd <- build_decon_command(model=fmri_mod, dataset,working_dir=working_dir,opts=defopts)
  
  structure(
    list(
      model=fmri_mod,
      dataset=dataset,
      working_dir=working_dir,
      options=defopts,
      cmd=cmd),
    class="afni_lm_spec")

}


#' @export
#' @rdname print
print.afni_lm_spec <- function(x,...) {
  cat("AFNI linear model via 3dDeconvolve \n")
  cat("  working_dir: ", x$working_dir, "\n")
  cat("  number of GLTs:", length(x$cmd$glts), "\n")
  cat("  command line: \n", x$cmd$cmd, "\n")
  
}


#' @keywords internal
#' @noRd
afni_stim_file <- function(label, file_name, values) {
  structure(
    list(label=label, file_name=file_name, values=values),
    class=c("afni_stim_file", "afni_stim")
  )
}

#' @keywords internal
#' @noRd
afni_stim_times <- function(label, file_name, hrf, onsets, blockids, iresp=FALSE, sresp=FALSE, tr_times=1) {
  structure(
    list(label=label, file_name=file_name, hrf=hrf, onsets=onsets, blockids, iresp=iresp, sresp=sresp, tr_times=tr_times),
    class=c("afni_stim_times","afni_stim")
  )
}


#' @keywords internal
#' @noRd
afni_stim_im_times <- function(label, file_name, hrf, onsets, blockids) {
  structure(
    list(label=label, file_name=file_name, hrf=hrf, onsets=onsets, blockids=blockids),
    class=c("afni_stim_im_times", "afni_stim_times","afni_stim")
  )
}


#' @keywords internal
#' @noRd
afni_command_switch <- function(x, ...) UseMethod("afni_command_switch")

#' @keywords internal
#' @noRd
write_afni_stim <- function(x, ...) UseMethod("write_afni_stim")


#' @keywords internal
#' @noRd
afni_command_switch.afni_stim_file <- function(x, k, type) {
  switch(type,
    label = paste(k, x$label, collapse = " "),
    file = paste(k, x$file_name, collapse = " "),
    ortvec = paste(x$file_name, x$label),
    times = NULL,
    iresp = NULL,
    sresp=NULL
  )
}


#' @keywords internal
#' @noRd
afni_command_switch.afni_stim_im_times <- function(x, k, type) {
  switch(type,
         label = paste(k, x$label, collapse = " "),
         times_IM = paste(k, x$file_name, as.character(x$hrf), collapse = " "),
         NULL
  )
}

#' @keywords internal
#' @noRd
afni_command_switch.afni_stim_times <- function(x, k, type) {
  switch(type,
         label = paste(k, x$label, collapse = " "),
         times=paste(k, x$file_name, as.character(x$hrf), collapse=" "),
         file=NULL,
         ortvec=NULL,
         iresp = if (x$iresp) paste(k, paste(paste(x$label, "_iresp", sep=""), collapse=" ")) else NULL,
         sresp = if (x$sresp) paste(k, paste(paste(x$label, "_sresp", sep=""), collapse=" ")) else NULL
  )
}

#' @keywords internal
#' @noRd
next_dir_name <- function(wd) {
  nd <- paste(wd, "+", sep="")
  if (!file.exists(nd)) {
    nd
  } else {
    Recall(nd)
  }
}


#' @noRd
#' @keywords internal
#' @importFrom utils write.table
write_baseline_mat <- function(stim, dir) {
  write.table(stim$mat, paste0(dir, "/", stim$file_name, sep=""), col.names=FALSE, row.names=FALSE)
}

#' @keywords internal
#' @noRd
write_baseline_mats <- function(blist) {
  purrr::walk(blist, ~ write_baseline_mat(., "."))
}
  
#' @keywords internal
#' @noRd
write_stim_files <- function(afni_stims) {
  purrr::walk(afni_stims, ~ write_afni_stim(., "."))
}

#' @keywords internal
#' @noRd
write_afni_stim.afni_stim_file <- function(stim, dir) {
  .write_values <- function(outname, vals) {
    hfile <- file(outname, "w")
    write(vals, file=hfile, ncolumns=1)
    close(hfile)
  }
  
  ## TODO stim$values is a data.frame sometimes (trialwise?), hence 'unlist' hack. Ensure uniformity.
  .write_values(paste0(dir, "/", stim$file_name, sep=""), unlist(stim$values))
}

#' @keywords internal
#' @noRd
write_afni_stim.afni_stim_times <- function(stim, dir) {
  ## TODO handle runs with no stimuli
  .write_onsets <- function(outname, vals) {
    hfile <- file(outname, "w")
    for (i in 1:length(vals)) {
      ons <- paste(stim$onsets[[i]], collapse=" ")
      writeLines(ons, con=hfile)
    }
    
    close(hfile)
  }
  
  .write_onsets(paste0(dir, "/", stim$file_name, sep=""), stim$onsets)
}



#' @keywords internal
#' @noRd
write_censor_file <- function(dir, censor) {
    outname <- paste0(dir, "/censor.1D")
    hfile <- file(outname, "w")
    write(censor, file=hfile, ncolumns=1)
    close(hfile)
}
  

#' @keywords internal
#' @noRd
#' @importFrom purrr imap
write_glts <- function(glts, gltfiles) {
  imap(glts, function(glt, i) {
    fout <- file(gltfiles[i], "w")
    write(glt, file=fout, sep="\n")
    close(fout)
  })
}

#' @keywords internal
#' @noRd
afni_baseline_matrix <- function(label, file_name, mat) {
  structure(
    list(label=label, file_name=file_name, mat=mat),
    class="afni_baseline_matrix"
  )
}

#' @keywords internal
#' @noRd
build_baseline_stims <- function(x) {
  blens <- fmrihrf::blocklens(x)
  nblocks <- length(blens)
  
  bterms <- terms(x$baseline_model)
  
  ret <- lapply(bterms, function(bt) {
    lapply(1:nblocks, function(i) {
      mat <- design_matrix(bt, i, allrows=TRUE)
      afni_baseline_matrix(paste0(bt$varname, "_", i), paste0(bt$varname, "_", i, ".1D"), mat)
    })
  })
  
  unlist(ret, recursive=FALSE)
}

#' Build AFNI stimulus definitions
#'
#' Generate stimulus files for AFNI's 3dDeconvolve.
#' @param x Term object
#' @param ... Additional arguments
#' @param iresp Logical; generate -iresp files when TRUE
#' @param tr_times Numeric; TR timing for stimuli
#' @return List of stimulus descriptors or `NULL`
#' @export
#' @family AFNI
#' @rdname build_afni_stims
build_afni_stims <- function(x, ...) UseMethod("build_afni_stims")

#' @keywords internal
#' @noRd
build_afni_stims.convolved_term <- function(x, iresp=FALSE, tr_times=1) {
  stimlabels <- longnames(x)
  stimfiles <- paste(stimlabels, "_reg.1D", sep = "")
  desmat <- design_matrix(x)
  
  lapply(1:length(stimlabels), function(i) {
    afni_stim_file(stimlabels[i], stimfiles[i], desmat[, i])
  })
}

#' @keywords internal
#' @noRd
build_afni_stims.afni_hrf_convolved_term <- function(x, iresp=FALSE, tr_times=1) {

  stimlabels <- longnames(x)
  stimfiles <- paste(stimlabels, "_times.1D", sep = "")
  dmat <- design_matrix(x$evterm)
  
  hrf_name <- as.character(x$hrfspec$hrf)
  blids <- unique(fmrihrf::blockids(x$evterm))
  split_ons <- split_onsets(x$evterm, x$sampling_frame, global=FALSE, blocksplit = TRUE)
  names(split_ons) <- stimlabels
  
  ret <- lapply(1:length(stimlabels), function(i) {
    afni_stim_times(stimlabels[i], stimfiles[i], hrf_name, 
                    split_ons[[stimlabels[[i]]]], blockids=blids, 
                    iresp=iresp, tr_times=tr_times)
  })
  
}

#' @keywords internal
#' @noRd
build_afni_stims.afni_trialwise_convolved_term <- function(x, iresp=FALSE, tr_times=1) {
  eterm <- x$evterm
  sf <- x$sampling_frame
  hrf_name <- as.character(x$hrfspec$hrf) # Get HRF name

  # Extract the trial variable name specified in afni_trialwise(TRIAL_VAR, ...)
  trial_var_name <- eterm$varname # This should be 'trial' or whatever was passed

  # Get the full event table associated with this term
  full_event_table <- event_table(eterm)

  # Ensure the trial variable and any modulators exist in the event table
  if (!(trial_var_name %in% names(full_event_table))) {
    stop(paste("Trial variable '", trial_var_name, "' not found in event table for term '", eterm$name, "'."))
  }
  
  # Get modulators if any. This part might need adjustment based on how modulators are stored.
  # For afni_trialwise, the main 'modulator' is the trial variable itself.
  # Additional parametric modulators would be columns in the event_table.
  # The -stim_times_IM format is onset*modulator1:modulator2 etc.
  # For basic afni_trialwise("trial"), the modulator is just the trial number (or a unique trial ID).
  # If the user supplied `afni_trialwise("trial", modulate_by="RT")`, then RT would be a modulator.
  # For now, let's assume the 'trial' column itself is the modulator if no other is specified.

  # Check if `x$hrfspec` contains explicit modulator names
  modulator_vars <- x$hrfspec$modulate_by
  if (is.null(modulator_vars)) {
      # If no explicit modulators, use the trial_var_name as the modulator for amplitude (common for -IM)
      # However, AFNI's -stim_times_IM expects onsets and *optional* modulators.
      # If only trial-wise HRFs are desired without amplitude modulation per trial,
      # we still need a 'modulator' for AFNI, often it's implicitly 1 for each trial.
      # For now, we will format it as onset*trial_value
      modulator_vars <- trial_var_name 
  }
  
  if (!all(modulator_vars %in% names(full_event_table))) {
      stop(paste("One or more modulator variables not found in event table:", paste(modulator_vars, collapse=", ")))
  }

  stim_label <- longnames(eterm) # Should be a unique label for this IM regressor
  stim_file_name <- paste0(stim_label, "_IM_times.1D")

  # Prepare onsets and modulators per run, in AFNI's asterisk format
  # e.g., run 1: "10*mod1a:mod2a 20*mod1b:mod2b ...", run 2: "5*mod1c:mod2c ..."
  onsets_IM_per_run <- lapply(unique(sf$blockids), function(run_id) {
    run_events <- full_event_table[full_event_table[[sf$block_name]] == run_id, , drop = FALSE]
    if (nrow(run_events) == 0) {
      return("*") # AFNI placeholder for empty run
    }
    
    sapply(1:nrow(run_events), function(i) {
      onset_val <- run_events$onset[i]
      mod_values <- sapply(modulator_vars, function(mvar) run_events[[mvar]][i])
      paste0(onset_val, "*", paste(mod_values, collapse=":"))
    }) |> paste(collapse = " ")
  })

  # Create the afni_stim_im_times object
  afni_stim_obj <- afni_stim_im_times(
    label = stim_label,
    file_name = stim_file_name,
    hrf = hrf_name, # HRF model string like 'SPMG1'
    onsets = onsets_IM_per_run, # This should be a list of character strings, one per run
    blockids = unique(sf$blockids) # The unique run/block identifiers
  )
  
  # build_afni_stims can return a single stim object or a list of them.
  # For afni_trialwise, it's typically one -stim_times_IM directive.
  return(list(afni_stim_obj)) # Ensure it's returned in a list as expected by build_decon_command processing
}

#' @export
#' @rdname build_afni_stims
build_afni_stims.event_term <- function(x, iresp=FALSE, tr_times=1,...) {
  # This method handles plain event_terms created from standard hrf() calls.
  # These are not directly translated into AFNI -stim_times or similar directives.
  # They might be relevant if we were generating -stim_regressor from an existing matrix,
  # but that's not the purpose of this function family.
  # We return NULL to indicate no AFNI-specific stimulus file for this term.
  # The calling function (build_decon_command) should filter out NULLs.
  
  # Optionally, issue a message if verbose debugging is needed.
  # message(paste("Skipping AFNI stimulus generation for standard event term:", x$name,"(class:", paste(class(x), collapse=", "), ") whose hrfspec is of class:", paste(class(hrfspec(x)), collapse=", ")))
  return(NULL)
}

# Helper function to extract HRF spec, safely returning NULL if not found
#' @keywords internal
if (!exists("hrfspec", mode = "function")) {
  hrfspec <- function(x) attr(x, "hrfspec")
}

#' @keywords internal
#' @noRd
.make_decon_command_str <- function(cmdlines) {
  cmdstr <- lapply(names(cmdlines), function(optname) {
    entry <- cmdlines[[optname]]
    if (is.list(entry)) {
      switchnames <- rep(paste("-", optname, sep=""), length(entry))
      paste(switchnames, entry, collapse=" ")
    } else if (is.null(entry) || length(entry) == 0) {
      ""
    } else if (is.numeric(entry[[1]])) {
      paste(paste("-", optname, sep=""), entry[[1]])
    } else if (entry[[1]] == TRUE) {
      paste("-", optname, sep="")
    } else if (entry[[1]] == FALSE) {
      paste("")
    } else {
      paste(paste("-", optname, sep=""), paste(entry, collapse=" "))
    }
  })

  cmdstr <- Filter(function(x) !is.null(x) & x != "", cmdstr)
  cmdstr <- paste(cmdstr, collapse=" ")
  cmdstr <- paste("3dDeconvolve", cmdstr)
  cmdstr
}



#' Build 3dDeconvolve command for an fMRI model
#'
#' This function constructs the command string and associated options required
#' to run 3dDeconvolve using the specified fMRI model, dataset, working directory,
#' and other options. This command string can then be used to perform the actual
#' fMRI analysis using the AFNI software.
#'
#' @param model The fMRI model, usually created using the fmri_model function
#' @param dataset The fMRI dataset, usually created using the fmri_dataset function
#' @param working_dir The working directory
#' @param opts A list of options for the 3dDeconvolve command
#'
#' @return A list containing:
#'         - cmd: The 3dDeconvolve command string
#'         - cmdlines: The command lines for the 3dDeconvolve command
#'         - afni_stims: A list of AFNI stimulus objects
#'         - afni_baseline_mats: A list of AFNI baseline matrices
#'         - gltfiles: A list of GLT (general linear test) filenames
#'         - gltnames: A list of GLT names
#'         - glts: A list of GLT objects
#'         - gltstr: A list of GLT strings
#'         - censor: The censoring vector
#'
#' @keywords internal
#' @noRd
build_decon_command <- function(model, dataset, working_dir, opts) {
  func_terms <- terms(model$event_model)
  message("number of functional terms: ", length(func_terms))
  
  # First, generate AFNI stims and filter out NULLs
  afni_stims <- lapply(func_terms, function(term) { build_afni_stims(term, iresp=opts[["iresp"]], tr_times=opts[["TR_times"]]) })
  afni_stims <- Filter(Negate(is.null), afni_stims) # Filter out NULLs
  afni_stims <- unlist(afni_stims, recursive = FALSE) # Unlist one level
  
  # Now get stimlabels only for terms that generated AFNI stims
  # We need to identify which terms generated stims
  terms_with_stims <- func_terms[!sapply(lapply(func_terms, function(term) { build_afni_stims(term, iresp=opts[["iresp"]], tr_times=opts[["TR_times"]]) }), is.null)]
  stimlabels <- unlist(lapply(terms_with_stims, longnames))
  
  ## all stims must be unique
  assert_that(length(unique(stimlabels)) == length(stimlabels))

  # Note: We can't assert stimlabels == conditions because conditions includes ALL terms,
  # but stimlabels only includes terms that generate AFNI stims
  # assert_that(length(stimlabels) == length(conditions(model$event_model)))
  
  ## extract all contrast matrices
  cons <- contrast_weights(model)
  ## ensure contrasts are valid objects
  cons <- Filter(function(x) inherits(x, "contrast"), cons)
  assert_that(is.list(cons), all(sapply(cons, inherits, "contrast")))
  
  ## convert to 'glt's
  glts <- lapply(cons, to_glt)
  
  gltfiles <- unlist(lapply(glts, function(x) paste0(x$name, ".txt")))
  gltnames <- unlist(lapply(glts, function(x) x$name))
  gltstr <- unlist(lapply(glts, function(x) x$glt_str))
  
  assert_that(sum(duplicated(gltnames))  == 0, msg="Cannot have two GLTs with the same name")
  
  afni_baseline_mats <- build_baseline_stims(model)
  
  purge_nulls <- function(A) {
    A[!sapply(A, is.null)]
  }
  
  opt_stim_labels <-  purge_nulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "label")))
  opt_stim_files  <-  purge_nulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "file")))
  opt_stim_times  <-  purge_nulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "times")))
  opt_stim_times_IM  <-  purge_nulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "times_IM")))
  opt_stim_ortvecs <- purge_nulls(lapply(seq_along(afni_baseline_mats), function(i) afni_command_switch.afni_stim_file(afni_baseline_mats[[i]], i, "ortvec")))
  opt_stim_iresp  <-  purge_nulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "iresp")))
  
  #browser()
  
  if ( (length(opt_stim_times) + length(opt_stim_times_IM)) > 0) {
    ## if we use `afni_hrf` that use -stim_times, then we use local times
    global_times <- FALSE
  } else {
    ## otherwise global_times is irrelevant, since values rather than times are provided.
    global_times <- TRUE
  }
  
  cmdlines <- list(input=paste0(dataset$scans),
                   mask=paste0(dataset$mask_file),
                   polort=if (opts[["polort"]] > 0) opts[["polort"]] else -1,
                   global_times=if (global_times) TRUE else NULL,
                   num_stimts=length(afni_stims),
                   num_glt=length(gltfiles),
                   stim_file=opt_stim_files,
                   stim_label=opt_stim_labels,
                   ortvec=opt_stim_ortvecs,
                   censor=if (!is.null(opts[["censor"]])) "censor.1D" else NULL,
                   stim_times=opt_stim_times,
                   stim_times_IM=opt_stim_times_IM,
                   TR_times=opts[["TR_times"]],
                   iresp=opt_stim_iresp,
                   gltsym=lapply(seq_along(gltfiles), function(i) paste(gltfiles[i], collapse=" ")),
                   glt_label=lapply(seq_along(gltnames), function(i) paste(i, gltnames[i], collapse=" ")),
                   nofullf_atall=opts[["nofullf_atall"]],
                   fout=opts[["fout"]],
                   rout=opts[["rout"]],
                   tout=opts[["tout"]],
                   bout=opts[["bout"]],
                   noFDR=opts[["noFDR"]],
                   cbucket=opts[["cbucket"]],
                   bucket=opts[["bucket"]],
                   nocond=opts[["nocond"]],
                   x1D_stop=opts[["x1D_stop"]],
                   jobs=opts[["jobs"]],
                   errts=if (!is.null(opts[["errts"]])) opts[["errts"]] else NULL,
                   float=TRUE)
  
  cmd <- .make_decon_command_str(cmdlines)
  
  list(cmd=cmd, cmdlines=cmdlines, afni_stims=afni_stims, afni_baseline_mats=afni_baseline_mats,
       gltfiles=gltfiles, gltnames=gltnames, glts=glts, gltstr=gltstr, censor=opts$censor)
}



#' Run an afni_lm_spec object
#'
#' This function runs the 3dDeconvolve command for the specified afni_lm_spec object.
#' It outputs the results to a specified directory and can either execute the command
#' or only output the shell '3dDeconvolve.sh' script.
#'
#' @param x An afni_lm_spec object containing the model, dataset, and options
#' @param outdir The output folder
#' @param execute Whether to execute the command or only output the shell '3dDeconvolve.sh' script (default is TRUE)
#' @param execfun Function used to execute external system command (default is system)
#' @param prepend Prepend string to command (default is an empty string)
#' @param ... Additional arguments passed to execfun
#'
#' @return NULL. The function is used for its side effects, such as writing output files.
#'
#' @examples
#' # Assuming you have created an afni_lm_spec object called alm
#' #run.afni_lm_spec(alm, outdir="results")
#' @export
run.afni_lm_spec <- function(x, outdir, execute=TRUE, execfun=system, prepend="",...) {
  start_dir <- getwd()
  res <- try({
    if (!file.exists(outdir)) {
      dir.create(outdir)
    } else {
      warning(paste("glm output directory: ", outdir, " already exists"))
      outdir <- next_dir_name(outdir)
      dir.create(outdir)
      warning(paste("outputting to: ", outdir))
    }
    print(paste("setting directory:", outdir))
    setwd(outdir)
    
    write_stim_files(x$cmd$afni_stims)
    
    if (!is.null(x$cmd$gltstr)) {
      write_glts(x$cmd$gltstr, x$cmd$gltfiles)
    }
    
    if (!is.null(x$cmd$afni_baseline_mats)) {
      write_baseline_mats(x$cmd$afni_baseline_mats)
    }
    
    if (!is.null(x$cmd$censor)) {
      write_censor_file(".", x$cmd$censor)
    }
    
    #if (reml) {
    #  x$cmd$cmd <- paste(x$cmd$cmd, "-x1D_stop")
    #}
    
    write(x$cmd$cmd, "3ddeconvolve.sh")
    
    if (execute) {
      execfun(paste(prepend, x$cmd$cmd))
      
      #if (reml) {
      #  execfun(paste0("./", x$options$bucket, ".REML_cmd"))
      #}
    }
  })
  
  setwd(start_dir)
}

#' convert a contrast to an AFNI 'GLT' 
#' 
#' @param x the contrast to convert
#' @param ... extra args
#' @export
to_glt <- function(x, ...) UseMethod("to_glt")


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


#' @importFrom assertthat assert_that
#' @importFrom rlang enexpr
NULL

#' AFNI HRF Constructor Function
#'
#' @description
#' The `AFNI_HRF` function creates an object representing an AFNI-specific hemodynamic response function (HRF). It is a class constructor for AFNI HRFs.
#'
#' @param name A string specifying the name of the AFNI HRF.
#' @param nbasis An integer representing the number of basis functions for the AFNI HRF.
#' @param params A list containing the parameter values for the AFNI HRF.
#' @param span A numeric value representing the span in seconds of the HRF. Default is 24.
#'
#' @return An AFNI_HRF object with the specified properties.
#'
#' @seealso HRF
#'
#' @export
#' @rdname AFNI_HRF-class
AFNI_HRF <- function(name, nbasis, params, span = 24) {
  structure(name,
            nbasis = as.integer(nbasis),
            params = params,
            span = span,
            class = c("AFNI_HRF", "HRF"))
  
}


#' @export
as.character.AFNI_HRF <- function(x,...) {
  paste(x, "\\(", paste(attr(x, "params"), collapse=","), "\\)", sep="")
}

#' construct an native AFNI hrf specification for '3dDeconvolve' with the 'stim_times' argument.
#' 
#' @param ... Variables to include in the HRF specification
#' @param basis Character string specifying the basis function type
#' @param onsets Numeric vector of event onset times
#' @param durations Numeric vector of event durations
#' @param prefix Character string prefix for the term
#' @param subset Expression for subsetting events
#' @param nbasis Number of basis functions
#' @param contrasts Contrast specifications
#' @param id Character string identifier for the term
#' @param lag Numeric lag in seconds
#' @param precision Numeric precision for convolution
#' @param summate Logical whether to summate overlapping responses
#' @param start the start of the window for sin/poly/csplin models
#' @param stop the stop time for sin/poly/csplin models
#' @return an \code{afni_hrfspec} instance with class "afni_hrfspec"
#' @examples
#' # Create SPM canonical HRF specification
#' hrf1 <- afni_hrf(onsets = c(10, 30, 50), basis = "spmg1")
#' 
#' # Create block HRF with duration
#' hrf2 <- afni_hrf(onsets = c(10, 30, 50), durations = 5, basis = "block")
#' 
#' # Create tent basis HRF
#' hrf3 <- afni_hrf(onsets = c(10, 30, 50), basis = "tent", 
#'                  start = 0, stop = 20, nbasis = 10)
#' @export
afni_hrf <- function(..., basis=c("spmg1", "block", "dmblock",
                                  "tent",   "csplin", "poly",  "sin", "sine",
                                  "gam", "gamma", "spmg2", "spmg3", "wav"),
                                  onsets=NULL, durations=NULL, prefix=NULL, subset=NULL,
                                  nbasis=1, contrasts=NULL, id=NULL,
                                  lag=0, precision = 0.3, summate = TRUE,
                                  start=NULL, stop=NULL) {
  
  ## TODO cryptic error message when argument is mispelled and is then added to ...
  basis <- tolower(basis)
  if (basis == "sin") basis <- "sine"
  if (basis == "gam") basis <- "gamma"
  basis <- match.arg(basis)
  
  vars <- as.list(base::substitute(list(...)))[-1]
  # Convert raw expressions to quosures for compatibility with construct_event_term
  vars_quos <- lapply(vars, function(expr) rlang::new_quosure(expr, env = rlang::caller_env()))
  # Requires parse_term (assuming it exists elsewhere now or is removed)
  # parsed <- parse_term(vars, "afni_hrf") 
  # term <- parsed$term
  # label <- parsed$label
  # --- Need to replace parse_term dependency --- 
  # Simplified naming based on input expressions/symbols for now
  var_labels <- sapply(vars, rlang::as_label)
  varnames <- sapply(vars, function(v) {
       if (rlang::is_symbol(v)) as.character(v) else make.names(rlang::as_label(v))
  })
  term <- vars_quos # Store quosures instead of raw expressions
  label <- paste0("afni_hrf(", paste0(var_labels, collapse=","), ")")
  # --- End replacement --- 
  
  hrf <- if (!is.null(durations)) {
    assert_that(length(durations) == 1, msg="afni_hrf does not currently accept variable durations")
    get_AFNI_HRF(basis, nbasis=nbasis, duration=durations[1], b=start, c=stop)
  } else {
    get_AFNI_HRF(basis, nbasis=nbasis, b=start, c=stop)
  }
  
  
  # varnames <- if (!is.null(prefix)) {
  #   paste0(prefix, "_", term)
  # } else {
  #   term
  # } 
  # Use the new varnames logic
  if (!is.null(prefix)) {
      varnames <- paste0(prefix, "_", varnames)
  }
  
  termname <- paste0(varnames, collapse="::")
  
  if (is.null(id)) {
    id <- termname
  }  
  
  # contrast_set function not yet implemented
  cset <- NULL
  # cset <- if (inherits(contrasts, "contrast_spec")) {
  #   contrast_set(con1=contrasts)
  # } else if (inherits(contrasts, "contrast_set")) {
  #   contrasts
  # } else { NULL } # Default to NULL if not correct type
  
  ret <- list(
    name=termname,
    id=id,
    varnames=varnames,
    vars=term, # Store original expressions/symbols
    label=label,
    hrf=hrf,
    onsets=onsets,
    durations=durations,
    prefix=prefix,
    subset=rlang::enexpr(subset), # Capture subset expression
    lag=lag,
    precision = precision,
    summate = summate,
    contrasts=cset)
  
  class(ret) <- c("afni_hrfspec", "hrfspec", "list")
  ret
  
}

#' construct a native AFNI hrf specification for '3dDeconvolve' and individually modulated events using the 'stim_times_IM' argument.
#' 
#' @param label name of regressor
#' @param basis Character string specifying the basis function type
#' @param onsets Numeric vector of event onset times
#' @param durations Numeric vector of event durations (default 0)
#' @param subset Expression for subsetting events
#' @param id Character string identifier for the term
#' @param precision Numeric precision for convolution (default 0.3)
#' @param summate Logical whether to summate overlapping responses (default TRUE)
#' @param start start of hrf (for multiple basis hrfs)
#' @param stop end of hrf (for multiple basis hrfs)
#' 
#' @examples
#' 
#' tw <- afni_trialwise("trialwise", basis="gamma", onsets=seq(1,100,by=5))
#' 
#' @export
#' @return an \code{afni_trialwise_hrfspec} instance
afni_trialwise <- function(label, basis=c("spmg1", "block", "dmblock", "gamma", "wav"),
                     onsets=NULL, durations=0, subset=NULL,
                      id=NULL, start=0, stop=22,
                      precision = 0.3, summate = TRUE) {
  
  ## TODO cryptic error message when argument is mispelled and is then added to ...
  basis <- match.arg(basis)
  basis <- tolower(basis)
  if (basis == "gam") basis <- "gamma"
  
  hrf <- if (!is.null(durations)) {
    assert_that(length(durations) == 1, msg="afni_trialwise does not currently accept variable durations")
    get_AFNI_HRF(basis, nbasis=1, duration=durations[1], b=start, c=stop)
  } else {
    get_AFNI_HRF(basis, nbasis=1, b=start, c=stop)
  }
  
  
  if (is.null(id)) {
    id <- label
  }  
  
  ret <- list(
    name=label,
    varname=label,
    id=id,
    hrf=hrf,
    onsets=onsets,
    durations=durations,
    subset=rlang::enexpr(subset),
    precision = precision,
    summate = summate)
  
  class(ret) <- c("afni_trialwise_hrfspec", "hrfspec", "list")
  ret
  
}

#' Construct AFNI HRF specification
#' 
#' Internal method for constructing AFNI HRF specifications.
#' 
#' @param x An afni_hrfspec object
#' @param model_spec Model specification
#' @param ... Additional arguments
#' @keywords internal
#' @export
construct.afni_hrfspec <- function(x, model_spec, ...) {
  
  # Assuming construct_event_term is defined in event_model_helpers.R
  # Note: construct_event_term might need adaptation if afni_hrfspec structure differs slightly
  et <- construct_event_term(x, model_spec)
  
  ## do not convolve an afni term
  ##cterm <- convolve(et, x$hrf, model_spec$sampling_frame, summate=x$summate)
  
  ret <- list(
    varname=et$varname,
    evterm=et,
    sampling_frame=model_spec$sampling_frame,
    hrfspec=x,
    precision = x$precision,
    summate = x$summate,
    lag = x$lag,
    contrasts=x$contrasts,
    id=if(!is.null(x$id)) x$id else et$varname
  )
  
  class(ret) <- c("afni_hrf_convolved_term", "convolved_term", "fmri_term", "list") 
  ret
}


#' Construct AFNI trialwise HRF specification
#' 
#' Internal method for constructing AFNI trialwise HRF specifications.
#' 
#' @param x An afni_trialwise_hrfspec object
#' @param model_spec Model specification
#' @param ... Additional arguments
#' @keywords internal
#' @export
construct.afni_trialwise_hrfspec <- function(x, model_spec, ...) {
  
  ## copied almost verbatim from construct.hrfspec
  onsets <- if (!is.null(x$onsets)) x$onsets else model_spec$onsets
  durations <- if (!is.null(x$durations)) x$durations else model_spec$durations
  
  trial_index <- factor(seq(1, length(onsets)))
  
  # Use the passed hrfspec varname (label)
  varlist <- list(trial_index)
  names(varlist) <- x$varname
  
  # Evaluate subset using model_spec$data and formula_env
  subs <- if (!is.null(x$subset) && !rlang::is_null(x$subset)) {
            eval_env <- rlang::env_bury(model_spec$formula_env %||% rlang::empty_env(), !!!model_spec$data)
            tryCatch(rlang::eval_tidy(x$subset, env = eval_env), 
                     error = function(e) stop(sprintf("Failed to evaluate subset expression for afni_trialwise term '%s': %s", x$name, e$message)))
          } else {
            rep(TRUE, length(onsets))
          }
  
  et <- event_term(varlist, onsets, model_spec$blockids, durations, subs)
  #cterm <- convolve(et, x$hrf, model_spec$sampling_frame)
  
  ret <- list(
    varname=et$varname,
    evterm=et,
    sampling_frame=model_spec$sampling_frame,
    hrfspec=x,
    precision = x$precision,
    summate = x$summate,
    id=x$id
  )
  
  class(ret) <- c("afni_trialwise_convolved_term", "convolved_term", "fmri_term", "list") 
  ret
}

#' @keywords internal
#' @noRd
AFNI_SPMG1 <- function(d=1) AFNI_HRF(name="SPMG1", nbasis=as.integer(1), params=list(d=d)) 

#' @keywords internal
#' @noRd
AFNI_SPMG2 <- function(d=1) AFNI_HRF(name="SPMG2", nbasis=as.integer(2), params=list(d=d))

#' @keywords internal
#' @noRd
AFNI_SPMG3 <- function(d=1) AFNI_HRF(name="SPMG3", nbasis=as.integer(3), params=list(d=d))

#' @keywords internal
#' @noRd
AFNI_BLOCK <- function(d=1,p=1) AFNI_HRF(name="BLOCK", nbasis=as.integer(1), params=list(d=d,p=p))

#' @keywords internal
#' @noRd
AFNI_dmBLOCK <- function(d=1,p=1) AFNI_HRF(name="dmBLOCK", nbasis=as.integer(1), params=list(d=d,p=p))

#' @keywords internal
#' @noRd
AFNI_TENT <- function(b=0,c=18, n=10) AFNI_HRF(name="TENT", nbasis=as.integer(n), params=list(b=b,c=c,n=n))

#' @keywords internal
#' @noRd
AFNI_CSPLIN <- function(b=0,c=18, n=6) AFNI_HRF(name="CSPLIN", nbasis=as.integer(n), params=list(b=b,c=c,n=n))

#' @keywords internal
#' @noRd
AFNI_POLY <- function(b=0,c=18, n=10) AFNI_HRF(name="POLY", nbasis=as.integer(n), params=list(b=b,c=c,n=n))

#' @keywords internal
#' @noRd
AFNI_SIN <- function(b=0,c=18, n=10) AFNI_HRF(name="SIN", nbasis=as.integer(n), params=list(b=b,c=c,n=n))

#' @keywords internal
#' @noRd
AFNI_GAM <- function(p=8.6,q=.547) AFNI_HRF(name="GAM", nbasis=as.integer(1), params=list(p=p,q=q))

#' @keywords internal
#' @noRd
AFNI_WAV <- function(d=1) AFNI_HRF(name="WAV", nbasis=as.integer(1), params=list(d=1))


#' @keywords internal
#' @noRd
get_AFNI_HRF <- function(name, nbasis=1, duration=1, b=0, c=18) {
  name <- tolower(name)
  if (name == "sin") name <- "sine"
  if (name == "gam") name <- "gamma"
  hrf <- switch(name,
                gamma=AFNI_GAM(),
                spmg1=AFNI_SPMG1(d=duration),
                spmg2=AFNI_SPMG2(d=duration),
                spmg3=AFNI_SPMG3(d=duration),
                csplin=AFNI_CSPLIN(b=b,c=c, n=nbasis),
                poly=AFNI_POLY(b=b,c=c, n=nbasis),
                sine=AFNI_SIN(b=b,c=c,n=nbasis),
                wav=AFNI_WAV(),
                block=AFNI_BLOCK(d=duration),
                dmblock=AFNI_dmBLOCK(),
                # Add tent as an alias or explicit entry
                tent=AFNI_TENT(b=b, c=c, n=nbasis))
  
  if (is.null(hrf)) {
    stop("could not find afni hrf named: ", name)
  }
  
  hrf
  
}







#           
#   afni.stims <- unlist(lapply(funcTerms, function(term) { buildAFNIStims(term, opts$iresp, opts$TR_times ) }))
# 
#   purgeNulls <- function(A) {
#     A[!sapply(A, is.null)]
#   }
# 
# 
#   opt_stim_labels <-  purgeNulls(lapply(seq_along(afni.stims), function(i) buildCommandSwitch(afni.stims[[i]], i, "label")))
#   opt_stim_files  <-  purgeNulls(lapply(seq_along(afni.stims), function(i) buildCommandSwitch(afni.stims[[i]], i, "file")))
#   opt_stim_times  <-  purgeNulls(lapply(seq_along(afni.stims), function(i) buildCommandSwitch(afni.stims[[i]], i, "times")))
#   opt_stim_iresp  <-  purgeNulls(lapply(seq_along(afni.stims), function(i) buildCommandSwitch(afni.stims[[i]], i, "iresp")))
# 
# 
#   cmdlines <- list(input=filelist(x@design, full.names=T),
#                              mask=x@mask,
#                              polort=opts[["polort"]],
#                              num_stimts=length(afni.stims),
#                              num_glt=length(gltlist),
#                              stim_file=opt_stim_files,
#                              stim_label=opt_stim_labels,
#                              stim_times=opt_stim_times,
#                              TR_times=opts[["TR_times"]],
#                              iresp=opt_stim_iresp,
#                              gltsym=lapply(seq_along(gltfiles), function(i) paste(gltfiles[i], collapse=" ")),
#                              glt_label=lapply(seq_along(gltnames), function(i) paste(i, gltnames[i], collapse=" ")),
#                              nofullf_atall=opts[["nofullf_atall"]],
#                              fout=opts[["fout"]],
#                              rout=opts[["rout"]],
#                              tout=opts[["tout"]],
#                              bout=opts[["bout"]],
#                              noFDR=opts[["noFDR"]],
#                              cbucket=opts[["cbucket"]],
#                              bucket=opts[["bucket"]],
#                              jobs=opts[["jobs"]],
#                              float=TRUE)
# 
# 
#             cmdstr <- .makeCommandStr(cmdlines)
# 
#             ret <- list()
#             wd <- workingDir(x)
# 
#             nextDirName <- function(wd) {
#               nd <- paste(wd, "+", sep="")
#               if (!file.exists(nd)) {
#                 nd
#               } else {
#                 Recall(nd)
#               }
#             }
# 
#             writeStimFiles <- function() {
#               sapply(afni.stims, function(stim) {
#                 writeAFNIStim(stim, ".")
#               })
#             }
# 
#
#             writeGLTs <- function() {
#               lapply(seq_along(gltlist), function(i) {
#                 fout <- file(gltfiles[i], "w")
#                 .glt <- gltlist[[i]]
# 
#                 write(unlist(.glt), file=fout, sep="\n")
# 
#                 close(fout)
#               })
#             }
# 
# 



   
