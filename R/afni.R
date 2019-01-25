
#' @keywords internal
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
  
  fmodel <- fmri_model(emodel, bmodel)
  alm <- afni_lm(fmodel, dset, censor=dset$censor)
  
}

  

#' afni_lm
#' 
#' 
#' @param fmri_mod
#' @param dataset
#' @param working_dir
#' @param polort
#' @param jobs
#' @param options
#' @param censor
#' 
#' @examples 
#' etab <- data.frame(onset=c(1,30,15,25), fac=factor(c("A", "B", "A", "B")), run=c(1,1,2,2))
#' dset <- fmri_dataset(scans=c("s1.nii", "s2.nii"), mask="mask.nii", TR=1, run_length=c(50,50),event_table=etab)
#' 
#' emodel <- event_model(onset ~ hrf(fac), block = ~ run, data=etab, sampling_frame=dset$sampling_frame)
#' bmodel <- baseline_model("bs", degree=4, sframe=dset$sampling_frame)
#' fmod <- fmri_model(emodel, bmodel)
#' alm <- afni_lm(fmod, dset)
#' @export
afni_lm <- function(fmri_mod, dataset, working_dir=".", polort=-1, jobs=1, censor=NULL, options=list()) {
  
  if (!is.null(censor)) {
    censor <- unlist(censor)
    assert_that(length(censor) == sum(dataset$sampling_frame$blocklens))
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
      options=options,
      cmd=cmd),
    class="afni_lm_spec")

}


#' @export
print.afni_lm_spec <- function(x) {
  cat("AFNI linear model via 3dDeconvolve \n")
  cat("  working_dir: ", x$working_dir, "\n")
  cat("  number of GLTs:", length(x$glts), "\n")
  cat("  command line: \n", x$cmd$cmd, "\n")
  
}



#' @keywords internal
afni_stim_file <- function(label, file_name, values) {
  structure(
    list(label=label, file_name=file_name, values=values),
    class="afni_stim_file"
  )
}

#' @keywords internal
afni_stim_times <- function(label, file_name, hrf, onsets, iresp=FALSE, sresp=FALSE, tr_times=1) {
  structure(
    list(label=label, file_name=file_name, hrf=hrf, onsets=onsets, iresp=iresp, sresp=sresp, tr_times=tr_times),
    class="afni_stim_times"
  )
}

# afni_command_switch.afni_stim_times <- function(x, k, type) {
#   switch(type,
#     label = paste(k, x$label, collapse = " "),
#     times = paste(k, x$file_name, as.character(x@hrf), collapse = " "),
#     file = NULL,
#     iresp = if (x$iresp) 
#       paste(k, paste(paste(x$label, "_iresp", sep = ""), collapse = " "))
#     else
#       NULL
#   )
# }



#' @keywords internal
afni_command_switch <- function(x, ...) UseMethod("afni_command_switch")

#' @keywords internal
write_afni_stim <- function(x, ...) UseMethod("write_afni_stim")


#' @keywords internal
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
next_dir_name <- function(wd) {
  nd <- paste(wd, "+", sep="")
  if (!file.exists(nd)) {
    nd
  } else {
    Recall(nd)
  }
}

#' @keywords internal
write_baseline_mat <- function(stim, dir) {
  write.table(stim$mat, paste0(dir, "/", stim$file_name, sep=""), col.names=FALSE, row.names=FALSE)
}

#' @keywords internal
write_baseline_mats <- function(blist) {
  purrr::walk(blist, ~ write_baseline_mat(., "."))
}
  
#' @keywords internal
write_stim_files <- function(afni_stims) {
  purrr::walk(afni_stims, ~ write_afni_stim(., "."))
}

#' @keywords internal
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
write_afni_stim.afni_stim_times <- function(stim, dir) {
  .write_onsets <- function(outname, vals) {
    hfile <- file(outname, "w")
    write(vals, file=hfile, ncolumns=1)
    close(hfile)
  }
  
  .write_onsets(paste0(dir, "/", stim$file_name, sep=""), unlist(stim$onsets))
}

#' @keywords internal
write_censor_file <- function(dir, censor) {
    outname <- paste0(dir, "/censor.1D")
    hfile <- file(outname, "w")
    write(censor, file=hfile, ncolumns=1)
    close(hfile)
}
  

#' @keywords internal
#' @importFrom purrr imap
write_glts <- function(glts, gltfiles) {
  imap(glts, function(glt, i) {
    fout <- file(gltfiles[i], "w")
    write(glt, file=fout, sep="\n")
    close(fout)
  })
}

#' @keywords internal
afni_baseline_matrix <- function(label, file_name, mat) {
  structure(
    list(label=label, file_name=file_name, mat=mat),
    class="afni_baseline_matrix"
  )
}

#' @keywords internal
build_baseline_stims <- function(x) {
  blens <- blocklens(x)
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

#' @keywords internal
build_afni_stims.convolved_term <- function(x) {
  stimlabels <- longnames(x)
  stimfiles <- paste(stimlabels, "_reg.1D", sep = "")
  desmat <- design_matrix(x)
  
  lapply(1:length(stimlabels), function(i) {
    afni_stim_file(stimlabels[i], stimfiles[i], desmat[, i])
  })
}

#' @keywords internal
build_afni_stims.afni_hrf_convolved_term <- function(x, iresp=FALSE, tr_times=1) {
  #browser()
  stimlabels <- longnames(x)
  stimfiles <- paste(stimlabels, "_times.1D", sep = "")
  dmat <- design_matrix(x$evterm)
  
  hrf_name <- as.character(x$hrf$hrf)
  
  split_ons <- split_onsets(x$evterm, x$sampling_frame, global=TRUE)
  names(split_ons) <- stimlabels
  
  lapply(1:length(stimlabels), function(i) {
    afni_stim_times(stimlabels[i], stimfiles[i], hrf_name, split_ons[[stimlabels[[i]]]], iresp, tr_times)
  })
  
}

#' @keywords internal
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

#' @keywords internal
build_decon_command <- function(model, dataset, working_dir, opts) {
  ## get the set of stimulus regressors
  stimlabels <- unlist(lapply(terms(model$event_model), longnames))
  
  ## all stims must be unique
  assert_that(length(unique(stimlabels)) == length(stimlabels))
  #browser()
  assert_that(length(stimlabels) == length(conditions(model$event_model)))
  
  ## extract all contrast matrices
  cons <- contrast_weights(model)
  cons <- unlist(cons, recursive=FALSE)
  
  ## convert to 'glt's
  glts <- lapply(cons, to_glt)
  
  gltfiles <- unlist(lapply(glts, function(x) paste0(x$name, ".txt")))
  gltnames <- unlist(lapply(glts, function(x) x$name))
  gltstr <- unlist(lapply(glts, function(x) x$glt_str))
  
  assert_that(sum(duplicated(gltnames))  == 0, msg="Cannot have two GLTs with the same name")
  
  func_terms <- terms(model$event_model)
  
  
  ## construct list of afni stims
  afni_stims <- unlist(lapply(func_terms, function(term) { build_afni_stims(term) }), recursive=FALSE)
  afni_baseline_mats <- build_baseline_stims(model)
  
  
  purge_nulls <- function(A) {
    A[!sapply(A, is.null)]
  }
  
  opt_stim_labels <-  purge_nulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "label")))
  opt_stim_files  <-  purge_nulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "file")))
  opt_stim_times  <-  purge_nulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "times")))
  opt_stim_ortvecs <- purge_nulls(lapply(seq_along(afni_baseline_mats), function(i) afni_command_switch.afni_stim_file(afni_baseline_mats[[i]], i, "ortvec")))
  opt_stim_iresp  <-  purge_nulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "iresp")))
  
  #browser()
  
  if (length(opt_stim_times) > 0) {
    global_times=TRUE
  } else {
    global_times <- FALSE
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
                   TR_times=opts[["TR_times"]],
                   iresp=opt_stim_iresp,
                   gltsym=lapply(seq_along(gltfiles), function(i) paste(gltfiles[i], collapse=" ")),
                   glt_label=lapply(seq_along(gltnames), function(i) paste(i, gltnames[i], collapse=" ")),
                   nofullf_atall=opts[["nofullf_atall"]],
                   fout=opts[["fout"]],
                   rout=opts[["rout"]],
                   tout=opts[["tout"]],
                   bout=opts[["bout"]],
                   nofullf_atall=opts[["nofullf_atall"]],
                   noFDR=opts[["noFDR"]],
                   cbucket=opts[["cbucket"]],
                   bucket=opts[["bucket"]],
                   nocond=opts[["nocond"]],
                   x1D_stop=opts[["x1D_stop"]],
                   jobs=opts[["jobs"]],
                   float=TRUE)
  
  cmd <- .make_decon_command_str(cmdlines)
  
  list(cmd=cmd, cmdlines=cmdlines, afni_stims=afni_stims, afni_baseline_mats=afni_baseline_mats,
       gltfiles=gltfiles, gltnames=gltnames, glts=glts, gltstr=gltstr, censor=opts$censor)
}



#' @export
#' @param outdir the output folder
#' @param execute whether to execute the command or only output shell '3dDeconvolve.sh' script
#' @rdname run
run.afni_lm_spec <- function(x, outdir, execute=TRUE) {
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
    
    write(x$cmd$cmd, "3ddeconvolve.sh")
    
    if (execute) {
      system(x$cmd$cmd)
    }
  })
  
  setwd(start_dir)
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



   
