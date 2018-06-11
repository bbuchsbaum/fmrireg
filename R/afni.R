

#' afni_lm
#' 
#' 
#' @param fmri_mod
#' @param dataset
#' @param working_dir
#' @param polort
#' @param jobs
#' @param options
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
afni_lm <- function(fmri_mod, dataset, working_dir=".", polort=-1, jobs=1, options=list()) {
  
  defopts <- list(noFDR=FALSE, 
                       fout=TRUE, 
                       rout=TRUE, 
                       tout=TRUE, 
                       float=TRUE, 
                       cbucket="coefout", 
                       bucket="statout", 
                       jobs=jobs, 
                       polort=polort, 
                       iresp=FALSE, 
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



#' @keywords internal
afni_stim_file <- function(label, file_name, values) {
  structure(
    list(label=label, file_name=file_name, values=values),
    class="afni_stim_file"
  )
}

#' @keywords internal
afni_stim_times <- function(label, file_name, hrf, onsets, iresp=FALSE, tr_times=1) {
  structure(
    list(label=label, file_name=file_name, hrf=hrf, onsets=onsets, iresp=iresp, tr_times=tr_times),
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
afni_command_switch <- function(x, k, type) {
  switch(type,
    label = paste(k, x$label, collapse = " "),
    file = paste(k, x$file_name, collapse = " "),
    ortvec = paste(x$file_name, x$label),
    times = NULL,
    iresp = NULL
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
  sapply(blist, function(bm) {
    write_baseline_mat(bm, ".")
  })
}
  
#' @keywords internal
write_stim_files <- function(afni_stims) {
  sapply(afni_stims, function(stim) {
    write_afni_stim(stim, ".")
  })
}

#' @keywords internal
write_afni_stim <- function(stim, dir) {
  .write_values <- function(outname, vals) {
    hfile <- file(outname, "w")
    write(vals, file=hfile, ncolumns=1)
    close(hfile)
  }
  
  .write_values(paste0(dir, "/", stim$file_name, sep=""), stim$values)
}

#' @keywords internal
write_glts <- function(gltstr, gltfiles) {
  lapply(seq_along(gltstr), function(i) {
    fout <- file(gltfiles[i], "w")
    .glt <- gltstr[[i]]
    write(.glt, file=fout, sep="\n")
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
build_afni_stims <- function(x) {
  stimlabels <- longnames(x)
  stimfiles <- paste(stimlabels, "_reg.1D", sep = "")
  desmat <- design_matrix(x)
  
  lapply(1:length(stimlabels), function(i) {
    afni_stim_file(stimlabels[i], stimfiles[i], desmat[, i])
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
  assert_that(sum(duplicated(gltfiles))  == 0, msg="Cannot have two GLTs with the same file name")
  
  func_terms <- terms(model$event_model)
  
  afni_stims <- unlist(lapply(func_terms, function(term) { build_afni_stims(term) }), recursive=FALSE)
  afni_baseline_mats <- build_baseline_stims(model)
  
  
  purgeNulls <- function(A) {
    A[!sapply(A, is.null)]
  }
  
  opt_stim_labels <-  purgeNulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "label")))
  opt_stim_files  <-  purgeNulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "file")))
  opt_stim_ortvecs <- purgeNulls(lapply(seq_along(afni_baseline_mats), function(i) afni_command_switch(afni_baseline_mats[[i]], i, "ortvec")))

  
  cmdlines <- list(input=paste0(dataset$base_path, "/", dataset$scans),
                   mask=paste0(dataset$base_path, "/", dataset$mask_file),
                   polort=opts[["polort"]],
                   num_stimts=length(afni_stims),
                   num_glt=length(gltfiles),
                   stim_file=opt_stim_files,
                   stim_label=opt_stim_labels,
                   ortvec=opt_stim_ortvecs,
                   #stim_times=opt_stim_times,
                   TR_times=opts[["TR_times"]],
                   #iresp=opt_stim_iresp,
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
                   jobs=opts[["jobs"]],
                   float=TRUE)
  
  cmd <- .make_decon_command_str(cmdlines)
  
  list(cmd=cmd, afni_stims=afni_stims, afni_baseline_mats=afni_baseline_mats,
       gltfiles=gltfiles, gltnames=gltnames, glts=glts, gltstr=gltstr)
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



   
