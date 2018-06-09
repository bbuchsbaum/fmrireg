
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
  switch(
    type,
    label = paste(k, x$label, collapse = " "),
    file = paste(k, x$file_name, collapse = " "),
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

#' @export
write_baseline_model <- function(bmod, dir, fname) {
  mat <- design_matrix(bmod)
  
  oname <- paste0(dir, "/", fname)
  write.table(mat, file=oname, col.names=FALSE)
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
write_glts <- function(glts, gltfiles) {
  lapply(seq_along(glts), function(i) {
    fout <- file(gltfiles[i], "w")
    .glt <- glts[[i]]
    write(.glt$glt_str, file=fout, sep="\n")
    close(fout)
  })
}

build_baseline_stims <- function(x) {
  bmat <- design_matrix(x)
}

#' @keywords internal
build_afni_stims <- function(x) {
  stimlabels <- longnames(x$event_model)
  stimfiles <- paste(stimlabels, "_reg.1D", sep = "")
  desmat <- design_matrix(x$event_model)
  
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
  
  cons <- contrast_weights(model)
  cons <- unlist(cons, recursive=FALSE)
  
  glts <- lapply(cons, to_glt)
  gltfiles <- sapply(glts, function(x) paste0(x$name, ".txt"))
  gltnames <- sapply(glts, function(x) x$con$name)
  
  assert_that(sum(duplicated(gltnames))  == 0, msg="Cannot have two GLTs with the same name")
  assert_that(sum(duplicated(gltfiles))  == 0, msg="Cannot have two GLTs with the same file name")
  
  func_terms <- terms(model$event_model)
  
  afni_stims <- unlist(lapply(func_terms, function(term) { build_afni_stims(term) }), recursive=FALSE)
  
  purgeNulls <- function(A) {
    A[!sapply(A, is.null)]
  }
  
  opt_stim_labels <-  purgeNulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "label")))
  opt_stim_files  <-  purgeNulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "file")))
  
  cmdlines <- list(input=paste0(dataset$base_path, "/", dataset$scans),
                   mask=paste0(dataset$base_path, "/", dataset$mask_file),
                   polort=opts[["polort"]],
                   num_stimts=length(afni_stims),
                   num_glt=length(glts),
                   stim_file=opt_stim_files,
                   stim_label=opt_stim_labels,
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
  list(cmd=cmd, afni_stims=afni_stims, gltfiles=gltfiles, gltnames=gltnames, glts=glts)
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
    write_glts(x$cmd$glts, x$cmd$gltfiles)
    
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



   
