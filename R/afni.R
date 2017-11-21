afni_lm <- function(fmri_mod, working_dir=".", options=list()) {
  
  default_opts <- list(noFDR=FALSE, 
                       fout=TRUE, 
                       rout=TRUE, 
                       tout=TRUE, 
                       float=TRUE, 
                       cbucket="coefout", 
                       bucket="statout", 
                       jobs=2, 
                       polort=2, 
                       iresp=FALSE, 
                       TR_times=1)
  
  for (optname in names(options)) {
    defopts[[optname]] <- options[[optname]]
  }
  
  structure(
    list(
      model=fmri_mod,
      working_dir=working_dir,
      options=options),
    class="afni_lm")
}

afni_stim_file <- function(label, file_name, values) {
  structure(
    list(label=label, file_name=file_name, values=values),
    class="afni_stim_file"
  )
}

afni_stim_times <- function(label, file_name, hrf, onsets, iresp=FALSE, tr_times=1) {
  structure(
    list(label=label, file_name=file_name, hrf=hrf, onsets=onsets, iresp=iresp, tr_times=tr_times),
    class="afni_stim_times"
  )
}

afni_command_switch.afni_stim_times <- function(x, k, type) {
  switch(type,
    label = paste(k, x$label, collapse = " "),
    times = paste(k, x$file_name, as.character(x@hrf), collapse = " "),
    file = NULL,
    iresp = if (x$iresp) 
      paste(k, paste(paste(x$label, "_iresp", sep = ""), collapse = " "))
    else
      NULL
  )
}

afni_command_switch.afni_stim_file <- function(x, k, type) {
  switch(
    type,
    label = paste(k, x$label, collapse = " "),
    file = paste(k, x$file_name, collapse = " "),
    times = NULL,
    iresp = NULL
  )
}



build_afni_stims <- function(x) {
  stimlabels <- longnames(x)
  stimfiles <- paste(stimlabels, "_reg.1D", sep = "")
  desmat <- convolve(x)
  
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


build_decon_command <- function(x, opts) {
   stimlabels <- unlist(lapply(terms(x$event_model), longnames))
   
   assert_that(length(unique(stimlabels)) == length(stimlabels))
   assert_that(length(stimlabels) == length(conditions(x$design)))
#   
#   gltlist <- glts(x)
#   gltnames <- if (length(gltlist) > 0) names(gltlist) else NULL
#   gltfiles <- if (length(gltnames) > 0) paste("glt_", gltnames, ".txt", sep="") else NULL
# 
   func_terms <- terms(x$event_model)
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
#             ret$run <- function() {
#               startDir <- getwd()
# 
#               res <- try({
# 
#                 if (!file.exists(wd)) {
#                   dir.create(wd)
#                 } else {
#                   warning(paste("glm output directory: ", wd, " already exists"))
#                   wd <- nextDirName(wd)
#                   dir.create(wd)
#                   warning(paste("outputting to: ", wd))
# 
#                 }
# 
#                 print(paste("setting directory:", wd))
#                 setwd(wd)
# 
#                 writeStimFiles()
#                 writeGLTs()
# 
# 
# 
#                 write(cmdstr, "3ddeconvolve.sh")
#                 system(cmdstr)
#               })
# 
#               setwd(startDir)
# 
# 
#             }
# 
#             ret$command <- cmdstr
#             #ret$designMat <- desmat
#             ret
# 
#           })
# 
# setMethod("workingDir", signature(x="AFNICommand"),
#           function(x) {
#             x@workingDir
#             
#           })
# 
# setMethod("getOptions", signature(x="AFNICommand"),
#           function(x) {
#             x@options           
#           })
# 
# setMethod("glts", signature(x="AFNICommand"),
#           function(x) {
#             x@glts
#           })

   # build_afni_stims_native <- function(x, iresp=FALSE, tr_times=1) {
   #           function(x, iresp=FALSE, tr.times=1) {
   #             
   #             #stimlabels <- shortnames(x, exclude.basis=TRUE)
   #             split.ons <- splitOnsets(x, global=TRUE)
   #             stimlabels <- names(split.ons)
   #             stimfiles <- paste(stimlabels, "_reg.1D", sep="")		
   #             
   #             
   #             lapply(1:length(stimlabels), function(i) {
   #               AFNIStimTimes(stimlabels[i], stimfiles[i], x@hrf, split.ons[[stimlabels[[i]]]], iresp, tr.times)
   #             })
   #           })
   
   
   # 
   # .buildBLOCK <- function(stimlabel, stimfile, stimnum, duration) {
   #   paste(paste("-stim_times", stimnum, stimfile,
   #               paste("\'BLOCK(", duration, ")\'"),
   #               paste("-stim_label", stimnum, stimlabel)))
   # }
   # 
   # .buildGAM <- function(stimlabel, stimfile, stimnum) {
   #   paste(paste("-stim_times", stimnum, stimfile,
   #               paste("\"GAM\""),
   #               paste("-stim_label", stimnum, stimlabel)))
   # }
   # 
   # .buildSPMG <- function(stimlabel, stimfile, stimnum) {
   #   paste(paste("-stim_times", stimnum, stimfile,
   #               paste("\"SPMG\""),
   #               paste("-stim_label", stimnum, stimlabel)))
   # }
   # .buildSPMG3 <- function(stimlabel, stimfile, stimnum) {
   #   paste(paste("-stim_times", stimnum, stimfile,
   #               paste("\"SPMG3\""),
   #               paste("-stim_label", stimnum, stimlabel)))
   # }
   # 
   # .buildBasis <- function(basisname, stimlabel, stimfile, stimnum, start, end, nparams) {
   #   paste("-stim_times", stimnum, stimfile,
   #         paste("\'", basisname, "(", start, ",", end, ",", nparams, ")\'", sep=""),
   #         paste("-stim_label", stimnum, stimlabel),
   #         paste("-iresp", stimnum, paste(stimlabel, "_iresp", sep="")))
   #   
   # }
   # 
   # setMethod("writeAFNIStim", signature(x="AFNIStimFile", dir="character"),
   #           function(x, dir) {
   #             
   #             .writeOnsets <- function(outname, vals) {
   #               hfile = file(outname, "w")
   #               write(x@values, file=hfile, ncolumns=1)
   #               close(hfile)
   #             }
   #             
   #             .writeOnsets(paste(dir, "/", x@fileName, sep=""), x@values)
   #           })
   # 
   # 
   # setMethod("writeAFNIStim", signature(x="AFNIStimTimes", dir="character"),
   #           function(x, dir) {
   #             
   #             .writeOnsets <- function(outname, onsets) {
   #               hfile = file(outname, "w")
   #               write(onsets, file=hfile, ncolumns=1)
   #               close(hfile)
   #             }
   #             
   #             .writeOnsets(paste(dir, "/", x@fileName, sep=""), x@onsets)
   #           })
   # 
   # 
   
