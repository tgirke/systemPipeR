###########################################################
## Additional utilities for SYSargs and SYSargs2 objects ##
###########################################################

######################################################
## Convenience write function for targetsout(args) ##
######################################################
writeTargetsout <- function (x, file = "default", silent = FALSE, overwrite = FALSE, step = NULL, new_col=NULL,
                             new_col_output_index=NULL, remove=FALSE, ...) {
  if(all(class(x) != "SYSargs" & class(x) != "SYSargs2")) stop("Argument 'x' needs to be assigned an object of class 'SYSargs' OR 'SYSargs2")
  ## SYSargs class
  if(class(x) == "SYSargs") {
    targets <- targetsout(x)
    software <- software(x)
    if(file == "default") {
      file <- paste("targets_", software, ".txt", sep = "")
      file <- gsub(" {1,}", "_", file)
    } else {
      file <- file
    }
    headerlines <- targetsheader(x)
    ## SYSargs2 class
  } else if(class(x) == "SYSargs2") {
    if(is.null(step))
      stop(paste("Argument 'step' needs to be assigned one of the following values:",
                 paste(names(x$clt), collapse = ", "), "OR the corresponding position"))
    if(all(!is.null(step) & is.character(step) & !any(names(x$clt) %in% step)))
      stop(paste("Argument 'step' can only be assigned one of the following values:",
                 paste(names(x$clt), collapse = ", "), "OR the corresponding position"))
    if(all(!is.null(step) & is.numeric(step) & !any(seq_along(names(x$clt)) %in% step)))
      stop(paste("Argument 'step' can only be assigned one of the following position:",
                 paste(seq_along(names(x$clt)), collapse = ", "), "OR the corresponding names"))
    targets <- targets.as.df(targets(x))
    if(remove==TRUE){
      targets <- targets[,-c(which(grepl("FileName", colnames(targets))))]
    }
    ## Adding the columns
    if ((!is.null(new_col) & is.null(new_col_output_index)) |
        (is.null(new_col) & !is.null(new_col_output_index)) |
        (is.null(new_col) & is.null(new_col_output_index))){
      cat("One of 'new_col' and 'new_col_output_index' is null. It is using default column naming and adding all the output files expected, and each one will be written in a different column. \n")
      for (i in seq_len(length(output(x)[[1]][[step]]))){
        pout <- sapply(names(output(x)), function(y) output(x)[[y]][[step]][[i]], simplify = FALSE)
        targets[[paste0(files(x)$steps, "_", i)]] = as.character(pout)
      }
    } else if(!is.null(new_col) & !is.null(new_col_output_index)){
      if(any(length(output(x)[[1]][[step]]) < new_col_output_index) | any(new_col_output_index < 1)) {
        stop(paste0("'new_col_output_index' argument needs to be equal or bigger than 1 and smaller than ", length(output(x)[[1]][[1]]), ", the maximum number of outputs files." ))
      }
      if(length(new_col) != length(new_col_output_index)){
        stop("'new_col' should have the same length as 'new_col_output_index'")
      }
      for(i in seq_along(new_col)){
        pout <- sapply(names(output(x)), function(y) output(x)[[y]][[step]][[new_col_output_index[i]]], simplify = FALSE)
        targets[[as.character(new_col[i])]] = as.character(pout)
      }
    }
    ## Workflow and Step Name
    software <- strsplit(basename(files(x)$cwl), split = "\\.")[[1]][1]
    if(is.character(step)) {
      step <- strsplit(step, split = "\\.")[[1]][1]
    } else {
      step <- strsplit(names(x$clt)[step], split = "\\.")[[1]][1]
    }
    if(file == "default") {
      file <- paste("targets_", step, ".txt", sep = "")
      file <- gsub(" {1,}", "_", file)
    } else {
      file <- file
    }
    headerlines <- targetsheader(x)[[1]]
  }
  if(file.exists(file) & overwrite == FALSE) stop(paste("I am not allowed to overwrite files; please delete existing file:", file, "or set 'overwrite=TRUE'"))
  names <- c(new_col,  colnames(targets[, -c(which(grepl(paste(new_col,collapse="|"), colnames(targets))))]))
  targets <- cbind(targets[, new_col], targets[, -c(which(grepl(paste(new_col,collapse="|"), colnames(targets))))])
  colnames(targets) <- names
  targetslines <- c(paste(colnames(targets), collapse = "\t"), apply(targets, 1, paste, collapse = "\t"))
  writeLines(c(headerlines, targetslines), file, ...)
  if(silent != TRUE) cat("\t", "Written content of 'targetsout(x)' to file:", file, "\n")
}
## Usage:
# writeTargetsout(x=args, file="default") ## SYSargs class
# writeTargetsout(x=WF, file="default", step=1, new_col = "FileName1", new_col_output_index = 1) ## SYSargs2 class

##############################################################################
## Function to run NGS aligners including sorting and indexing of BAM files ##
##############################################################################
runCommandline <- function(args, runid="01", make_bam=FALSE, del_sam=TRUE, dir=TRUE, dir.name=NULL, force=FALSE, ...) {
  ## Validation for 'args'
  if(any(!inherits(args, "SYSargs") & !inherits(args, "SYSargs2"))) stop("Argument 'args' needs to be assigned an object of class 'SYSargs' OR 'SYSargs2'")
    ## Environment Modules section
  .moduleload(args)
  ## SYSargs class ##
  if(class(args)=="SYSargs") {
    .sysargsrunCommandline(args=args, runid=runid, make_bam=make_bam, del_sam=del_sam)
  } else if(class(args)=="SYSargs2") {
    ## SYSargs2 class ##
    ## Check if the command is in the PATH
    if(!baseCommand(args) == c("bash")){
      tryCatch(system(baseCommand(args), ignore.stdout = TRUE, ignore.stderr = TRUE), warning=function(w) message("\n", paste0("ERROR: ", "\n", baseCommand(args), ": command not found. ", '\n', "Please make sure to configure your PATH environment variable according to the software in use."), "\n"))
    }
    ## Workflow Name (Workflow OR CommandLineTool class)
    if(length(wf(args)$steps)==0){
      cwl.wf <- gsub( "[[:space:]]", "_", paste(baseCommand(args), collapse = "_"), perl=TRUE)
    } else {
      cwl.wf <- strsplit(basename(files(args)$cwl), split="\\.")[[1]][1]
    }
    ## Check if "results" folders exists...
    if(!dir.exists(normalizePath(file.path(yamlinput(args)$results_path$path))))
      stop("Please check our current directory ('getwd()').
           The PATH defined at `yamlinput(args)` was not found and it is required.")
    if(is.null(dir.name)) {
      logdir <- normalizePath(yamlinput(args)$results_path$path)
      dir.name <- cwl.wf
    } else {
      logdir <- normalizePath(yamlinput(args)$results_path$path)
      dir.name <- dir.name
    }
    ## Check if expected files exists or not
    ## Important, this happen in the new location in the case of dir=TRUE
    return <- .checkOutArgs2(args, make_bam=make_bam, dir=dir, dir.name=dir.name, force=force)
    args.return <- return$args_complete
    completed <- return$completed
    args <- return$args
    ## Create log files
    file_log <- file.path(logdir, paste0("submitargs", runid, "_", dir.name, "_log_", format(Sys.time(), "%b%d%Y_%H%Ms%S")))
    ## Sample status and Time
    sample_status <- sapply(names(cmdlist(args)), function(x) list(NULL))
    time_status <- data.frame(Targets=names(cmdlist(args)), time_start=NA, time_end=NA)
    ## Progress bar
    pb <- txtProgressBar(min = 0, max = length(cmdlist(args)), style = 3)
    ## Check input
    if(length(args$inputvars) >= 1){
      inpVar <- args$inputvars
      check.inp <- colSums(sapply(inpVar, function(y) sapply(yamlinput(args), function(x) x["path"] %in% y)))
      check.inp[check.inp > 0]
      df.targets <- targets.as.df(args$targets)[check.inp[check.inp > 0]]
      inp_targets2 <- FALSE
    } else {
      inp_targets2 <- TRUE
    }
    ## LOOP
    for(i in seq_along(cmdlist(args))){
      setTxtProgressBar(pb, i)
      cat("## ", names(cmdlist(args)[i]), "\n", file=file_log, fill=TRUE, append=TRUE)
      ## Time
      time_status$time_start[i] <- Sys.time()
      for(j in seq_along(cmdlist(args)[[i]])){
        ## Run the commandline only for samples for which no output file is available.
        if(all(force==FALSE && all(as.logical(completed[[i]][[j]])))) {
          cat("The expected output file(s) already exist", file=file_log, fill=TRUE, append=TRUE, "/n")
          sample_status[[i]][[args$files$steps[j]]] <- "Success"
          next()
        } else {
          ## Add_command file at log file
          cat(c(
                paste0("Time: ", paste0(format(Sys.time(), "%b%d%Y_%H%Ms%S"))), "\n",
                "### Code: ",
                "```{r, eval=FALSE} ",
                cmdlist(args)[[i]][[j]][[1]],
                "```", "\n",
                "### Stdout: ",
                "```{r, eval=FALSE}" ), file = file_log, sep = "\n", append = TRUE)
          ## Create an object for executable
          command <- gsub(" .*", "", as.character(cmdlist(args)[[i]][[j]]))
          commandargs <- gsub("^.*? ", "",as.character(cmdlist(args)[[i]][[j]]))
          ## Return an error if input file doesn't exist
          if(all(inp_targets2) || all(inp_targets <- file.exists(as.character(df.targets[i,])))){
            stdout <- .tryRunC(command, commandargs)
          } else {
            stdout <- list(stdout = paste(paste0(as.character(df.targets[i,])[!inp_targets], collapse = ", "), "\n are missing"),
                           warning= "", error= "We have an error because target(s) file doesn't exit" )
          }
          cat(stdout$stdout, file=file_log, sep = "\n", append=TRUE)
          if(!is.null(stdout$error)) {
            cat("## Error", file=file_log, sep = "\n", append=TRUE)
            cat(stdout$error, file=file_log, sep = "\n", append=TRUE)
            sample_status[[i]][[args$files$steps[j]]] <- "Error"
          } else if(!is.null(stdout$warning)) {
            cat("## Warning", file=file_log, sep = "\n", append=TRUE)
            cat(stdout$warning, file=file_log, sep = "\n", append=TRUE)
            sample_status[[i]][[args$files$steps[j]]] <- "Warning"
          } else if(all(is.null(stdout$error) && is.null(stdout$warning))){
            sample_status[[i]][[args$files$steps[j]]] <- "Success"
          }
          ## close R chunk
          cat("```", file=file_log, sep = "\n", append=TRUE)
        }
        ## converting sam to bam using Rsamtools package
        .makeBam(output(args)[[i]][[j]], make_bam=make_bam, del_sam=del_sam)
      }
      time_status$time_end[i] <- Sys.time()
    }
    ## Status and log.files
    df.status <- data.frame(matrix(do.call("c", sample_status), nrow=length(sample_status), byrow=TRUE))
    colnames(df.status) <- files(args.return)$steps
    check <- check.output(args)
    df.status.f <- cbind(check, df.status)
    df.status.f[c(2:4)] <- sapply(df.status.f[c(2:4)],as.numeric)
    ## time
    time_status$time_start <- as.POSIXct(time_status$time_end, origin="1970-01-01")
    time_status$time_end <- as.POSIXct(time_status$time_end, origin="1970-01-01")
    ## updating object
    args.return[["status"]]$status.summary <- .statusSummary(df.status.f)
    args.return[["status"]]$status.completed <- df.status.f
    args.return[["status"]]$status.time <- time_status
    args.return[["files"]][["log"]] <- file_log
    ## Create recursive the subfolders
    if(dir==TRUE){
      for(i in seq_along(names(cmdlist(args)))){
        full_path <- file.path(logdir, dir.name)
        if(dir.exists(full_path)==FALSE){
          dir.create(full_path, recursive = TRUE) }
      }
      ## copy log files
      if(dir.exists(file.path(logdir, dir.name, "_logs"))==FALSE){
        dir.create(file.path(logdir, dir.name, "_logs"), recursive = TRUE)
      }
      file.rename(from=file_log, to=file.path(logdir, dir.name, "_logs", basename(file_log)))
      args.return[["files"]][["log"]] <- file.path(logdir, dir.name, "_logs", basename(file_log))
      ## output FILES
      if(make_bam==TRUE) args.return <- .checkOutArgs2(args, make_bam=make_bam, dir=FALSE, dir.name=dir.name)$args
      for(i in seq_along(output(args))){
        if(length(output(args)[[i]]) > 1){
            for(j in seq_along(output(args)[[i]])){
            for(k in seq_along(output(args)[[i]][[j]])){
              if(file.exists(output(args)[[i]][[j]])){
                name <- strsplit(output(args)[[i]][[j]][[k]], split="\\/")[[1]]
                name <- name[length(name)]
                file.rename(from=output(args)[[i]][[j]][[k]], to=file.path(logdir, dir.name, name))
                #outputList_new <- c(outputList_new, file.path(logdir, dir.name, name))
                # file.rename(from=output(args.return)[[i]][[j]][[k]], to=file.path(logdir, dir.name, names(output(args.return)[i]), name))
                # outputList_new <- c(outputList_new, file.path(logdir, dir.name, names(output(args.return)[i]), name))
              } else if(!file.exists(output(args.return)[[i]][[j]][[k]])){
                dump <- "No such file or directory"
              }
            }
          }
        } else if(length(output(args)[[i]]) == 1){
          for(j in seq_along(output(args)[[i]][[1]])){
            if(file.exists(output(args)[[i]][[1]][[j]])){
              name <- strsplit(output(args)[[i]][[1]][[j]], split="\\/")[[1]]
              name <- name[length(name)]
              file.rename(from=output(args)[[i]][[1]][[j]], to=file.path(logdir, dir.name, name))
              # outputList_new <- c(outputList_new, file.path(logdir, dir.name, name))
              #file.rename(from=output(args.return)[[i]][[1]][[j]], to=file.path(logdir, dir.name, names(output(args.return)[i]), name))
              # outputList_new <- c(outputList_new, file.path(logdir, dir.name, names(output(args.return)[i]), name))
            } else if(!file.exists(output(args)[[i]][[1]][[j]])){
              dump <- "No such file or directory"
            }
          }
        }
      }
      #args.return <- output_update(args.return, dir=TRUE, dir.name=dir.name, replace=FALSE)
    }
    ## double check output file
    cat("\n")
    cat(crayon::blue("---- Summary ----"), "\n")
    print(df.status.f)
    #print(S4Vectors::DataFrame(df.status.f))
    close(pb)
    return(args.return)
  }
}

## Usage:
# WF <- runCommandline(WF, make_bam=TRUE) # creates the files in the ./results folder
# WF <- runCommandline(WF, dir=TRUE) # creates the files in the ./results/workflowName/Samplename folder
# WF <- runCommandline(WF, make_bam = FALSE, dir=TRUE) ## For hisat2-mapping.cwl template
# runid="01"; make_bam=FALSE; del_sam=FALSE; dir=TRUE; dir.name=NULL; force=FALSE
# runid="01"; make_bam=FALSE; del_sam=FALSE; dir=TRUE; dir.name=NULL; force=TRUE
# runid="01"; make_bam=FALSE; del_sam=FALSE; dir=TRUE; dir.name=NULL; force=TRUE

.tryRunC <- function(command, commandargs){
  warning <- error <- NULL
  value <- withCallingHandlers(
    tryCatch(
      if(command %in% "bwa") {
        bwa_err <- tempfile()
        system2(command, args=commandargs, stdout = T, stderr = bwa_err)
        readLines(bwa_err)
      } else if(command %in% c("bash")) {
        system(paste(command, commandargs))
      } else if(isTRUE(grep('\\$', command)==1)) {
        system(paste(command, commandargs))
      } else {
        system2(command, args=commandargs, stdout=TRUE, stderr=TRUE)
      },
      error = function(e) {
        error <<- conditionMessage(e)
        NULL
        }),
    warning = function(w) {
      warning <<- append(warning, conditionMessage(w))
      invokeRestart("muffleWarning")
      })
  list(stdout = value, warning = warning, error = error)
}

## Usage:
# .tryRunC("ls", "la")


###########################################################################
## .moduleload function: internal function to module load <modules> ##
###########################################################################
.moduleload <- function(args){
  ## Environment Modules section
  if(any(nchar(gsub(" {1,}", "", modules(args))) > 0)) {
    ## Check if "Environment Modules" is in the PATH
    try(suppressWarnings(modulecmd_path <- system("which modulecmd", intern=TRUE, ignore.stderr=TRUE)),
        silent=TRUE)
    ## "Environment Modules" is not available
    if(length(modulecmd_path) == 0 ) {
      dump <- "do nothing"
      ## "Environment Modules" is available and proceed the module load
    } else if (length(modulecmd_path) > 0) {
      for(j in modules(args)) moduleload(j) # loads specified software from module system
    }
  }
}


###########################################################################
## .makeBam function: internal function to convert *.sam to *.bam file ##
###########################################################################
.makeBam <- function(output_args, make_bam=TRUE, del_sam=TRUE){
  if(all(!is.logical(c(make_bam, del_sam)))) stop("Arguments needs to be assigned 'TRUE' and 'FALSE'")
  if(make_bam==TRUE) {
    sam_files <- grepl(".sam$", output_args)
    others_files <- grepl("vcf$|bcf$|xls$|bed$", output_args)
    completed.bam <- grepl(".bam$", output_args)
    if(any(sam_files)){
      for(k in which(sam_files)){
        if(file.exists(output_args[k])){
          Rsamtools::asBam(file=output_args[k], destination=gsub("\\.sam$", "", output_args[k]), overwrite=TRUE, indexDestination=TRUE)
          if(del_sam==TRUE){
            unlink(output_args[k])
          } else if(del_sam==FALSE){
            dump <- "do nothing"
          }
        } else next ()
      } } else if(any(others_files)){
        dump <- "do nothing"
      }
    if(any(completed.bam)){ # If output is unindexed *.bam file (e.g. Tophat2)
      for(k in which(completed.bam)){
        Rsamtools::sortBam(file=output_args[k], destination=gsub("\\.bam$", "", output_args[k]))
        Rsamtools::indexBam(output_args[k])
      }
    }
  } else if(make_bam==FALSE){
    dump <- "do nothing"
  }
}
## Usage:
# .makeBam(output(args)[[1]][[1]], make_bam=TRUE, del_sam=FALSE)

########################################################
## .checkOutArgs2 function: internal function to check
## if the expected output has been created  ##
########################################################
.checkOutArgs2 <- function(args, make_bam, dir=dir, dir.name, force){
  args_complete <- args
  if(make_bam==TRUE) {
    ## Validation for Hisat2
    if(any(grepl("samtools", names(clt(args_complete))))){ stop("argument 'make_bam' should be 'FALSE' when using the workflow with 'SAMtools'")}
    args_complete <- output_update(args_complete, dir=dir, dir.name=dir.name, replace=TRUE, extension=c(".sam", ".bam"), make_bam=make_bam)
  }
  if(dir==TRUE) {
    args_complete <- output_update(args_complete, dir=dir, dir.name=dir.name)
  }
  completed <- output(args_complete)
  for(i in seq_along(output(args_complete))){
    for(j in seq_along(output(args_complete)[[i]])){
      completed[[i]][[j]] <- file.exists(output(args_complete)[[i]][[j]])
      names(completed[[i]][[j]]) <- output(args_complete)[[i]][[j]]
    }
  }
  check <- check.output(args_complete, "list")
  if(force) check <- replace(check, check==TRUE, FALSE)
  for(i in seq_along(check)){
    if(check[i]==FALSE){
      args[["output"]][[i]] <- args$internal_outfiles[[i]]
    }
  }
  return <- list(args= args, args_complete=args_complete, completed=completed)
  return(return)
}
## Usage:
# return <- .checkOutArgs2(args, make_bam=TRUE)
# args <- return$args
# completed <- return$completed

##########################################################################################
## .sysargsrunCommandline function: Old version of runCommandline accepts SYSargs class ##
##########################################################################################
.sysargsrunCommandline <- function(args, runid="01", make_bam=TRUE, del_sam=TRUE, ...) {
  message("This method is using the previus version of SYSargs workflow control modules.
          Please check the new version 'SYSargs2'.")
  commands <- sysargs(args)
  completed <- file.exists(outpaths(args))
  names(completed) <- outpaths(args)
  logdir <- results(args)
  for(i in seq(along=commands)) {
    ## Run alignmets only for samples for which no BAM file is available.
    if(as.logical(completed)[i]) {
      next()
    } else {
      ## Create soubmitargsID_command file
      cat(commands[i], file=paste(logdir, "submitargs", runid, sep=""), sep = "\n", append=TRUE)
      ## Run executable
      command <- gsub(" .*", "", as.character(commands[i]))
      commandargs <- gsub("^.*? ", "",as.character(commands[i]))
      ## Execute system command; note: BWA needs special treatment in stderr handling since it writes
      ## some stderr messages to sam file if used with system2()
      if(software(args) %in% c("bwa aln", "bwa mem")) {
        stdout <- system2(command, args=commandargs, stdout=TRUE, stderr=FALSE)
      } else if(software(args) %in% c("bash_commands")) {
        stdout <- system(paste(command, commandargs))
      } else {
        stdout <- system2(command, args=commandargs, stdout=TRUE, stderr=TRUE)
      }
      ## Create submitargsID_stdout file
      cat(commands[i], file=paste(logdir, "submitargs", runid, "_log", sep=""), sep = "\n", append=TRUE)
      cat(unlist(stdout), file=paste(logdir, "submitargs", runid, "_log", sep=""), sep = "\n", append=TRUE)
      ## Conditional postprocessing of results
      if(make_bam==TRUE) {
        if(grepl(".sam$", outfile1(args)[i])) { # If output is *.sam file (e.g. Bowtie2)
          asBam(file=outfile1(args)[i], destination=gsub("\\.sam$", "", outfile1(args)[i]), overwrite=TRUE, indexDestination=TRUE)
          if(del_sam==TRUE){
            unlink(outfile1(args)[i])
          } else if(del_sam==FALSE){
            dump <- "do nothing"
          }
        } else if(grepl("vcf$|bcf$|xls$|bed$", outpaths(args)[i])) {
          dump <- "do nothing"
        } else { # If output is unindexed *.bam file (e.g. Tophat2)
          sortBam(file=names(completed[i]), destination=gsub("\\.bam$", "", names(completed[i])))
          indexBam(names(completed[i]))
        }
      }
    }
  }
  bamcompleted <- gsub("sam$", "bam$", file.exists(outpaths(args)))
  names(bamcompleted) <- SampleName(args)
  cat("Missing alignment results (bam files):", sum(!as.logical(bamcompleted)), "\n"); cat("Existing alignment results (bam files):", sum(as.logical(bamcompleted)), "\n")
  return(bamcompleted)
}
## Usage:
# args <- systemArgs(sysma="param/hisat2.param", mytargets="targets.txt")
# sysargs(args)[1] # Command-line parameters for first FASTQ file
# system("hisat2-build ./data/tair10.fasta ./data/tair10.fasta")
# .sysargsrunCommandline(args=args)

############################################################################################
## batchtools-based function to submit runCommandline jobs to queuing system of a cluster ##
############################################################################################
## The advantage of this function is that it should work with most queuing/scheduling systems such as SLURM, Troque, SGE, ...
clusterRun <- function(args, FUN = runCommandline, more.args = list(args = args, make_bam = TRUE), conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl", Njobs, runid = "01", resourceList) {
  ## Validity checks of inputs
  if (any(class(args) != "SYSargs" & class(args) != "SYSargs2")) stop("Argument 'args' needs to be assigned an object of class 'SYSargs' OR 'SYSargs2'")
  if (class(FUN) != "function") stop("Value assigned to 'FUN' argument is not an object of class function.")
  if (!file.exists(conffile)) stop("Need to point under 'conffile' argument to proper config file. See more information here: https://mllg.github.io/batchtools/reference/makeRegistry.html. Note: in this file *.tmpl needs to point to a valid template file.")
  if (!file.exists(template)) stop("Need to point under 'template' argument to proper template file. Sample template files for different schedulers are available here: https://github.com/mllg/batchtools/blob/master/inst/templates/")
  if (!class(more.args) == "list") stop("'more.args' needs to be object of class 'list'.")
  if (any(!names(more.args) %in% names(as.list(formals(FUN))))) stop(paste("The list of arguments assigned to 'more.args' can only be the following arguments defined in the function 'FUN':", paste(names(as.list(formals(FUN))), collapse = ", ")))
  ## SYSargs class
  if (class(args) == "SYSargs") {
    path <- normalizePath(results(args))
    args.f <- seq(along = args)
    ## SYSargs2 class
  } else if (class(args) == "SYSargs2") {
    path <- normalizePath(args$yamlinput$results_path$path)
    args.f <- seq(along = cmdlist(args))
  }
  ## batchtools routines
  f <- function(i, args, ...) FUN(args = args[i], ...)
  logdir1 <- paste0(path, "/submitargs", runid, "_btdb_", paste(sample(0:9, 4), collapse = ""))
  reg <- makeRegistry(file.dir = logdir1, conf.file = conffile, packages = "systemPipeR")
  ids <- batchMap(fun = f, args.f, more.args = more.args, reg = reg)
  chunk <- chunk(ids$job.id, n.chunks = Njobs, shuffle = FALSE)
  ids$chunk <- chunk
  done <- submitJobs(ids = ids, reg = reg, resources = resourceList)
  return(reg)
}
## Usage:
# resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024)
# reg <- clusterRun(args, conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl", Njobs=18, runid="01", resourceList=resources)
# getStatus(reg=reg)
# waitForJobs(reg=reg)

########################
## Read preprocessing ##
########################
preprocessReads <- function(args, Fct, batchsize = 100000, overwrite = TRUE, ...) {
  if (all(class(args) != "SYSargs" & class(args) != "SYSargs2")) stop("Argument 'args' needs to be assigned an object of class 'SYSargs' OR 'SYSargs2")
  if (class(Fct) != "character") stop("Argument 'Fct' needs to be of class character")
  if (class(args) == "SYSargs") {
    colnames_args <- colnames(targetsout(args)) # SYSargs
    outpaths <- outpaths(args) # SYSargs
    targets_in <- targetsin(args)
  } else if (class(args) == "SYSargs2") {
    colnames_args <- colnames(targets.as.df(args$targets)) # SYSargs2
    outpaths <- subsetWF(args = args, slot = "output", subset = 1, index = 1)
    targets_in <- targets.as.df(args$targets)
  }
  ## Run function in loop over all fastq files
  ## Single end fastq files
  if (!all(c("FileName1", "FileName2") %in% colnames_args)) {
    for (i in seq(along = args)) {
      outfile <- outpaths[i]
      ## Delete existing fastq files with same names, since writeFastq will append to them
      if (overwrite == TRUE) {
        if (any(file.exists(outfile))) unlink(outfile)
      } else {
        if (any(file.exists(outfile))) stop(paste("File", outfile, "exists. Please delete file first or set overwrite=TRUE."))
      }
      ## Run preprocessor function with FastqStreamer
      counter <- 0
      f <- FastqStreamer(infile1(args)[i], batchsize)
      while (length(fq <- yield(f))) {
        fqtrim <- eval(parse(text = Fct))
        writeFastq(fqtrim, outfile, mode = "a", ...)
        counter <- counter + length(fqtrim)
        cat(counter, "processed reads written to file:", outfile, "\n")
      }
      close(f)
    }
  }
  ## Paired end fastq files
  if (all(c("FileName1", "FileName2") %in% colnames_args)) {
    for (i in seq(along = args)) {
      p1 <- as.character(targets_in$FileName1[i])
      p2 <- as.character(targets_in$FileName2[i])
      if (class(args) == "SYSargs") {
        p1out <- as.character(targetsout(args)$FileName1[i])
        p2out <- as.character(targetsout(args)$FileName2[i])
      } else if (class(args) == "SYSargs2") {
        p1out <- args$output[[i]][[1]][[1]]
        p2out <- args$output[[i]][[1]][[2]]
      }
      ## Delete existing fastq files with same names, since writeFastq will append to them
      if (overwrite == TRUE) {
        if (any(file.exists(p1out))) unlink(p1out)
        if (any(file.exists(p2out))) unlink(p2out)
      } else {
        if (any(file.exists(p1out))) stop(paste("File", p1out, "exists. Please delete file first or set overwrite=TRUE."))
        if (any(file.exists(p2out))) stop(paste("File", p2out, "exists. Please delete file first or set overwrite=TRUE."))
      }
      ## Run preprocessor function with FastqStreamer
      counter1 <- 0
      counter2 <- 0
      f1 <- FastqStreamer(p1, batchsize)
      f2 <- FastqStreamer(p2, batchsize)
      while (length(fq1 <- yield(f1))) {
        fq2 <- yield(f2)
        if (length(fq1) != length(fq2)) stop("Paired end files cannot have different read numbers.")
        ## Process p1
        fq <- fq1 # for simplicity in eval
        fq1trim <- eval(parse(text = Fct))
        ## Index for p1
        index1 <- as.character(id(fq1)) %in% as.character(id(fq1trim))
        names(index1) <- seq(along = index1)
        index1 <- names(index1[index1])
        ## Process p2
        fq <- fq2 # for simplicity in eval
        fq2trim <- eval(parse(text = Fct))
        ## Index for p1
        index2 <- as.character(id(fq2)) %in% as.character(id(fq2trim))
        names(index2) <- seq(along = index2)
        index2 <- names(index2[index2])
        ## Export to processed paired files
        indexpair1 <- index1 %in% index2
        writeFastq(fq1trim[indexpair1], p1out, mode = "a", ...)
        indexpair2 <- index2 %in% index1
        writeFastq(fq2trim[indexpair2], p2out, mode = "a", ...)
        counter1 <- counter1 + sum(indexpair1)
        cat(counter1, "processed reads written to file:", p1out, "\n")
        counter2 <- counter2 + sum(indexpair2)
        cat(counter2, "processed reads written to file:", p2out, "\n")
      }
      close(f1)
      close(f2)
    }
  }
}
## Usage:
# preprocessReads(args=args, Fct="trimLRPatterns(Rpattern="GCCCGGGTAA", subject=fq)", batchsize=100000, overwrite=TRUE, compress=TRUE)

##################################################################
## Function to create sym links to bam files for viewing in IGV ##
##################################################################
symLink2bam <- function(sysargs, command = "ln -s", htmldir, ext = c(".bam", ".bai"), urlbase, urlfile) {
  ## Create URL file
  if(all(is(sysargs, "SYSargs") & is(sysargs, "SYSargs2"))) stop("Argument 'sysargs' needs to be assigned an object of class 'SYSargs' OR 'SYSargs2")
  ## SYSargs class
  if(is(sysargs, "SYSargs")) {
    bampaths <- outpaths(sysargs)
    symname <- SampleName(sysargs)
    ## SYSargs2 class ##
  } else if (is(sysargs, "SYSargs2")) {
    bampaths <- normalizePath(subsetWF(args = sysargs, slot = "output", subset = 1, index = 1))
    symname <- names(targets(sysargs))
  }
  urls <- paste(urlbase, htmldir[2], symname, ext[1], "\t", symname, sep = "")
  writeLines(urls, urlfile)
  ## Creat correspoding sym links
  dir.create(paste(htmldir, collapse = ""))
  symname <- rep(symname, each = 2)
  symname <- paste(symname, c(ext[1], paste(ext, collapse = "")), sep = "")
  bampaths2 <- as.vector(t(cbind(bampaths, paste(bampaths, ext[2], sep = ""))))
  symcommands <- paste(command, " ", bampaths2, " ", paste(htmldir, collapse = ""), symname, sep = "")
  for (i in symcommands) system(i)
}
## Usage:
# symLink2bam(sysargs=args, command="ln -s", htmldir=c("~/.html/", "somedir/"), ext=c(".bam", ".bai"), urlbase="http://cluster.hpcc.ucr.edu/~tgirke/", urlfile="IGVurl.txt")

#####################
## Alignment Stats ##
#####################
#' @param files named character vector of indexed bam files. Names of this vector will
#' be used as sample names
#' @param out_dir character, directory path to store individual sample stats files
alignStats <- function(files, out_dir="results") {
  if (class(files) %in% c("SYSargs", "SYSargs2")) {
    return(warning('alignStats: New version of SPR no longer accept "SYSargs", "SYSargs2" objects as inputs.\n',
                   'Use `getColumn` to get a vector of paths instead.'))
  }
  stopifnot(is.character(files))
  stopifnot(is.character(out_dir) && length(out_dir) == 1)
  if(is.null(names(files))) stop("Files must be a named vector, names will be used as sample names")
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(out_dir)) stop("Cannot create output directory", out_dir)
  if(!all(check_files <- file.exists(files))) stop("Some files are missing:\n", paste0(files[!check_files], collapse = ",\n"))
  if(!all(check_ext <- stringr::str_detect(files, "\\.bam$"))) stop("alignStats: All files need to end with .bam\n", paste0(files[!check_ext], collapse = ",\n"))
  file_base <- basename(files)
  out_files <- file.path(out_dir, gsub("\\.bam$", "_bam_stats.csv", file_base))
  samplenames <- names(files)
  bam_df <- data.frame(sample= samplenames, files= file_base, mapped = NA, unmapped = NA)
  for (i in seq_along(files)) {
    bam_stats <- Rsamtools::idxstatsBam(files[i])
    message("Write bam stats for ", samplenames[i], " to ", basename(out_files[i]))
    write.csv(bam_stats, out_files[i], row.names = FALSE)
    map_info <- colSums(bam_stats[, c(3, 4)])
    bam_df[i, c(3, 4)] <- map_info
  }
  bam_df$map_rate <- round(bam_df$mapped/(bam_df$mapped + bam_df$unmapped), 3)
  bam_df
}
## Usage:
# read_statsDF <- alignStats(files=files)

########################
## RPKM Normalization ##
########################
returnRPKM <- function(counts, ranges) {
  geneLengthsInKB <- sum(width(reduce(ranges))) / 1000 # Length of exon union per gene in kbp
  millionsMapped <- sum(counts) / 1e+06 # Factor for converting to million of mapped reads.
  rpm <- counts / millionsMapped # RPK: reads per kilobase of exon model.
  rpkm <- rpm / geneLengthsInKB # RPKM: reads per kilobase of exon model per million mapped reads.
  return(rpkm)
}
## Usage:
# countDFrpkm <- apply(countDF, 2, function(x) returnRPKM(counts=x, ranges=eByg))

###############################################
## Read Sample Comparisons from Targets File ##
###############################################
## Parses sample comparisons from <CMP> line(s) in targets.txt file or SYSars object.
## All possible comparisons can be specified with 'CMPset: ALL'.
readComp <- function(file, format = "vector", delim = "-") {
  if (!format %in% c("vector", "matrix")) stop("Argument format can only be assigned: vector or matrix!")
  if (class(file) == "SYSargs") {
    if (length(targetsheader(file)) == 0) stop("Input has no targets header lines.")
    comp <- targetsheader(file)
    ## SYSargs2 class
  } else if (class(file) == "SYSargs2") {
    if (length(targetsheader(file))[[1]] == 0) stop("Input has no targets header lines.")
    comp <- targetsheader(file)[[1]]
  } else {
    comp <- readLines(file)
  }
  comp <- comp[grepl("<CMP>", comp)]
  comp <- gsub("#.*<CMP>| {1,}", "", comp)
  comp <- gsub("\t", "", comp)
  comp <- gsub("^\"|\"$", "", comp) # Often required if Excel is used for editing targets file
  comp <- strsplit(comp, ":|,")
  names(comp) <- lapply(seq(along = comp), function(x) comp[[x]][1])
  comp <- sapply(names(comp), function(x) comp[[x]][-1], simplify = FALSE)
  ## Check whether all samples are present in Factor column of targets file
  checkvalues <- unique(unlist(strsplit(unlist(comp), "-")))
  checkvalues <- checkvalues[checkvalues != "ALL"]
  if (class(file) == "SYSargs") {
    all <- unique(as.character(targetsin(file)$Factor))
  } else if (class(file) == "SYSargs2") {
    all <- unique(as.character(targets.as.df(targets(file))$Factor))
  } else {
    all <- unique(as.character(read.delim(file, comment.char = "#")$Factor))
  }
  if (any(!checkvalues %in% all)) stop(paste("The following samples are not present in Factor column of targets file:", paste(checkvalues[!checkvalues %in% all], collapse = ", ")))
  ## Generate outputs
  allindex <- sapply(names(comp), function(x) any(grepl("ALL", comp[[x]])))
  if (any(allindex)) for (i in which(allindex)) comp[[i]] <- combn(all, m = 2, FUN = paste, collapse = delim)
  if (format == "vector" & delim != "-") comp <- sapply(names(comp), function(x) gsub("-", delim, comp[[x]]), simplify = FALSE)
  if (format == "vector") {
    return(comp)
  }
  if (format == "matrix") {
    return(sapply(names(comp), function(x) do.call("rbind", strsplit(comp[[x]], "-")), simplify = FALSE))
  }
}
## Usage:
# targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
# cmp <- readComp(targetspath, format="vector", delim="-")
# cmp <- readComp(WF, format="vector", delim="-")

#######################################################################
## Run edgeR GLM with entire count matrix or subsetted by comparison ##
#######################################################################
## If independent=TRUE then countDF will be subsetted for each comparison
run_edgeR <- function(countDF, targets, cmp, independent = TRUE, paired = NULL, mdsplot = "") {
  ## if(class(cmp) != "matrix" & length(cmp)==2) cmp <- t(as.matrix(cmp)) # If cmp is vector of length 2, convert it to matrix.
  ## fix for _R_CHECK_LENGTH_1_LOGIC2_ error: " --- failure: the condition has length > 1 ---"
  if (all(class(cmp) != "matrix" & length(cmp) == 2)) cmp <- t(as.matrix(cmp))
  samples <- as.character(targets$Factor)
  names(samples) <- paste(as.character(targets$SampleName), "", sep = "")
  countDF <- countDF[, names(samples)]
  countDF[is.na(countDF)] <- 0
  edgeDF <- data.frame(row.names = rownames(countDF))
  group <- as.character(samples)
  if (independent == TRUE) {
    loopv <- seq(along = cmp[, 1])
  } else {
    loopv <- 1
  }
  for (j in loopv) {
    ## Filtering and normalization
    y <- DGEList(counts = countDF, group = group) # Constructs DGEList object
    if (independent == TRUE) {
      subset <- samples[samples %in% cmp[j, ]]
      y <- y[, names(subset)]
      y$samples$group <- factor(as.character(y$samples$group))
    }
    keep <- rowSums(edgeR::cpm(y) > 1) >= 2
    y <- y[keep, ]
    y <- calcNormFactors(y)
    ## Design matrix
    if (length(paired) == 0) {
      design <- model.matrix(~ 0 + y$samples$group, data = y$samples)
      colnames(design) <- levels(y$samples$group)
    } else {
      if (length(paired) > 0 & independent == FALSE) stop("When providing values under 'paired' also set independent=TRUE")
      Subject <- factor(paired[samples %in% cmp[j, ]]) # corrected Jun 2014 (won't change results)
      Treat <- y$samples$group
      design <- model.matrix(~ Subject + Treat)
      levels(design) <- levels(y$samples$group)
    }
    ## Estimate dispersion
    y <- estimateGLMCommonDisp(y, design, verbose = TRUE) # Estimates common dispersions
    y <- estimateGLMTrendedDisp(y, design) # Estimates trended dispersions
    y <- estimateGLMTagwiseDisp(y, design) # Estimates tagwise dispersions
    fit <- glmFit(y, design) # Fits the negative binomial GLM for each tag and produces an object of class DGEGLM with some new components.
    ## Contrast matrix is optional but makes anlysis more transparent
    if (independent == TRUE) {
      mycomp <- paste(cmp[j, 1], cmp[j, 2], sep = "-")
    } else {
      mycomp <- paste(cmp[, 1], cmp[, 2], sep = "-")
    }
    if (length(paired) == 0) contrasts <- makeContrasts(contrasts = mycomp, levels = design)
    for (i in seq(along = mycomp)) {
      if (length(paired) == 0) {
        lrt <- glmLRT(fit, contrast = contrasts[, i]) # Takes DGEGLM object and carries out the likelihood ratio test.
      } else {
        lrt <- glmLRT(fit) # No contrast matrix with paired design
      }
      deg <- as.data.frame(topTags(lrt, n = length(rownames(y))))
      colnames(deg) <- paste(paste(mycomp[i], collapse = "_"), colnames(deg), sep = "_")
      edgeDF <- cbind(edgeDF, deg[rownames(edgeDF), ])
    }
    if (nchar(mdsplot) > 0) {
      pdf(paste("./results/sample_MDS_", paste(unique(subset), collapse = "-"), ".pdf", sep = ""))
      plotMDS(y)
      dev.off()
    }
  }
  return(edgeDF)
}
## Usage:
# cmp <- readComp(file=targetspath, format="matrix", delim="-")
# edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=TRUE, mdsplot="")

####################################################################
## Run DESeq2 with entire count matrix or subsetted by comparison ##
####################################################################
## If independent=TRUE then countDF will be subsetted for each comparison
run_DESeq2 <- function(countDF, targets, cmp, independent = FALSE, lfcShrink=FALSE, type="normal") {
  ## if(class(cmp) != "matrix" & length(cmp)==2) cmp <- t(as.matrix(cmp)) # If cmp is vector of length 2, convert it to matrix.
  ## fix for _R_CHECK_LENGTH_1_LOGIC2_ error: " --- failure: the condition has length > 1 ---"
  if (all(class(cmp) != "matrix" & length(cmp) == 2)) cmp <- t(as.matrix(cmp))
  samples <- as.character(targets$Factor)
  names(samples) <- paste(as.character(targets$SampleName), "", sep = "")
  countDF <- countDF[, names(samples)]
  countDF[is.na(countDF)] <- 0
  deseqDF <- data.frame(row.names = rownames(countDF))
  if (independent == TRUE) {
    loopv <- seq(along = cmp[, 1])
  } else {
    loopv <- 1
  }
  for (j in loopv) {
    if (independent == TRUE) {
      ## Create subsetted DESeqDataSet object
      subset <- samples[samples %in% cmp[j, ]]
      countDFsub <- countDF[, names(subset)]
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(countDFsub), colData = data.frame(condition = subset), design = ~condition)
      mycmp <- cmp[j, , drop = FALSE]
    } else {
      ## Create full DESeqDataSet object
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(countDF), colData = data.frame(condition = samples), design = ~condition)
      mycmp <- cmp
    }
    ## Estimate of (i) size factors, (ii) dispersion, (iii) negative binomial GLM fitting and (iv) Wald statistics
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    for (i in seq(along = mycmp[, 1])) {
      ## Extracts DEG results for specific contrasts from DESeqDataSet object
      res <- DESeq2::results(dds, contrast = c("condition", mycmp[i, ]))

      ## lfcShrink
      if (lfcShrink == FALSE) {
        res <- DESeq2::results(dds, contrast = c("condition", mycmp[i, ]))
      } else if (lfcShrink == TRUE) {
        suppressMessages(
        res <- DESeq2::lfcShrink(dds, contrast = c("condition", mycmp[i, ]), type=type)
        )
      }
      ## Set NAs to reasonable values to avoid errors in downstream filtering steps
      res[is.na(res[, "padj"]), "padj"] <- 1
      res[is.na(res[, "log2FoldChange"]), "log2FoldChange"] <- 0
      deg <- as.data.frame(res)
      colnames(deg)[colnames(deg) %in% c("log2FoldChange", "padj")] <- c("logFC", "FDR")
      colnames(deg) <- paste(paste(mycmp[i, ], collapse = "-"), colnames(deg), sep = "_")
      deseqDF <- cbind(deseqDF, deg[rownames(deseqDF), ])
    }
  }
  return(deseqDF)
}
## Usage:
# cmp <- readComp(file=targetspath, format="matrix", delim="-")
# degseqDF <- run_DESeq2(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=FALSE)

###########################################
## Filter DEGs by p-value and fold change ##
############################################
filterDEGs <- function(degDF, filter, plot = TRUE) {
  ## global functions or variables
  Comparisons <- Counts <- Type <- NULL
  pval <- degDF[, grep("_FDR$", colnames(degDF)), drop = FALSE]
  log2FC <- degDF[, grep("_logFC$", colnames(degDF)), drop = FALSE]
  ## DEGs that are up or down regulated
  pf <- pval <= filter["FDR"] / 100 & (log2FC >= log2(filter["Fold"]) | log2FC <= -log2(filter["Fold"]))
  colnames(pf) <- gsub("_FDR", "", colnames(pf))
  pf[is.na(pf)] <- FALSE
  DEGlistUPorDOWN <- sapply(colnames(pf), function(x) rownames(pf[pf[, x, drop = FALSE], , drop = FALSE]), simplify = FALSE)
  ## DEGs that are up regulated
  pf <- pval <= filter["FDR"] / 100 & log2FC >= log2(filter["Fold"])
  colnames(pf) <- gsub("_FDR", "", colnames(pf))
  pf[is.na(pf)] <- FALSE
  DEGlistUP <- sapply(colnames(pf), function(x) rownames(pf[pf[, x, drop = FALSE], , drop = FALSE]), simplify = FALSE)
  ## DEGs that are down regulated
  pf <- pval <= filter["FDR"] / 100 & log2FC <= -log2(filter["Fold"])
  colnames(pf) <- gsub("_FDR", "", colnames(pf))
  pf[is.na(pf)] <- FALSE
  DEGlistDOWN <- sapply(colnames(pf), function(x) rownames(pf[pf[, x, drop = FALSE], , drop = FALSE]), simplify = FALSE)
  df <- data.frame(Comparisons = names(DEGlistUPorDOWN), Counts_Up_or_Down = sapply(DEGlistUPorDOWN, length), Counts_Up = sapply(DEGlistUP, length), Counts_Down = sapply(DEGlistDOWN, length))
  resultlist <- list(UporDown = DEGlistUPorDOWN, Up = DEGlistUP, Down = DEGlistDOWN, Summary = df)
  if (plot == TRUE) {
    mytitle <- paste("DEG Counts (", names(filter)[1], ": ", filter[1], " & ", names(filter)[2], ": ", filter[2], "%)", sep = "")
    df_plot <- data.frame(Comparisons = rep(as.character(df$Comparisons), 2), Counts = c(df$Counts_Up, df$Counts_Down), Type = rep(c("Up", "Down"), each = length(df[, 1])))
    p <- ggplot(df_plot, aes(Comparisons, Counts, fill = Type)) +
      geom_bar(position = "stack", stat = "identity") +
      coord_flip() +
      theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
      ggtitle(mytitle)
    print(p)
  }
  return(resultlist)
}
## Usage:
# DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=1))
