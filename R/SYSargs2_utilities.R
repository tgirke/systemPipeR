#################################################################
## Functions to construct SYSargs2 objects and other utilities ##
################################################################

###################
## Load Workflow ##
###################
loadWorkflow <- function(targets=NULL, wf_file, input_file, dir_path=".") {
    if(any(is(wf_file)=="list")) {
        wf <- wf_file
    } else {
        if(!file.exists(file.path(dir_path, wf_file))==TRUE) stop("Provide valid '.cwl' file. Check the file PATH.")
        wf <- yaml::read_yaml(file.path(dir_path, wf_file))
    }
    if(any(is(input_file)=="list")) {
        input <- input_file
    } else {
        if(!file.exists(file.path(dir_path, input_file))==TRUE) stop("Provide valid 'files.'.yml' file. Check the file PATH.")
        input <- yaml::read_yaml(file.path(dir_path, input_file))
    }
  modules <- input$ModulesToLoad
  if(is.null(modules)) modules <- list()
  if(is.null(dir_path)){
      cwlfiles <- list(cwl=NA, yml=NA)
  } else {
      cwlfiles <- list(cwl=normalizePath(file.path(dir_path, wf_file)), yml=normalizePath(file.path(dir_path, input_file)))
  }
  inputvars <- list()
  if(tolower(wf$class) == "workflow") { 
    steps <- names(wf$steps)
    cwlfiles$steps <- steps
    cltpaths <- sapply(seq_along(steps), function(x) normalizePath(file.path(dir_path, wf$steps[[steps[x]]]$run)))
    names(cltpaths) <- strsplit(basename(cltpaths), ".cwl")
    cwlfiles$cltpaths <- cltpaths
    cltlist <- sapply(cltpaths, function(x) yaml::read_yaml(file.path(x)), simplify = FALSE) 
    names(cltlist) <- sapply(seq_along(steps), function(x) wf$steps[[steps[x]]]$run)
    #names(cltlist) <- steps
    cmdlist <- sapply(names(cltlist), function(x) list(NULL))
    myinput <- sapply(names(cltlist), function(x) list(NULL))
    myoutput <- sapply(names(cltlist), function(x) list(NULL))
    cwlfiles$output_names <- sapply(names(cltlist), function(x) names(cltlist[[x]]$outputs))
    WF <- list(modules=modules, wf=wf, clt=cltlist, yamlinput=input, cmdlist=cmdlist,
               input=myinput, output=myoutput, files=cwlfiles, inputvars=inputvars,
               cmdToCwl=list(), status=list())
  } else if(tolower(wf$class) == "commandlinetool") {
    cltlist <- list(wf)
    names(cltlist) <- basename(wf_file)
    cmdlist <- sapply(names(cltlist), function(x) list(NULL))
    myinput <- sapply(names(cltlist), function(x) list(NULL))
    myoutput <- sapply(names(cltlist), function(x) list(NULL))
    cwlfiles$steps <- strsplit(basename(wf_file), ".cwl")[[1]]
    cwlfiles$output_names <-  names(cltlist[[1]]$outputs)
    WF <- list(modules=modules, wf=list(), clt=cltlist, yamlinput=input, cmdlist=cmdlist,
               input=myinput, output=myoutput, files=cwlfiles, inputvars=inputvars, 
               cmdToCwl=list(), status=list())
  } else {
    stop("Class slot in '<wf_file>.cwl' needs to be 'Workflow' or 'CommandLineTool'.")
  }
  if(!is.null(targets)) {
      if(inherits(targets, "SummarizedExperiment")){
          mytargets <- as.data.frame(colData(targets))
          targetsheader <- metadata(targets)
          WF <- c(list(targets=mytargets, targetsheader=targetsheader), WF)
      } else if(class(targets)=="SYSargs2"){
      mytargets <- targets(targets)
      targetsheader <- targetsheader(targets)[[1]]
    } else {
      if(!file.exists(file.path(targets))==TRUE) stop("Provide valid 'targets' file. Check the file PATH.")
      ext <- strsplit(basename(targets), split="\\.")[[1]]
      ext <- ext[[-1]]
      if("txt" %in% ext){
        mytargets <- read.delim(normalizePath(file.path(targets)), comment.char = "#")
        mytargets <- targets.as.list(mytargets)
      } else if( any(c("yml", "yaml") %in% ext)){
        mytargets <- yaml::read_yaml(targets)
      }
      targetsheader <- readLines(normalizePath(file.path(targets)))
      targetsheader <- targetsheader[grepl("^#", targetsheader)] 
      WF$files["targets"] <- file.path(targets)
      WF <- c(list(targets=mytargets, targetsheader=list(targetsheader=targetsheader)), WF)
    }
  } else {
    WF$files["targets"] <- NA
    WF <- c(list(targets=data.frame(), targetsheader=list()), WF)
  }
  return(as(WF, "SYSargs2"))
}

## Wrapper for loadWorkflow: Short and consistent name for the function
loadWF <- loadWorkflow

## Usage:
# targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
# dir_path <- system.file("extdata/cwl/hisat2", package="systemPipeR")
# WF <- loadWF(targets=targetspath, wf_file="hisat2-se/hisat2-mapping-se.cwl",
#              input_file="hisat2-se/hisat2-mapping-se.yml", dir_path=dir_path)

###################################################
##   Create CommandLineTools from Command-line   ##
###################################################
createWF <- function(targets=NULL, commandLine, results_path="./results", module_load="baseCommand", file = "default", 
                     overwrite = FALSE, cwlVersion = "v1.0", class = "CommandLineTool", writeout=FALSE, silent=FALSE){
  ## TODO if is not default, module load and file
  if(all(!is(commandLine)=="list")) stop("'commandLine' needs to be object of class 'list'.")  
  if(any(!c("baseCommand", "inputs", "outputs") %in% names(commandLine))) stop("Argument 'commandLine' needs to be assigned at least to: 'baseCommand', 'input' or 'output'.")
  if(all(!c("Workflow", "CommandLineTool") %in% class)) stop("Class slot in '<wf_file>.cwl' needs to be 'Workflow' or 'CommandLineTool'.")
  if(dir.exists(results_path)==FALSE) dir.create(path=results_path)
  ## module_load 
  ## 1.  module_load="baseCommand" will use the name of the software 
  ## 2.  module_load=c("ncbi-blast/2.2.30+", "hisat2/2.1.0") will use the specific version and names
  if(module_load == "baseCommand"){
    module_load <- commandLine$baseCommand[[1]]
  } else {
    module_load <- module_load
  }
  ## File Path
  ## 1. file="default"
  ## 2. file = c("test.cwl", "test.yml")
  if("default" %in% file){
    if(dir.exists(paste("param/cwl/", commandLine$baseCommand, sep=""))==FALSE) dir.create(path=paste("param/cwl/", commandLine$baseCommand, sep=""), recursive = TRUE)
    file.cwl <- paste("param/cwl/", commandLine$baseCommand, "/", commandLine$baseCommand, ".cwl", sep="")
    file.yml <- paste("param/cwl/", commandLine$baseCommand, "/", commandLine$baseCommand, ".yml", sep="")
  } else {
    for(i in seq_along(file)){
      extension <- sub('.*\\.', '', file[[i]])
      if(!c("cwl") %in% extension & !c("yml") %in% extension) stop ("Argument 'file' needs to be assigned as a character vector with the names of the two param file. For example, 'test.cwl' and 'test.yml'.")
      if(c("yml") %in% extension){
        file.yml <- file[[i]]
      } else if(c("cwl") %in% extension){
        file.cwl <- file[[i]]
      }
    }
  }
  if(file.exists(file.cwl) & overwrite == FALSE) 
    stop(paste("I am not allowed to overwrite files; please delete existing file:", 
               file, "or set 'overwrite=TRUE', or provide a different name in the 'file' argument"))
  if(file.exists(file.yml) & overwrite == FALSE) 
    stop(paste("I am not allowed to overwrite files; please delete existing file:", 
               file, "or set 'overwrite=TRUE', or provide a different name in the 'file' argument"))
  ## class("CommandLineTool", "Workflow")
  # WF.temp <- SYScreate("SYSargs2")
  WF.temp <- as(SYScreate("SYSargs2"), "list")
  WF.temp$wf <- list()
  WF.temp$clt <- write.clt(commandLine, cwlVersion, class, file.cwl, writeout= writeout, silent = silent) 
  WF.temp$yamlinput <- write.yml(commandLine, file.yml, results_path, module_load, writeout=writeout, silent = silent) 
  WF.temp$modules <- WF.temp$yamlinput$ModulesToLoad
  WF.temp$cmdlist <- sapply(names(WF.temp$clt), function(x) list(NULL))
  WF.temp$input <- sapply(names(WF.temp$clt), function(x) list(NULL))
  WF.temp$output <- sapply(names(WF.temp$clt), function(x) list(NULL))
  WF.temp$files <- list(cwl=file.path(file.cwl), 
                           yml=file.path(file.yml), 
                           steps=names(WF.temp$clt))
  WF.temp$cmdToCwl <- commandLine
  ## targets
  if(!is.null(targets)) {
    mytargets <- read.delim(normalizePath(file.path(targets)), comment.char = "#")
    mytargets <- targets.as.list(mytargets)
    targetsheader <- readLines(normalizePath(file.path(targets)))
    targetsheader <- targetsheader[grepl("^#", targetsheader)]
    WF.temp$files["targets"] <- normalizePath(file.path(targets))
    WF.temp <- c(list(targets=mytargets, targetsheader=list(targetsheader=targetsheader)), WF.temp)
  } else {
    WF.temp$files["targets"] <- NA
    WF.temp <- c(list(targets=data.frame(), targetsheader=list()), WF.temp)
  }
  WF.temp <- as(WF.temp, "SYSargs2")
 # WF.temp[["status"]] <- .statusPending(WF.temp)
  return(WF.temp)
}

## Usage: 
## Example commandLine
# "hisat2 -S ./results/_SampleName_.sam  -x ./data/tair10.fasta  -k 1  --min-intronlen 30  --max-intronlen 3000 --threads 4 -U _FASTQ_PATH1_"
## Provide a list 
# baseCommand <- list(baseCommand="hisat2")
# inputs <- list(
#   "S"=list(type="File", preF="-S", yml="./results/_SampleName_.sam"),
#   "x"=list(type="File", preF="-x", yml="./data/tair10.fasta"),
#   "k"= list(type="int", preF="-k", yml= 1L),
#   "threads"= list(type="int", preF="-threads", yml=4L),
#   "min-intronlen"= list(type="int", preF="--min-intronlen", yml= 30L),
#   "max-intronlen"= list(type="int", preF="--max-intronlen", yml=3000L),
#   "U"=list(type="File", preF="-U", yml="./data/_FASTQ_PATH1_") )
# outputs <- list("hisat2_sam"=list(type="File", "./results/_SampleName_.sam"))
# commandLine <- list(baseCommand=baseCommand, inputs=inputs, outputs=outputs)
# 
## Creates a SYSargs2 object, populate all the command-line, and creates *.cwl and *.yml files
# targets <- system.file("extdata", "targets.txt", package="systemPipeR")
# WF <- createWF(targets=targets, commandLine, results_path="./results", module_load="baseCommand",
#                file = "default", overwrite = FALSE, cwlVersion = "v1.0",
#                class = "CommandLineTool")
# WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))

###############################################
## Render WF for all samples in targets slot ##
###############################################
renderWF <- function(WF, inputvars=NULL) {
  if(any(length(cmdlist(WF)[[1]])!=0)) stop("Argument 'WF' needs to be assigned an object of class 'SYSargs2' and an object created by the 'loadWorkflow' function")  
  ids <- names(targets(WF))
  if(length(ids)==0) ids <- "defaultid"
  bucket <- sapply(ids, function(x) "", simplify=FALSE)
  bucketlist <- list(cmd=bucket, input=bucket, output=bucket)
  if(length(targets(WF))==0){
    if(!is.null(inputvars)){
      
    }
  }
  for(i in ids) {
    tmplist <- .renderWFsingle(WF=WF, id=i, inputvars=inputvars) 
    bucketlist[["cmd"]][[i]] <- tmplist[["cmd"]]
    bucketlist[["input"]][[i]] <- tmplist[["input"]]
    bucketlist[["output"]][[i]] <- tmplist[["output"]]
  }
  WF <- as(WF, "list") 
  WF$cmdlist <- bucketlist$cmd
  WF$input <- bucketlist$input
  WF$output <- bucketlist$output
  WF$inputvars <- tmplist$inputvars
  WF <- as(WF, "SYSargs2")
  WF[["status"]] <- .statusPending(WF)
  return(WF)
}

## Usage:
# WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))

###############################################################
## Update  WF container for all samples in targets slot  ##
###############################################################
updateWF <- function(args, write.yaml=FALSE, name.yaml="default", new_targets=NULL,
                     new_targetsheader=NULL, inputvars=NULL, silent=FALSE){
    if(!inherits(args, "SYSargs2")) stop("args needs to be object of class 'SYSargs2'.")  
    args <- sysargs2(args)
    if(is.null(inputvars)){
        args$inputvars <- inputvars <- args$inputvars
    } else {
        args$inputvars <- inputvars
    }
    if(length(args$cmdToCwl)>1){
        cwlVersion <- args$clt[[1]]$cwlVersion
        class <- args$clt[[1]]$class
        module_load <- args$yamlinput$ModulesToLoad[[1]]
        results_path <- args$yamlinput$results_path$path
        args$clt <- write.clt(args$cmdToCwl, cwlVersion=cwlVersion, class=class, writeout= FALSE, silent = silent) 
        args$yamlinput <- write.yml(args$cmdToCwl, results_path=results_path, module_load=module_load, writeout=FALSE, silent = silent) 
    } else if (length(args$cmdToCwl)==0){
        ## targets
        if(!is.null(new_targets)){
            args$targets <- new_targets
        } else{
            args$targets <- args$targets
        }
        ## targetsheader
        if(!is.null(new_targetsheader)){
            args$targetsheader <- new_targetsheader
        } else{
            args$targetsheader <- args$targetsheader
        }
        
        args$yamlinput <- args$yamlinput
        args$clt <- args$clt
        results_path <- args$yamlinput$results_path$path
    }
    args$input <- args$output <- args$cmdlist <- sapply(names(args$clt), function(x) list(NULL))
    ## write the new yaml
    if(write.yaml==TRUE){
        if(name.yaml=="default"){
            path <- .getPath(args$files$yml)
            name <- paste0(.getFileName(args$files$yml), 
                           format(Sys.time(), "%b%d%Y_%H%M"), ".yml")
            name.yaml <- file.path(path, name)
        } else{
            name.yaml <- name.yaml
        }
        yaml::write_yaml(args$yamlinput, name.yaml)
        if (silent != TRUE) 
            cat("\t", "Written content of 'yamlinput(x)' to file:",  "\n",
                name.yaml, "\n")
        args$files$yml <- name.yaml
    }
    args <- as(args, "SYSargs2")
    args <- renderWF(args, inputvars=inputvars)
    ## Update status slot
    args[["status"]] <- .statusPending(args)
    return(args)
}

## Usage:
# yamlinput(WF)$thread <- 6
# WF <- updateWF(args, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))


###############################################################
## Subsetting the input and output slots by name or position ##
###############################################################
subsetWF <- function(args, slot, subset=NULL, index=NULL, delete=FALSE){
  ## Check the class and slot
  if(!class(args)=="SYSargs2") stop("args needs to be object of class 'SYSargs2'.")  
  if(all(!c("input", "output", "step") %in% slot)) stop("Argument slot can only be assigned one of: 'input' or 'output' or 'step'.")
  
  ## slot input
  if(slot %in% "input"){
    ## Check the subset
    if(all(!is.null(subset) & is.character(subset) & !any(names(inputvars(args)) %in% subset))) stop(paste("For the", slot, "slot, can only be assigned one of the following values in the subset argument:", paste(names(inputvars(args)), collapse=", "), "OR the corresponding position OR NULL")) 
    if(all(!is.null(subset) & is.numeric(subset) & !any(seq_along(names(inputvars(args))) %in% subset))) stop(paste("For the", slot, "slot, can only be assigned one of the following position in the subset argument:", paste(seq_along(names(inputvars(args))), collapse=", "), "OR the names OR NULL")) 
    subset_input <- input(args)
    subset_sample <- sapply(names(subset_input), function(x) list(NULL))
    if(!is.null(subset)) {
      for(i in seq_along(names(subset_input)))
        subset_sample[[i]] <- subset_input[[i]][[subset]]
    } else {
      subset_sample <- subset_input
    }
  } 
  ## slot output
  if(slot %in% "output"){
    ## Check the subset
    if(all(!is.null(subset) & is.character(subset) & !any(names(args$clt) %in% subset))) stop(paste("For the", slot, "slot, can only be assigned one of the following values in the subset argument:", paste(names(args$clt), collapse=", "), "OR the corresponding position OR NULL")) 
    if(all(!is.null(subset) & is.numeric(subset) & !any(seq_along(names(args$clt)) %in% subset))) stop(paste("For the", slot, "slot, can only be assigned one of the following position in the subset argument:", paste(seq_along(names(args$clt)), collapse=", "), "OR the names OR NULL")) 
    if(!is.null(subset)){
      if(!any(seq_along(output(args)[[1]][[subset]]) %in% index)) stop(paste("For the 'index' argument, can only be assigned one of the following position:", paste(seq_along(output(args)[[1]][[subset]]), collapse=", ")))
    }
    
    subset_output <- output(args)
    subset_sample <- as.character()
    if(all(!is.null(subset) & !is.null(index))) {
      for(i in seq_along(names(subset_output)))
        # subset_sample[[i]] <- subset_output[[i]][[subset]]
        subset_sample <- c(subset_sample, subset_output[[i]][[subset]][index])
    } else {
      subset_sample <- subset_output
    }
    names(subset_sample) <- rep(names(subset_output), each=length(subset_sample[1]))
  }
  ## slot step
  if(slot %in% "step"){
    ## Check the subset
    if(all(!is.null(subset) & is.character(subset) & !any(names(args$clt) %in% subset))) stop(paste("For the", slot, "slot, can only be assigned one of the following values in the subset argument:", paste(names(args$clt), collapse=", "), "OR the corresponding position OR NULL")) 
    if(all(!is.null(subset) & is.numeric(subset) & !any(seq_along(names(args$clt)) %in% subset))) stop(paste("For the", slot, "slot, can only be assigned one of the following position in the subset argument:", paste(seq_along(names(args$clt)), collapse=", "), "OR the names OR NULL")) 
    subset_step <- cmdlist(args)
    subset_sample <- sapply(names(subset_step), function(x) list(NULL))
    if(!is.null(subset)) {
      for(i in seq_along(names(subset_sample)))
        subset_sample[[i]] <- subset_step[[i]][[subset]]
    } else {
      subset_sample <- subset_step
    }
  }
  ## IF subset=NULL returns a list 
  if(!is.null(subset)){
    names <- names(subset_sample)
    subset_sample <- as.character(subset_sample)
    names(subset_sample) <- names 
  }
  ## delete == TRUE
  if(delete==TRUE){
    ## delete option only works if the subset is define
    if(!is.character(subset_sample)) { 
      stop(paste("Please define the 'subset' to be deleted in the subset argument")) 
    }
    if(all(file.exists(subset_sample))){
      del <- file.remove(subset_sample)
      cat("The following files has been deleted:", paste(subset_sample, collapse = ", "), "\n")
    } else if(all(!file.exists(subset_sample))) {
      cat("The subset cannot be deleted: no such file ", "\n")
    }
  }
  return(subset_sample)
}

## Usage:
# subsetWF(WF, slot="input", subset='FileName')
# subsetWF(WF, slot="step", subset=1)
# subsetWF(WF, slot="output", subset=1, index=1)
# subsetWF(WF, slot="output", subset=1, index=1, delete=TRUE) ## in order to delete the subset files list

###############################################################
## Subsetting the input and output slots by name or position ##
###############################################################
check.output <- function(args, type="data.frame"){
    ## Check the class and slot
    if(inherits(args, c("SYSargs2"))){
        return(.check.output.sysargs2(args, type=type))
    } else if(inherits(args, c("SYSargsList"))){
        steps <- sapply(names(stepsWF(args)), function(x) list(NULL))
        for(i in seq_along(stepsWF(args))){
            steps[[i]] <- .check.output.sysargs2(stepsWF(args)[[i]], type=type)
        }
        return(steps)
    } else {
        stop("args needs to be object of class 'SYSargs2' or 'SYSargsList'.")   
    }
}

## Usage:
# check.output(WF)

## .check.output.sysargs2
.check.output.sysargs2 <- function(args, type){
    if(type=="data.frame"){
        targets <- sapply(names(output(args)), function(x) list(NULL))
        for(i in seq_along(output(args))){
            targets[[i]][['Total']] <- length(unlist(output(args)[[i]]))
            targets[[i]][['Existing']] <- sum(file.exists(unlist(output(args)[[i]])))
            targets[[i]][['Missing']] <- length(unlist(output(args)[[i]])) - sum(file.exists(unlist(output(args)[[i]])))
            #targets[[i]][['Status']] <- ifelse(targets[[i]]$Missing > 0, "Missing", "Completed")
        }
        targets <- data.frame(matrix(unlist(targets), nrow=length(targets), byrow=T))
        targets <- as.data.frame(cbind(Targets=names(output(args)),  Total_Files = targets$X1,
                                       Existing_Files = targets$X2, Missing_Files=targets$X3 #, status=targets$X4
                                       ))
        
        return(targets)
    } else if(type=="list"){
        targets <- sapply(names(output(args)), function(x) list(NULL))
        for(i in seq_along(output(args))){
            targets[[i]] <- all(file.exists(unlist(output(args)[[i]])))
        }
        return(unlist(targets))
    }
}

#########################################################
## Update the output location after run runCommandline ##
#########################################################
output_update <- function(args, dir=FALSE, dir.name=NULL, replace=FALSE, extension=NULL, make_bam=FALSE){
  ## Folder name provide in the yml file
  ## this file will exists, because its create on the definition of the project or when runCommandline is used.
  # if(is.null(args$yamlinput$results_path$path)) {
  #   if(is.null(dir.name)) {
  #     stop("argument 'dir.name' missing. The argument can only be assigned 'NULL' when directory name is provided in the yml template. The argument should be assigned as a character vector of length 1")
  #   }}
  ## Validation for 'args'
  if(any(!inherits(args, "SYSargs") & !inherits(args, "SYSargs2"))) stop("Argument 'args' needs to be assigned an object of class 'SYSargs' OR 'SYSargs2'")
  ## If the argument 'replace' is TRUE, it is required to specify the 'extension' argument
  if(replace!=FALSE){
    if(is.null(extension)) {
      stop("argument 'extension' missing. The argument can only be assigned 'NULL' when no replacement is required. The argument should be assigned as a character vector with the current extension and the new one.")  
    }}
  if(replace==TRUE){
    args <- as(args, "list")
    for(i in seq_along(args$output)){
      for(j in seq_along(args$output[[i]])){
        for(k in seq_along(args$output[[i]][[j]])){
          name <- basename(args$output[[i]][[j]][k])
          dirRep <- sub("/([^/]*)$", "", args$output[[i]][[j]][k])
          if(grepl(extension[1], name)){
            sam <- .getExt(name)
            args$output[[i]][[j]][k] <- file.path(dirRep, paste0(.getFileName(name), extension[2]))
          } else {
            # message if the extension are not matching, return the same object
            args$output[[i]][[j]][k] <- args$output[[i]][[j]][k]
          }
          if(make_bam==TRUE){
            if(grepl("bam", args$output[[i]][[j]][k])){
              args$output[[i]][[j]][length(args$output[[i]][[j]])+1] <- paste0(gsub("\\.bam$", "", args$output[[i]][[j]][k]), ".bam.bai")
            }
          }
        }
      }
    }
    args <- as(args, "SYSargs2")
  }
  if(dir==TRUE){
    args <- as(args, "list")
    ## Results path
    logdir <- normalizePath(args$yamlinput$results_path$path)
    ## Workflow Name: Detail: if the folder was not created during 'runCommandline,' it will return a warning message pointing 'no such directory'
    if(is.null(dir.name)) {
        cwl.wf <- file.path(logdir, .getFileName(args$files$cwl))
    } else if(!is.null(dir.name)){
        cwl.wf <- file.path(logdir, dir.name)
    }
    ## New path
    for(i in seq_along(args$output)){
      for(j in seq_along(args$output[[i]])){
        for(k in seq_along(args$output[[i]][[j]])){
          name <- basename(args$output[[i]][[j]][k])
          args$output[[i]][[j]][k] <- file.path(cwl.wf, name)
        }
      }
    }
    args <- as(args, "SYSargs2")
  }
  return(args)
}
## Usage:
# WF <- output_update(WF, dir=FALSE, replace=TRUE, extension=c(".sam", ".bam"))
# WF <- output_update(WF, dir=TRUE, replace=TRUE, extension=c(".sam", ".bam"))

##################################
## Unexported helper functions ##
##################################

##############################################
## Resolve CWL variable-based path instance ##
##############################################
pathInstance <- function(pathvar, input, altinput) {
  pathvar <- unlist(strsplit(pathvar[[1]], "/"))
  extension <- gsub("^.*\\)", "", pathvar)
  pathvar <- gsub("(^.*\\)).*", "\\1", pathvar)
  pathvar <- gsub("\\$|\\(|\\)", "", pathvar)
  pathvarlist <- strsplit(pathvar, "\\.")
  filenametype <- unlist(lapply(seq_along(pathvarlist), function(x) pathvarlist[[x]][pathvarlist[[x]] %in% c("basename", "nameroot")]))
  filenametype <- sapply(seq_along(pathvarlist), function(x) filenametype[x]) # In case of empty filenamelist NA are returned instead
  filenametype <- ifelse(is.na(filenametype), "NA", filenametype)
  myvalue_list <- sapply((pathvarlist), function(x) list(NULL), simplify=FALSE)
  for(i in seq_along(pathvarlist)){
    myvalue <- NULL
    for(j in seq_along(pathvarlist[[i]])){
      myvalue_input <- input[[pathvarlist[[i]][[j]]]]
      if("class" %in% names(myvalue_input)) myvalue_input$class <- NULL
      if (!is.null(myvalue_input)) {
        myvalue <- c(myvalue, list(myvalue_input))
      }
    }
    myvalue_list[i] <- myvalue
  }
  ## Use altinput if no values were assigned by input
  if(all(sapply(myvalue_list, length)==0)) {
    myvalue_list <- sapply(seq_along(pathvarlist), function(x) altinput[[pathvarlist[[x]][[2]]]]$default, simplify=FALSE)
  }
  enforce_list <- length(myvalue_list) == 1
  if(any(enforce_list)) for(i in which(enforce_list)) myvalue_list[[i]] <- list(class=NULL, path=myvalue_list[[i]]) 
  mypathvec <- list()
  for(i in seq_along(myvalue_list)){
    for(j in seq_along(myvalue_list[[i]])){
      vec_input <- myvalue_list[[i]][[j]]
      if (any(names(myvalue_list[[i]][[j]]) %in% c("path"))){
        vec_input <- myvalue_list[[i]][[j]]$path
      }
      if (!is.null(vec_input )) {
        vec <- sapply(seq_along(filenametype), function(x) pathUtils(vec_input [x], type=filenametype[x]))
        mypathvec <- c(mypathvec, list(vec))}
    }
    mypathvec <- lapply(mypathvec, function(x) x[!is.na(x)])
  }
  ## Assign '.' to 'runtime.outdir'
  pwdindex <- any(pathvar %in% c("runtime.outdir"))
  nullindex <- any(sapply(mypathvec, length) == 0)
  if(any(pwdindex & nullindex)) {
    mypathvec[pwdindex & nullindex] <- "."
  }
  mypathvec <- sapply(seq_along(mypathvec), function(x) if(is.null(mypathvec[[x]])) {mypathvec[[x]] <- "" } else {mypathvec[[x]] <- mypathvec[[x]]})
  ## Generate output
  if(length(mypathvec)==1){
    extension <- extension[!is.na(extension)]
    if(all("" %in% extension & length(extension)==1)){
      mypath <- paste(mypathvec)
    } else {
      extension <- extension[extension != ""]
      mypath <- paste(mypathvec, extension, sep="/") }
    } else if(any(is.na(extension))){
      mypath <- sapply(seq_along(mypathvec), function(x) paste0(mypathvec[x]))
    } else {
      mypath <- sapply(seq_along(mypathvec), function(x) paste0(mypathvec[x], extension[x]))
      if(length(mypathvec) < length(extension)){
        extension <- extension[extension != ""]
        mypath <- c(mypath, extension)
      }
    } 
  returnpath <- file.path(paste(mypath, collapse="/"))
  return(returnpath)
}

#################################################################
## Helper function to construct/manipulate path and file names ##
#################################################################
pathUtils <- function(x, type, dropdir=TRUE) {
    if(type=="dirname") {
        mypath <- dirname(x)
    } else if(type=="basename") {
        mypath <- basename(x)    
    } else if(type=="nameroot") {
        mypath <- gsub("(^.*)\\..*$", "\\1", basename(x))    
    } else {
        mypath <- x # Return unchanged input if 'type' is not one the above three values
        # warning("Argument 'type' needs to be assigned one of: 'dirname', 'basename', 'nameroot'")
    }
    ## Construct output
    if(dropdir==TRUE) {
        return(mypath) 
    } else if(dropdir==FALSE & type!="dirname") {
        return(file.path(dirname(x), mypath))
    } else {
        stop("Invalid path request, e.g. drop and appenrequest, e.g. drop and append dir.")
    }
}

#############################################
## Assemble Commandline Components in List ##
#############################################
assembleCommandlineList <- function(clt=WF$clt[[1]]) {
    ## global functions or variables
    WF <- NULL
    ## Base command and arguments
    basecommand <- clt$baseCommand
    arguments <- clt$arguments
    ## Special case for Rscript --  get the absolute path to the Rscript command
    if("Rscript" %in% basecommand){
        basecommand <- file.path(R.home('bin'), 'Rscript')
    }
    if(!is.null(arguments)){
      for(i in seq_along(arguments)) arguments[[i]][["position"]] <- ""
    }
    # ## Handling of special cases (here mkdir)
    # if(basecommand[1]=="mkdir") {
    #     clt$inputs <- ""
    #     clt$outputs[[1]]$type <- NULL
    # }
    ## Inputs
    inputargs <- renderInputs(x=clt$inputs, returntags=list(inputBinding=c("prefix"), type="any"))
    ## Outputs
    outputargs <- renderOutputs(x=clt$outputs, stdout=clt$stdout, returntags=list(outputBinding=c("prefix", "glob"), type="any"))
    ## Assemble command-line
    arglist <- list(baseCommand=basecommand, arguments=arguments, inputs=inputargs, outputs=outputargs)
    return(arglist)
}

########################################
## Populate Variables in Command-line ##
########################################
populateCommandline <- function(WF, cltid, exclude=exclude, mmp) {
    ## Populate inputs
    cmdlist <- assembleCommandlineList(nameUnnamed(WF$clt[[cltid]]))
    # yamlinput <- list(inputs=WF$yamlinput) # fix 29Dec18
    yamlinput <- list(inputs=nameUnnamed(WF$yamlinput))
    cmdnames <- nestedNames(cmdlist)
    yamlnames <- nestedNames(yamlinput)
    ## Create yamlnames list without 'unnamed' tags if present
    yamlnamestmp <- sapply(names(yamlnames), function(x) yamlnames[[x]][!grepl("_\\d{1,}_", yamlnames[[x]])], simplify=FALSE)
    matchmap <- findBestMatch(x=cmdnames, y=yamlnamestmp, mmp, minmatch=2)
    excludepos <- sapply(names(cmdnames), function(x) tail(cmdnames[[x]], 1) %in% exclude)
    removetype <- sapply(cmdnames, function(x) tail(x %in% "type", 1))
    for(i in seq_along(cmdnames)) {
        inputslice <- yamlnames[[matchmap[[i]]]]
        ## Inject yamlinput in corresponding cmdlist components
        if((length(inputslice) > 0) & (excludepos[[i]])!=TRUE) {
            cmdlist <- injectListbyName(cmdlist, cmdnames[[i]], value=subsetListbyName(yamlinput, inputslice))
        }
        ## Eliminate type values (e.g. boolean/string) that are not associated with yamlinput
        if((length(inputslice) == 0) & (excludepos[[i]]!=TRUE & removetype[[i]]==TRUE)) {
            cmdlist <- injectListbyName(cmdlist, cmdnames[[i]], value="")
        }
    }
    ## Populate variable path instances
    pos <- sapply(names(cmdnames), function(x) grepl("^\\$\\(", subsetListbyName(cmdlist, cmdnames[[x]])))
    if(any(pos)) {
        for(i in which(pos)) {
            mypath <- pathInstance(pathvar=subsetListbyName(cmdlist, cmdnames[[i]]), input=WF$yamlinput, altinput=nameUnnamed(WF$clt[[cltid]])$inputs)
            cmdlist <- injectListbyName(cmdlist, cmdnames[[i]], value=mypath)
        }
    }
    return(cmdlist)
}

####################################################################################
## Inject command-line lists from populateCommandline into cmdlist slot WF object ##
####################################################################################
injectCommandlinelist <- function(WF) {
    if(class(WF)=="SYSargs2") {
        WF <- as(WF, "list") 
    } else {
        stop("WF needs to be object of class 'SYSargs2'.")
    }
    cmdsteps <- names(WF$cmdlist)
    mmp <- list(class=c("class", "type"), path=c("path", "type"))
    exclude <- c("baseCommand", "prefix")
    for(i in seq_along(cmdsteps)) { 
        WF$cmdlist[[i]] <- populateCommandline(WF, cltid=cmdsteps[i], exclude=exclude, mmp)
    } 
    ## Update inputs that depend on output of one of the previous steps
    .connectInout <- function(WF) {
        steps <- WF$wf$steps
        stepnames <- sapply(names(steps), function(x) steps[[x]]$run)
        connectedlist <- sapply(seq_along(steps), function(x) which(grepl("/", steps[[x]]$`in`)))
        names(connectedlist) <- names(steps)
        connectedlist <- connectedlist[sapply(connectedlist, length) > 0]
        if(length(connectedlist) > 0) {
            for(j in seq_along(connectedlist)) {
                connected <- connectedlist[j]
                ## Inner loop to support multiple connections
                for(i in seq_along(connected[[1]])) {
                    ## Generate output string from corresponding upstream step
                    connectednames <- strsplit(unlist(steps[[names(connected)]])[connected[[1]][i]], "/")
                    connectedstep <- stepnames[[connectednames[[1]][1]]]
                    connectedoutput <- WF$cmdlist[[connectedstep]]$outputs
                    connectedoutput <- unlist(connectedoutput[connectednames[[1]][[2]]])
                    ## Inject output string into input of corresponding downstream step
                    termname <- as.character(sapply(nestedNames(steps[[names(connected)]]), tail, 1)[connected[[1]][i]])
                    namepath <- c("cmdlist", stepnames[[names(connected)]], "inputs", termname)
                    #fix 02Jan19 # namepath <- c("cmdlist", stepnames[[names(connected[i])]], "inputs", 
                    #                            names(steps[[names(connected[i])]][connected[[i]]][["in"]]))
                    WF <- injectListbyName(l=WF, namepath, value=as.character(connectedoutput), type="value") 
                }
            }
        }
        names(WF$cmdlist) <- WF$files$steps
        names(WF$clt) <- WF$files$steps
        return(WF)
    }
    WF <- .connectInout(WF)
    return(as(WF, "SYSargs2"))
}

################################
## Render Commandline Strings ##
################################
renderCommandline <- function(x, dropoutput=TRUE, redirect=">") {
    if(class(x)=="SYSargs2") {
        x <- as(x, "list") 
    } else {
        stop("x needs to be object of class 'SYSargs2'.")
    }
    ## Fct to render single command-line string
    .renderCommandline <- function(x=x, redirect=redirect) {
        cmdnames <- nestedNames(x)
        ## Check for stdout instance, if present prepend value of redirect
        stdoutpos <- sapply(cmdnames, function(y) {
                    tmp <- y %in% c("outputs", "stdout")
                    tmp[1] == TRUE & tail(tmp, 1) == TRUE
                    })
        if(any(stdoutpos)) {
            cmd <- x
            cmdnamessub <- cmdnames[stdoutpos]
            for(i in seq_along(cmdnamessub)) {
                stdoutstr <- unlist(subsetListbyName(l=cmd, cmdnamessub[[i]]))
                stdoutstr <- paste(redirect, stdoutstr, sep=" ")
                cmd <- injectListbyName(l=cmd, cmdnamessub[[i]], value=stdoutstr, type="value") 
            }
        } else {
          #if(dropoutput==TRUE & x$baseCommand[1]!="mkdir") x <- x[names(x) != "outputs"] 
          if(dropoutput==TRUE) x <- x[names(x) != "outputs"] 
            cmd <- x
        }
        ## Collapse list to single string
        cmd <- paste(unlist(cmd), collapse=" ")
        return(cmd)
    }
    ## Check here for WF class instead of names once implemented
    if(all(names(x) == c("targets", "targetsheader", "modules", "wf", "clt", "yamlinput", "cmdlist", "input", "output", "files", "inputvars", "cmdToCwl", "status"))) {
        cmdstring <- sapply(names(x$cmdlist), 
                            function(i) .renderCommandline(x=x$cmdlist[[i]], redirect=redirect),
                            simplify=FALSE)
    } else {
        cmdstring <- .renderCommandline(x=x, redirect=redirect)
    }
    return(cmdstring)
}

########################################################### 
## Helper Function to render inputs of single clt object ##
###########################################################
renderInputs <- function(x=WF$clt[[1]]$inputs, returntags=list(inputBinding=c("prefix"), type="any")) {
    ## global functions or variables
    WF <- NULL
    inputnames <- names(x)
    ## Remove entries where inputBinding is missing since those parameters do not appear on the command-line
    inputnames <- inputnames[!sapply(inputnames, function(i) is.null(x[[i]]$inputBinding))]
    x <- x[inputnames]
    inputlist <- sapply(inputnames, function(i) {
            tmp <- c(x[[i]][[names(returntags["inputBinding"])]][returntags$inputBinding],
                  x[[i]][names(returntags["type"])]) 
                  # input=x[[i]][[names(returntags["type"])]]) 
            tmp <- tmp[!sapply(tmp, length) < 1]
            tmp}, simplify=FALSE)
    inputposition <- sapply(inputnames, function(i) x[[i]]$inputBinding$position)
    if(class(inputposition)=="list") inputposition <- seq_along(inputposition)
    inputlist <- inputlist[order(inputposition)]
    return(inputlist)
}

############################################################
## Helper Function to render outputs of single clt object ##
############################################################
renderOutputs <- function(x=WF$clt[[1]]$outputs, stdout=WF$clt[[1]]$stdout, returntags=list(outputBinding=c("prefix", "glob"), type="any")) {
    ## global functions or variables
    WF <- NULL
    outputnames <- names(x)
    outputlist <- sapply(outputnames, function(i) {
            tmp <- c(x[[i]][[names(returntags["outputBinding"])]][returntags$outputBinding[1]],
                    x[[i]][[names(returntags["outputBinding"])]][returntags$outputBinding[2]],
                    type=x[[i]][[names(returntags["type"])]]) 
                    # x[[i]][[names(returntags["type"])]]) 
            tmp <- tmp[!sapply(tmp, length) < 1]
            tmp}, simplify=FALSE)
    ## Remove entries 'type: Directory' since they are only relevant for CWL internals
    types <- sapply(seq_along(outputlist), function(i) outputlist[[i]][["type"]]) 
    outputlist <- outputlist[tolower(types) != "directory"]
    for(i in seq_along(outputlist)) {
        check <- outputlist[[i]][["type"]]
        if(is.null(check)) check <- "absent"
        if(check == "stdout") {
            # outputlist[[i]] <- c(outputlist[[i]]["type"], stdout=stdout) 
            outputlist[[i]] <- list(stdout=stdout)
        }
    }
    return(outputlist)
}

##################
## Name unnamed ##
##################
## Yaml format may have unnamed components; not having names creates
## problems in some cases with name-based subsetting of lists; thus the 
## following function assigns names to unnamed components. A tagging syntax
## can be used to ignore them in subsetting routines if necessary.
nameUnnamed <- function(l) {
    for(a in seq_along(l)) {
        for(b in seq_along(l[[a]])) {
            for(c in seq_along(l[[a]][[b]])) {
                for(d in seq_along(l[[a]][[b]][[c]])) {
                    if(is.null(names(l[[a]][[b]][[c]][[d]])) & class(l[[a]][[b]][[c]][[d]])=="list") names(l[[a]][[b]][[c]][[d]]) <- paste0("_", seq_along(l[[a]][[b]][[c]][[d]]), "_") 
                    if(is.null(names(l[[a]][[b]][[c]])) & class(l[[a]][[b]][[c]])=="list") names(l[[a]][[b]][[c]]) <- paste0("_", seq_along(l[[a]][[b]][[c]]), "_") 
                    if(is.null(names(l[[a]][[b]])) & class(l[[a]][[b]])=="list") names(l[[a]][[b]]) <- paste0("_", seq_along(l[[a]][[b]]), "_") 
                    if(is.null(names(l[[a]])) & class(l[[a]])=="list") names(l[[a]]) <- paste0("_", seq_along(l[[a]]), "_") 
                }
            }
        }
    }
    return(l)
}

###########################################################################
## Return names of nested lists (with up to 5 nesting levels) as vectors ##
###########################################################################
nestedNames <- function(l, sep="_") {
    nestednames <- lapply(seq_along(l), function(a) 
                      lapply(seq_along(l[[a]]), function(b) 
                             lapply(seq_along(l[[a]][[b]]), function(c) 
                                lapply(seq_along(l[[a]][[b]][[c]]), function(d) 
                                    lapply(seq_along(l[[a]][[b]][[c]][[d]]), function(e) 
                                    c(names(l[a]), 
                                      names(l[[a]][b]), 
                                      names(l[[a]][[b]][c]), 
                                      names(l[[a]][[b]][[c]][d]),
                                      names(l[[a]][[b]][[c]][[d]][e])))))))
    while(any(sapply(nestednames, is.list))) nestednames <- unlist(nestednames, recursive=FALSE)
    names(nestednames) <- sapply(nestednames, paste, collapse=sep)
    return(nestednames)
}

########################################################
## Find best matches among name vectors for two lists ##
########################################################
findBestMatch <- function(x, y, mmp, minmatch=2) {
	matchindex <- setNames(rep(NA, length(x)), names(x))	
	## Perfect matches
	matchpos <- sapply(names(y), function(j) { 
					sapply(names(x), function(i) {
						ident <- all(y[[j]] %in% x[[i]])
						lengthcheck <- length(y[[j]]) == length(x[[i]])
						ident & lengthcheck
						})
					})
	## Note, if there are several input matches to one command-line slot then only the last one is used (see tail)
	matchpos <- sapply(rownames(matchpos), function(i) tail(which(matchpos[i,]),1), simplify=FALSE)
	for(i in seq_along(matchindex)) {
		if(is.na(matchindex[i])==TRUE) matchindex[i] <- matchpos[[i]][1]
	}
	## Matches with tag specific terminal mismatches
	matchpos <- sapply(names(y), function(j) { 
					sapply(names(x), function(i) {
						all(termMMatch(x=x[[i]], y=y[[j]], mmp=mmp, minmatch=minmatch, returntype="logical"))
						})
					})
	## Note, if there are several input matches to one command-line slot then only the last one is used (see tail)
	matchpos <- sapply(rownames(matchpos), function(i) tail(which(matchpos[i,]),1), simplify=FALSE)
	for(i in seq_along(matchpos)) {
		if(is.na(matchindex[i])==TRUE) matchindex[i] <- matchpos[[i]][1]
	}
	## Matches with terminal (max 1) length differences 
	matchpos <- sapply(names(y), function(j) { 
					sapply(names(x), function(i) {
						ident <- all(y[[j]] %in% x[[i]][1:length(y[[j]])])
						lengthcheck <- length(y[[j]]) >= minmatch
						ident & lengthcheck
						})
					})
	## Note, if there are several input matches to one command-line slot then only the last one is used (see tail)
	matchpos <- sapply(rownames(matchpos), function(i) tail(which(matchpos[i,]),1), simplify=FALSE)
	for(i in seq_along(matchpos)) {
		if(is.na(matchindex[i])==TRUE) matchindex[i] <- matchpos[[i]][1]
	}
	return(matchindex)
}

###############################################
## Inject into List Subsetted by Name Vector ##
###############################################
injectListbyName <- function(l, name_index, value, type="value") {
    ## Input validity check for l and name_index
    if(all(!c("name", "value") %in% type[1])) stop("Argument type can only be assigned one of: 'name' or 'value'.")
    checknames <- nestedNames(l) 
    checknames <- checknames[!sapply(checknames, length) < length(name_index)]
    checknames <- sapply(checknames, function(x) which(name_index == x[seq_along(name_index)]), simplify=FALSE)
    checknames <- checknames[sapply(checknames, length) == length(name_index)]    
    checkindex <- sapply(seq_along(checknames), function(x) all(checknames[[x]] == seq_along(checknames[[x]])))
    if(length(checkindex)==0) checkindex <- FALSE
    if(any(checkindex==FALSE)) stop("Invalid 'name_index' lacking consecutive name matches in any list component.") 
    if(length(name_index)==1) {
        if(type=="name") names(l)[which(names(l) %in% name_index[1])] <- value
        if(type=="value") l[[name_index[1]]] <- value
    } else if(length(name_index)==2) {
        if(type=="name") names(l[[name_index[1]]])[1] <- value
        if(type=="value") l[[name_index[1]]][[name_index[2]]] <- value
    } else if(length(name_index)==3) {
        if(type=="name") names(l[[name_index[1]]][[name_index[2]]]) <- value
        if(type=="value") l[[name_index[1]]][[name_index[2]]][[name_index[3]]] <- value
    } else if(length(name_index)==4) {
        if(type=="name") names(l[[name_index[1]]][[name_index[2]]][[name_index[3]]]) <- value
        if(type=="value") l[[name_index[1]]][[name_index[2]]][[name_index[3]]][[name_index[4]]] <- value
    } else if(length(name_index)==5) {
        if(type=="name") names(l[[name_index[1]]][[name_index[2]]][[name_index[3]]][[name_index[4]]]) <- value 
        if(type=="value") l[[name_index[1]]][[name_index[2]]][[name_index[3]]][[name_index[4]]][[name_index[5]]] <- value 
    } else {
        stop("Nesting level (length of name_index) cannot exceed 5.")
    }
    return(l)
}

################################
## Subset List by Name vector ##
################################
subsetListbyName <- function(l, name_index) {
    ## Input validity check for l and name_index
    checknames <- nestedNames(l) 
    checknames <- checknames[!sapply(checknames, length) < length(name_index)]
    checknames <- sapply(checknames, function(x) which(name_index == x[seq_along(name_index)]), simplify=FALSE)
    checknames <- checknames[sapply(checknames, length) == length(name_index)]    
    checkindex <- sapply(seq_along(checknames), function(x) all(checknames[[x]] == seq_along(checknames[[x]])))
    if(length(checkindex)==0) checkindex <- FALSE
    if(any(checkindex==FALSE)) stop("Invalid 'name_index' lacking consecutive name matches in any list component.") 
    if(length(name_index)==1) {
        lsub <- l[name_index[1]] 
    } else if(length(name_index)==2) {
        lsub <- l[[name_index[1]]][name_index[2]]
    } else if(length(name_index)==3) {
        lsub <- l[[name_index[1]]][[name_index[2]]][name_index[3]] 
    } else if(length(name_index)==4) {
        lsub <- l[[name_index[1]]][[name_index[2]]][[name_index[3]]][name_index[4]] 
    } else if(length(name_index)==5) {
        lsub <- l[[name_index[1]]][[name_index[2]]][[name_index[3]]][[name_index[4]]][name_index[5]] 
    } else {
        stop("Nesting level (length of name_index) cannot exceed 5.")
    }
    return(lsub)
}

############################################################
## Matching vectors with tag specific terminal mismatches ##
############################################################
termMMatch <- function(x, y, mmp, minmatch=2, returntype="values") {
    ## Input validity checks
    if(!all(sapply(mmp, length) == 2)) stop("List components of 'mmp' need to be vectors of length = 2.")
    if(!all(sapply(names(mmp), function(x) x %in% mmp[[x]]))) stop("Names of 'mmp' list components need to match one of the entries in the corresponding list.")
    ## Check for candidate matches
    lengthcheck <- length(x) == length(y) 
    term <- unique(c(tail(x, 1), tail(y, 1)))
    l <- length(term)==2 
    minmatchcheck <- sum(x[-length(x)] == y[-length(y)]) >= minmatch
    index <- which(sapply(names(mmp), function(x) all(term %in% mmp[[x]])))
    p <- length(index) > 0
	if(lengthcheck & l & p & minmatchcheck) {
        for(i in index) {
            x[length(x)] <- names(mmp)[i]
            y[length(y)] <- names(mmp)[i]
        }
		if(returntype=="values") {
			return(list(x=x, y=y))
		} else if(returntype=="logical") {
			return(rep(TRUE, length(x)))
		} else {
			stop("Argument 'returntype' can only be assigned one of: 'values' or 'logical'.")
		}
    ## Return input if there are not matches
	} else {
		if(returntype=="values") {
			return(list(x=x, y=y))
		} else if(returntype=="logical") {
			x <- x == y[1:length(x)]
			x[is.na(x)] <- FALSE
			return(x) 
		} else {
			stop("Argument 'returntype' can only be assigned one of: 'values' or 'logical'.")
		}
	}
}

#########################
##   .renderWFsingle   ##
#########################
.renderWFsingle <- function(WF, id, inputvars) { 
  inputvarslist <- sapply(names(inputvars), function(x) "", simplify=FALSE)
  if(!length(names(targets(WF)))==0) (if(any(!names(inputvars) %in% colnames(targets.as.df(WF$targets)))) stop("Please note that the 'inputvars' variables need to be defined in the 'input_file'; as well it needs to match the column names defined in the 'targets' file."))
  input <- yamlinput(WF)
  for(i in seq_along(inputvars)) {
    subvalue <- targets(WF)[[id]][[names(inputvars)[i]]]
    if(length(subvalue)!=0) {
      input <- rapply(input, function(x) gsub(inputvars[[i]], subvalue, x), how = "replace")
      inputvarslist[[i]] <- subvalue
    }
  }
  WF <- as(WF, "list")
  WF$yamlinput <- input
  WF <- as(WF, "SYSargs2")
  WF <- injectCommandlinelist(WF)
  ## Fix for cases like IDX/tophat/STAR, more than one output file per cmdlist
  outfilelist <- sapply(names(cmdlist(WF)), function(x) list(NULL))
  for(i in seq_along(outfilelist)) {
    for (j in seq_along(cmdlist(WF)[[names(outfilelist[i])]]$output))
        if("stdout" %in% names(cmdlist(WF)[[names(outfilelist[i])]]$output[[j]])){
            outfilelist[[i]][j] <- cmdlist(WF)[[names(outfilelist[i])]]$output[[j]]$stdout
        } else {
            outfilelist[[i]][j] <- cmdlist(WF)[[names(outfilelist[i])]]$output[[j]]$glob
        }
  }
  cmdlist <- renderCommandline(WF, redirect=">")
  inputvars <- as.list(inputvars)
  return(list(cmd=cmdlist, input=inputvarslist, output=outfilelist, inputvars=inputvars))
}

###################
##   write.cwl   ##
###################
write.clt <- function(commandLine, cwlVersion, class, file.cwl, writeout=TRUE, silent=FALSE) {
  cwlVersion <- cwlVersion 
  class <- class
  baseCommand <- commandLine$baseCommand[[1]]
  ## File
  clt <- list(cwlVersion=cwlVersion, class=class)
  ##requirements
  if(is.null(commandLine$requeriments)){
    dump <- "do nothing" 
  } else if(!is.null(commandLine$requeriments)) {
    requeriments <- list() ##TODO
    clt <- c(clt, list(requeriments=requeriments))
  }
  ##ARGUMENTS
  if(is.null(commandLine$arguments)){
      dump <- "do nothing"
  } else if(!is.null(commandLine$arguments)) {
    arguments <- sapply(seq_along(commandLine$arguments), function(x) list())
    for(i in seq_along(commandLine$arguments)){
      arguments[[i]]["prefix"] <- commandLine$arguments[[i]]$preF
      arguments[[i]]["valueFrom"] <- commandLine$arguments[[i]]$valueFrom
    }
    clt <- c(clt, list(arguments=arguments))
  }
  ##INPUTS
  if(any(names(commandLine$inputs)=="")) stop("Each element of the list 'commandLine' needs to be assigned a name")
  if(is.null(names(commandLine$inputs))) stop("Each element of the list 'commandLine' needs to be assigned a name")
  inputs <- sapply(names(commandLine$inputs), function(x) list())
  for(i in seq_along(commandLine$inputs)){
    if("type" %in% names(commandLine$inputs[[i]])){
      if(!c("type") %in% names(commandLine$inputs[[i]])) stop("Each element of the sublist 'inputs' in 'commandLine' needs to be defined the type of the argument, for example: type='Directory' or type='File' or type='int' or type='string'")
      inputs[[i]]["type"] <-commandLine$inputs[[i]]$type 
    } 
    if("preF" %in% names(commandLine$inputs[[i]])){
      if(commandLine$inputs[[i]]$preF==""){
        inputs[[i]]["inputBinding"] <- list(list(prefix=NULL))
      } else{
        inputs[[i]]["inputBinding"] <- list(list(prefix=commandLine$inputs[[i]]$preF))
      }
    } 
    if(any(c("label", "secondaryFiles", "doc", "default", "format", "streamable") %in% names(commandLine$inputs[[i]]))){
      for(j in which(c("label", "secondaryFiles", "doc", "default", "format", "streamable") %in% names(commandLine$inputs[[i]]))){
        inputs[[i]][names(commandLine$inputs[[i]][j])] <- commandLine$inputs[[i]][names(commandLine$inputs[[i]])][[j]]
      }
    }
  }
  ##OUTPUTS
  outputs <- sapply(names(commandLine$outputs), function(x) list())
  if(!c("type") %in% names(commandLine$inputs[[i]])) stop("Each element of the sublist 'outputs' in 'commandLine' needs to be defined the type of the argument, for example: type='Directory' or type='File'.")
  if(all(!c("File", "Directory") %in% commandLine$outputs[[1]]$type)) stop("Each element of the sublist 'outputs' in 'commandLine' needs to be defined the type = 'Directory', 'File'.")
  for(i in seq_along(commandLine$outputs)){
    outputs[[i]]["type"] <- commandLine$outputs[[i]]$type
    outputs[[i]]["outputBinding"] <- list(list(glob=commandLine$outputs[[i]][[2]]))
  }
  ## FILE
    clt <- c(clt, list(baseCommand=baseCommand, inputs=inputs, outputs=outputs))
  ## writing file
    if (writeout == TRUE) {
        yaml::write_yaml(x=clt, file = file.cwl)
        ## print message 
        if (silent != TRUE) cat("\t", "Written content of 'commandLine' to file:", "\n", file.cwl, "\n")
    }
  clt <- list(clt)
  names(clt) <- baseCommand
  return(clt)
}

## Usage:
# clt_cwl <- write.clt(commandLine, cwlVersion, class, file.cwl) 

###################
##   write.yml   ##
###################
## Write the yaml file  
write.yml <- function(commandLine, file.yml, results_path, module_load, writeout=TRUE, silent=FALSE){
  inputs <- commandLine$inputs
  if(any(names(inputs)=="")) stop("Each element of the list 'commandLine' needs to be assigned a name")
  if(is.null(names(inputs))) stop("Each element of the list 'commandLine' needs to be assigned a name")
  ##yamlinput_yml 
  yamlinput_yml <- sapply(names(inputs), function(x) list())
  for(i in seq_along(inputs)){
      ##TODO!!!
    if(!c("type") %in% names(inputs[[i]])) stop("Each element of the sublist 'inputs' in 'commandLine' needs to be defined the type of the argument, for example: type='Directory' or type='File' or type='int' or type='string'")
    if("type" %in% names(inputs[[i]])){
      if(any(c("File", "Directory") %in% inputs[[i]])){
        yamlinput_yml[[i]]["class"] <- inputs[[i]]$type
        yamlinput_yml[[i]]["path"] <- inputs[[i]]$yml
      } else if (any(c("int", "string") %in% inputs[[i]])){
        yamlinput_yml[[i]]  <- inputs[[i]]$yml
      }
    } else {
      print("do something")
    }
  }
  ## results_path
  yamlinput_yml[["results_path"]]["class"] <- list("Directory")
  yamlinput_yml[["results_path"]]["path"] <- list(results_path)
  ## moduleload
  for(i in seq_along(module_load)){
    yamlinput_yml[["ModulesToLoad"]][paste0("module", i)] <- list(module_load[[i]])
  }
  ## write out the '.yml' file
  if (writeout == TRUE) {
      yaml::write_yaml(x=yamlinput_yml, file = file.yml)
      ## print message 
    if (silent != TRUE) cat("\t", "Written content of 'commandLine' to file:", "\n", file.yml, "\n")
  }
  return(yamlinput_yml)
}

## Usage: 
# yamlinput_yml <- write.yml(commandLine, file.yml, results_path, module_load) 

###################
##  cmdTool2wf   ##
###################
## CommandlineTool --> Workflow class
cmdTool2wf <- function(cmdTool_path, file.cwl, writeout=TRUE, silent=FALSE){
    cmdTools <- yaml::read_yaml(file.path(cmdTool_path))
    cwlVersion <- cmdTools$cwlVersion 
    class <- "Workflow"
    ## Input
    inputs_names <- names(cmdTools$inputs)
    inputs <- sapply(inputs_names, function(x) list(cmdTools$inputs[[x]]$type))
    ## output
    outputs <- sapply(names(cmdTools$outputs), function(x) list(NULL)) 
    for( i in seq_along(outputs)){
        outputs[[i]][["outputSource"]] <- paste0(cmdTools$baseCommand, "/", names(outputs)[i])
        outputs[[i]][["type"]] <- cmdTools$outputs[[i]]$type
    }
    ## Steps
    step.in <- sapply(names(inputs), function(x) list(x)) 
    step.out <- paste0("[", paste0(names(outputs), collapse = ", "), "]")
    steps <- list(list(`in`=step.in, `out`=step.out, run=cmdTool_path))
    if(length( cmdTools$baseCommand) >1){
        cmdTools$baseCommand <- paste0(cmdTools$baseCommand, collapse = "_")
    } else if(is.null(cmdTools$baseCommand)){
        cmdTools$baseCommand <-.getFileName(cmdTool_path)
    }
    names(steps) <- cmdTools$baseCommand
    ## Combine
    wf2 <- list(class= class, cwlVersion= cwlVersion, inputs=inputs, 
                outputs=outputs, steps=steps)
    ## write out the '.cwl' file
    if (writeout == TRUE) {
        yaml::write_yaml(x=wf2, file = file.cwl)
        ## print message 
        if (silent != TRUE) cat("\t", "Written content of 'Workflow' to file:", "\n", file.cwl, "\n")
    }
    ## Return
    return(wf2)
}

# cmdTools_path <- "param/cwl/hisat2/hisat2-pe/hisat2-mapping-pe.cwl"
# cmdTools_path <- "param/cwl/hisat2/hisat2-idx/hisat2-index.cwl"
# 
# ## Usage: 
# wf <- cmdTool2wf(cmdTools_path, file.cwl = "test.cwl") 
