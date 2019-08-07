##############################################
## Class and Method Definitions for SYSargs ##
##############################################
## Define SYSargs class
setClass("SYSargs", representation(targetsin="data.frame",	
					targetsout="data.frame",
		            targetsheader="character",
					modules="character", 
					software="character",
					cores="numeric", 
					other="character",
					reference="character",
					results="character",
					infile1="character",
					infile2="character",
					outfile1="character",
					sysargs="character", 
					outpaths="character")
)

## Methods to return SYSargs components 
setGeneric(name="targetsin", def=function(x) standardGeneric("targetsin"))
setMethod(f="targetsin", signature="SYSargs", definition=function(x) {return(x@targetsin)})
setGeneric(name="targetsout", def=function(x) standardGeneric("targetsout"))
setMethod(f="targetsout", signature="SYSargs", definition=function(x) {return(x@targetsout)})
setGeneric(name="targetsheader", def=function(x) standardGeneric("targetsheader"))
setMethod(f="targetsheader", signature="SYSargs", definition=function(x) {return(x@targetsheader)})
setGeneric(name="modules", def=function(x) standardGeneric("modules"))
setMethod(f="modules", signature="SYSargs", definition=function(x) {return(as.character(x@modules))})
setGeneric(name="software", def=function(x) standardGeneric("software"))
setMethod(f="software", signature="SYSargs", definition=function(x) {return(as.character(x@software))})
setGeneric(name="cores", def=function(x) standardGeneric("cores"))
setMethod(f="cores", signature="SYSargs", definition=function(x) {return(x@cores)})
setGeneric(name="other", def=function(x) standardGeneric("other"))
setMethod(f="other", signature="SYSargs", definition=function(x) {return(x@other)})
setGeneric(name="reference", def=function(x) standardGeneric("reference"))
setMethod(f="reference", signature="SYSargs", definition=function(x) {return(x@reference)})
setGeneric(name="results", def=function(x) standardGeneric("results"))
setMethod(f="results", signature="SYSargs", definition=function(x) {return(as.character(x@results))})
setGeneric(name="infile1", def=function(x) standardGeneric("infile1"))
setMethod(f="infile1", signature="SYSargs", definition=function(x) {return(x@infile1)})
setGeneric(name="infile2", def=function(x) standardGeneric("infile2"))
setMethod(f="infile2", signature="SYSargs", definition=function(x) {return(x@infile2)})
setGeneric(name="outfile1", def=function(x) standardGeneric("outfile1"))
setMethod(f="outfile1", signature="SYSargs", definition=function(x) {return(x@outfile1)})
setGeneric(name="SampleName", def=function(x) standardGeneric("SampleName"))
setMethod(f="SampleName", signature="SYSargs", definition=function(x) {return(names(x@sysargs))})
setGeneric(name="sysargs", def=function(x) standardGeneric("sysargs"))
setMethod(f="sysargs", signature="SYSargs", definition=function(x) {return(x@sysargs)})
setGeneric(name="outpaths", def=function(x) standardGeneric("outpaths"))
setMethod(f="outpaths", signature="SYSargs", definition=function(x) {return(x@outpaths)})

## Constructor methods
## List to SYSargs with: as(mylist, "SYSargs")
setAs(from="list", to="SYSargs",  
        def=function(from) {
		new("SYSargs", targetsin=from$targetsin,
					targetsout=from$targetsout,
					targetsheader=from$targetsheader,
					modules=from$modules, 
					software=from$software,
					cores=from$cores, 
					other=from$other,
					reference=from$reference,
					results=from$results,
					infile1=from$infile1, 
					infile2=from$infile2,
					outfile1=from$outfile1,
					sysargs=from$sysargs, 
					outpaths=from$outpaths)
})

## Define print behavior for SYSargs
setMethod(f="show", signature="SYSargs", 
	definition=function(object) {    
	cat("An instance of '", class(object), "' for running '", object@software, "' on ", length(object@sysargs), " samples ", "\n", sep="")
})

## Extend names() method
setMethod(f="names", signature="SYSargs",
		definition=function(x) {
		return(slotNames(x))
})

## Extend length() method
setMethod(f="length", signature="SYSargs",
	definition=function(x) {
	return(length(x@infile1))
})

## Behavior of "[" operator for SYSargs
setMethod(f="[", signature="SYSargs", definition=function(x, i, ..., drop) {
        if(is.logical(i)) {
                i <- which(i)
        }
        x@targetsin <- x@targetsin[i,]
        x@targetsout <- x@targetsout[i,]
        x@infile1 <- x@infile1[i]
        x@infile2 <- x@infile2[i]
        x@outfile1 <- x@outfile1[i]
        x@sysargs <- x@sysargs[i]
        x@outpaths <- x@outpaths[i]
        return(x)
})

## Construct SYSargs object from param and targets files
systemArgs <- function(sysma, mytargets, type="SYSargs") {
	## Read sysma and convert to arglist; if NULL is assigned to sysma then a dummy version is generated
	sysmapath <- sysma
	if(length(sysmapath)!=0) {
		sysma <- as.matrix(read.delim(sysma, comment.char = "#"))
		sysma[is.na(sysma)] <- ""
	} else {
		sysma <- cbind(PairSet=c("software", "cores", "other", "outfile1", "reference", "infile1", "infile1", "infile2", "infile2"), 
					Name=c("", "", "", "", "", "", "path", "", "path"), 
					Value=c("", "1", "", "<FileName1>", "", "<FileName1>", "", "<FileName2>", ""))
	}
	if(any(sysma[,1] %in% "type")) { # Detects software type: commandline or R
		iscommandline <- sysma[sysma[,1] %in% "type",, drop = FALSE]
		iscommandline <- as.logical(iscommandline[1, "Value"])
		sysma <- sysma[!sysma[,1] %in% "type",] # removes type row(s)
	} else {
		iscommandline <- TRUE # If type line not present then 'commandline' will be assumed
	}
	arglist <- sapply(as.character(unique(sysma[,"PairSet"])), function(x) as.vector(t(as.matrix(sysma[sysma[,"PairSet"]==x, 2:3]))), simplify=FALSE)
	for(i in seq(along=arglist)) names(arglist[[i]]) <- paste(rep(c("n", "v"), length(arglist[[i]])/2), rep(1:(length(arglist[[i]])/2), 2), sep="")
	if(type=="json") return(toJSON(arglist))
	## Read comment/header lines from targets file
	targetsheader <- readLines(mytargets)
        targetsheader <- targetsheader[grepl("^#", targetsheader)]
	## Validity checks
	mytargets <- read.delim(mytargets, comment.char = "#")
	mytargetsorig <- mytargets
	if(any(duplicated(mytargets$SampleName))) stop("SampleName column of mytargets cannot contain duplicated entries!")
	## Preprocessing of targets input
	colnames(mytargets)[1] <- "FileName1" # To support FileName column for SE data
	## Insert empty FileName2 column if not present
	if(length(mytargets$FileName2)==0) mytargets <- data.frame(FileName1=mytargets$FileName1, FileName2="", mytargets[,!colnames(mytargets) %in% "FileName1"])
	## Check name:value violations in arglist
	check <- sapply(names(arglist), function(x) sum(grepl("^n", names(arglist[[x]]))) == sum(grepl("^n", names(arglist[[x]]))))
	if(any(!check)) stop(paste("Name:Value violation in arglist component(s):", paste(names(check[check]), collapse=", ")))
	
	## Modify arglist object as specified in arglist and mytargets
	## Remove module component and store values in separate container
	modules <- as.character(arglist$modules[grepl("v", names(arglist$modules))])
	arglist <- arglist[!names(arglist) %in% "modules"]
	
	## Extract single value components
	software <- as.character(arglist$software[grepl("v", names(arglist$software))])
	other <- as.character(arglist$other[grepl("v", names(arglist$other))])
	if(!(is.na(arglist[["reference"]][["v1"]]) | nchar(arglist[["reference"]][["v1"]])==0)) {
		reference <- as.character(arglist$reference[grepl("v", names(arglist$reference))])
		if(!grepl("^/", reference)) reference <- paste0(getwd(), gsub("^\\.", "", reference)) # Turn relative into absolute path.
		arglist[["reference"]]["v1"] <- reference
	} else {
		reference <- ""
	}
	cores <- as.numeric(arglist$cores[grepl("v", names(arglist$cores))])	

	## Populate arglist$infile1
	if(any(grepl("^<.*>$", arglist$infile1))) { 
		infile1 <- gsub("<|>", "", arglist$infile1[grepl("^<.*>$", arglist$infile1)][[1]])
		infile1 <- as.character(mytargets[,infile1])	
		infile1 <- normalizePath(infile1)
		argname <- arglist$infile1[grep("<.*>", arglist$infile1)[1] -1]
		path <- arglist$infile1[grep("path", arglist$infile1)[1] +1]
		infile1back <- paste(path, infile1, sep="")
		names(infile1back) <- as.character(mytargets$SampleName)	
		infile1 <- paste(argname, " ", path, infile1, sep="")
		arglist[["infile1"]] <- gsub("(^ {1,})|( ${1,})", "", infile1)
	} else {
		infile1back <- rep("", length(mytargets[,1])) 
		infile1 <- infile1back
		names(infile1back) <- as.character(mytargets$SampleName)	
		arglist[["infile1"]] <- infile1back
	}
	
	## Populate arglist$infile2
	if(any(grepl("^<.*>$", arglist$infile2))) { 
		infile2 <- gsub("<|>", "", arglist$infile2[grepl("^<.*>$", arglist$infile2)][[1]])
		infile2 <- as.character(mytargets[,infile2])	
		if(nchar(infile2[1]) > 0) infile2 <- normalizePath(infile2)
		argname <- arglist$infile2[grep("<.*>", arglist$infile2)[1] -1]
		path <- arglist$infile2[grep("path", arglist$infile2)[1] +1]
		infile2back <- paste(path, infile2, sep="")
		names(infile2back) <- as.character(mytargets$SampleName)	
		infile2 <- paste(argname, " ", path, infile2, sep="")
		arglist[["infile2"]] <- gsub("(^ {1,})|( ${1,})", "", infile2)
	} else {
		infile2back <- rep("", length(mytargets[,1])) 
		infile2 <- infile2back
		names(infile2back) <- as.character(mytargets$SampleName)	
		arglist[["infile2"]] <- infile2back
	}
		
	## Populate arglist$outfile1
	outfile1 <- gsub("<|>", "", arglist$outfile1[grepl("^<.*>$", arglist$outfile1)][[1]])
	outfile1 <- as.character(mytargets[,outfile1])
	outfile1 <- gsub("^.*/", "", outfile1) 	
	remove <- arglist$outfile1[grep("remove", arglist$outfile1)[1] +1]
	outfile1 <- gsub(as.character(remove), "", outfile1)
	outfile1back <- outfile1
	outpaths <- outfile1
	outextension <- as.character(arglist$outfile1[grep("outextension", arglist$outfile1)+1])
	append <- arglist$outfile1[grep("append", arglist$outfile1)[1] +1]
	outfile1 <- paste(outfile1, append, sep="")
	argname <- arglist$outfile1[grep("<.*>", arglist$outfile1)[1] -1]
	path <- arglist$outfile1[grep("path", arglist$outfile1)[1] +1]
	path <- gsub("^\\./|^/|/$", "", path)
	resultpath <- paste(getwd(), "/", path, "/", sep="")	
	outfile1back <- paste(getwd(), "/", path, "/", outfile1, sep="")
	names(outfile1back) <- as.character(mytargets$SampleName)	
	outfile1 <- paste(argname, " ", getwd(), "/", path, "/", outfile1, sep="")
	arglist[["outfile1"]] <- gsub("(^ {1,})|( ${1,})", "", outfile1)

	## Populate arglist$outfile2 if it exists (usually only required for PE trimming)
	if("outfile2" %in% names(arglist)) {
		outfile2 <- gsub("<|>", "", arglist$outfile2[grepl("^<.*>$", arglist$outfile2)][[1]])
		outfile2 <- as.character(mytargets[,outfile2])
		outfile2 <- gsub("^.*/", "", outfile2) 	
		remove2 <- arglist$outfile2[grep("remove", arglist$outfile2)[1] +1]
		outfile2 <- gsub(as.character(remove2), "", outfile2)
		outfile2back <- outfile2
		outpaths2 <- outfile2
		outextension2 <- as.character(arglist$outfile2[grep("outextension", arglist$outfile2)+1])
		append2 <- arglist$outfile2[grep("append", arglist$outfile2)[1] +1]
		outfile2 <- paste(outfile2, append2, sep="")
		argname2 <- arglist$outfile2[grep("<.*>", arglist$outfile2)[1] -1]
		path2 <- arglist$outfile2[grep("path", arglist$outfile2)[1] +1]
		path2 <- gsub("^\\./|^/|/$", "", path2)
		resultpath2 <- paste(getwd(), "/", path2, "/", sep="")	
		outfile2back <- paste(getwd(), "/", path2, "/", outfile2, sep="")
		names(outfile2back) <- as.character(mytargets$SampleName)	
		outfile2 <- paste(argname2, " ", getwd(), "/", path2, "/", outfile2, sep="")
		arglist[["outfile2"]] <- gsub("(^ {1,})|( ${1,})", "", outfile2)
	}

    ## Generate arglist$outpaths
	outpaths <- paste(getwd(), "/", path, "/", outpaths, outextension, sep="")
	names(outpaths) <- as.character(mytargets$SampleName)	

	## Generate targetsout
	targetsout <- mytargetsorig
	targetsout[,1] <- outpaths[as.character(targetsout$SampleName)]
	if("outfile2" %in% names(arglist)) { 
		outpaths2 <- paste(getwd(), "/", path2, "/", outpaths2, outextension2, sep="")
		names(outpaths2) <- as.character(mytargets$SampleName)	
		targetsout[,2] <- outpaths2[as.character(targetsout$SampleName)]
		arglist <- arglist[!names(arglist) %in% "outfile2"]
	} else {
		colnames(targetsout)[1] <- "FileName"
		targetsout <- targetsout[ ,!colnames(targetsout) %in% "FileName2"]
	}

	## Collapse remaining components to single string vectors
	remaining <- names(arglist)[!names(arglist) %in% c("outfile1", "infile1", "infile2", "outpaths")]
	for(i in remaining) arglist[[i]] <- rep(gsub("(^ {1,})|( ${1,})", "", paste(arglist[[i]], collapse=" ")), length(arglist$infile1))
	args <- do.call("cbind", arglist)	
	rownames(args) <- as.character(mytargets$SampleName)
	args <- apply(args, 1, paste, collapse=" ")
	if(software=="bash_commands") { # If command-line is series of bash commands
		args <- gsub("' {1,}| {1,}'", "'", args)
		args <- gsub("bash_commands {1,}", "", args)
	}
	
	## When software is R-based then system commands make no sense and NA is used instead
	if(iscommandline==FALSE) args[] <- "" 

	## If sysma=NULL then make adjustments that are most reasonable
	if(length(sysmapath)==0) {
		targetsout <- mytargetsorig
		modules <- ""; software <- "R functions"; other=""; reference=""; resultpath=""
		outfile1back <- infile1back; args[] <- ""; outpaths[] <- infile1back
	}
		
	## Construct SYSargs object from components
	syslist <- list(targetsin=mytargetsorig,
					targetsout=targetsout,
					targetsheader=targetsheader,
					modules=modules, 
					software=software, 
					cores=cores,
					other=other,
					reference=reference,
					results=resultpath,
					infile1=infile1back,
					infile2=infile2back,
					outfile1=outfile1back,
					sysargs=args, 
					outpaths=outpaths)
	sysargs <- as(syslist, "SYSargs")
	if(type=="SYSargs") return(sysargs)
}

## Usage:
# args <- systemArgs(sysma="../inst/extdata/tophat.param", mytargets="../inst/extdata/targets.txt")
# names(args); modules(args); cores(args); outpaths(args); sysargs(args)

###########################################################
## Additional utilities for SYSargs and SYSargs2 objects ##
###########################################################
## Convenience write function for targetsout(args)
writeTargetsout <- function (x, file = "default", silent = FALSE, overwrite = FALSE, step = NULL, ...) {
  if(all(class(x)!="SYSargs" & class(x)!="SYSargs2")) stop("Argument 'x' needs to be assigned an object of class 'SYSargs' OR 'SYSargs2")
  ## SYSargs class
  if(class(x)=="SYSargs") {
    targets <- targetsout(x)
    software <- software(x)
    if (file == "default") {
      file <- paste("targets_", software, ".txt", sep = "")
      file <- gsub(" {1,}", "_", file)
    }
    else {
      file <- file
    }
    headerlines <- targetsheader(x)
    ## SYSargs2 class
  } else if (class(x)=="SYSargs2") {
    if(is.null(step)) stop(paste("Argument 'step' needs to be assigned one of the following values:", paste(names(x$clt), collapse=", "), "OR the corresponding position"))
    if(all(!is.null(step) & is.character(step) & !any(names(x$clt) %in% step))) stop(paste("Argument 'step' can only be assigned one of the following values:", paste(names(x$clt), collapse=", "), "OR the corresponding position")) 
    if(all(!is.null(step) & is.numeric(step) & !any(seq_along(names(x$clt)) %in% step))) stop(paste("Argument 'step' can only be assigned one of the following position:", paste(seq_along(names(x$clt)), collapse=", "), "OR the corresponding names")) 
    targets <- targets.as.df(targets(x))
    if(any(colnames(targets) %in% "FileName")){
      for(i in which(colnames(targets) %in% "FileName")){
        pout <- sapply(names(output(x)), function(y) paste(getwd(), "/", output(x)[[y]][[1]][[1]], sep=""), simplify = F) 
        targets[,i] <- as.character(pout)
      }
    } 
    if(any(colnames(targets) %in% "FileName1")){
      for(i in which(colnames(targets) %in% "FileName1")){
        p1out <- sapply(names(output(x)), function(y) paste(getwd(), "/", output(x)[[y]][[1]][[1]], sep=""), simplify = F) 
        targets[,i] <- as.character(p1out)
      }
    } 
    if(any(colnames(targets) %in% "FileName2")){
      for(i in which(colnames(targets) %in% "FileName2")){
        p2out <- sapply(names(output(x)), function(y) paste(getwd(), "/", output(x)[[y]][[1]][[2]], sep=""), simplify = F) 
        targets[,i] <- as.character(p2out)
      }
    }
    ## Workflow and Step Name
    software <- strsplit(basename(cwlfiles(x)$cwl), split="\\.")[[1]][1]
    if(is.character(step)) {
      step <- strsplit(step, split="\\.")[[1]][1]
    } else {
      step <- strsplit(names(x$clt)[step], split="\\.")[[1]][1]
    }
    if (file == "default") {
      file <- paste("targets_", step, ".txt", sep = "")
      file <- gsub(" {1,}", "_", file)
    } else {
      file <- file
    }
    headerlines <- targetsheader(x)[[1]]
  }
  if(file.exists(file) & overwrite==FALSE) stop(paste("I am not allowed to overwrite files; please delete existing file:", file, "or set 'overwrite=TRUE'"))
  targetslines <- c(paste(colnames(targets), collapse="\t"), apply(targets, 1, paste, collapse="\t"))
  writeLines(c(headerlines, targetslines), file, ...)
  if(silent!=TRUE) cat("\t", "Written content of 'targetsout(x)' to file:", file, "\n")
}

## Usage:
# writeTargetsout(x=args, file="default") ## SYSargs class
# writeTargetsout(x=WF, file="default", step=1) ## SYSargs2 class

##############################################################################
## Function to run NGS aligners including sorting and indexing of BAM files ##
##############################################################################
runCommandline <- function(args, runid="01", make_bam=TRUE, dir=FALSE, dir.name=NULL, ...) {
  if(any(nchar(gsub(" {1,}", "", modules(args))) > 0)) {
    ## Check if "Environment Modules" is installed in the system
    ## "Environment Modules" is not available
    if(suppressWarnings(system("modulecmd bash -V", ignore.stderr = TRUE, ignore.stdout = TRUE))!=1) {
      message("Environment Modules is not available. Please make sure to configure your PATH environment variable according to the software in use.", "\n")
    } else {
      ## "Environment Modules" is available and proceed the module load
      if(suppressWarnings(system("modulecmd bash -V", ignore.stderr = TRUE, ignore.stdout = TRUE))==1) { # Returns TRUE if module system is present.
        for(j in modules(args)) moduleload(j) # loads specified software from module system
      }
    }
  }
  ## SYSargs class
  if(class(args)=="SYSargs") { 
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
            unlink(outfile1(args)[i])
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
    ## SYSargs2 class
  } else if(class(args)=="SYSargs2") {
    ## Workflow Name
    cwl.wf <- strsplit(basename(cwlfiles(args)$cwl), split="\\.")[[1]][1]
    ## Folder name provide in the yml file or in the dir.name
    if(is.null(args$yamlinput$results_path$path)) {
      if(is.null(dir.name)) {
        stop("argument 'dir.name' missing. The argument can only be assigned 'NULL' when directory name is provided in the yml template. The argument should be assigned as a character vector of length 1")
      } else { logdir <- dir.name }
    } else { logdir <- normalizePath(args$yamlinput$results_path$path) 
    }
    args.return <- args
    ## Check what expected outputs have been generated
    if(make_bam==FALSE){
      completed <- output(args)
      outputList <- as.character()
      for(i in seq_along(output(args))){
        for(j in seq_along(output(args)[[i]])){
          completed[[i]][[j]] <- file.exists(output(args)[[i]][[j]])
          names(completed[[i]][[j]]) <- output(args)[[i]][[j]]
          outputList <- c(outputList, output(args)[[i]][[j]])
        }
      }
      names(outputList) <- rep(names(output(args)), each=length(output(args)[[1]][[1]]))
    } else if(make_bam==TRUE) {
      if(any(grepl("samtools", names(clt(args))))){ stop("argument 'make_bam' should be 'FALSE' when using the workflow with 'SAMtools'")} 
      args1 <- output_update(args, dir=FALSE, replace=TRUE, extension=c(".sam", ".bam"))
      completed <- output(args1)
      outputList <- as.character()
      for(i in seq_along(output(args1))){
        for(j in seq_along(output(args1)[[i]])){
          completed[[i]][[j]] <- file.exists(output(args1)[[i]][[j]])
          names(completed[[i]][[j]]) <- output(args1)[[i]][[j]]
          outputList <- c(outputList, output(args1)[[i]][[j]])
          if(any(grepl(".bam", output(args1)[[i]][[j]]))){
            for(k in which(grepl(".bam", output(args1)[[i]][[j]]))){
              outputList <- c(outputList, paste0(gsub("\\.bam$", "", output(args1)[[i]][[j]][k]), ".bam.bai"))
            }
          }
        }
      }
      names(outputList) <- rep(names(output(args)), each=length(output(args)[[1]][[1]])+1)
      args.return <- output_update(args.return, dir=FALSE, replace=TRUE, extension=c(".sam", ".bam"))
    }
    for(i in seq_along(cmdlist(args))){
      for(j in seq_along(cmdlist(args)[[i]])){
        ## Run the commandline only for samples for which no outpu file is available.
        if(all(as.logical(completed[[i]][[j]]))) {
          next()
        } else {
          # Create soubmitargsID_command file
          cat(cmdlist(args)[[i]][[j]], file=paste(logdir, "/submitargs", runid, "_", cwl.wf, sep=""), fill=TRUE, labels=paste0(names(cmdlist(args))[[i]], ":"), append=TRUE)
          ## Run executable 
          command <- gsub(" .*", "", as.character(cmdlist(args)[[i]][[j]]))
          commandargs <- gsub("^.*? ", "",as.character(cmdlist(args)[[i]][[j]]))
          if(command %in% "bwa") {
            stdout <- system2(command, args=commandargs, stdout=TRUE, stderr=FALSE)
          } else {
            stdout <- system2(command, args=commandargs, stdout=TRUE, stderr=TRUE)
          }
          ## Create submitargsID_stdout file
          cat(cmdlist(args)[[i]][[j]], file=paste(logdir, "/submitargs", runid, "_", cwl.wf, "_log", sep=""), fill=TRUE, labels=paste0(names(cmdlist(args))[[i]], ":"), sep = "\n", append=TRUE)
          cat(unlist(stdout), file=paste(logdir, "/submitargs", runid, "_", cwl.wf, "_log", sep=""), sep = "\n", append=TRUE)
        }
        cat("################", file=paste(logdir, "/submitargs", runid, "_", cwl.wf, "_log", sep=""), sep = "\n", append=TRUE)
        if(make_bam==TRUE) {
          sam_files <- grepl(".sam$", output(args)[[i]][[j]])
          others_files <- grepl("vcf$|bcf$|xls$|bed$", output(args)[[i]][[j]])
          completed.bam <- grepl(".bam$", output(args)[[i]][[j]])
          if(any(sam_files)){
            for(k in which(sam_files)){
              Rsamtools::asBam(file=output(args)[[i]][[j]], destination=gsub("\\.sam$", "", output(args)[[i]][[j]]), overwrite=TRUE, indexDestination=TRUE)
              unlink(output(args)[[i]][[k]])
            } } else if(any(others_files)){
              dump <- "do nothing"
            }
          if(any(completed.bam)){ # If output is unindexed *.bam file (e.g. Tophat2)
            for(k in which(completed.bam)){
              Rsamtools::sortBam(file=output(args)[[i]][[j]][k], destination=gsub("\\.bam$", "", output(args)[[i]][[j]][k]))
              Rsamtools::indexBam(output(args)[[i]][[j]][k]) 
            }
          }
        }
      } 
    }
    ## Create recursive the subfolders
    if(dir==TRUE){
      for(i in seq_along(names(cmdlist(args)))){
        full_path <- paste0(logdir, "/", cwl.wf, "/", names(cmdlist(args)[i]))
        if(dir.exists(full_path)==FALSE){
          dir.create(full_path, recursive = TRUE) }
      }
      if(dir.exists(paste0(logdir, "/", cwl.wf, "/_logs/"))==FALSE){
        dir.create(paste0(logdir, "/", cwl.wf, "/_logs/"), recursive = TRUE) }
      # 
      files_log <- list.files(path=logdir, pattern = "submitargs")
      for(i in seq_along(files_log)){
        file.rename(from=paste0(logdir, "/", files_log[i]), to=paste0(logdir, "/", cwl.wf, "/_logs/", files_log[i]))
      }
      outputList_new <- as.character()
      for(i in seq_along(outputList)){
        if(file.exists(outputList[i])){
          name <- strsplit(outputList[i], split="\\/")[[1]]
          name <- name[length(name)]
          file.rename(from=outputList[i], to=paste0(logdir, "/", cwl.wf, "/", names(outputList[i]), "/", name))
          outputList_new <- c(outputList_new, paste0(logdir, "/", cwl.wf, "/", names(outputList[i]), "/", name))
        } else if(file.exists(outputList[i])==FALSE){
          dump <- "No such file or directory"
        }
      }
      outputList <- outputList_new
      args.return <- output_update(args.return, dir=TRUE, replace=FALSE)
    }
    output_completed <- as.character()
    for(i in seq_along(outputList)){
      output_completed[i] <- file.exists(outputList[i])
    }
    names(output_completed) <- outputList
    cat("Missing expected outputs files:", sum(!as.logical(output_completed)), "\n"); cat("Existing expected outputs files:", sum(as.logical(output_completed)), "\n")
    print(output_completed)
    return(args.return)
    #return(output_completed)
  }
}

## Usage:
# WF <- runCommandline(WF) # creates the files in the ./results folder
# WF <- runCommandline(WF, dir=TRUE) # creates the files in the ./results/workflowName/Samplename folder
# WF <- runCommandline(WF, make_bam = FALSE, dir=TRUE) ## For hisat2-mapping.cwl template

#########################
## Old: qsub Arguments ##
#########################
getQsubargs <- function(software="qsub", queue="batch", Nnodes="nodes=1", cores=as.numeric(gsub("^.* ", "", tophatargs$args["p"])), memory="mem=10gb", time="walltime=20:00:00") {
	.Deprecated("clusterRun")
	qsubargs <- list(software=software, 
			queue=queue, 
                	Nnodes=Nnodes, 
                 	cores=cores, 
                 	memory=memory, 
                 	time=time)
	return(qsubargs)
}
## Usage:
# qsubargs <- getQsubargs(queue="batch", Nnodes="nodes=1", cores=cores(tophat), memory="mem=10gb", time="walltime=20:00:00")

#######################################################################################
## Old: custom function to submit runCommandline jobs to queuing system of a cluster ##
#######################################################################################
qsubRun <- function(appfct="runCommandline(args=args, runid='01')", args, qsubargs, Nqsubs=1, package="systemPipeR", shebang="#!/bin/bash") {
	.Deprecated("clusterRun")
	args2 <- sysargs(args)
	mydir <- getwd()
	setwd(results(args))
	splitvector <- sort(rep_len(1:Nqsubs, length.out=length(args2)))
	commands <- split(args2, splitvector)
	qsub_command <- paste(qsubargs$software, " -q ", qsubargs$queue, " -l ", qsubargs$Nnodes, ":ppn=", qsubargs$cores, ",", qsubargs$memory, ",", qsubargs$time, sep="")	
	jobids <- NULL
	for(i in 1:Nqsubs) {
		args2 <- commands[[i]]
		counter <- formatC(i, width = 2, format = "d", flag = "0")
		appfct <- gsub("runid.*)", paste("runid=", "'", counter, "'", ")", sep=""), appfct) # Passes on proper runid
		splitargs <- args[splitvector==i]
		save(splitargs, file=paste("submitargs", counter, sep="")) 
		rscript <- c(paste("library('", package, "')", sep=""), paste("load('submitargs", counter, "')", sep=""), "args <- splitargs", appfct)
		writeLines(rscript, paste("submitargs", counter, ".R", sep=""))
		writeLines(c(shebang, "cd $PBS_O_WORKDIR", paste("Rscript --verbose submitargs", counter, ".R", sep="")), paste("submitargs", counter, ".sh", sep=""))
		myqsub <- paste(qsub_command, paste("submitargs", counter, ".sh", sep=""))
		(jobids <- c(jobids, system(myqsub, intern=TRUE)))
	}
	setwd(mydir)
	names(commands) <- jobids
	return(commands)
}
## Usage:
# qsubRun(args=args, qsubargs=qsubargs, Nqsubs=1, package="systemPipeR")

###########################################################################################
## batchtools-based function to submit runCommandline jobs to queuing system of a cluster ##
###########################################################################################
## The advantage of this function is that it should work with most queuing/scheduling systems such as SLURM, Troque, SGE, ...
clusterRun <- function(args, FUN = runCommandline, more.args=list(args=args, make_bam=TRUE), conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl", Njobs, runid = "01", resourceList) {
  ## Validity checks of inputs
  if(any(class(args)!="SYSargs" & class(args)!="SYSargs2")) stop("Argument 'args' needs to be assigned an object of class 'SYSargs' OR 'SYSargs2'")
  if(class(FUN)!="function") stop("Value assigned to 'FUN' argument is not an object of class function.")
  if(!file.exists(conffile)) stop("Need to point under 'conffile' argument to proper config file. See more information here: https://mllg.github.io/batchtools/reference/makeRegistry.html. Note: in this file *.tmpl needs to point to a valid template file.")
  if(!file.exists(template)) stop("Need to point under 'template' argument to proper template file. Sample template files for different schedulers are available here: https://github.com/mllg/batchtools/blob/master/inst/templates/")
  if(!class(more.args)=="list") stop("'more.args' needs to be object of class 'list'.")    
  if(any(!names(more.args) %in% names(as.list(formals(FUN))))) stop(paste("The list of arguments assigned to 'more.args' can only be the following arguments defined in the function 'FUN':", paste(names(as.list(formals(FUN))), collapse=", "))) 
  ## SYSargs class
  if(class(args)=="SYSargs") {
    path <- normalizePath(results(args))
    args.f <- seq(along = args)
    ## SYSargs2 class
  } else if (class(args)=="SYSargs2") {
    path <- normalizePath(args$yamlinput$results_path$path)
    args.f <- seq(along=cmdlist(args))
  }
  ## batchtools routines
  f <- function(i, args, ...) FUN(args=args[i], ...)
  logdir1 <- paste0(path, "/submitargs", runid, "_btdb_", paste(sample(0:9, 4), collapse = ""))
  reg <- makeRegistry(file.dir = logdir1, conf.file = conffile, packages = "systemPipeR")
  ids <- batchMap(fun = f, args.f, more.args = more.args, reg=reg)
  chunk <- chunk(ids$job.id, n.chunks = Njobs, shuffle = FALSE)
  ids$chunk <- chunk
  done <- submitJobs(ids=ids, reg=reg, resources = resourceList)
  return(reg)
}
## Usage: 
# resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024) 
# reg <- clusterRun(args, conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl", Njobs=18, runid="01", resourceList=resources)
# getStatus(reg=reg)  
# waitForJobs(reg=reg)