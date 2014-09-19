##############################################
## Class and Method Definitions for SYSargs ##
##############################################
## Define SYSargs class
setClass("SYSargs", representation(modules="character", 
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

## Methods to return SYSargs as list
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
		new("SYSargs", modules=from$modules, 
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
        x@infile1 <- x@infile1[i]
        x@infile2 <- x@infile2[i]
        x@outfile1 <- x@outfile1[i]
        x@sysargs <- x@sysargs[i]
        x@outpaths <- x@outpaths[i]
        return(x)
})

## Construct SYSargs object from param and targets files
systemArgs <- function(sysma, mytargets, type="SYSargs") {
	## Read sysma and convert to arglist
	sysma <- as.matrix(read.delim(sysma, comment.char = "#"))
	sysma[is.na(sysma)] <- ""
	arglist <- sapply(as.character(unique(sysma[,"PairSet"])), function(x) as.vector(t(as.matrix(sysma[sysma[,"PairSet"]==x, 2:3]))))
	for(i in seq(along=arglist)) names(arglist[[i]]) <- paste(rep(c("n", "v"), length(arglist[[i]])/2), rep(1:(length(arglist[[i]])/2), 2), sep="")
	if(type=="json") return(toJSON(arglist))
	## Validity checks
	mytargets <- read.delim(mytargets, comment.char = "#")
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
	reference <- as.character(arglist$reference[grepl("v", names(arglist$reference))])
	if(!grepl("^/", reference)) reference <- paste0(getwd(), gsub("^\\.", "", reference)) # Turn relative into absolute path.
	arglist[["reference"]]["v1"] <- reference
	cores <- as.numeric(arglist$cores[grepl("v", names(arglist$cores))])	

	## Populate arglist$infile1
	infile1 <- gsub("<|>", "", arglist$infile1[grepl("^<.*>$", arglist$infile1)][[1]])
	infile1 <- as.character(mytargets[,infile1])	
	infile1 <- normalizePath(infile1)
	argname <- arglist$infile1[grep("<.*>", arglist$infile1)[1] -1]
	path <- arglist$infile1[grep("path", arglist$infile1)[1] +1]
	infile1back <- paste(path, infile1, sep="")
	names(infile1back) <- as.character(mytargets$SampleName)	
	infile1 <- paste(argname, " ", path, infile1, sep="")
	arglist[["infile1"]] <- gsub("(^ {1,})|( ${1,})", "", infile1)
	
	## Populate arglist$infile2
	infile2 <- gsub("<|>", "", arglist$infile2[grepl("^<.*>$", arglist$infile2)][[1]])
	infile2 <- as.character(mytargets[,infile2])	
	if(nchar(infile2[1]) > 0) infile2 <- normalizePath(infile2)
	argname <- arglist$infile2[grep("<.*>", arglist$infile2)[1] -1]
	path <- arglist$infile2[grep("path", arglist$infile2)[1] +1]
	infile2back <- paste(path, infile2, sep="")
	names(infile2back) <- as.character(mytargets$SampleName)	
	infile2 <- paste(argname, " ", path, infile2, sep="")
	arglist[["infile2"]] <- gsub("(^ {1,})|( ${1,})", "", infile2)
	
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
	path <- arglist$outfile1[grep("path", arglist$outfile)[1] +1]
	path <- gsub("^\\./|^/|/$", "", path)
	resultpath <- paste(getwd(), "/", path, "/", sep="")	
	outfile1back <- paste(getwd(), "/", path, "/", outfile1, sep="")
	names(outfile1back) <- as.character(mytargets$SampleName)	
	outfile1 <- paste(argname, " ", getwd(), "/", path, "/", outfile1, sep="")
	arglist[["outfile1"]] <- gsub("(^ {1,})|( ${1,})", "", outfile1)

	## Generate arglist$outpaths
	outpaths <- paste(getwd(), "/", path, "/", outpaths, outextension, sep="")
	names(outpaths) <- as.character(mytargets$SampleName)	

	## Collapse remaining components to single string vectors
	remaining <- names(arglist)[!names(arglist) %in% c("outfile1", "infile1", "infile2", "outpaths")]
	for(i in remaining) arglist[[i]] <- rep(gsub("(^ {1,})|( ${1,})", "", paste(arglist[[i]], collapse=" ")), length(arglist$infile1))
	args <- do.call("cbind", arglist)	
	rownames(args) <- as.character(mytargets$SampleName)
	args <- apply(args, 1, paste, collapse=" ")
	
	## Construct SYSargs object from components
	syslist <- list(modules=modules, 
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
# args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
# names(args); modules(args); cores(args); outpaths(args); sysargs(args)

##############################################################################
## Function to run NGS aligners including sorting and indexing of BAM files ##
##############################################################################
runCommandline <- function(args, runid="01", usemodule=TRUE) {
	if(usemodule==TRUE) {
		for(j in modules(args)) moduleload(j) # loads specified software from module system
	}	
	commands <- sysargs(args)
	completed <- file.exists(outpaths(args))
	names(completed) <- outpaths(args)
	resultdir <- results(args)
	for(i in seq(along=commands)) {
		## Run alignmets only for samples for which no BAM file is available.
		if(as.logical(completed)[i]) {
			next()
		} else {
			## Create submitrunID_log file
			cat(commands[i], file=paste(resultdir, "submitargs", runid, "_log", sep=""), sep = "\n", append=TRUE)
        		## Run executable  
			system(as.character(commands[i]))
			if(grepl(".sam$", outfile1(args)[i])) { # If output is *.sam file (e.g. Bowtie2)
				asBam(file=outfile1(args)[i], destination=gsub("\\.sam$", "", outfile1(args)[i]), overwrite=TRUE, indexDestination=TRUE)
				unlink(outfile1(args)[i])
			} else { # If output is unindexed *.bam file (e.g. Tophat2)
				sortBam(file=names(completed[i]), destination=gsub("\\.bam$", "", names(completed[i])))
        			indexBam(names(completed[i]))
			}
		}
	}
	bamcompleted <- gsub("sam$", "bam$", file.exists(outpaths(args)))
	names(bamcompleted) <- SampleName(args)
	cat("Missing alignment results (bam files):", sum(!as.logical(bamcompleted)), "\n"); cat("Existing alignment results (bam files):", sum(as.logical(bamcompleted)), "\n")
	return(bamcompleted)
}

## Usage: 
# runCommandline(args=args)

####################
## qsub Arguments ##
####################
getQsubargs <- function(software="qsub", queue="batch", Nnodes="nodes=1", cores=as.numeric(gsub("^.* ", "", tophatargs$args["p"])), memory="mem=10gb", time="walltime=20:00:00") {
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

###########################################################################
## Function to submit runCommandline jobs to queuing system of a cluster ##
###########################################################################
qsubRun <- function(appfct="runCommandline(args=args, runid='01')", args, qsubargs, Nqsubs=1, package="systemPipeR", shebang="#!/bin/bash") {
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



