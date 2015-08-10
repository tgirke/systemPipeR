#######################################################################
## Merge BAM files based on factor and return updated SYSargs object ##
#######################################################################
## Useful prior to peak calling in ChIP-Seq or miRNA gene prediction experiments
## where pooling of replicates maximizes depth of read coverage.
## Note: default factor is targetsin(args)$Factor
mergeBamByFactor <- function(args, mergefactor=targetsin(args)$Factor, overwrite=FALSE, silent=FALSE, ...) {
    ## Check validity of input
    allbam <- !grepl("\\.bam$|\\.BAM", infile1(args))
    if(any(allbam)) stop("The following files lack the extension '.bam': ", paste(basename(infile1(args)[allbam]), collapse=", "))
    if(colnames(targetsin(args)[, 1, drop=FALSE]) != "FileName") stop("Name of first column in 'targetsin(arg)' is expected to be 'FileName'.")

    ## Unique values in Factor column of targetsin(args) 
    sample_index <- !duplicated(as.character(mergefactor)); names(sample_index) <- names(infile1(args))
    uniorder <- unique(mergefactor) 
    unifreq <- table(mergefactor)[uniorder] # Important: in order of targetsin(args)!
    if(!any(unifreq >= 2)) warning("Values in Factor column are all unique. Thus, there are no BAM files to merge.")

    ## Store BAM file paths in a list that have >=2 identical values (replicates) in Factor column of targets file
    filelist <- tapply(targetsin(args)$FileName, factor(mergefactor), as.character)
    filelist <- filelist[names(unifreq)] # Assures original order

    ## Create vector containing paths of output BAM files 
    filelist_merge <- filelist[names(unifreq[unifreq>=2])] # Merge only those with >= files
    outname_vec <- character(length(filelist_merge)); names(outname_vec) <- names(filelist_merge)
    for(i in seq(along=outname_vec)) {
        outname <- gsub("\\.bam$|\\.BAM$", "", filelist_merge[[i]][1])
        outname <- paste(outname, "_", names(filelist_merge[i]), ".bam", sep="")
        outname_vec[i] <- outname
    }
    
    ## If any output BAM file exists and 'overwrite=FALSE' then stop
    file_exists <- file.exists(outname_vec); names(file_exists) <- names(outname_vec) 
    if(any(file_exists) & overwrite==FALSE) stop("The following files exist: ", paste(names(file_exists)[file_exists], collapse=", ") , ". Delete/rename them or set 'overwrite=TRUE'")
    
    ## Generate collapsed BAM files
    for(i in seq(along=filelist_merge)) {
        mergeBam(filelist_merge[[i]], outname_vec[i], indexDestination=TRUE, overwrite=overwrite, ...)
        if(silent!=TRUE) {
            cat("Merged BAM files:", basename(filelist_merge[[i]]), "and saved to file", basename(outname_vec[i]), "\n\n")
        }
    }
    
    ## Generate updated SYSargs object 
    filelist[names(outname_vec)] <- outname_vec # Assign new file names to proper slots
    outfile_names <- unlist(filelist) 
	args_sub <- args[sample_index]
    targets_out <- targetsout(args_sub)
    targets_out[,"FileName"] <- outfile_names
    rownames(targets_out) <- NULL
    syslist <- list(targetsin=targetsin(args_sub),
					targetsout=targets_out,
					targetsheader=targetsheader(args_sub),
					modules=modules(args_sub), 
					software="mergeBamByFactor", 
					cores=cores(args_sub),
					other=other(args_sub),
					reference=reference(args_sub),
					results=results(args_sub),
					infile1=infile1(args_sub),
					infile2=infile2(args_sub),
					outfile1=outfile_names,
					sysargs=sysargs(args_sub), 
					outpaths=outfile_names)
	args_sub_out <- as(syslist, "SYSargs")
	return(args_sub_out)
}
## Usage:
# args <- systemArgs(sysma=NULL, mytargets="targets_bam.txt")
# args_merge <- mergeBamByFactor(args, overwrite=TRUE, silent=FALSE)
# writeTargetsout(x=args_merge, file="targets_mergeBamByFactor.txt", overwrite=TRUE) 

############################################################################
## Creat targets file with refence sample, e.g. input sample for ChIP-Seq ##
############################################################################
writeTargetsRef <- function(infile, outfile, silent=FALSE, overwrite=FALSE, ...) {
    ## Import
    headerlines <- readLines(infile)
    targets <- read.delim(infile, comment="#")
    
    ## Check for expected input
    if(!c("SampleReference") %in% colnames(targets)) stop("Targets file lacks SampleReference column")
    if(!c("FileName") %in% colnames(targets)) stop("Targets file lacks FileName column")
    if(all(c("FileName1", "FileName2") %in% colnames(targets))) stop("Targets file is expected to have only one FileName column")
    if(file.exists(outfile) & overwrite==FALSE) stop(paste("I am not allowed to overwrite files; please delete existing file:", outfile, "or set 'overwrite=TRUE'")) 
    testv <- as.character(targets$SampleReference); testv <- testv[!is.na(testv)]; testv <- testv[testv!=""] 
    myfiles <- as.character(targets$FileName); names(myfiles) <- as.character(targets$SampleName)
    if(!all(testv %in% names(myfiles))) stop(paste("Value(s)", paste(testv[!testv %in% names(myfiles)], collapse=", "), "from SampleReference column have no matches in SampleName column!"))
    
    ## Rearrange targets file
    targets <- data.frame(FileName1=targets$FileName, FileName2=NA, targets[,2:length(targets[1,])])
    targets[,"FileName2"] <- myfiles[as.character(targets$SampleReference)]
    targets <- targets[!is.na(as.character(targets$SampleReference)), , drop=FALSE] 
    targets <- targets[targets$SampleReference!="", , drop=FALSE] 
    
    ## Export targets file including header lines
    headerlines <- headerlines[grepl("^#", headerlines)]
    targetslines <- c(paste(colnames(targets), collapse="\t"), apply(targets, 1, paste, collapse="\t"))
    writeLines(c(headerlines, targetslines), outfile, ...)
    if(silent!=TRUE) cat("\t", "Modified", infile, "file with sample-wise reference has been written to outfile", outfile, "\n") 
}
## Usage:
# writeTargetsRef(infile="~/targets.txt", outfile="~/targets_refsample.txt", silent=FALSE, overwrite=FALSE)

