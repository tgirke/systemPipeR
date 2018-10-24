########################
## Read preprocessing ##
########################
preprocessReads <- function(args, Fct, batchsize=100000, overwrite=TRUE, ...) {
	if(class(args)!="SYSargs") stop("Argument 'args' needs to be of class SYSargs")
	if(class(Fct)!="character") stop("Argument 'Fct' needs to be of class character")
	## Run function in loop over all fastq files
	## Single end fastq files
	if(!all(c("FileName1", "FileName2") %in% colnames(targetsout(args)))) {
    	for(i in seq(along=args)) {
			outfile <- outpaths(args)[i]
			## Delete existing fastq files with same names, since writeFastq will append to them
			if(overwrite==TRUE) {
				if(any(file.exists(outfile))) unlink(outfile)
			} else {
				if(any(file.exists(outfile))) stop(paste("File", outfile , "exists. Please delete file first or set overwrite=TRUE."))
			}
			## Run preprocessor function with FastqStreamer
			counter <- 0
			f <- FastqStreamer(infile1(args)[i], batchsize)
			while(length(fq <- yield(f))) {
				fqtrim <- eval(parse(text=Fct))
				writeFastq(fqtrim, outfile, mode="a", ...)
				counter <- counter + length(fqtrim)
				cat(counter, "processed reads written to file:", outfile, "\n")
			}
			close(f)
		}
	}
	## Paired end fastq files
	if(all(c("FileName1", "FileName2") %in% colnames(targetsout(args)))) {
    	for(i in seq(along=args)) {
			p1 <- as.character(targetsin(args)$FileName1[i])
			p2 <- as.character(targetsin(args)$FileName2[i])
			p1out <- as.character(targetsout(args)$FileName1[i])
			p2out <- as.character(targetsout(args)$FileName2[i])
			## Delete existing fastq files with same names, since writeFastq will append to them
			if(overwrite==TRUE) {
				if(any(file.exists(p1out))) unlink(p1out)
				if(any(file.exists(p2out))) unlink(p2out)
			} else {
				if(any(file.exists(p1out))) stop(paste("File", p1out , "exists. Please delete file first or set overwrite=TRUE."))
				if(any(file.exists(p2out))) stop(paste("File", p2out , "exists. Please delete file first or set overwrite=TRUE."))
			}
			## Run preprocessor function with FastqStreamer
			counter1 <- 0
			counter2 <- 0
			f1 <- FastqStreamer(p1, batchsize)
			f2 <- FastqStreamer(p2, batchsize)
			while(length(fq1 <- yield(f1))) {
				fq2 <- yield(f2)
				if(length(fq1)!=length(fq2)) stop("Paired end files cannot have different read numbers.")
				## Process p1
				fq <- fq1 # for simplicity in eval
				fq1trim <- eval(parse(text=Fct))
				## Index for p1
				index1 <- as.character(id(fq1)) %in% as.character(id(fq1trim))
				names(index1) <- seq(along=index1)
				index1 <- names(index1[index1])
				## Process p2
				fq <- fq2 # for simplicity in eval
				fq2trim <- eval(parse(text=Fct))
				## Index for p1
				index2 <- as.character(id(fq2)) %in% as.character(id(fq2trim))
				names(index2) <- seq(along=index2)
				index2 <- names(index2[index2])
				## Export to processed paired files
				indexpair1 <- index1 %in% index2
				writeFastq(fq1trim[indexpair1], p1out, mode="a", ...)
				indexpair2 <- index2 %in% index1
				writeFastq(fq2trim[indexpair2], p2out, mode="a", ...)
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
symLink2bam <- function(sysargs, command="ln -s", htmldir, ext=c(".bam", ".bai"), urlbase, urlfile) {
	## Create URL file 
	bampaths <- outpaths(sysargs)
	symname <- SampleName(sysargs)
	urls <- paste(urlbase, htmldir[2], symname, ext[1], "\t", symname, sep="")
	writeLines(urls, urlfile)
	## Creat correspoding sym links
	dir.create(paste(htmldir, collapse=""))
	symname <- rep(symname, each=2)
	symname <- paste(symname, c(ext[1], paste(ext, collapse="")), sep="")
	bampaths2 <- as.vector(t(cbind(bampaths, paste(bampaths, ext[2], sep=""))))	
	symcommands <- paste(command, " ", bampaths2, " ", paste(htmldir, collapse=""), symname, sep="") 
	for(i in symcommands) system(i)
}
## Usage: 
# symLink2bam(sysargs=args, command="ln -s", htmldir=c("~/.html/", "somedir/"), ext=c(".bam", ".bai"), urlbase="http://biocluster.ucr.edu/~tgirke/", urlfile="IGVurl.txt") 

#####################
## Alignment Stats ##
#####################
alignStats <- function(args) {
	fqpaths <- infile1(args)
	bampaths <- outpaths(args)
	bamexists <- file.exists(bampaths)
	fqpaths <- fqpaths[bamexists]
	bampaths <- bampaths[bamexists]
	## Obtain total read number from FASTQ files
	Nreads <- countLines(fqpaths)/4
	names(Nreads) <- names(fqpaths)		
	## If reads are PE multiply by 2 as a rough approximation
	if(nchar(infile2(args))[1] > 0) Nreads <- Nreads * 2	
	## Obtain total number of alignments from BAM files
	bfl <- BamFileList(bampaths, yieldSize=50000, index=character())
	param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE))
	Nalign <- countBam(bfl, param=param)
	## Obtain number of primary alignments from BAM files
	param <- ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE, isUnmappedQuery=FALSE))
	Nalignprim <- countBam(bfl, param=param)
	statsDF <- data.frame(FileName=names(Nreads), 
                              Nreads=Nreads, 
                              Nalign=Nalign$records, 
                              Perc_Aligned=Nalign$records/Nreads*100, 
                              Nalign_Primary=Nalignprim$records, 
                              Perc_Aligned_Primary=Nalignprim$records/Nreads*100
	)
	if(nchar(infile2(args))[1] > 0) colnames(statsDF)[which(colnames(statsDF)=="Nreads")] <- "Nreads2x"
	return(statsDF)
}
## Usage:
# read_statsDF <- alignStats(args=args)

########################
## RPKM Normalization ##
########################
returnRPKM <- function(counts, ranges) {
        geneLengthsInKB <- sum(width(reduce(ranges)))/1000 # Length of exon union per gene in kbp
        millionsMapped <- sum(counts)/1e+06 # Factor for converting to million of mapped reads.
        rpm <- counts/millionsMapped # RPK: reads per kilobase of exon model.
        rpkm <- rpm/geneLengthsInKB # RPKM: reads per kilobase of exon model per million mapped reads.
        return(rpkm)
}

## Usage:
# countDFrpkm <- apply(countDF, 2, function(x) returnRPKM(counts=x, ranges=eByg))

###############################################
## Read Sample Comparisons from Targets File ##
###############################################
## Parses sample comparisons from <CMP> line(s) in targets.txt file or SYSars object. 
## All possible comparisons can be specified with 'CMPset: ALL'.
readComp <- function(file, format="vector", delim="-") {
	if(!format %in% c("vector", "matrix")) stop("Argument format can only be assigned: vector or matrix!")
	## Parse <CMP> line
	if(class(file)=="SYSargs") {
		if(length(targetsheader(file))==0) stop("Input has no targets header lines.")
		comp <- targetsheader(file)
	} else {
		comp <- readLines(file)
	}
	comp <- comp[grepl("<CMP>", comp)]
	comp <- gsub("#.*<CMP>| {1,}", "", comp)
	comp <- gsub("\t", "", comp); comp <- gsub("^\"|\"$", "", comp) # Often required if Excel is used for editing targets file
	comp <- strsplit(comp, ":|,")
	names(comp) <- lapply(seq(along=comp), function(x) comp[[x]][1])	
	comp <- sapply(names(comp), function(x) comp[[x]][-1], simplify=FALSE)	
	
	## Check whether all samples are present in Factor column of targets file
	checkvalues <- unique(unlist(strsplit(unlist(comp), "-")))
	checkvalues <- checkvalues[checkvalues!="ALL"]
	if(class(file)=="SYSargs") {
		all <- unique(as.character(targetsin(file)$Factor))
	} else {
		all <- unique(as.character(read.delim(file, comment.char = "#")$Factor))
	}
	if(any(!checkvalues %in% all)) stop(paste("The following samples are not present in Factor column of targets file:", paste(checkvalues[!checkvalues %in% all], collapse=", ")))	

	## Generate outputs 
	allindex <- sapply(names(comp), function(x) any(grepl("ALL", comp[[x]])))
	if(any(allindex)) for(i in which(allindex)) comp[[i]] <- combn(all, m=2, FUN=paste, collapse=delim)
	if(format == "vector" & delim != "-") comp <- sapply(names(comp), function(x) gsub("-", delim, comp[[x]]), simplify=FALSE)
	if(format == "vector") return(comp)
	if(format == "matrix") return(sapply(names(comp), function(x) do.call("rbind", strsplit(comp[[x]], "-")), simplify=FALSE))
}
## Usage:
# cmp <- readComp("targets.txt", format="vector", delim="-")
# cmp <- readComp(args, format="vector", delim="-")

#################################
## Access module system from R ##
#################################

## New module functions 

# S3 Class for handling function calls
myEnvModules <- structure(list(), class="EnvModules")

## Main function to allow avail, list and list
myEnvModules$init <- function(){
  # Module function assumes MODULEPATH and MODULEDIR are set in login profile
  # Get base environment from login profile
  base_env <- strsplit(system('bash -l -c "env"',intern = TRUE),'\n')
  base_env <- strsplit(as.character(base_env),'=')
  
  # Iterate through base environment
  for (x in seq(1,length(base_env))) {
    
    # Set environment based on login profile
    if (base_env[[x]][1]=="LOADEDMODULES" || base_env[[x]][1]=="MODULESHOME" || base_env[[x]][1]=="MODULEPATH" || base_env[[x]][1]=="MODULES_DIR" || base_env[[x]][1]=="IIGB_MODULES"){
      if (base_env[[x]][1]=="LOADEDMODULES"){
        default_modules <- strsplit(base_env[[x]][2],":")
      }
      else{
        l <- list(base_env[[x]][2])
        names(l) <- base_env[[x]][1]
        do.call(Sys.setenv, l)
      }
    }
  }
  
  # Make sure to process default modules after the environment is set with the above loop
  for (x in seq(1,length(default_modules[[1]]))){
    module_name <- default_modules[[1]][x]
    print(paste("Loading module",module_name))
    try(myEnvModules$load_unload("load",module_name))
  }
}

# Print available modules or currently loaded modules on stderr
myEnvModules$avail_list <- function(action_type){
  try(module_vars <- system(paste('modulecmd bash',action_type),intern = TRUE))
}

# Eventually we could support clearing all modules, however for now just print
myEnvModules$clear <- function(action_type){
  loaded_modules <-  strsplit(Sys.getenv("LOADEDMODULES"),":")
  if (length(loaded_modules[[1]]) > 0) {
    for (x in seq(1,length(loaded_modules[[1]]))){
      module_name <- loaded_modules[[1]][x]
      print(paste("Unloading module",module_name))
      try(myEnvModules$load_unload("unload",module_name))
    }
  }
}

# Load and unload actions are basically the same, set environment variables given by modulecmd
myEnvModules$load_unload <- function(action_type, module_name=""){
  # Use the low level C binary for generating module environment variables
  try(module_vars <- system(paste('modulecmd bash',action_type, module_name),intern = TRUE))
  
  if (length(module_vars) > 0){
    # Separate environment variables
    module_var <- strsplit(module_vars,";")
    
    # Iterate through all environment variables
    for (x in seq(1,length(module_var[[1]]))) {
      # Isolate key, value pair
      evar <- module_var[[1]][x]
      
      # Filter export commands
      if (length(grep('^export',evar)) == 0 && length(evar) > 0) {
        
        # Seprate key and value
        evar <- strsplit(as.character(evar),'=')
        # Stip spaces at the end of the value
        evar_val <- gsub('[[:space:]]','',evar[[1]][2])
        # Remove extra backslashes
        l <- list(gsub('\\$','',evar_val))
        
        # Unset variables that need to be unset
        if(length(grep("^unset ",evar[[1]][1])) > 0){
          evar <- gsub("^unset (.*)$","\\1",evar[[1]][1])
          Sys.unsetenv(evar)
        } 
        else {
          # Assign names to each value in list
          names(l) <- evar[[1]][1]
          # Set environment variable in current environment
          do.call(Sys.setenv, l)
        }
      }
    }
  }
}
# Define what happens bases on action
module <- function(action_type,module_name=""){
  # Check to see if modulecmd is in current PATH
  try(
    suppressWarnings(modulecmd_path <- system("which modulecmd",intern=TRUE,ignore.stderr=TRUE)),
    silent=TRUE
  )
  
  # Only initialize module system if it has not yet been initialized and the modulecmd exisits
  if ( Sys.getenv('MODULEPATH') == "" && length(modulecmd_path) > 0) {
    myEnvModules$init()
  } else if (Sys.getenv('MODULEPATH') == "" && length(modulecmd_path) == 0) {
    stop("Cound not find the installation of Environment Modules: \"modulecmd\"")
  }
  
  switch(action_type,
    "load"   = myEnvModules$load_unload(action_type,module_name),
    "unload" = myEnvModules$load_unload(action_type,module_name),
    "list"   = myEnvModules$avail_list(action_type),
    "avail"  = myEnvModules$avail_list(action_type),
    "clear"  = myEnvModules$clear(action_type),
    "init"   = myEnvModules$init(),
    stop("That action is not supported.")
  )
}
## Usage: 
# module("load","tophat")
# module("load","tophat/2.1.1")
# module("list")
# module("avail")
# module("init")
# module("clear")
# module("unload", "tophat")
# module("unload", "tophat/2.1.1")

## Old module functions 

## List software available in module system
modulelist <- function() {
	system("bash -lc \"module avail\"")
}

## New version of load software from module system
moduleload <- function(module, envir="PATH") { 
    for(i in seq_along(envir)) {
        modpath <- system(paste0("bash -lc \"module load ", module, "; export | grep '^declare -x ", envir[i], "='\""), intern = TRUE) 
        modpath <- gsub(paste0("^declare -x ", envir[i], "=\"(.*)\"$"), "\\1", modpath, perl = TRUE)
        l <- list(modpath); names(l) <- envir[i]
        do.call(Sys.setenv, l)
    }
}
## Usage: 
# moduleload(module="hisat2/2.0.1", envir="PATH)
# moduleload(module="python", envir=c("PATH", "LD_LIBRARY_PATH", "PYTHONPATH"))

## Old version of load software from module system
# moduleload <- function(module) { 
#	modpath <- system(paste("bash -c \"module load", module, "; export | grep '^declare -x PATH='\""), intern = TRUE) 
#	modpath <- gsub("^declare -x PATH=\"(.*)\"$", "\\1", modpath, perl = TRUE)
#	Sys.setenv(PATH=modpath) 
# }

#######################################################################
## Run edgeR GLM with entire count matrix or subsetted by comparison ##
#######################################################################
## If independent=TRUE then countDF will be subsetted for each comparison
run_edgeR <- function(countDF, targets, cmp, independent=TRUE, paired=NULL, mdsplot="") {
    if(class(cmp) != "matrix" & length(cmp)==2) cmp <- t(as.matrix(cmp)) # If cmp is vector of length 2, convert it to matrix.
    samples <- as.character(targets$Factor); names(samples) <- paste(as.character(targets$SampleName), "", sep="")
    countDF <- countDF[, names(samples)]
    countDF[is.na(countDF)] <- 0
    edgeDF <- data.frame(row.names=rownames(countDF))
    group <- as.character(samples)
    if(independent==TRUE) {
        loopv <- seq(along=cmp[,1])
    } else {
	loopv <- 1
    }
    for(j in loopv) {
	## Filtering and normalization
	y <- DGEList(counts=countDF, group=group) # Constructs DGEList object
	if(independent == TRUE) {
	    subset <- samples[samples %in% cmp[j,]]
	    y <- y[, names(subset)]
	    y$samples$group <- factor(as.character(y$samples$group))
        }
	keep <- rowSums(cpm(y)>1) >= 2; y <- y[keep, ]
	y <- calcNormFactors(y)
	## Design matrix
	if(length(paired)==0) {
		design <- model.matrix(~0+y$samples$group, data=y$samples)
		colnames(design) <- levels(y$samples$group)
	} else {
        	if(length(paired)>0 & independent==FALSE) stop("When providing values under 'paired' also set independent=TRUE")
		Subject <- factor(paired[samples %in% cmp[j,]]) # corrected Jun 2014 (won't change results)
		Treat <- y$samples$group
        	design <- model.matrix(~Subject+Treat)
		levels(design) <- levels(y$samples$group)
	}
        ## Estimate dispersion
	y <- estimateGLMCommonDisp(y, design, verbose=TRUE) # Estimates common dispersions
	y <- estimateGLMTrendedDisp(y, design) # Estimates trended dispersions
	y <- estimateGLMTagwiseDisp(y, design) # Estimates tagwise dispersions 
	fit <- glmFit(y, design) # Fits the negative binomial GLM for each tag and produces an object of class DGEGLM with some new components.
	## Contrast matrix is optional but makes anlysis more transparent
	if(independent == TRUE) {
		mycomp <- paste(cmp[j,1], cmp[j,2], sep="-")
	} else {
		mycomp <- paste(cmp[,1], cmp[,2], sep="-")
	}
	if(length(paired)==0) contrasts <- makeContrasts(contrasts=mycomp, levels=design)
	for(i in seq(along=mycomp)) {
	    if(length(paired)==0) {
            	lrt <- glmLRT(fit, contrast=contrasts[,i]) # Takes DGEGLM object and carries out the likelihood ratio test. 
	    } else {
                lrt <- glmLRT(fit) # No contrast matrix with paired design
            }
            deg <- as.data.frame(topTags(lrt, n=length(rownames(y))))
	    colnames(deg) <- paste(paste(mycomp[i], collapse="_"), colnames(deg), sep="_")
	    edgeDF <- cbind(edgeDF, deg[rownames(edgeDF),]) 
	}
	if(nchar(mdsplot)>0) {
		pdf(paste("./results/sample_MDS_", paste(unique(subset), collapse="-"), ".pdf", sep=""))
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
run_DESeq2 <- function(countDF, targets, cmp, independent=FALSE) {
    if(class(cmp) != "matrix" & length(cmp)==2) cmp <- t(as.matrix(cmp)) # If cmp is vector of length 2, convert it to matrix.
    samples <- as.character(targets$Factor); names(samples) <- paste(as.character(targets$SampleName), "", sep="")
    countDF <- countDF[, names(samples)]
    countDF[is.na(countDF)] <- 0
    deseqDF <- data.frame(row.names=rownames(countDF))
    if(independent==TRUE) {
	loopv <- seq(along=cmp[,1])
    } else {
	loopv <- 1
    }
    for(j in loopv) {
	if(independent==TRUE) {
	    ## Create subsetted DESeqDataSet object
	    subset <- samples[samples %in% cmp[j,]]
	    countDFsub <- countDF[, names(subset)]
	    dds <- DESeq2::DESeqDataSetFromMatrix(countData=as.matrix(countDFsub), colData=data.frame(condition=subset), design = ~ condition)
	    mycmp <- cmp[j, , drop=FALSE]	
	} else {
	    ## Create full DESeqDataSet object
	    dds <- DESeq2::DESeqDataSetFromMatrix(countData=as.matrix(countDF), colData=data.frame(condition=samples), design = ~ condition)
	    mycmp <- cmp
	}
	## Estimate of (i) size factors, (ii) dispersion, (iii) negative binomial GLM fitting and (iv) Wald statistics
	dds <- DESeq2::DESeq(dds, quiet=TRUE)
	for(i in seq(along=mycmp[,1])) { 
	    ## Extracts DEG results for specific contrasts from DESeqDataSet object 
	    res <- DESeq2::results(dds, contrast=c("condition", mycmp[i,])) 
	    ## Set NAs to reasonable values to avoid errors in downstream filtering steps
	    res[is.na(res[,"padj"]), "padj"] <- 1
	    res[is.na(res[,"log2FoldChange"]), "log2FoldChange"] <- 0
	    deg <- as.data.frame(res)	
	    colnames(deg)[colnames(deg) %in% c("log2FoldChange", "padj")] <- c("logFC", "FDR")
	    colnames(deg) <- paste(paste(mycmp[i,], collapse="-"), colnames(deg), sep="_")
	    deseqDF <- cbind(deseqDF, deg[rownames(deseqDF),]) 
	}
    }
    return(deseqDF)
} 
## Usage:
# cmp <- readComp(file=targetspath, format="matrix", delim="-")
# degseqDF <- run_DESeq2(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=TRUE)

############################################
## Filter DEGs by p-value and fold change ##
############################################
filterDEGs <- function(degDF, filter, plot=TRUE) {
	pval <- degDF[, grep("_FDR$", colnames(degDF)), drop=FALSE]
	log2FC <- degDF[, grep("_logFC$", colnames(degDF)), drop=FALSE]
	## DEGs that are up or down regulated 
	pf <- pval <= filter["FDR"]/100 & (log2FC >= log2(filter["Fold"]) | log2FC <= -log2(filter["Fold"]))
	colnames(pf) <- gsub("_FDR", "", colnames(pf))
	pf[is.na(pf)] <- FALSE
	DEGlistUPorDOWN <- sapply(colnames(pf), function(x) rownames(pf[pf[,x,drop=FALSE],,drop=FALSE]), simplify=FALSE)
	## DEGs that are up regulated 
	pf <- pval <= filter["FDR"]/100 & log2FC >= log2(filter["Fold"])
	colnames(pf) <- gsub("_FDR", "", colnames(pf))
	pf[is.na(pf)] <- FALSE
	DEGlistUP <- sapply(colnames(pf), function(x) rownames(pf[pf[,x,drop=FALSE],,drop=FALSE]), simplify=FALSE)
	## DEGs that are down regulated 
	pf <- pval <= filter["FDR"]/100 & log2FC <= -log2(filter["Fold"])
	colnames(pf) <- gsub("_FDR", "", colnames(pf))
	pf[is.na(pf)] <- FALSE
	DEGlistDOWN <- sapply(colnames(pf), function(x) rownames(pf[pf[,x,drop=FALSE],,drop=FALSE]), simplify=FALSE)
	df <- data.frame(Comparisons=names(DEGlistUPorDOWN), Counts_Up_or_Down=sapply(DEGlistUPorDOWN, length), Counts_Up=sapply(DEGlistUP, length), Counts_Down=sapply(DEGlistDOWN, length))
	resultlist <- list(UporDown=DEGlistUPorDOWN, Up=DEGlistUP, Down=DEGlistDOWN, Summary=df)
	if(plot==TRUE) {
		mytitle <- paste("DEG Counts (", names(filter)[1], ": ", filter[1], " & " , names(filter)[2], ": ", filter[2], "%)", sep="")
		df_plot <- data.frame(Comparisons=rep(as.character(df$Comparisons), 2), Counts=c(df$Counts_Up, df$Counts_Down), Type=rep(c("Up", "Down"), each=length(df[,1])))
		p <- ggplot(df_plot, aes(Comparisons, Counts, fill = Type)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.y=element_text(angle=0, hjust=1)) + ggtitle(mytitle)
		print(p)
	}
	return(resultlist)
}
## Usage:
# DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=1))
