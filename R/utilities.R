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
	Nalign <- countBam(bfl)
	## Obtain number of primary alignments from BAM files
	param <- ScanBamParam(flag=scanBamFlag(isNotPrimaryRead=FALSE, isUnmappedQuery=FALSE))
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
## Parses sample comparisons from <CMP> line(s) in targets.txt file. All possible
## comparisons can be specified with 'CMPset: ALL'.
readComp <- function(file, format="vector", delim="-") {
	if(!format %in% c("vector", "matrix")) stop("Argument format can only be assigned: vector or matrix!")
	## Parse <CMP> line
	comp <- readLines(file)
	comp <- comp[grepl("<CMP>", comp)]
	comp <- gsub("#.*<CMP>| {1,}", "", comp)
	comp <- strsplit(comp, ":|,")
	names(comp) <- lapply(seq(along=comp), function(x) comp[[x]][1])	
	comp <- sapply(names(comp), function(x) comp[[x]][-1], simplify=FALSE)	
	
	## Check whether all samples are present in Factor column of targets file
	checkvalues <- unique(unlist(strsplit(unlist(comp), "-")))
	checkvalues <- checkvalues[checkvalues!="ALL"]
	all <- unique(as.character(read.delim(file, comment.char = "#")$Factor))
	if(any(!checkvalues %in% all)) stop(paste("The following samples are not present in Factor column of targets file:", paste(checkvalues[!checkvalues %in% all], collapse=", ")))	

	## Generate outputs 
	allindex <- sapply(names(comp), function(x) any(grepl("ALL", comp[[x]])))
	if(any(allindex)) for(i in which(allindex)) comp[[i]] <- combn(all, m=2, FUN=paste, collapse=delim)
	if(format == "vector" & delim != "-") comp <- sapply(names(comp), function(x) gsub("-", delim, comp[[x]]), simplify=FALSE)
	if(format == "vector") return(comp)
	if(format == "matrix") return(sapply(names(comp), function(x) do.call("rbind", strsplit(comp[[x]], "-")), simplify=FALSE))
}
## Usage:
# cmp <- readComp(file="targets.txt", format="vector", delim="-")

#################################
## Access module system from R ##
#################################
## List software available in module system
modulelist <- function() {
	system("bash -c \"module avail\"")
}

## Load software from module system
moduleload <- function(module) { 
	modpath <- system(paste("bash -c \"module load", module, "; export | grep '^declare -x PATH='\""), intern = TRUE) 
	modpath <- gsub("^declare -x PATH=\"(.*)\"$", "\\1", modpath, perl = TRUE)
	Sys.setenv(PATH=modpath) 
}

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

## Filter DEGs by p-value and fold change
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
