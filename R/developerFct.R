##########################################
## Helper Functions not Exposed to User ##
##########################################

####################################################################
## Subset Reads by Mapping Regions to Create Mini FASTQ/BAM Files ##
####################################################################

## (A) Download FASTQ files of study SRP010938 from SRA at NCBI
.getSRAfastq <- function(sraid, maxreads) {
    moduleload("sratoolkit/2.5.0")
    system(paste("fastq-dump --split-files --gzip --maxSpotId", maxreads, sraid))
}
## Usage:
# library(systemPipeR)
# sraidv <- paste("SRR4460", 27:44, sep="")
# bplapply(sraidv, .getSRAfastq, maxreads="1000000000", BPPARAM = MulticoreParam(workers=4))
## Non-parallized download
# for(i in sraidv) .getSRAfastq(sraid=i, maxreads = "1000000000")

## (B) Align reads with systemPipeR/tophat against truncated TAIR10 reference (each chr truncated to 100kbp)
## (C) Extract 100,000 reads from each fastq file mapping to the truncated regions in reference
## Step (C) is performed with the following function
.subsetReadsByMappingRegion <- function(args) {
    pkg <- c("GenomicAlignments", "GenomeInfoDb", "IRanges")
    checkPkg(pkg, quietly = FALSE)
	chromosomelength <- 100000
	mydir <- getwd()
	setwd("./data/SRP010938_sub") # outdir
	for(i in seq(along=outpaths(args))) {
		fl <- outpaths(args)[i]
		si <- GenomeInfoDb::seqinfo(BamFile(fl))                                                                                                                                                                                                                 
		gr <- GenomicRanges::GRanges(GenomeInfoDb::seqnames(si), IRanges::IRanges(100, GenomeInfoDb::seqlengths(si)-100))                                                                                                                                                                                
		aligns <- GenomicAlignments::readGAlignments(fl, param=Rsamtools::ScanBamParam(which=gr), use.names=TRUE)
		keepids <- names(aligns[start(aligns) < chromosomelength]) # Return read ids mapping in first 100000 nucleotides of chromosomes
		myN <- sample(90000:100000, 1) # Keep number of random sampled reads between 90-100K
		keepids <- sample(unique(keepids), myN) # random sample x reads
		reads1 <- ShortRead::readFastq(infile1(args)[i]) # Reads in a FASTQ file from the current folder
		index <- gsub(" .*", "", as.vector(ShortRead::id(reads1))) %in% keepids
		reads1 <- reads1[index] # subset by keepids
		ShortRead::writeFastq(reads1, gsub("^.*/", "", infile1(args)[i]), full=TRUE) # writes ShortReadQ object to output file
		rm(reads1); gc() # Clean up memory
		reads2 <- ShortRead::readFastq(infile2(args)[i]) 
		index <- gsub(" .*", "", as.vector(ShortRead::id(reads2))) %in% keepids
		reads2 <- reads2[index] 
		ShortRead::writeFastq(reads2, gsub("^.*/", "", infile2(args)[i]), full=TRUE) 
		rm(reads2); gc() # Clean up memory
		print(i)
	}
	setwd(mydir)
}
## Usage:
# args <- systemArgs(sysma="tophat.param", mytargets="./data/SRP010938/targets.txt")
# .subsetReadsByMappingRegion(args=args)
