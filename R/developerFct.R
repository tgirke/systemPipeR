##########################################
## Helper Functions not Exposed to User ##
##########################################

####################################################################
## Subset Reads by Mapping Regions to Create Mini FASTQ/BAM Files ##
####################################################################

## (A) Download FASTQ files of study SRP010938 from SRA at NCBI
.getSRAfastq <- function(sraid, maxreads) {
    system(paste("fastq-dump --defline-seq '@$sn[_$rn]/$ri' --gzip --split-files --maxSpotId", maxreads, sraid))
}
## Usage:
# library(systemPipeR)
# moduleload("sratoolkit/2.9.2")
# sraidv <- paste("SRR4460", 27:44, sep="")
# bplapply(sraidv, .getSRAfastq, maxreads="1000000000", BPPARAM = MulticoreParam(workers=18))
## Non-parallized download
# for(i in sraidv) .getSRAfastq(sraid=i, maxreads = "1000000000")

## (B) Align reads with systemPipeR/tophat against truncated TAIR10 reference (each chr truncated to 100kbp)
## (C) Extract 100,000 reads from each fastq file mapping to the truncated regions in reference
## Step (C) is performed with the following function
.subsetReadsByMappingRegion <- function(args, mydir = getwd(), 
                                        outdir = "./data/SRP010938_sub",
                                        chromosomelength =100000) {
    pkg <- c("GenomicAlignments", "GenomeInfoDb", "IRanges")
    checkPkg(pkg, quietly = FALSE)
    if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    setwd(outdir) # outdir
    for(i in seq(along=outpaths(args))) {
        fl <- outpaths(args)[i]
        si <- GenomeInfoDb::seqinfo(BamFile(fl))
        gr <- GenomicRanges::GRanges(GenomeInfoDb::seqnames(si), IRanges::IRanges(100, GenomeInfoDb::seqlengths(si)-100))
        aligns <- GenomicAlignments::readGAlignments(fl, param=Rsamtools::ScanBamParam(which=gr), use.names=TRUE)
        keepids <- names(aligns[start(aligns) < chromosomelength]) # Return read ids mapping in first 100000 nucleotides of chromosomes
        myN <- sample(90000:100000, 1) # Keep number of random sampled reads between 90-100K
        keepids <- sample(unique(keepids), myN, replace=TRUE) # random sample x reads
        keepids1 <- paste(keepids, "/1",sep="") 
        reads1 <- ShortRead::readFastq(infile1(args)[i]) # Reads in a FASTQ file from the current folder
        index <- gsub(" .*", "", as.vector(id(reads1))) %in% keepids1 
        length(index[index ==TRUE]) 
        reads1 <- reads1[index] # subset by keepids
        ShortRead::writeFastq(reads1, gsub("^.*/", "", infile1(args)[i]), full=TRUE) # writes ShortReadQ object to output file
        rm(reads1); gc() # Clean up memory
        reads2 <- ShortRead::readFastq(infile2(args)[i])
        keepids2 <- paste(keepids, "/2",sep="") 
        index <- gsub(" .*", "", as.vector(ShortRead::id(reads2))) %in% keepids2
        length(index[index ==TRUE]) 
        reads2 <- reads2[index] 
        ShortRead::writeFastq(reads2, gsub("^.*/", "", infile2(args)[i]), full=TRUE) 
        rm(reads2); gc() # Clean up memory
        print(i)
	}
	setwd(mydir)
}
## Usage:
# getwd()
# targets <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
# dir_path <- system.file("extdata/cwl", package="systemPipeR")
# args <- loadWorkflow(targets=targets, wf_file="hisat2/hisat2-mapping-pe.cwl", 
#                    input_file="hisat2/hisat2-mapping-pe.yml", dir_path=dir_path)
# args <- renderWF(args, inputvars = c(FileName1 = "_FASTQ_PATH1_", FileName2 = "_FASTQ_PATH2_", 
#                                  SampleName = "_SampleName_"))
# args
# 
# cmdlist(args)[1]
# 
# resources <- list(walltime = 1200, ntasks = 1, ncpus = 14, memory = 20480)
# reg <- clusterRun(args, FUN = runCommandline, more.args = list(dir =FALSE, make_bam = TRUE),
#    conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl",
#    Njobs = 18, runid = "01", resourceList = resources)
# batchtools::getStatus(reg = reg)
# 
# args <- output_update(args, replace=".bam", dir=FALSE)
# .subsetReadsByMappingRegion(args = args)  