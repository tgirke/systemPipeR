## pre code {

## white-space: pre !important;

## overflow-x: scroll !important;

## word-break: keep-all !important;

## word-wrap: initial !important;

## }


## ----style, echo = FALSE, results = 'asis'----------------
BiocStyle::markdown()
options(width=60, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")), 
    tidy.opts=list(width.cutoff=60), tidy=TRUE)


## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE, eval=FALSE----
## suppressPackageStartupMessages({
##     library(systemPipeR)
##     library(BiocParallel)
##     library(Biostrings)
##     library(Rsamtools)
##     library(GenomicRanges)
##     library(ggplot2)
##     library(GenomicAlignments)
##     library(ShortRead)
##     library(ape)
##     library(batchtools)
## })


## ----genChip_workflow, eval=FALSE-------------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="chipseq")
## setwd("chipseq")


## Rscript -e "systemPipeRdata::genWorkenvir(workflow='chipseq')"


## ----closeR, eval=FALSE-----------------------------------
## q("no") # closes R session on head node


## srun --x11 --partition=short --mem=2gb --cpus-per-task 4 --ntasks 1 --time 2:00:00 --pty bash -l

## module load R/3.6.0

## R


## ----r_environment, eval=FALSE----------------------------
## system("hostname") # should return name of a compute node starting with i or c
## getwd() # checks current working directory of R session
## dir() # returns content of current working directory


## ----load_systempiper, eval=TRUE--------------------------
library(systemPipeR)


## ----load_custom_fct, eval=FALSE--------------------------
## source("systemPipeChIPseq_Fct.R")


## ----load_targets_file, eval=TRUE-------------------------
targetspath <- system.file("extdata", "targets_chip.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")
targets[1:4,-c(5,6)]


## ----construct_SYSargs2_trim-se, eval=FALSE---------------
## dir_path <- system.file("extdata/cwl/preprocessReads/trim-se", package="systemPipeR")
## trim <- loadWF(targets=targetspath, wf_file="trim-se.cwl", input_file="trim-se.yml", dir_path=dir_path)
## trim <- renderWF(trim, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
## trim
## output(trim)[1:2]


## ----proprocess_reads, eval=FALSE, messages=FALSE, warning=FALSE, cache=TRUE----
## # args <- systemArgs(sysma="param/trim.param", mytargets="targets_chip.txt")
## filterFct <- function(fq, cutoff=20, Nexceptions=0) {
##     qcount <- rowSums(as(quality(fq), "matrix") <= cutoff, na.rm=TRUE)
##     fq[qcount <= Nexceptions]
##     # Retains reads where Phred scores are >= cutoff with N exceptions
## }
## preprocessReads(args=trim, Fct="filterFct(fq, cutoff=20, Nexceptions=0)",
##                 batchsize=100000)
## writeTargetsout(x=trim, file="targets_chip_trim.txt", step = 1,
##                 new_col = "FileName", new_col_output_index = 1, overwrite = TRUE)


## ----fastq_report, eval=FALSE-----------------------------
## library(BiocParallel); library(batchtools)
## f <- function(x) {
##   library(systemPipeR)
##   targets <- system.file("extdata", "targets_chip.txt", package="systemPipeR")
##   dir_path <- system.file("extdata/cwl/preprocessReads/trim-se", package="systemPipeR")
##   trim <- loadWorkflow(targets=targets, wf_file="trim-se.cwl", input_file="trim-se.yml", dir_path=dir_path)
##   trim <- renderWF(trim, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
##   seeFastq(fastq=infile1(trim)[x], batchsize=100000, klength=8)
## }
## 
## resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024)
## param <- BatchtoolsParam(workers = 4, cluster = "slurm", template = "batchtools.slurm.tmpl", resources = resources)
## fqlist <- bplapply(seq(along=trim), f, BPPARAM = param)
## 
## pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
## seeFastqPlot(unlist(fqlist, recursive=FALSE))
## dev.off()


## ----bowtie2_index, eval=FALSE----------------------------
## dir_path <- system.file("extdata/cwl/bowtie2/bowtie2-idx", package="systemPipeR")
## idx <- loadWorkflow(targets=NULL, wf_file="bowtie2-index.cwl", input_file="bowtie2-index.yml", dir_path=dir_path)
## idx <- renderWF(idx)
## idx
## cmdlist(idx)
## 
## ## Run in single machine
## runCommandline(idx, make_bam = FALSE)


## ----bowtie2_align, eval=FALSE----------------------------
## targets <- system.file("extdata", "targets_chip.txt", package="systemPipeR")
## dir_path <- system.file("extdata/cwl/bowtie2/bowtie2-se", package="systemPipeR")
## args <- loadWF(targets = targets, wf_file = "bowtie2-mapping-se.cwl",
##     input_file = "bowtie2-mapping-se.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", SampleName = "_SampleName_"))
## args
## cmdlist(args)[1:2]
## output(args)[1:2]


## ----bowtie2_align_cluster, eval=FALSE--------------------
## moduleload(modules(args)) # Skip if a module system is not used
## resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024)
## reg <- clusterRun(args, FUN = runCommandline, more.args = list(args=args, dir = FALSE),
##     conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl", Njobs = 18, runid = "01", resourceList = resources)
## getStatus(reg=reg)
## waitForJobs(reg=reg)
## args <- output_update(args, dir=FALSE, replace=TRUE, extension=c(".sam", ".bam")) ## Updates the output(args) to the right location in the subfolders
## output(args)


## ----bowtie2_align_single, eval=FALSE---------------------
## args <- runCommandline(args, force=F)


## ----check_files_exist, eval=FALSE------------------------
## writeTargetsout(x=args, file="targets_bam.txt", step = 1,
##                 new_col = "FileName", new_col_output_index = 1, overwrite = TRUE)
## outpaths <- subsetWF(args , slot="output", subset=1, index=1)
## file.exists(outpaths)


## ----align_stats, eval=FALSE------------------------------
## read_statsDF <- alignStats(args=args)
## write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
## read.delim("results/alignStats.xls")


## ----symbol_links, eval=FALSE-----------------------------
## symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"),
##             urlbase="http://cluster.hpcc.ucr.edu/~tgirke/",
##             urlfile="./results/IGVurl.txt")


## ----rle_object, eval=FALSE-------------------------------
## library(rtracklayer); library(GenomicRanges)
## library(Rsamtools); library(GenomicAlignments)
## outpaths <- subsetWF(args , slot="output", subset=1, index=1)
## aligns <- readGAlignments(outpaths[1])
## cov <- coverage(aligns)
## cov


## ----resize_align, eval=FALSE-----------------------------
## trim(resize(as(aligns, "GRanges"), width = 200))


## ----rle_slice, eval=FALSE--------------------------------
## islands <- slice(cov, lower = 15)
## islands[[1]]


## ----plot_coverage, eval=FALSE----------------------------
## library(ggbio)
## myloc <- c("Chr1", 1, 100000)
## ga <- readGAlignments(outpaths[1], use.names=TRUE,
##                       param=ScanBamParam(which=GRanges(myloc[1],
##                         IRanges(as.numeric(myloc[2]), as.numeric(myloc[3])))))
## autoplot(ga, aes(color = strand, fill = strand), facets = strand ~ seqnames, stat = "coverage")


## ----merge_bams, eval=FALSE-------------------------------
## dir_path <- system.file("extdata/cwl/mergeBamByFactor", package="systemPipeR")
## args <- loadWF(targets = "targets_bam.txt", wf_file = "merge-bam.cwl",
##     input_file = "merge-bam.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName = "_BAM_PATH_", SampleName = "_SampleName_"))
## 
## args_merge <- mergeBamByFactor(args=args, overwrite=TRUE)
## writeTargetsout(x=args_merge, file="targets_mergeBamByFactor.txt", step = 1,
##                 new_col = "FileName", new_col_output_index = 1, overwrite = TRUE )


## ----call_peaks_macs_envVar_settings, eval=FALSE----------
## # Skip if a module system is not used
## module("list")
## module("unload", "miniconda2")
## module("load", "python/2.7.14") # Make sure to set up your enviroment variable for MACS2


## ----call_peaks_macs_noref, eval=FALSE--------------------
## dir_path <- system.file("extdata/cwl/MACS2/MACS2-noinput/", package="systemPipeR")
## args <- loadWF(targets = "targets_mergeBamByFactor.txt", wf_file = "macs2.cwl",
##     input_file = "macs2.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", SampleName = "_SampleName_"))
## 
## runCommandline(args, make_bam = FALSE, force=T)
## outpaths <- subsetWF(args, slot="output", subset=1, index=1)
## file.exists(outpaths)
## writeTargetsout(x=args, file="targets_macs.txt", step = 1,
##                 new_col = "FileName", new_col_output_index = 1, overwrite = TRUE )


## ----call_peaks_macs_withref, eval=FALSE------------------
## writeTargetsRef(infile="targets_mergeBamByFactor.txt",
##                 outfile="targets_bam_ref.txt", silent=FALSE, overwrite=TRUE)
## dir_path <- system.file("extdata/cwl/MACS2/MACS2-input/", package="systemPipeR")
## args_input <- loadWF(targets = "targets_bam_ref.txt", wf_file = "macs2-input.cwl",
##     input_file = "macs2.yml", dir_path = dir_path)
## args_input <- renderWF(args_input, inputvars = c(FileName1 = "_FASTQ_PATH1_", FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_"))
## cmdlist(args_input)[1]
## ## Run
## args_input <- runCommandline(args_input, make_bam = FALSE, force=T)
## outpaths_input <- subsetWF(args_input , slot="output", subset=1, index=1)
## file.exists(outpaths_input)
## writeTargetsout(x=args_input, file="targets_macs_input.txt", step = 1,
##                 new_col = "FileName", new_col_output_index = 1, overwrite = TRUE )


## ----consensus_peaks, eval=FALSE--------------------------
## # source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/rangeoverlapper.R")
## outpaths <- subsetWF(args , slot="output", subset=1, index=1) ## escolher um dos outputs index
## peak_M1A <- outpaths["M1A"]
## peak_M1A <- as(read.delim(peak_M1A, comment="#")[,1:3], "GRanges")
## peak_A1A <- outpaths["A1A"]
## peak_A1A <- as(read.delim(peak_A1A, comment="#")[,1:3], "GRanges")
## (myol1 <- subsetByOverlaps(peak_M1A, peak_A1A, minoverlap=1))
##             # Returns any overlap
## myol2 <- olRanges(query=peak_M1A, subject=peak_A1A, output="gr")
##             # Returns any overlap with OL length information
## myol2[values(myol2)["OLpercQ"][,1]>=50]
##             # Returns only query peaks with a minimum overlap of 50%


## ----chip_peak_anno, eval=FALSE---------------------------
## library(ChIPpeakAnno); library(GenomicFeatures)
## dir_path <- system.file("extdata/cwl/annotate_peaks", package="systemPipeR")
## args <- loadWF(targets = "targets_macs.txt", wf_file = "annotate-peaks.cwl",
##     input_file = "annotate-peaks.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", SampleName = "_SampleName_"))
## 
## txdb <- makeTxDbFromGFF(file="data/tair10.gff", format="gff", dataSource="TAIR",
##                         organism="Arabidopsis thaliana")
## ge <- genes(txdb, columns=c("tx_name", "gene_id", "tx_type"))
## for(i in seq(along=args)) {
##     peaksGR <- as(read.delim(infile1(args)[i], comment="#"), "GRanges")
##     annotatedPeak <- annotatePeakInBatch(peaksGR, AnnotationData=genes(txdb))
##     df <- data.frame(as.data.frame(annotatedPeak),
##                      as.data.frame(values(ge[values(annotatedPeak)$feature,])))
##     outpaths <- subsetWF(args , slot="output", subset=1, index=1)
##     write.table(df, outpaths[i], quote=FALSE, row.names=FALSE, sep="\t")
## }
## writeTargetsout(x=args, file="targets_peakanno.txt", step = 1,
##                 new_col = "FileName", new_col_output_index = 1, overwrite = TRUE )


## ----chip_peak_anno_full_annotation, include=FALSE, eval=FALSE----
## ## Perform previous step with full genome annotation from Biomart
## # txdb <- makeTxDbFromBiomart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host="plants.ensembl.org")
## # tx <- transcripts(txdb, columns=c("tx_name", "gene_id", "tx_type"))
## # ge <- genes(txdb, columns=c("tx_name", "gene_id", "tx_type")) # works as well
## # seqlevels(ge) <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrC", "ChrM")
## # table(mcols(tx)$tx_type)
## # tx <- tx[!duplicated(unstrsplit(values(tx)$gene_id, sep=","))] # Keeps only first transcript model for each gene]
## # annotatedPeak <- annotatePeakInBatch(peaksGR, AnnotationData = tx)


## ----chip_peak_seeker, eval=FALSE-------------------------
## library(ChIPseeker)
## for(i in seq(along=args)) {
##     peakAnno <- annotatePeak(infile1(args)[i], TxDb=txdb, verbose=FALSE)
##     df <- as.data.frame(peakAnno)
##     outpaths <- subsetWF(args , slot="output", subset=1, index=1)
##     write.table(df, outpaths[i], quote=FALSE, row.names=FALSE, sep="\t")
## }
## writeTargetsout(x=args, file="targets_peakanno.txt", step = 1,
##                 new_col = "FileName", new_col_output_index = 1, overwrite = TRUE )


## ----chip_peak_seeker_plots, eval=FALSE-------------------
## peak <- readPeakFile(infile1(args)[1])
## covplot(peak, weightCol="X.log10.pvalue.")
## outpaths <- subsetWF(args , slot="output", subset=1, index=1)
## peakHeatmap(outpaths[1], TxDb=txdb, upstream=1000, downstream=1000,
##             color="red")
## plotAvgProf2(outpaths[1], TxDb=txdb, upstream=1000, downstream=1000,
##              xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")


## ----count_peak_ranges, eval=FALSE------------------------
## library(GenomicRanges)
## dir_path <- system.file("extdata/cwl/count_rangesets", package="systemPipeR")
## args <- loadWF(targets = "targets_macs.txt", wf_file = "count_rangesets.cwl",
##     input_file = "count_rangesets.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", SampleName = "_SampleName_"))
## 
## ## Bam Files
## targets <- system.file("extdata", "targets_chip.txt", package="systemPipeR")
## dir_path <- system.file("extdata/cwl/bowtie2/bowtie2-se", package="systemPipeR")
## args_bam <- loadWF(targets = targets, wf_file = "bowtie2-mapping-se.cwl",
##     input_file = "bowtie2-mapping-se.yml", dir_path = dir_path)
## args_bam <- renderWF(args_bam, inputvars = c(FileName = "_FASTQ_PATH1_", SampleName = "_SampleName_"))
## args_bam <- output_update(args_bam, dir=FALSE, replace=TRUE, extension=c(".sam", ".bam"))
## outpaths <- subsetWF(args_bam, slot="output", subset=1, index=1)
## 
## bfl <- BamFileList(outpaths, yieldSize=50000, index=character())
## countDFnames <- countRangeset(bfl, args, mode="Union", ignore.strand=TRUE)
## writeTargetsout(x=args, file="targets_countDF.txt", step = 1,
##                 new_col = "FileName", new_col_output_index = 1, overwrite = TRUE )


## ----diff_bind_analysis, eval=FALSE-----------------------
## dir_path <- system.file("extdata/cwl/rundiff", package="systemPipeR")
## args_diff <- loadWF(targets = "targets_countDF.txt", wf_file = "rundiff.cwl",
##     input_file = "rundiff.yml", dir_path = dir_path)
## args_diff <- renderWF(args_diff, inputvars = c(FileName = "_FASTQ_PATH1_", SampleName = "_SampleName_"))
## 
## cmp <- readComp(file=args_bam, format="matrix")
## dbrlist <- runDiff(args=args_diff, diffFct=run_edgeR,
##                    targets=targets.as.df(targets(args_bam)), cmp=cmp[[1]],
##                    independent=TRUE, dbrfilter=c(Fold=2, FDR=1))
## writeTargetsout(x=args_diff, file="targets_rundiff.txt", step = 1,
##                 new_col = "FileName", new_col_output_index = 1, overwrite = TRUE )


## ----go_enrich, eval=FALSE--------------------------------
## dir_path <- system.file("extdata/cwl/annotate_peaks", package="systemPipeR")
## args <- loadWF(targets = "targets_bam_ref.txt", wf_file = "annotate-peaks.cwl",
##     input_file = "annotate-peaks.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName1 = "_FASTQ_PATH1_", FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_"))
## 
## args_anno <- loadWF(targets = "targets_macs.txt", wf_file = "annotate-peaks.cwl",
##     input_file = "annotate-peaks.yml", dir_path = dir_path)
## args_anno <- renderWF(args_anno, inputvars = c(FileName = "_FASTQ_PATH1_", SampleName = "_SampleName_"))
## annofiles <- subsetWF(args_anno, slot="output", subset=1, index=1)
## gene_ids <- sapply(names(annofiles),
##                    function(x) unique(as.character
##                     (read.delim(annofiles[x])[,"geneId"])), simplify=FALSE)
## load("data/GO/catdb.RData")
## BatchResult <- GOCluster_Report(catdb=catdb, setlist=gene_ids, method="all",
##                                 id_type="gene", CLSZ=2, cutoff=0.9,
##                                 gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)


## ----parse_peak_sequences, eval=FALSE---------------------
## library(Biostrings); library(seqLogo); library(BCRANK)
## dir_path <- system.file("extdata/cwl/annotate_peaks", package="systemPipeR")
## args <- loadWF(targets = "targets_macs.txt", wf_file = "annotate-peaks.cwl",
##     input_file = "annotate-peaks.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", SampleName = "_SampleName_"))
## 
## rangefiles <- infile1(args)
## for(i in seq(along=rangefiles)) {
##     df <- read.delim(rangefiles[i], comment="#")
##     peaks <- as(df, "GRanges")
##     names(peaks) <- paste0(as.character(seqnames(peaks)), "_", start(peaks),
##                            "-", end(peaks))
##     peaks <- peaks[order(values(peaks)$X.log10.pvalue., decreasing=TRUE)]
##     pseq <- getSeq(FaFile("./data/tair10.fasta"), peaks)
##     names(pseq) <- names(peaks)
##     writeXStringSet(pseq, paste0(rangefiles[i], ".fasta"))
## }


## ----bcrank_enrich, eval=FALSE----------------------------
## set.seed(0)
## BCRANKout <- bcrank(paste0(rangefiles[1], ".fasta"), restarts=25,
##                     use.P1=TRUE, use.P2=TRUE)
## toptable(BCRANKout)
## topMotif <- toptable(BCRANKout, 1)
## weightMatrix <- pwm(topMotif, normalize = FALSE)
## weightMatrixNormalized <- pwm(topMotif, normalize = TRUE)
## pdf("results/seqlogo.pdf")
## seqLogo(weightMatrixNormalized)
## dev.off()


## ----sessionInfo------------------------------------------
sessionInfo()

