## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")))

## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE-------------------
suppressPackageStartupMessages({
    library(systemPipeR)
    library(BiocParallel)
    library(Biostrings)
    library(Rsamtools)
    library(GenomicRanges)
    library(ggplot2)
    library(GenomicAlignments)
    library(ShortRead)
    library(ape)
})

## ----genChip_workflow, eval=FALSE----------------------------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="chipseq")
## setwd("chipseq")

## Rscript -e "systemPipeRdata::genWorkenvir(workflow='chipseq')"

## ----node_environment, eval=FALSE----------------------------------------
## q("no") # closes R session on head node
## srun --x11 --partition=short --mem=2gb --cpus-per-task 4 --ntasks 1 --time 2:00:00 --pty bash -l
## module load R/3.3.0
## R

## ----r_environment, eval=FALSE-------------------------------------------
## system("hostname") # should return name of a compute node starting with i or c
## getwd() # checks current working directory of R session
## dir() # returns content of current working directory

## ----load_systempiper, eval=TRUE-----------------------------------------
library(systemPipeR)

## ----load_custom_fct, eval=FALSE-----------------------------------------
## source("systemPipeChIPseq_Fct.R")

## ----load_targets_file, eval=TRUE----------------------------------------
targetspath <- system.file("extdata", "targets_chip.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")
targets[1:4,-c(5,6)]

## ----proprocess_reads, eval=FALSE, messages=FALSE, warning=FALSE, cache=TRUE----
## args <- systemArgs(sysma="param/trim.param", mytargets="targets_chip.txt")
## filterFct <- function(fq, cutoff=20, Nexceptions=0) {
##     qcount <- rowSums(as(quality(fq), "matrix") <= cutoff)
##     fq[qcount <= Nexceptions] # Retains reads where Phred scores are >= cutoff with N exceptions
## }
## preprocessReads(args=args, Fct="filterFct(fq, cutoff=20, Nexceptions=0)", batchsize=100000)
## writeTargetsout(x=args, file="targets_chip_trim.txt", overwrite=TRUE)

## ----fastq_report, eval=FALSE--------------------------------------------
## args <- systemArgs(sysma="param/tophat.param", mytargets="targets_chip.txt")
## library(BiocParallel); library(BatchJobs)
## f <- function(x) {
##     library(systemPipeR)
##     args <- systemArgs(sysma="param/tophat.param", mytargets="targets_chip.txt")
##     seeFastq(fastq=infile1(args)[x], batchsize=100000, klength=8)
## }
## funs <- makeClusterFunctionsSLURM("slurm.tmpl")
## param <- BatchJobsParam(length(args), resources=list(walltime="00:20:00", ntasks=1, ncpus=1, memory="2G"), cluster.functions=funs)
## register(param)
## fqlist <- bplapply(seq(along=args), f)
## pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
## seeFastqPlot(unlist(fqlist, recursive=FALSE))
## dev.off()

## ----bowtie2_align, eval=FALSE-------------------------------------------
## args <- systemArgs(sysma="param/bowtieSE.param", mytargets="targets_chip_trim.txt")
## sysargs(args)[1] # Command-line parameters for first FASTQ file
## moduleload(modules(args)) # Skip if a module system is not used
## system("bowtie2-build ./data/tair10.fasta ./data/tair10.fasta") # Indexes reference genome
## resources <- list(walltime="1:00:00", ntasks=1, ncpus=cores(args), memory="10G")
## reg <- clusterRun(args, conffile=".BatchJobs.R", template="slurm.tmpl", Njobs=18, runid="01",
##                   resourceList=resources)
## waitForJobs(reg)
## writeTargetsout(x=args, file="targets_bam.txt", overwrite=TRUE)

## ----bowtie2_align_seq, eval=FALSE---------------------------------------
## runCommandline(args)

## ----check_files_exist, eval=FALSE---------------------------------------
## file.exists(outpaths(args))

## ----align_stats, eval=FALSE---------------------------------------------
## read_statsDF <- alignStats(args=args)
## write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
## read.delim("results/alignStats.xls")

## ----symbol_links, eval=FALSE--------------------------------------------
## symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"),
##             urlbase="http://biocluster.ucr.edu/~tgirke/",
##             urlfile="./results/IGVurl.txt")

## ----rle_object, eval=FALSE----------------------------------------------
## library(rtracklayer); library(GenomicRanges); library(Rsamtools); library(GenomicAlignments)
## aligns <- readGAlignments(outpaths(args)[1])
## cov <- coverage(aligns)
## cov

## ----resize_align, eval=FALSE--------------------------------------------
## trim(resize(as(aligns, "GRanges"), width = 200))

## ----rle_slice, eval=FALSE-----------------------------------------------
## islands <- slice(cov, lower = 15)
## islands[[1]]

## ----plot_coverage, eval=FALSE-------------------------------------------
## library(ggbio)
## myloc <- c("Chr1", 1, 100000)
## ga <- readGAlignments(outpaths(args)[1], use.names=TRUE, param=ScanBamParam(which=GRanges(myloc[1], IRanges(as.numeric(myloc[2]), as.numeric(myloc[3])))))
## autoplot(ga, aes(color = strand, fill = strand), facets = strand ~ seqnames, stat = "coverage")

## ----merge_bams, eval=FALSE----------------------------------------------
## args <- systemArgs(sysma=NULL, mytargets="targets_bam.txt")
## args_merge <- mergeBamByFactor(args, overwrite=TRUE)
## writeTargetsout(x=args_merge, file="targets_mergeBamByFactor.txt", overwrite=TRUE)

## ----call_peaks_macs_noref, eval=FALSE-----------------------------------
## args <- systemArgs(sysma="param/macs2_noinput.param", mytargets="targets_mergeBamByFactor.txt")
## sysargs(args)[1] # Command-line parameters for first FASTQ file
## runCommandline(args)
## file.exists(outpaths(args))
## writeTargetsout(x=args, file="targets_macs.txt", overwrite=TRUE)

## ----call_peaks_macs_withref, eval=FALSE---------------------------------
## writeTargetsRef(infile="targets_mergeBamByFactor.txt", outfile="targets_bam_ref.txt", silent=FALSE, overwrite=TRUE)
## args_input <- systemArgs(sysma="param/macs2.param", mytargets="targets_bam_ref.txt")
## sysargs(args_input)[1] # Command-line parameters for first FASTQ file
## runCommandline(args_input)
## file.exists(outpaths(args_input))
## writeTargetsout(x=args_input, file="targets_macs_input.txt", overwrite=TRUE)

## ----consensus_peaks, eval=FALSE-----------------------------------------
## source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/rangeoverlapper.R")
## peak_M1A <- outpaths(args)["M1A"]
## peak_M1A <- as(read.delim(peak_M1A, comment="#")[,1:3], "GRanges")
## peak_A1A <- outpaths(args)["A1A"]
## peak_A1A <- as(read.delim(peak_A1A, comment="#")[,1:3], "GRanges")
## (myol1 <- subsetByOverlaps(peak_M1A, peak_A1A, minoverlap=1)) # Returns any overlap
## myol2 <- olRanges(query=peak_M1A, subject=peak_A1A, output="gr") # Returns any overlap with OL length information
## myol2[values(myol2)["OLpercQ"][,1]>=50] # Returns only query peaks with a minimum overlap of 50%

## ----chip_peak_anno, eval=FALSE------------------------------------------
## library(ChIPpeakAnno); library(GenomicFeatures)
## args <- systemArgs(sysma="param/annotate_peaks.param", mytargets="targets_macs.txt")
## # txdb <- loadDb("./data/tair10.sqlite")
## txdb <- makeTxDbFromGFF(file="data/tair10.gff", format="gff", dataSource="TAIR", organism="Arabidopsis thaliana")
## ge <- genes(txdb, columns=c("tx_name", "gene_id", "tx_type"))
## for(i in seq(along=args)) {
##     peaksGR <- as(read.delim(infile1(args)[i], comment="#"), "GRanges")
##     annotatedPeak <- annotatePeakInBatch(peaksGR, AnnotationData=genes(txdb))
##     df <- data.frame(as.data.frame(annotatedPeak), as.data.frame(values(ge[values(annotatedPeak)$feature,])))
##     write.table(df, outpaths(args[i]), quote=FALSE, row.names=FALSE, sep="\t")
## }
## writeTargetsout(x=args, file="targets_peakanno.txt", overwrite=TRUE)

## ----chip_peak_anno_full_annotation, include=FALSE, eval=FALSE-----------
## ## Perform previous step with full genome annotation from Biomart
## # txdb <- makeTxDbFromBiomart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host="plants.ensembl.org")
## # tx <- transcripts(txdb, columns=c("tx_name", "gene_id", "tx_type"))
## # ge <- genes(txdb, columns=c("tx_name", "gene_id", "tx_type")) # works as well
## # seqlevels(ge) <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrC", "ChrM")
## # table(mcols(tx)$tx_type)
## # tx <- tx[!duplicated(unstrsplit(values(tx)$gene_id, sep=","))] # Keeps only first transcript model for each gene]
## # annotatedPeak <- annotatePeakInBatch(peaksGR, AnnotationData = tx)

## ----chip_peak_seeker, eval=FALSE----------------------------------------
## library(ChIPseeker)
## for(i in seq(along=args)) {
##     peakAnno <- annotatePeak(infile1(args)[i], TxDb=txdb, verbose=FALSE)
##     df <- as.data.frame(peakAnno)
##     write.table(df, outpaths(args[i]), quote=FALSE, row.names=FALSE, sep="\t")
## }
## writeTargetsout(x=args, file="targets_peakanno.txt", overwrite=TRUE)

## ----chip_peak_seeker_plots, eval=FALSE----------------------------------
## peak <- readPeakFile(infile1(args)[1])
## covplot(peak, weightCol="X.log10.pvalue.")
## peakHeatmap(outpaths(args)[1], TxDb=txdb, upstream=1000, downstream=1000, color="red")
## plotAvgProf2(outpaths(args)[1], TxDb=txdb, upstream=1000, downstream=1000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

## ----count_peak_ranges, eval=FALSE---------------------------------------
## library(GenomicRanges)
## args <- systemArgs(sysma="param/count_rangesets.param", mytargets="targets_macs.txt")
## args_bam <- systemArgs(sysma=NULL, mytargets="targets_bam.txt")
## bfl <- BamFileList(outpaths(args_bam), yieldSize=50000, index=character())
## countDFnames <- countRangeset(bfl, args, mode="Union", ignore.strand=TRUE)
## writeTargetsout(x=args, file="targets_countDF.txt", overwrite=TRUE)

## ----diff_bind_analysis, eval=FALSE--------------------------------------
## args_diff <- systemArgs(sysma="param/rundiff.param", mytargets="targets_countDF.txt")
## cmp <- readComp(file=args_bam, format="matrix")
## dbrlist <- runDiff(args=args_diff, diffFct=run_edgeR, targets=targetsin(args_bam),
##                     cmp=cmp[[1]], independent=TRUE, dbrfilter=c(Fold=2, FDR=1))
## writeTargetsout(x=args_diff, file="targets_rundiff.txt", overwrite=TRUE)

## ----go_enrich, eval=FALSE-----------------------------------------------
## args <- systemArgs(sysma="param/macs2.param", mytargets="targets_bam_ref.txt")
## args_anno <- systemArgs(sysma="param/annotate_peaks.param", mytargets="targets_macs.txt")
## annofiles <- outpaths(args_anno)
## gene_ids <- sapply(names(annofiles), function(x) unique(as.character(read.delim(annofiles[x])[,"geneId"])), simplify=FALSE)
## load("data/GO/catdb.RData")
## BatchResult <- GOCluster_Report(catdb=catdb, setlist=gene_ids, method="all", id_type="gene", CLSZ=2, cutoff=0.9, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)

## ----parse_peak_sequences, eval=FALSE------------------------------------
## library(Biostrings); library(seqLogo); library(BCRANK)
## args <- systemArgs(sysma="param/annotate_peaks.param", mytargets="targets_macs.txt")
## rangefiles <- infile1(args)
## for(i in seq(along=rangefiles)) {
##     df <- read.delim(rangefiles[i], comment="#")
##     peaks <- as(df, "GRanges")
##     names(peaks) <- paste0(as.character(seqnames(peaks)), "_", start(peaks), "-", end(peaks))
##     peaks <- peaks[order(values(peaks)$X.log10.pvalue, decreasing=TRUE)]
##     pseq <- getSeq(FaFile("./data/tair10.fasta"), peaks)
##     names(pseq) <- names(peaks)
##     writeXStringSet(pseq, paste0(rangefiles[i], ".fasta"))
## }

## ----bcrank_enrich, eval=FALSE-------------------------------------------
## set.seed(0)
## BCRANKout <- bcrank(paste0(rangefiles[1], ".fasta"), restarts=25, use.P1=TRUE, use.P2=TRUE)
## toptable(BCRANKout)
## topMotif <- toptable(BCRANKout, 1)
## weightMatrix <- pwm(topMotif, normalize = FALSE)
## weightMatrixNormalized <- pwm(topMotif, normalize = TRUE)
## pdf("results/seqlogo.pdf")
## seqLogo(weightMatrixNormalized)
## dev.off()

## ----sessionInfo---------------------------------------------------------
sessionInfo()

