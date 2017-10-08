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

## ----genVAR_workflow, eval=FALSE-----------------------------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="varseq")
## setwd("varseq")

## ----genVar_workflow_command_line, eval=FALSE, engine="sh"---------------
## Rscript -e "systemPipeRdata::genWorkenvir(workflow='varseq')"

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
## source("systemPipeVARseq_Fct.R")

## ----load_targets_file, eval=TRUE----------------------------------------
targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")
targets[,-c(5,6)]

## ----preprocess_reads, eval=FALSE----------------------------------------
## args <- systemArgs(sysma="param/trimPE.param", mytargets="targetsPE.txt")[1:4] # Note: subsetting!
## filterFct <- function(fq, cutoff=20, Nexceptions=0) {
##     qcount <- rowSums(as(quality(fq), "matrix") <= cutoff)
##     fq[qcount <= Nexceptions] # Retains reads where Phred scores are >= cutoff with N exceptions
## }
## preprocessReads(args=args, Fct="filterFct(fq, cutoff=20, Nexceptions=0)", batchsize=100000)
## writeTargetsout(x=args, file="targets_PEtrim.txt", overwrite=TRUE)

## ----fastq_report, eval=FALSE--------------------------------------------
## args <- systemArgs(sysma="param/tophat.param", mytargets="targets.txt")
## fqlist <- seeFastq(fastq=infile1(args), batchsize=100000, klength=8)
## pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
## seeFastqPlot(fqlist)
## dev.off()

## ----load_sysargs, eval=FALSE--------------------------------------------
## args <- systemArgs(sysma="param/bwa.param", mytargets="targets.txt")
## sysargs(args)[1] # Command-line parameters for first FASTQ file

## ----bwa_serial, eval=FALSE----------------------------------------------
## moduleload(modules(args))
## system("bwa index -a bwtsw ./data/tair10.fasta")
## bampaths <- runCommandline(args=args)
## writeTargetsout(x=args, file="targets_bam.txt", overwrite=TRUE)

## ----bwa_parallel, eval=FALSE--------------------------------------------
## moduleload(modules(args))
## system("bwa index -a bwtsw ./data/tair10.fasta")
## resources <- list(walltime="1:00:00", ntasks=1, ncpus=cores(args), memory="10G")
## reg <- clusterRun(args, conffile=".BatchJobs.R", template="slurm.tmpl", Njobs=18, runid="01",
##                   resourceList=resources)
## waitForJobs(reg)
## writeTargetsout(x=args, file="targets_bam.txt", overwrite=TRUE)

## ----check_file_presence, eval=FALSE-------------------------------------
## file.exists(outpaths(args))

## ----gsnap_parallel, eval=FALSE------------------------------------------
## library(gmapR); library(BiocParallel); library(BatchJobs)
## args <- systemArgs(sysma="param/gsnap.param", mytargets="targetsPE.txt")
## gmapGenome <- GmapGenome(systemPipeR::reference(args), directory="data", name="gmap_tair10chr", create=TRUE)
## f <- function(x) {
##     library(gmapR); library(systemPipeR)
##     args <- systemArgs(sysma="param/gsnap.param", mytargets="targetsPE.txt")
##     gmapGenome <- GmapGenome(reference(args), directory="data", name="gmap_tair10chr", create=FALSE)
##     p <- GsnapParam(genome=gmapGenome, unique_only=TRUE, molecule="DNA", max_mismatches=3)
##     o <- gsnap(input_a=infile1(args)[x], input_b=infile2(args)[x], params=p, output=outfile1(args)[x])
## }
## funs <- makeClusterFunctionsSLURM("slurm.tmpl")
## param <- BatchJobsParam(length(args), resources=list(walltime="00:20:00", ntasks=1, ncpus=1, memory="6G"), cluster.functions=funs)
## register(param)
## d <- bplapply(seq(along=args), f)
## writeTargetsout(x=args, file="targets_gsnap_bam.txt", overwrite=TRUE)

## ----align_stats, eval=FALSE---------------------------------------------
## read_statsDF <- alignStats(args=args)
## write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")

## ----symbolic_links, eval=FALSE------------------------------------------
## symLink2bam(sysargs=args, htmldir=c("~/.html/", "projects/gen242/"),
##             urlbase="http://biocluster.ucr.edu/~tgirke/",
##             urlfile="./results/IGVurl.txt")

## ----run_gatk, eval=FALSE------------------------------------------------
## moduleload("picard/1.130"); moduleload("samtools/1.3")
## system("picard CreateSequenceDictionary R=./data/tair10.fasta O=./data/tair10.dict")
## system("samtools faidx data/tair10.fasta")
## args <- systemArgs(sysma="param/gatk.param", mytargets="targets_bam.txt")
## resources <- list(walltime="00:20:00", ntasks=1, ncpus=1, memory="6G")
## reg <- clusterRun(args, conffile=".BatchJobs.R", template="slurm.tmpl", Njobs=18, runid="01",
##                   resourceList=resources)
## waitForJobs(reg)
## # unlink(outfile1(args), recursive = TRUE, force = TRUE)
## writeTargetsout(x=args, file="targets_gatk.txt", overwrite=TRUE)

## ----run_bcftools, eval=FALSE--------------------------------------------
## args <- systemArgs(sysma="param/sambcf.param", mytargets="targets_bam.txt")
## resources <- list(walltime="00:20:00", ntasks=1, ncpus=1, memory="6G")
## reg <- clusterRun(args, conffile=".BatchJobs.R", template="slurm.tmpl", Njobs=18, runid="01",
##                   resourceList=resources)
## waitForJobs(reg)
## # unlink(outfile1(args), recursive = TRUE, force = TRUE)
## writeTargetsout(x=args, file="targets_sambcf.txt", overwrite=TRUE)

## ----run_varianttools, eval=FALSE----------------------------------------
## library(gmapR); library(BiocParallel); library(BatchJobs)
## args <- systemArgs(sysma="param/vartools.param", mytargets="targets_gsnap_bam.txt")
## f <- function(x) {
##     library(VariantTools); library(gmapR); library(systemPipeR)
##     args <- systemArgs(sysma="param/vartools.param", mytargets="targets_gsnap_bam.txt")
##     gmapGenome <- GmapGenome(systemPipeR::reference(args), directory="data", name="gmap_tair10chr", create=FALSE)
##     tally.param <- TallyVariantsParam(gmapGenome, high_base_quality = 23L, indels = TRUE)
##     bfl <- BamFileList(infile1(args)[x], index=character())
##     var <- callVariants(bfl[[1]], tally.param)
##     sampleNames(var) <- names(bfl)
##     writeVcf(asVCF(var), outfile1(args)[x], index = TRUE)
## }
## funs <- makeClusterFunctionsSLURM("slurm.tmpl")
## param <- BatchJobsParam(length(args), resources=list(walltime="00:20:00", ntasks=1, ncpus=1, memory="6G"), cluster.functions=funs)
## register(param)
## d <- bplapply(seq(along=args), f)
## writeTargetsout(x=args, file="targets_vartools.txt", overwrite=TRUE)

## ----inspect_vcf, eval=FALSE---------------------------------------------
## library(VariantAnnotation)
## args <- systemArgs(sysma="param/filter_gatk.param", mytargets="targets_gatk.txt")
## vcf <- readVcf(infile1(args)[1], "A. thaliana")
## vcf
## vr <- as(vcf, "VRanges")
## vr

## ----filter_gatk, eval=FALSE---------------------------------------------
## library(VariantAnnotation)
## library(BBmisc) # Defines suppressAll()
## args <- systemArgs(sysma="param/filter_gatk.param", mytargets="targets_gatk.txt")[1:4]
## filter <- "totalDepth(vr) >= 2 & (altDepth(vr) / totalDepth(vr) >= 0.8) & rowSums(softFilterMatrix(vr))>=1"
## # filter <- "totalDepth(vr) >= 20 & (altDepth(vr) / totalDepth(vr) >= 0.8) & rowSums(softFilterMatrix(vr))==6"
## suppressAll(filterVars(args, filter, varcaller="gatk", organism="A. thaliana"))
## writeTargetsout(x=args, file="targets_gatk_filtered.txt", overwrite=TRUE)

## ----filter_bcftools, eval=FALSE-----------------------------------------
## args <- systemArgs(sysma="param/filter_sambcf.param", mytargets="targets_sambcf.txt")[1:4]
## filter <- "rowSums(vr) >= 2 & (rowSums(vr[,3:4])/rowSums(vr[,1:4]) >= 0.8)"
## # filter <- "rowSums(vr) >= 20 & (rowSums(vr[,3:4])/rowSums(vr[,1:4]) >= 0.8)"
## suppressAll(filterVars(args, filter, varcaller="bcftools", organism="A. thaliana"))
## writeTargetsout(x=args, file="targets_sambcf_filtered.txt", overwrite=TRUE)

## ----filter_varianttools, eval=FALSE-------------------------------------
## library(VariantAnnotation)
## library(BBmisc) # Defines suppressAll()
## args <- systemArgs(sysma="param/filter_vartools.param", mytargets="targets_vartools.txt")[1:4]
## filter <- "(values(vr)$n.read.pos.ref + values(vr)$n.read.pos) >= 2 & (values(vr)$n.read.pos / (values(vr)$n.read.pos.ref + values(vr)$n.read.pos) >= 0.8)"
## # filter <- "(values(vr)$n.read.pos.ref + values(vr)$n.read.pos) >= 20 & (values(vr)$n.read.pos / (values(vr)$n.read.pos.ref + values(vr)$n.read.pos) >= 0.8)"
## filterVars(args, filter, varcaller="vartools", organism="A. thaliana")
## writeTargetsout(x=args, file="targets_vartools_filtered.txt", overwrite=TRUE)

## ----check_filter, eval=FALSE--------------------------------------------
## length(as(readVcf(infile1(args)[1], genome="Ath"), "VRanges")[,1])
## length(as(readVcf(outpaths(args)[1], genome="Ath"), "VRanges")[,1])

## ----annotate_basics, eval=FALSE-----------------------------------------
## library("GenomicFeatures")
## args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_gatk_filtered.txt")
## txdb <- loadDb("./data/tair10.sqlite")
## vcf <- readVcf(infile1(args)[1], "A. thaliana")
## locateVariants(vcf, txdb, CodingVariants())

## ----annotate_basics_non-synon, eval=FALSE-------------------------------
## fa <- FaFile(systemPipeR::reference(args))
## predictCoding(vcf, txdb, seqSource=fa)

## ----annotate_gatk, eval=FALSE-------------------------------------------
## library("GenomicFeatures")
## args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_gatk_filtered.txt")
## txdb <- loadDb("./data/tair10.sqlite")
## fa <- FaFile(systemPipeR::reference(args))
## suppressAll(variantReport(args=args, txdb=txdb, fa=fa, organism="A. thaliana"))

## ----annotate_bcftools, eval=FALSE---------------------------------------
## args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_sambcf_filtered.txt")
## txdb <- loadDb("./data/tair10.sqlite")
## fa <- FaFile(systemPipeR::reference(args))
## suppressAll(variantReport(args=args, txdb=txdb, fa=fa, organism="A. thaliana"))

## ----annotate_varianttools, eval=FALSE-----------------------------------
## args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_vartools_filtered.txt")
## txdb <- loadDb("./data/tair10.sqlite")
## fa <- FaFile(systemPipeR::reference(args))
## suppressAll(variantReport(args=args, txdb=txdb, fa=fa, organism="A. thaliana"))

## ----view_annotation, eval=FALSE-----------------------------------------
## read.delim(outpaths(args)[1])[38:40,]

## ----combine_gatk, eval=FALSE--------------------------------------------
## args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_gatk_filtered.txt")
## combineDF <- combineVarReports(args, filtercol=c(Consequence="nonsynonymous"))
## write.table(combineDF, "./results/combineDF_nonsyn_gatk.xls", quote=FALSE, row.names=FALSE, sep="\t")

## ----combine_bcftools, eval=FALSE----------------------------------------
## args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_sambcf_filtered.txt")
## combineDF <- combineVarReports(args, filtercol=c(Consequence="nonsynonymous"))
## write.table(combineDF, "./results/combineDF_nonsyn_sambcf.xls", quote=FALSE, row.names=FALSE, sep="\t")

## ----combine_varianttools, eval=FALSE------------------------------------
## args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_vartools_filtered.txt")
## combineDF <- combineVarReports(args, filtercol=c(Consequence="nonsynonymous"))
## write.table(combineDF, "./results/combineDF_nonsyn_vartools.xls", quote=FALSE, row.names=FALSE, sep="\t")
## combineDF[2:4,]

## ----summary_gatk, eval=FALSE--------------------------------------------
## args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_gatk_filtered.txt")
## varSummary(args)
## write.table(varSummary(args), "./results/variantStats_gatk.xls", quote=FALSE, col.names = NA, sep="\t")

## ----summary_bcftools, eval=FALSE----------------------------------------
## args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_sambcf_filtered.txt")
## varSummary(args)
## write.table(varSummary(args), "./results/variantStats_sambcf.xls", quote=FALSE, col.names = NA, sep="\t")

## ----summary_varianttools, eval=FALSE------------------------------------
## args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_vartools_filtered.txt")
## varSummary(args)
## write.table(varSummary(args), "./results/variantStats_vartools.xls", quote=FALSE, col.names = NA, sep="\t")

## ----venn_diagram, eval=FALSE--------------------------------------------
## args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_gatk_filtered.txt")
## varlist <- sapply(names(outpaths(args))[1:4], function(x) as.character(read.delim(outpaths(args)[x])$VARID))
## vennset_gatk <- overLapper(varlist, type="vennsets")
## args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_sambcf_filtered.txt")
## varlist <- sapply(names(outpaths(args))[1:4], function(x) as.character(read.delim(outpaths(args)[x])$VARID))
## vennset_bcf <- overLapper(varlist, type="vennsets")
## args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_vartools_filtered.txt")
## varlist <- sapply(names(outpaths(args))[1:4], function(x) as.character(read.delim(outpaths(args)[x])$VARID))
## vennset_vartools <- overLapper(varlist, type="vennsets")
## pdf("./results/vennplot_var.pdf")
## vennPlot(list(vennset_gatk, vennset_bcf, vennset_vartools), mymain="", mysub="GATK: red; BCFtools: blue; VariantTools: green", colmode=2, ccol=c("red", "blue", "green"))
## dev.off()

## ----plot_variant, eval=FALSE--------------------------------------------
## library(ggbio)
## mychr <- "ChrC"; mystart <- 11000; myend <- 13000
## args <- systemArgs(sysma="param/bwa.param", mytargets="targets.txt")
## ga <- readGAlignments(outpaths(args)[1], use.names=TRUE, param=ScanBamParam(which=GRanges(mychr, IRanges(mystart, myend))))
## p1 <- autoplot(ga, geom = "rect")
## p2 <- autoplot(ga, geom = "line", stat = "coverage")
## p3 <- autoplot(vcf[seqnames(vcf)==mychr], type = "fixed") + xlim(mystart, myend) + theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y=element_blank())
## p4 <- autoplot(txdb, which=GRanges(mychr, IRanges(mystart, myend)), names.expr = "gene_id")
## png("./results/plot_variant.png")
## tracks(Reads=p1, Coverage=p2, Variant=p3, Transcripts=p4, heights = c(0.3, 0.2, 0.1, 0.35)) + ylab("")
## dev.off()

## ----sessionInfo---------------------------------------------------------
sessionInfo()

