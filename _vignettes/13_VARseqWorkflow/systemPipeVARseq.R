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


## ----genVAR_workflow, eval=FALSE--------------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="varseq")
## setwd("varseq")


## Rscript -e "systemPipeRdata::genWorkenvir(workflow='varseq')"


## ----closeR, eval=FALSE-----------------------------------
## q("no") # closes R session on head node


## srun --x11 --partition=short --mem=2gb --cpus-per-task 4 --ntasks 1 --time 2:00:00 --pty bash -l

## module load R/3.6.0

## R


## ----r_environment, eval=FALSE----------------------------
## system("hostname") # should return the computer name or cluster name
## getwd() # checks current working directory of R session
## dir() # returns content of current working directory


## ----load_systempiper, eval=TRUE, messages=FALSE, warnings=FALSE----
library(systemPipeR)


## ----load_custom_fct, eval=FALSE--------------------------
## source("systemPipeVARseq_Fct.R")


## ----load_targets_file, eval=TRUE-------------------------
targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")
targets[1:4, 1:4]


## ----construct_SYSargs2_trim-pe, eval=FALSE---------------
## targetsPE <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
## dir_path <- system.file("extdata/cwl/preprocessReads/trim-pe", package="systemPipeR")
## trim <- loadWorkflow(targets=targetsPE, wf_file="trim-pe.cwl", input_file="trim-pe.yml", dir_path=dir_path)
## trim <- renderWF(trim, inputvars=c(FileName1="_FASTQ_PATH1_", FileName2="_FASTQ_PATH2_", SampleName="_SampleName_"))
## trim
## output(trim)[1:2]
## 
## preprocessReads(args=trim, Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)",
##                 batchsize=100000, overwrite=TRUE, compress=TRUE)
## writeTargetsout(x=trim, file="targets_trimPE.txt", step=1, new_col = c("FileName1", "FileName2"),
##                 new_col_output_index = c(1,2), overwrite = TRUE)


## ----fastq_report, eval=FALSE-----------------------------
## fqlist <- seeFastq(fastq=infile1(trim), batchsize=100000, klength=8)
## pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
## seeFastqPlot(fqlist)
## dev.off()


## ----dir_path, eval=FALSE---------------------------------
## dir_path <- system.file("extdata/cwl/gatk", package="systemPipeR")


## ----index, eval=FALSE------------------------------------
## ## Index for BWA
## args <- loadWorkflow(targets = NULL, wf_file = "bwa-index.cwl", input_file = "gatk.yaml", dir_path = dir_path)
## args <- renderWF(args)
## cmdlist(args) # shows the command
## output(args) # shows the expected output files
## # Run single Machine
## runCommandline(args, make_bam = FALSE)
## 
## ## Index needed for gatk tools
## args <- loadWorkflow(wf_file = "fasta_dict.cwl", input_file = "gatk.yaml", dir_path = dir_path)
## args <- renderWF(args)
## args <- runCommandline(args, make_bam = FALSE)
## 
## ## Index
## args <- loadWorkflow(wf_file = "fasta_faidx.cwl", input_file = "gatk.yaml", dir_path = dir_path)
## args <- renderWF(args)
## args <- runCommandline(args, make_bam = FALSE)


## ----bwa-pe_alignment, eval=FALSE-------------------------
## targetsPE <- system.file("extdata", "targetsPE.txt", package = "systemPipeR")
## args <- loadWorkflow(targets = targetsPE, wf_file = "bwa-pe.cwl",
##                        input_file = "gatk.yaml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName1 = "_FASTQ_PATH1_", FileName2 = "_FASTQ_PATH2_",
##                                      SampleName = "_SampleName_"))
## cmdlist(args)[1:2]
## output(args)[1:2]


## ----start_BWA, eval=FALSE--------------------------------
## args <- runCommandline(args = args, make_bam=FALSE)
## writeTargetsout(x = args[1:2], file = "./results/targetsPE.txt",
##                 step = 1, new_col = "BWA_SAM", new_col_output_index = 1, overwrite = TRUE)


## ----bwa_parallel, eval=FALSE-----------------------------
## library(batchtools)
## resources <- list(walltime = 120, ntasks = 1, ncpus = 4, memory = 1024)
## reg <- clusterRun(args, FUN = runCommandline, more.args = list(args = args, dir = FALSE, make_bam=FALSE),
##     conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl", Njobs = 18,
##     runid = "01", resourceList = resources)
## getStatus(reg = reg)
## writeTargetsout(x = args, file = "./results/targetsPE.txt",
##                 step = 1, new_col = "BWA_SAM", new_col_output_index = 1, overwrite = TRUE)


## ----check_files_exist, eval=FALSE------------------------
## outpaths <- subsetWF(args , slot="output", subset=1, index=1)
## file.exists(outpaths)


## ----align_stats, eval=FALSE------------------------------
## read_statsDF <- alignStats(args=args)
## write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")


## ----symbolic_links, eval=FALSE---------------------------
## symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"),
##             urlbase="http://cluster.hpcc.ucr.edu/~tgirke/", urlfile="IGVurl.txt")


## ----fastq2ubam, eval=FALSE-------------------------------
## dir_path <- system.file("extdata/cwl/gatk", package="systemPipeR")
## targets.gatk <- "./results/targetsPE.txt" ## targets generated from BWA
## args <- loadWorkflow(targets = targets.gatk, wf_file = "gatk_fastq2ubam.cwl",
##                        input_file = "gatk.yaml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName1 = "_FASTQ_PATH1_", FileName2 = "_FASTQ_PATH2_",
##                                      SampleName = "_SampleName_"))
## cmdlist(args)[1:2]
## output(args)[1:2]
## args <- runCommandline(args= args[1:2],  make_bam=FALSE)
## writeTargetsout(x = args, file = "./results/targets_gatk.txt",
##                 step = 1, new_col = "GATK_UBAM", new_col_output_index = 1, overwrite = TRUE)


## ----merge_bam, eval=FALSE--------------------------------
## targets.gatk <- "./results/targets_gatk.txt"
## args <- loadWorkflow(targets = targets.gatk, wf_file = "gatk_mergebams.cwl",
##                        input_file = "gatk.yaml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(BWA_SAM = "_bwasam_", GATK_UBAM = "_ubam_",
##                                      SampleName = "_SampleName_"))
## cmdlist(args)[1:2]
## output(args)[1:2]
## args <- runCommandline(args= args,  make_bam=FALSE)
## writeTargetsout(x = args, file = "./results/targets_gatk.txt",
##                 step = 1, new_col = "GATK_MERGE", new_col_output_index = 1, overwrite = TRUE)


## ----sort, eval=FALSE-------------------------------------
## targets.gatk <- "./results/targets_gatk.txt"
## args <- loadWorkflow(targets = targets.gatk, wf_file = "gatk_sort.cwl",
##                      input_file = "gatk.yaml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(SampleName = "_SampleName_", GATK_MERGE = "_mergebam_"))
## cmdlist(args)[1:2]
## output(args)[1:2]
## args <- runCommandline(args= args,  make_bam=FALSE)
## writeTargetsout(x = args, file = "./results/targets_gatk.txt",
##                 step = 1, new_col = "GATK_SORT", new_col_output_index = 1, overwrite = TRUE)


## ----mark_dup, eval=FALSE---------------------------------
## targets.gatk <- "./results/targets_gatk.txt"
## args <- loadWorkflow(targets = targets.gatk, wf_file = "gatk_markduplicates.cwl",
##                      input_file = "gatk.yaml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(SampleName = "_SampleName_", GATK_SORT = "_sort_"))
## cmdlist(args)[1:2]
## output(args)[1:2]
## args <- runCommandline(args= args,  make_bam=FALSE)
## writeTargetsout(x = args, file = "./results/targets_gatk.txt",
##                 step = 1, new_col = c("GATK_MARK", "GATK_MARK_METRICS"),
##                 new_col_output_index = c(1,2), overwrite = TRUE)


## ----fix_tag, eval=FALSE----------------------------------
## targets.gatk <- "./results/targets_gatk.txt"
## args <- loadWorkflow(targets = targets.gatk, wf_file = "gatk_fixtag.cwl",
##                      input_file = "gatk.yaml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(SampleName = "_SampleName_", GATK_MARK = "_mark_"))
## cmdlist(args)[1:2]
##   output(args)[1:2]
##   args <- runCommandline(args= args,  make_bam=FALSE)
##   writeTargetsout(x = args, file = "./results/targets_gatk.txt",
##                   step = 1, new_col = "GATK_FIXED", new_col_output_index = 1, overwrite = TRUE)


## ----hc, eval=FALSE---------------------------------------
## targets.gatk <- "./results/targets_gatk.txt"
## args <- loadWorkflow(targets = targets.gatk, wf_file = "gatk_haplotypecaller.cwl",
##                      input_file = "gatk.yaml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(SampleName = "_SampleName_", GATK_FIXED = "_fixed_"))
## cmdlist(args)[1:2]
## output(args)[1:2]
## args <- runCommandline(args= args,  make_bam=FALSE)
## writeTargetsout(x = args, file = "./results/targets_gatk.txt",
##                 step = 1, new_col = "GVCF", new_col_output_index = 1, overwrite = TRUE)


## ----import, eval=FALSE-----------------------------------
## # drop all  *.g.vcf.gz files into results folder, make sure the tbi index is also there.
## args <- loadWorkflow(targets = NULL, wf_file = "gatk_genomicsDBImport.cwl",
##                      input_file = "gatk.yaml", dir_path = dir_path)
## args <- renderWF(args)
## cmdlist(args)
## output(args)
## args <- runCommandline(args= args,  make_bam=FALSE)


## ----call_variants, eval=FALSE----------------------------
## args <- loadWorkflow(targets = NULL, wf_file = "gatk_genotypeGVCFs.cwl",
##                        input_file = "gatk.yaml", dir_path = dir_path)
## args <- renderWF(args)
## cmdlist(args)
## output(args)
## args <- runCommandline(args= args,  make_bam=FALSE)


## ----filter, eval=FALSE-----------------------------------
## args <- loadWorkflow(wf_file = "gatk_variantFiltration.cwl",
##                        input_file = "gatk.yaml", dir_path = dir_path)
## args <- renderWF(args)
## cmdlist(args)
## output(args)
## args <- runCommandline(args= args,  make_bam=FALSE)


## ----extract_single_vcf, eval=F---------------------------
## targets.gatk <- "./results/targets_gatk.txt"
## args <- loadWorkflow(targets = targets.gatk, wf_file = "gatk_select_variant.cwl",
##                        input_file = "gatk.yaml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(SampleName = "_SampleName_"))
## cmdlist(args)[1:2]
## output(args)[1:2]
## args <- runCommandline(args= args, make_bam=FALSE)
## writeTargetsout(x = args, file = "./results/targets_gatk.txt",
##                 step = 1, new_col = "FileName1", new_col_output_index = 1, overwrite = TRUE)


## ----run_bcftools, eval=FALSE-----------------------------
## dir_path <- system.file("extdata/cwl/workflow-bcftools", package="systemPipeR")
## targetsPE <- './results/targetsPE.txt'
## args <- loadWorkflow(targets = targetsPE, wf_file = "workflow_bcftools.cwl",
##                      input_file = "bcftools.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(SampleName = "_SampleName_", BWA_SAM = "_SAM_"))
## cmdlist(args[1])
## output(args[1])
## args <- runCommandline(args= args, make_bam=FALSE)
## writeTargetsout(x = args, file = "./results/targets_bcf.txt",
##                 step = 5, new_col = "FileName1", new_col_output_index = 1, overwrite = TRUE)


## dir_path <- system.file("extdata/cwl/varseq", package="systemPipeR")


## ----inspect_vcf, eval=FALSE------------------------------
## library(VariantAnnotation)
## args <- loadWorkflow(targets = './results/targets_gatk.txt', wf_file = "filter.cwl",
##                      input_file = "varseq.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_"))
## vcf <- readVcf(infile1(args)[1], "A. thaliana")
## vcf
## vr <- as(vcf, "VRanges")
## vr


## ----filter_gatk, eval=FALSE------------------------------
## dir_path <- system.file("extdata/cwl/varseq", package="systemPipeR")
## library(VariantAnnotation)
## library(BBmisc) # Defines suppressAll()
## args <- loadWorkflow(targets = './results/targets_gatk.txt',
##                      wf_file = "filter.cwl",input_file = "varseq.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName="_SampleName_"))
## filter <- "totalDepth(vr) >= 2 & (altDepth(vr) / totalDepth(vr) >= 0.8)"
## # filter <- "totalDepth(vr) >= 20 & (altDepth(vr) / totalDepth(vr) >= 0.8)"
## suppressAll(filterVars(args, filter, varcaller="gatk", organism="A. thaliana"))
## writeTargetsout(x = args, file = "./results/targets_filter_gatk.txt",
##                 step = 1, new_col = "FileName1", new_col_output_index = 1, overwrite = TRUE)


## ----filter_bcftools, eval=FALSE--------------------------
## args <- loadWorkflow(targets = './results/targets_bcf.txt',
##                      wf_file = "filter.cwl",
##                      input_file = "varseq.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = '_SampleName_'))
## filter <- "rowSums(vr) >= 2 & (rowSums(vr[,3:4])/rowSums(vr[,1:4]) >= 0.8)"
## # filter <- "rowSums(vr) >= 20 & (rowSums(vr[,3:4])/rowSums(vr[,1:4]) >= 0.8)"
## suppressAll(filterVars(args, filter, varcaller="bcftools", organism="A. thaliana"))
## writeTargetsout(x = args, file = "./results/targets_filter_bcf.txt",
##                 step = 1, new_col = "FileName1", new_col_output_index = 1, overwrite = TRUE)


## ----check_filter, eval=FALSE-----------------------------
## length(as(readVcf(infile1(args)[1], genome="Ath"), "VRanges")[,1])
## length(as(readVcf(subsetWF(args, slot='output', subset = 1, index=1)[1], genome="Ath"), "VRanges")[,1])


## ----annotate_basics, eval=FALSE--------------------------
## dir_path <- system.file("extdata/cwl/varseq", package="systemPipeR")
## library("GenomicFeatures")
## args <- loadWorkflow(targets = './results/targets_filter_gatk.txt',
##                      wf_file = "annotate.cwl", input_file = "varseq.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = '_SampleName_'))
## txdb <- loadDb("./data/tair10.sqlite")
## vcf <- readVcf(infile1(args)[1], "A. thaliana")
## locateVariants(vcf, txdb, CodingVariants())


## ----annotate_basics_non-synon, eval=FALSE----------------
## fa <- FaFile(normalizePath(file.path(args$yamlinput$data_path$path,args$yamlinput$ref_name)))
## predictCoding(vcf, txdb, seqSource=fa)


## ----annotate_gatk, eval=FALSE----------------------------
## library("GenomicFeatures")
## args <- loadWorkflow(targets = './results/targets_filter_gatk.txt',
##                      wf_file = "annotate.cwl", input_file = "varseq.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = '_SampleName_'))
## txdb <- loadDb("./data/tair10.sqlite")
## fa <- FaFile(normalizePath(file.path(args$yamlinput$data_path$path,args$yamlinput$ref_name)))
## suppressAll(variantReport(args=args, txdb=txdb, fa=fa, organism="A. thaliana"))
## writeTargetsout(x = args, file = "./results/targets_report_gatk.txt",
##                 step = 1, new_col = "FileName1", new_col_output_index = 1, overwrite = TRUE)


## ----annotate_bcf, eval=FALSE-----------------------------
## library("GenomicFeatures")
## args <- loadWorkflow(targets = './results/targets_filter_bcf.txt',
##                      wf_file = "annotate.cwl",
##                      input_file = "varseq.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = '_SampleName_'))
## txdb <- loadDb("./data/tair10.sqlite")
## fa <- FaFile(normalizePath(file.path(args$yamlinput$data_path$path,args$yamlinput$ref_name)))
## suppressAll(variantReport(args=args, txdb=txdb, fa=fa, organism="A. thaliana"))
## writeTargetsout(x = args, file = "./results/targets_report_bcf.txt",
##                 step = 1, new_col = "FileName1", new_col_output_index = 1, overwrite = TRUE)


## ----view_annotation, eval=FALSE--------------------------
## read.delim(output(args)[[1]][[1]])[38:40,]


## ----combine_gatk, eval=FALSE-----------------------------
## dir_path <- system.file("extdata/cwl/varseq", package="systemPipeR")
## args <- loadWorkflow(targets = 'results/targets_report_gatk.txt',
##                      wf_file = "combine.cwl", input_file = "varseq.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = '_SampleName_'))
## combineDF <- combineVarReports(args, filtercol=c(Consequence="nonsynonymous"))
## write.table(combineDF, "./results/combineDF_nonsyn_gatk.xls", quote=FALSE, row.names=FALSE, sep="\t")


## ----combine_bcf, eval=FALSE------------------------------
## args <- loadWorkflow(targets = 'results/targets_report_bcf.txt',
##                      wf_file = "combine.cwl", input_file = "varseq.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = '_SampleName_'))
## combineDF <- combineVarReports(args, filtercol=c(Consequence="nonsynonymous"))
## write.table(combineDF, "./results/combineDF_nonsyn_bcf.xls", quote=FALSE, row.names=FALSE, sep="\t")


## ----summary_gatk, eval=FALSE-----------------------------
## args <- loadWorkflow(targets = './results/targets_report_gatk.txt',
##                      wf_file = "combine.cwl", input_file = "varseq.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = '_SampleName_'))
## varSummary(args)
## write.table(varSummary(args), "./results/variantStats_gatk.xls", quote=FALSE, col.names = NA, sep="\t")


## ----summary_bcf, eval=FALSE------------------------------
## args <- loadWorkflow(targets = './results/targets_report_bcf.txt',
##                      wf_file = "combine.cwl", input_file = "varseq.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = '_SampleName_'))
## varSummary(args)
## write.table(varSummary(args), "./results/variantStats_bcf.xls", quote=FALSE, col.names = NA, sep="\t")


## ----venn_diagram, eval=FALSE-----------------------------
## dir_path <- system.file("extdata/cwl/varseq", package="systemPipeR")
## ## gatk
## args <- loadWorkflow(targets = 'results/targets_report_gatk.txt',
##                      wf_file = "combine.cwl",  input_file = "varseq.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = '_SampleName_'))
## varlist <- sapply(names(subsetWF(args[1:2], slot='output', subset = 1, index=1)),
##                   function(x) as.character(read.delim(subsetWF(args[1:2], slot='output', subset = 1, index=1)[x])$VARID))
## vennset_gatk <- overLapper(varlist, type="vennsets")
## 
## ## bcf
## args <- loadWorkflow(targets = './results/targets_report_bcf.txt',
##                      wf_file = "combine.cwl", input_file = "varseq.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = '_SampleName_'))
## varlist <- sapply(names(subsetWF(args[1:2], slot='output', subset=1, index=1)),
##                   function(x) as.character(read.delim(subsetWF(args[1:2], slot='output', subset=1, index=1)[x])$VARID))
## vennset_bcf <- overLapper(varlist, type="vennsets")
## 
## pdf("./results/vennplot_var.pdf")
## vennPlot(list(vennset_gatk, vennset_bcf), mymain="", mysub="GATK: red; BCFtools: blue", colmode=2, ccol=c("red", "blue"))
## dev.off()


## ----plot_variant, eval=FALSE-----------------------------
## library(ggbio)
## mychr <- "ChrC"; mystart <- 11000; myend <- 13000
## args <- loadWorkflow(targets = 'results/targets_gatk.txt', wf_file = "combine.cwl",
##                        input_file = "varseq.yml", dir_path = 'param/cwl/varseq_downstream/')
## args <- renderWF(args, inputvars = c(GATK_FIXED = "_FILE1_", SampleName = "_SampleName_"))
## ga <- readGAlignments(subsetWF(args, slot='input', subset = 1)[1], use.names=TRUE,
##                       param=ScanBamParam(which=GRanges(mychr, IRanges(mystart, myend))))
## p1 <- autoplot(ga, geom = "rect")
## p2 <- autoplot(ga, geom = "line", stat = "coverage")
## p3 <- autoplot(vcf[seqnames(vcf)==mychr], type = "fixed") +
##                 xlim(mystart, myend) + theme(legend.position = "none",
##                     axis.text.y = element_blank(), axis.ticks.y=element_blank())
## p4 <- autoplot(loadDb("./data/tair10.sqlite"), which=GRanges(mychr, IRanges(mystart, myend)), names.expr = "gene_id")
## png("./results/plot_variant.png")
## tracks(Reads=p1, Coverage=p2, Variant=p3, Transcripts=p4, heights = c(0.3, 0.2, 0.1, 0.35)) + ylab("")
## dev.off()


## ----sessionInfo------------------------------------------
sessionInfo()

