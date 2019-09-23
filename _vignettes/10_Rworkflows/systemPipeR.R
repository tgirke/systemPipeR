## pre code {

## white-space: pre !important;

## overflow-x: scroll !important;

## word-break: keep-all !important;

## word-wrap: initial !important;

## }


## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()
options(width=80, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")), 
    tidy.opts=list(width.cutoff=80), tidy=TRUE)


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
    library(batchtools)
})


## ----install, eval=FALSE-------------------------------------------------
## if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
## BiocManager::install("systemPipeR")
## BiocManager::install("systemPipeRdata")


## ----documentation, eval=FALSE-------------------------------------------
## library("systemPipeR") # Loads the package
## library(help="systemPipeR") # Lists package info
## vignette("systemPipeR") # Opens vignette


## ----genRna_workflow, eval=FALSE-----------------------------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="rnaseq")
## setwd("rnaseq")


## ----targetsSE, eval=TRUE------------------------------------------------
library(systemPipeR)
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR") 
read.delim(targetspath, comment.char = "#")


## ----targetsPE, eval=TRUE------------------------------------------------
targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
read.delim(targetspath, comment.char = "#")[1:2,1:6]


## ----comment_lines, echo=TRUE--------------------------------------------
readLines(targetspath)[1:4]


## ----targetscomp, eval=TRUE----------------------------------------------
readComp(file=targetspath, format="vector", delim="-")


## ----CWL_structure, echo = FALSE, eval=FALSE-----------------------------
## hisat2.cwl <- system.file("extdata", "cwl/hisat2/hisat2-se/hisat2-mapping-se.cwl", package="systemPipeR")
## yaml::read_yaml(hisat2.cwl)


## ----yaml_structure, echo = FALSE, eval=FALSE----------------------------
## hisat2.yml <- system.file("extdata", "cwl/hisat2/hisat2-se/hisat2-mapping-se.yml", package="systemPipeR")
## yaml::read_yaml(hisat2.yml)


## ----SYSargs2_structure, eval=TRUE---------------------------------------
library(systemPipeR)
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl/hisat2/hisat2-se", package="systemPipeR")
WF <- loadWorkflow(targets=targets, wf_file="hisat2-mapping-se.cwl",
                   input_file="hisat2-mapping-se.yml",
                   dir_path=dir_path)

WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))


## ----names_WF, eval=TRUE-------------------------------------------------
names(WF)


## ----cmdlist, eval=TRUE--------------------------------------------------
cmdlist(WF)[1]


## ----output_WF, eval=TRUE------------------------------------------------
output(WF)[1]
modules(WF)
targets(WF)[1]
targets.as.df(targets(WF))[1:4,1:4]
output(WF)[1]
cwlfiles(WF)
inputvars(WF)


## ----param_structure, eval=TRUE------------------------------------------
parampath <- system.file("extdata", "tophat.param", package="systemPipeR")
read.delim(parampath, comment.char = "#")


## ----param_import, eval=TRUE---------------------------------------------
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR") 
args <- suppressWarnings(systemArgs(sysma=parampath, mytargets=targetspath))
args


## ----sysarg_access, eval=TRUE--------------------------------------------
names(args)


## ----sysarg_access2, eval=TRUE-------------------------------------------
sysargs(args)[1]
modules(args)
cores(args)
outpaths(args)[1]


## ----sysarg_json, eval=TRUE----------------------------------------------
systemArgs(sysma=parampath, mytargets=targetspath, type="json")


## ----load_package, eval=FALSE--------------------------------------------
## library(systemPipeR)
## library(systemPipeRdata)
## genWorkenvir(workflow="rnaseq", mydirname=NULL)
## setwd("rnaseq")


## ----construct_SYSargs2_trim-se, echo = FALSE, eval=FALSE----------------
## targets <- system.file("extdata", "targets.txt", package="systemPipeR")
## dir_path <- system.file("extdata/cwl/preprocessReads/trim-se", package="systemPipeR")
## trim <- loadWorkflow(targets=targets, wf_file="trim-se.cwl", input_file="trim-se.yml", dir_path=dir_path)
## trim <- renderWF(trim, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
## trim


## ----construct_SYSargs2_trim-pe, eval=FALSE------------------------------
## targetsPE <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
## dir_path <- system.file("extdata/cwl/preprocessReads/trim-pe", package="systemPipeR")
## trim <- loadWorkflow(targets=targetsPE, wf_file="trim-pe.cwl", input_file="trim-pe.yml", dir_path=dir_path)
## trim <- renderWF(trim, inputvars=c(FileName1="_FASTQ_PATH1_", FileName2="_FASTQ_PATH2_", SampleName="_SampleName_"))
## trim
## output(trim)[1:2]


## ----preprocessing, eval=FALSE-------------------------------------------
## preprocessReads(args=trim, Fct="trimLRPatterns(Rpattern='GCCCGGGTAA',
##                 subject=fq)",
##                 batchsize=100000, overwrite=TRUE, compress=TRUE)
## writeTargetsout(x=trim, file="targets_trimPE.txt", step=1)


## ----custom_preprocessing, eval=FALSE------------------------------------
## filterFct <- function(fq, cutoff=20, Nexceptions=0) {
##     qcount <- rowSums(as(quality(fq), "matrix") <= cutoff, na.rm=TRUE)
##     # Retains reads where Phred scores are >= cutoff with N exceptions
##     fq[qcount <= Nexceptions]
## }
## preprocessReads(args=trim, Fct="filterFct(fq, cutoff=20, Nexceptions=0)",
##                 batchsize=100000)
## writeTargetsout(x=trim, file="targets_trimPE.txt", step=1)


## ----trimGalore, eval=FALSE----------------------------------------------
## targets <- system.file("extdata", "targets.txt", package="systemPipeR")
## dir_path <- system.file("extdata/cwl/trim_galore/trim_galore-se", package="systemPipeR")
## trimG <- loadWorkflow(targets=targets, wf_file="trim_galore-se.cwl", input_file="trim_galore-se.yml", dir_path=dir_path)
## trimG <- renderWF(trimG, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
## trimG
## cmdlist(trimG)[1:2]
## output(trimG)[1:2]
## ## Run Single Machine Option
## runCommandline(trimG[1], make_bam = FALSE)
## writeTargetsout(x=trimG, file="targets_trimG.txt", step=1)


## ----trimmomatic, eval=FALSE---------------------------------------------
## targetsPE <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
## dir_path <- system.file("extdata/cwl/trimmomatic/trimmomatic-pe", package="systemPipeR")
## trimM <- loadWorkflow(targets=targetsPE, wf_file="trimmomatic-pe.cwl", input_file="trimmomatic-pe.yml", dir_path=dir_path)
## trimM <- renderWF(trimM, inputvars=c(FileName1="_FASTQ_PATH1_", FileName2="_FASTQ_PATH2_", SampleName="_SampleName_"))
## trimM
## cmdlist(trimM)[1:2]
## output(trimM)[1:2]
## ## Run Single Machine Option
## runCommandline(trimM[1], make_bam = FALSE)
## writeTargetsout(x=trimM, file="targets_trimM.txt", step=1)


## ----fastq_quality, eval=FALSE-------------------------------------------
## fqlist <- seeFastq(fastq=infile1(trim), batchsize=10000, klength=8)
## pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
## seeFastqPlot(fqlist)
## dev.off()


## ----fastq_quality_parallel_single, eval=FALSE---------------------------
## f <- function(x) seeFastq(fastq=infile1(trim)[x], batchsize=100000, klength=8)
## fqlist <- bplapply(seq(along=trim), f, BPPARAM = MulticoreParam(workers=4))
## seeFastqPlot(unlist(fqlist, recursive=FALSE))


## ----fastq_quality_parallel_cluster, eval=FALSE--------------------------
## library(BiocParallel); library(batchtools)
## f <- function(x) {
##   library(systemPipeR)
##   targetsPE <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
##   dir_path <- system.file("extdata/cwl/preprocessReads/trim-pe", package="systemPipeR")
##   trim <- loadWorkflow(targets=targetsPE, wf_file="trim-pe.cwl", input_file="trim-pe.yml", dir_path=dir_path)
##   trim <- renderWF(trim, inputvars=c(FileName1="_FASTQ_PATH1_", FileName2="_FASTQ_PATH2_", SampleName="_SampleName_"))
##   seeFastq(fastq=infile1(trim)[x], batchsize=100000, klength=8)
## }
## resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024)
## param <- BatchtoolsParam(workers = 4, cluster = "slurm", template = "batchtools.slurm.tmpl", resources = resources)
## fqlist <- bplapply(seq(along=trim), f, BPPARAM = param)
## seeFastqPlot(unlist(fqlist, recursive=FALSE))


## ----hisat_SYSargs2_object, eval=TRUE------------------------------------
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl/hisat2/hisat2-se", package="systemPipeR")
align <- loadWorkflow(targets=targets, wf_file="hisat2-mapping-se.cwl", input_file="hisat2-mapping-se.yml", dir_path=dir_path)
align <- renderWF(align, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
align
cmdlist(align)[1:2]
output(align)[1:2]


## ----subset, eval=TRUE---------------------------------------------------
subsetWF(align, slot="input", subset='FileName')[1:2] ## Subsetting the input files for this particular workflow 
subsetWF(align, slot="output", subset=1)[1:2]  ## Subsetting the output files for one particular step in the workflow 
subsetWF(align, slot="step", subset=1)[1] ## Subsetting the command-lines for one particular step in the workflow 
subsetWF(align, slot="output", subset=1, delete=TRUE)[1] ## DELETING specific output files


## ----hisat_index, eval=FALSE---------------------------------------------
## dir_path <- system.file("extdata/cwl/hisat2/hisat2-idx", package="systemPipeR")
## idx <- loadWorkflow(targets=NULL, wf_file="hisat2-index.cwl", input_file="hisat2-index.yml", dir_path=dir_path)
## idx <- renderWF(idx)
## idx
## cmdlist(idx)
## 
## ## Run
## runCommandline(idx, make_bam = FALSE)


## ----runCommandline_align, eval=FALSE------------------------------------
## runCommandline(align, make_bam = FALSE) ## generates alignments and writes *.sam files to ./results folder
## align <- runCommandline(align, make_bam = TRUE) ## same as above but writes files and converts *.sam files to sorted and indexed BAM files. Assigning the new extention of the output files to the object align.


## ----clusterRun_align, eval=FALSE----------------------------------------
## library(batchtools)
## resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024)
## reg <- clusterRun(align, FUN = runCommandline, more.args = list(args=align, make_bam=TRUE, dir=FALSE),
##                   conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl",
##                   Njobs=18, runid="01", resourceList=resources)
## getStatus(reg=reg)
## waitForJobs(reg=reg)


## ----output, eval=FALSE--------------------------------------------------
## align <- output_update(align, dir=FALSE, replace=TRUE, extension=c(".sam", ".bam")) ## Updates the output(align) to the right location in the subfolders
## output(align)


## ----writeTargetsout, eval=FALSE-----------------------------------------
## names(clt(align))
## writeTargetsout(x=align, file="default", step=1)


## ----hisat_alignment, eval=FALSE-----------------------------------------
## targets <- system.file("extdata", "targets.txt", package="systemPipeR")
## dir_path <- system.file("extdata/cwl/workflow-hisat2/workflow-hisat2-se", package="systemPipeR")
## WF <- loadWorkflow(targets=targets, wf_file="workflow_hisat2-se.cwl", input_file="workflow_hisat2-se.yml", dir_path=dir_path)
## WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
## WF
## cmdlist(WF)[1:2]
## output(WF)[1:2]


## ----bowtie_index, eval=FALSE--------------------------------------------
## dir_path <- system.file("extdata/cwl/bowtie2/bowtie2-idx", package="systemPipeR")
## idx <- loadWorkflow(targets=NULL, wf_file="bowtie2-index.cwl", input_file="bowtie2-index.yml", dir_path=dir_path)
## idx <- renderWF(idx)
## idx
## cmdlist(idx)
## 
## ## Run in single machine
## runCommandline(idx, make_bam = FALSE)


## ----tophat2-pe, eval=FALSE----------------------------------------------
## targetsPE <- system.file("extdata", "targetsPE.txt", package = "systemPipeR")
## dir_path <- system.file("extdata/cwl/tophat2/tophat2-pe", package="systemPipeR")
## tophat2PE <- loadWorkflow(targets = targetsPE, wf_file = "tophat2-mapping-pe.cwl",
##     input_file = "tophat2-mapping-pe.yml", dir_path = dir_path)
## tophat2PE <- renderWF(tophat2PE, inputvars = c(FileName1 = "_FASTQ_PATH1_", FileName2 = "_FASTQ_PATH2_",
##     SampleName = "_SampleName_"))
## tophat2PE
## cmdlist(tophat2PE)[1:2]
## output(tophat2PE)[1:2]
## 
## ## Run in single machine
## tophat2PE <- runCommandline(tophat2PE[1], make_bam = TRUE)


## ----tophat2-pe_parallel, eval=FALSE-------------------------------------
## resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024)
## reg <- clusterRun(tophat2PE, FUN = runCommandline, more.args = list(args=tophat2PE, make_bam=TRUE, dir=FALSE),
##                   conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl",
##                   Njobs=18, runid="01", resourceList=resources)
## waitForJobs(reg=reg)


## ----writeTargetsout_tophat2PE, eval=FALSE-------------------------------
## names(clt(tophat2PE))
## writeTargetsout(x=tophat2PE, file="default", step=1)


## ----bowtie2_index, eval=FALSE-------------------------------------------
## dir_path <- system.file("extdata/cwl/bowtie2/bowtie2-idx", package="systemPipeR")
## idx <- loadWorkflow(targets=NULL, wf_file="bowtie2-index.cwl", input_file="bowtie2-index.yml", dir_path=dir_path)
## idx <- renderWF(idx)
## idx
## cmdlist(idx)
## 
## ## Run in single machine
## runCommandline(idx, make_bam = FALSE)


## ----bowtie2_SYSargs2, eval=FALSE----------------------------------------
## targetsPE <- system.file("extdata", "targetsPE.txt", package = "systemPipeR")
## dir_path <- system.file("extdata/cwl/bowtie2/bowtie2-pe", package="systemPipeR")
## bowtiePE <- loadWorkflow(targets = targetsPE, wf_file = "bowtie2-mapping-pe.cwl",
##     input_file = "bowtie2-mapping-pe.yml", dir_path = dir_path)
## bowtiePE <- renderWF(bowtiePE, inputvars = c(FileName1 = "_FASTQ_PATH1_", FileName2 = "_FASTQ_PATH2_",
##     SampleName = "_SampleName_"))
## bowtiePE
## cmdlist(bowtiePE)[1:2]
## output(bowtiePE)[1:2]


## ----bowtie2_cluster, eval=FALSE-----------------------------------------
## resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024)
## reg <- clusterRun(bowtiePE, FUN = runCommandline, more.args = list(args=bowtiePE, dir = FALSE),
##     conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl",
##     Njobs = 18, runid = "01", resourceList = resources)
## getStatus(reg = reg)


## ----bowtie2_sm, eval=FALSE----------------------------------------------
## bowtiePE <- runCommandline(bowtiePE)


## ----writeTargetsout_bowtiePE, eval=FALSE--------------------------------
## names(clt(bowtiePE))
## writeTargetsout(x=bowtiePE, file="default", step=1)


## ----bwa_index, eval=FALSE-----------------------------------------------
## dir_path <- system.file("extdata/cwl/bwa/bwa-idx", package="systemPipeR")
## idx <- loadWorkflow(targets=NULL, wf_file="bwa-index.cwl", input_file="bwa-index.yml", dir_path=dir_path)
## idx <- renderWF(idx)
## idx
## cmdlist(idx) # Indexes reference genome
## 
## ## Run
## runCommandline(idx, make_bam = FALSE)


## ----bwa-pe_alignment, eval=FALSE----------------------------------------
## targetsPE <- system.file("extdata", "targetsPE.txt", package = "systemPipeR")
## dir_path <- system.file("extdata/cwl/bwa/bwa-pe", package="systemPipeR")
## bwaPE <- loadWorkflow(targets = targetsPE, wf_file = "bwa-pe.cwl",
##     input_file = "bwa-pe.yml", dir_path = dir_path)
## bwaPE <- renderWF(bwaPE, inputvars = c(FileName1 = "_FASTQ_PATH1_", FileName2 = "_FASTQ_PATH2_",
##     SampleName = "_SampleName_"))
## bwaPE
## cmdlist(bwaPE)[1:2]
## output(bwaPE)[1:2]
## ## Single Machine
## bwaPE <- runCommandline(args= bwaPE, make_bam=FALSE)
## 
## ## Cluster
## library(batchtools)
## resources <- list(walltime = 120, ntasks = 1, ncpus = 4, memory = 1024)
## reg <- clusterRun(bwaPE, FUN = runCommandline, more.args = list(args=bwaPE, dir = FALSE),
##     conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl",
##     Njobs = 18, runid = "01", resourceList = resources)
## getStatus(reg = reg)


## ----writeTargetsout_bwaPE, eval=FALSE-----------------------------------
## names(clt(bwaPE))
## writeTargetsout(x=bwaPE, file="default", step=1)


## ----rsubread, eval=FALSE------------------------------------------------
## ## Build the index:
## dir_path <- system.file("extdata/cwl/rsubread/rsubread-idx", package="systemPipeR")
## idx <- loadWorkflow(targets = NULL, wf_file = "rsubread-index.cwl",
##                      input_file = "rsubread-index.yml", dir_path = dir_path)
## idx <- renderWF(idx)
## idx
## cmdlist(idx)
## runCommandline(args= idx, make_bam = FALSE)
## 
## ## Running the alignment:
## targets <- system.file("extdata", "targets.txt", package = "systemPipeR")
## dir_path <- system.file("extdata/cwl/rsubread/rsubread-se", package="systemPipeR")
## rsubread <- loadWorkflow(targets = targets, wf_file = "rsubread-mapping-se.cwl",
##                      input_file = "rsubread-mapping-se.yml", dir_path = dir_path)
## rsubread <- renderWF(rsubread, inputvars = c(FileName = "_FASTQ_PATH1_", SampleName = "_SampleName_"))
## rsubread
## cmdlist(rsubread)[1]
## 
## ## Single Machine
## rsubread <- runCommandline(args=rsubread[1])


## ----writeTargetsout_rsubread, eval=FALSE--------------------------------
## names(clt(rsubread))
## writeTargetsout(x=rsubread, file="default", step=1)


## ----gsnap, eval=FALSE---------------------------------------------------
## ## Build the index:
## dir_path <- system.file("extdata/cwl/gsnap/gsnap-idx", package="systemPipeR")
## idx <- loadWorkflow(targets = NULL, wf_file = "gsnap-index.cwl",
##                      input_file = "gsnap-index.yml", dir_path = dir_path)
## idx <- renderWF(idx)
## idx
## cmdlist(idx)
## runCommandline(args= idx, make_bam = FALSE)
## 
## ## Running the alignment:
## targetsPE <- system.file("extdata", "targetsPE.txt", package = "systemPipeR")
## dir_path <- system.file("extdata/cwl/gsnap/gsnap-pe", package="systemPipeR")
## gsnap <- loadWorkflow(targets = targetsPE, wf_file = "gsnap-mapping-pe.cwl",
##                      input_file = "gsnap-mapping-pe.yml", dir_path = dir_path)
## gsnap <- renderWF(gsnap, inputvars = c(FileName1 = "_FASTQ_PATH1_", FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_"))
## gsnap
## cmdlist(gsnap)[1]
## output(gsnap)[1]
## 
## ## Cluster
## library(batchtools)
## resources <- list(walltime = 120, ntasks = 1, ncpus = 4, memory = 1024)
## reg <- clusterRun(gsnap, FUN = runCommandline, more.args = list(args=gsnap, make_bam=FALSE),
##     conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl",
##     Njobs = 18, runid = "01", resourceList = resources)
## getStatus(reg = reg)


## ----writeTargetsout_gsnap, eval=FALSE-----------------------------------
## names(clt(gsnap))
## writeTargetsout(x=gsnap, file="default", step=1)


## ----igv, eval=FALSE-----------------------------------------------------
## symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"),
##             urlbase="http://myserver.edu/~username/",
##         urlfile="IGVurl.txt")


## ----create_txdb, eval=FALSE---------------------------------------------
## library(GenomicFeatures)
## txdb <- makeTxDbFromGFF(file="data/tair10.gff", format="gff", dataSource="TAIR", organism="Arabidopsis thaliana")
## saveDb(txdb, file="./data/tair10.sqlite")


## ----read_counting_multicore, eval=FALSE---------------------------------
## library(BiocParallel)
## txdb <- loadDb("./data/tair10.sqlite")
## eByg <- exonsBy(txdb, by="gene")
## outpaths <- subsetWF(align, slot="output", subset=1)
## bfl <- BamFileList(outpaths, yieldSize=50000, index=character())
## multicoreParam <- MulticoreParam(workers=4); register(multicoreParam); registered()
## counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", ignore.strand=TRUE, inter.feature=TRUE, singleEnd=TRUE))
## 
## # Note: for strand-specific RNA-Seq set 'ignore.strand=FALSE' and for PE data set 'singleEnd=FALSE'
## countDFeByg <- sapply(seq(along=counteByg),
##                       function(x) assays(counteByg[[x]])$counts)
## rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
## rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
## write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
## write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")


## ----read_counting_multinode, eval=FALSE---------------------------------
## library(BiocParallel)
## f <- function(x) {
##     library(systemPipeR); library(BiocParallel); library(GenomicFeatures)
##     txdb <- loadDb("./data/tair10.sqlite")
##     eByg <- exonsBy(txdb, by="gene")
##     args <- systemArgs(sysma="param/tophat.param", mytargets="targets.txt")
##     outpaths <- subsetWF(align, slot="output", subset=1)
##     bfl <- BamFileList(outpaths, yieldSize=50000, index=character())
##     summarizeOverlaps(eByg, bfl[x], mode="Union", ignore.strand=TRUE, inter.feature=TRUE, singleEnd=TRUE)
## }
## resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024)
## param <- BatchtoolsParam(workers = 4, cluster = "slurm", template = "batchtools.slurm.tmpl", resources = resources)
## counteByg <- bplapply(seq(along=args), f, BPPARAM = param)
## countDFeByg <- sapply(seq(along=counteByg),
##                       function(x) assays(counteByg[[x]])$counts)
## rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(outpaths)


## ----process_monitoring, eval=FALSE--------------------------------------
## getStatus(reg=reg)
## file.exists(output(tophat2PE))
## sapply(1:length(tophat2PE), function(x) loadResult(reg, id=x)) # Works after job completion


## ----align_stats1, eval=FALSE--------------------------------------------
## read_statsDF <- alignStats(align)
## write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")


## ----align_stats2, eval=TRUE---------------------------------------------
read.table(system.file("extdata", "alignStats.xls", package="systemPipeR"), header=TRUE)[1:4,]


## ----align_stats_parallel, eval=FALSE------------------------------------
## f <- function(x) alignStats(align[x])
## read_statsList <- bplapply(seq(along=align), f,
##                            BPPARAM = MulticoreParam(workers=8))
## read_statsDF <- do.call("rbind", read_statsList)


## ----align_stats_parallel_cluster, eval=FALSE----------------------------
## library(BiocParallel); library(batchtools)
## f <- function(x) {
##     library(systemPipeR)
##     targets <- system.file("extdata", "targets.txt", package="systemPipeR")
##     dir_path <- "param/cwl/hisat2/hisat2-se" ## TODO: replace path to system.file
##     align <- loadWorkflow(targets=targets, wf_file="hisat2-mapping-se.cwl", input_file="hisat2-mapping-se.yml", dir_path=dir_path)
##     align <- renderWF(align, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
##     align <- output_update(align, dir=FALSE, replace=TRUE, extension=c(".sam", ".bam"))
##     alignStats(align[x])
## }
## resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024)
## param <- BatchtoolsParam(workers = 4, cluster = "slurm", template = "batchtools.slurm.tmpl", resources = resources)
## read_statsList <- bplapply(seq(along=align), f, BPPARAM = param)
## read_statsDF <- do.call("rbind", read_statsList)


## ----read_counting_mirna, eval=FALSE-------------------------------------
## system("wget ftp://mirbase.org/pub/mirbase/19/genomes/My_species.gff3 -P ./data/")
## gff <- import.gff("./data/My_species.gff3")
## gff <- split(gff, elementMetadata(gff)$ID)
## bams <- names(bampaths); names(bams) <- targets$SampleName
## bfl <- BamFileList(bams, yieldSize=50000, index=character())
## countDFmiR <- summarizeOverlaps(gff, bfl, mode="Union", ignore.strand=FALSE, inter.feature=FALSE) # Note: inter.feature=FALSE important since pre and mature miRNA ranges overlap
## rpkmDFmiR <- apply(countDFmiR, 2,
##                    function(x) returnRPKM(counts=x, gffsub=gff))
## write.table(assays(countDFmiR)$counts, "results/countDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")
## write.table(rpkmDFmiR, "results/rpkmDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")


## ----sample_tree_rlog, eval=TRUE-----------------------------------------
library(DESeq2, warn.conflicts=FALSE, quietly=TRUE); library(ape, warn.conflicts=FALSE)
countDFpath <- system.file("extdata", "countDFeByg.xls", package="systemPipeR")
countDF <- as.matrix(read.table(countDFpath))
colData <- data.frame(row.names=targetsin(args)$SampleName, condition=targetsin(args)$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
d <- cor(assay(rlog(dds)), method="spearman")
hc <- hclust(dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)


## ----sample_tree_rpkm, eval=FALSE----------------------------------------
## rpkmDFeBygpath <- system.file("extdata", "rpkmDFeByg.xls", package="systemPipeR")
## rpkmDFeByg <- read.table(rpkmDFeBygpath, check.names=FALSE)
## rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
## d <- cor(rpkmDFeByg, method="spearman")
## hc <- hclust(as.dist(1-d))
## plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)


## ----edger_wrapper, eval=TRUE--------------------------------------------
targets <- read.delim(targetspath, comment="#")
cmp <- readComp(file=targetspath, format="matrix", delim="-")
cmp[[1]]
countDFeBygpath <- system.file("extdata", "countDFeByg.xls", package="systemPipeR")
countDFeByg <- read.delim(countDFeBygpath, row.names=1)
edgeDF <- run_edgeR(countDF=countDFeByg, targets=targets, cmp=cmp[[1]], independent=FALSE, mdsplot="")


## ----edger_deg_counts, eval=TRUE-----------------------------------------
DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=10))


## ----edger_deg_stats, eval=TRUE------------------------------------------
names(DEG_list)
DEG_list$Summary[1:4,]


## ----deseq2_wrapper, eval=TRUE-------------------------------------------
degseqDF <- run_DESeq2(countDF=countDFeByg, targets=targets, cmp=cmp[[1]], independent=FALSE)


## ----deseq2_deg_counts, eval=TRUE----------------------------------------
DEG_list2 <- filterDEGs(degDF=degseqDF, filter=c(Fold=2, FDR=10))


## ----vennplot, eval=TRUE-------------------------------------------------
vennsetup <- overLapper(DEG_list$Up[6:9], type="vennsets")
vennsetdown <- overLapper(DEG_list$Down[6:9], type="vennsets")
vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))


## ----get_go_biomart, eval=FALSE------------------------------------------
## library("biomaRt")
## listMarts() # To choose BioMart database
## listMarts(host="plants.ensembl.org")
## m <- useMart("plants_mart", host="plants.ensembl.org")
## listDatasets(m)
## m <- useMart("plants_mart", dataset="athaliana_eg_gene", host="plants.ensembl.org")
## listAttributes(m) # Choose data types you want to download
## go <- getBM(attributes=c("go_id", "tair_locus", "namespace_1003"), mart=m)
## go <- go[go[,3]!="",]; go[,3] <- as.character(go[,3])
## go[go[,3]=="molecular_function", 3] <- "F"; go[go[,3]=="biological_process", 3] <- "P"; go[go[,3]=="cellular_component", 3] <- "C"
## go[1:4,]
## dir.create("./data/GO")
## write.table(go, "data/GO/GOannotationsBiomart_mod.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
## catdb <- makeCATdb(myfile="data/GO/GOannotationsBiomart_mod.txt", lib=NULL, org="", colno=c(1,2,3), idconv=NULL)
## save(catdb, file="data/GO/catdb.RData")


## ----go_enrichment, eval=FALSE-------------------------------------------
## load("data/GO/catdb.RData")
## DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=50), plot=FALSE)
## up_down <- DEG_list$UporDown; names(up_down) <- paste(names(up_down), "_up_down", sep="")
## up <- DEG_list$Up; names(up) <- paste(names(up), "_up", sep="")
## down <- DEG_list$Down; names(down) <- paste(names(down), "_down", sep="")
## DEGlist <- c(up_down, up, down)
## DEGlist <- DEGlist[sapply(DEGlist, length) > 0]
## BatchResult <- GOCluster_Report(catdb=catdb, setlist=DEGlist, method="all", id_type="gene", CLSZ=2, cutoff=0.9, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
## library("biomaRt")
## m <- useMart("plants_mart", dataset="athaliana_eg_gene", host="plants.ensembl.org")
## goslimvec <- as.character(getBM(attributes=c("goslim_goa_accession"), mart=m)[,1])
## BatchResultslim <- GOCluster_Report(catdb=catdb, setlist=DEGlist, method="slim", id_type="gene", myslimv=goslimvec, CLSZ=10, cutoff=0.01, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)


## ----plot_go_enrichment, eval=FALSE--------------------------------------
## gos <- BatchResultslim[grep("M6-V6_up_down", BatchResultslim$CLID), ]
## gos <- BatchResultslim
## pdf("GOslimbarplotMF.pdf", height=8, width=10); goBarplot(gos, gocat="MF"); dev.off()
## goBarplot(gos, gocat="BP")
## goBarplot(gos, gocat="CC")


## ----hierarchical_clustering, eval=FALSE---------------------------------
## library(pheatmap)
## geneids <- unique(as.character(unlist(DEG_list[[1]])))
## y <- assay(rlog(dds))[geneids, ]
## pdf("heatmap1.pdf")
## pheatmap(y, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
## dev.off()


## Rscript -e "rmarkdown::render('systemPipeRNAseq.Rmd')"


## ----genRna_workflow_single, eval=FALSE----------------------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="rnaseq")
## setwd("rnaseq")


## ----genChip_workflow_single, eval=FALSE---------------------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="chipseq")
## setwd("chipseq")


## ----genVar_workflow_single, eval=FALSE----------------------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="varseq")
## setwd("varseq")


## ----genRibo_workflow_single, eval=FALSE---------------------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="riboseq")
## setwd("riboseq")


## ----sessionInfo---------------------------------------------------------
sessionInfo()

