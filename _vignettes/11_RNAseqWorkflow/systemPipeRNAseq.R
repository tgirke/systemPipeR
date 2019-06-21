## pre code {

## white-space: pre !important;

## overflow-x: scroll !important;

## word-break: keep-all !important;

## word-wrap: initial !important;

## }


## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()
options(width=60, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")), 
    tidy.opts=list(width.cutoff=60), tidy=TRUE)


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


## ----genRna_workflow, eval=FALSE-----------------------------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="rnaseq")
## setwd("rnaseq")


## Rscript -e "systemPipeRdata::genWorkenvir(workflow='rnaseq')"


## ----closeR, eval=FALSE--------------------------------------------------
## q("no") # closes R session on head node


## srun --x11 --partition=short --mem=2gb --cpus-per-task 4 --ntasks 1 --time 2:00:00 --pty bash -l

## module load R/3.4.2

## R


## ----r_environment, eval=FALSE-------------------------------------------
## system("hostname") # should return name of a compute node starting with i or c
## getwd() # checks current working directory of R session
## dir() # returns content of current working directory


## ----load_systempiper, eval=TRUE-----------------------------------------
library(systemPipeR)


## ----source_helper_fcts, eval=FALSE--------------------------------------
## source("systemPipeRNAseq_Fct.R")


## ----load_targets, eval=TRUE---------------------------------------------
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")[,1:4]
targets


## ----fastq_filter, eval=FALSE--------------------------------------------
## args <- systemArgs(sysma="param/trim.param", mytargets="targets.txt")
## preprocessReads(args=args, Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)",
##                 batchsize=100000, overwrite=TRUE, compress=TRUE)
## writeTargetsout(x=args, file="targets_trim.txt", overwrite=TRUE)


## ----fastq_report, eval=FALSE--------------------------------------------
## args <- systemArgs(sysma="param/tophat.param", mytargets="targets.txt")
## fqlist <- seeFastq(fastq=infile1(args), batchsize=100000, klength=8)
## pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
## seeFastqPlot(fqlist)
## dev.off()


## ----tophat_alignment1, eval=FALSE---------------------------------------
## args <- systemArgs(sysma="param/tophat.param", mytargets="targets.txt")
## sysargs(args)[1] # Command-line parameters for first FASTQ file


## ----tophat_alignment2, eval=FALSE---------------------------------------
## moduleload(modules(args))
## system("bowtie2-build ./data/tair10.fasta ./data/tair10.fasta")
## resources <- list(walltime=120, ntasks=1, ncpus=cores(args), memory=1024)
## reg <- clusterRun(args, conffile = ".batchtools.conf.R", Njobs=18, template = "batchtools.slurm.tmpl", runid="01", resourceList=resources)
## getStatus(reg=reg)
## waitForJobs(reg=reg)


## ----hisat_alignment2, eval=FALSE----------------------------------------
## args <- systemArgs(sysma="param/hisat2.param", mytargets="targets.txt")
## sysargs(args)[1] # Command-line parameters for first FASTQ file
## moduleload(modules(args))
## system("hisat2-build ./data/tair10.fasta ./data/tair10.fasta")
## resources <- list(walltime=120, ntasks=1, ncpus=cores(args), memory=1024)
## reg <- clusterRun(args, conffile = ".batchtools.conf.R", Njobs=18, template = "batchtools.slurm.tmpl", runid="01", resourceList=resources)
## getStatus(reg=reg)
## waitForJobs(reg=reg)


## ----check_files_exist, eval=FALSE---------------------------------------
## file.exists(outpaths(args))


## ----align_stats, eval=FALSE---------------------------------------------
## read_statsDF <- alignStats(args=args)
## write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")


## ----align_stats_view, eval=TRUE-----------------------------------------
read.table(system.file("extdata", "alignStats.xls", package="systemPipeR"), header=TRUE)[1:4,]


## ----bam_urls, eval=FALSE------------------------------------------------
## symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"),
##             urlbase="http://biocluster.ucr.edu/~tgirke/",
## 	        urlfile="./results/IGVurl.txt")


## ----read_counting1, eval=FALSE------------------------------------------
## library("GenomicFeatures"); library(BiocParallel)
## txdb <- makeTxDbFromGFF(file="data/tair10.gff", format="gff", dataSource="TAIR", organism="Arabidopsis thaliana")
## saveDb(txdb, file="./data/tair10.sqlite")
## txdb <- loadDb("./data/tair10.sqlite")
## (align <- readGAlignments(outpaths(args)[1])) # Demonstrates how to read bam file into R
## eByg <- exonsBy(txdb, by=c("gene"))
## bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
## multicoreParam <- MulticoreParam(workers=2); register(multicoreParam); registered()
## counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union",
##                                                ignore.strand=TRUE,
##                                                inter.feature=FALSE,
##                                                singleEnd=TRUE))
## countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
## rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
## rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
## write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
## write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")


## ----view_counts, eval=FALSE---------------------------------------------
## read.delim("results/countDFeByg.xls", row.names=1, check.names=FALSE)[1:4,1:5]


## ----view_rpkm, eval=FALSE-----------------------------------------------
## read.delim("results/rpkmDFeByg.xls", row.names=1, check.names=FALSE)[1:4,1:4]


## ----sample_tree, eval=FALSE---------------------------------------------
## library(DESeq2, quietly=TRUE); library(ape,  warn.conflicts=FALSE)
## countDF <- as.matrix(read.table("./results/countDFeByg.xls"))
## colData <- data.frame(row.names=targetsin(args)$SampleName, condition=targetsin(args)$Factor)
## dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
## d <- cor(assay(rlog(dds)), method="spearman")
## hc <- hclust(dist(1-d))
## pdf("results/sample_tree.pdf")
## plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
## dev.off()


## ----run_edger, eval=FALSE-----------------------------------------------
## library(edgeR)
## countDF <- read.delim("results/countDFeByg.xls", row.names=1, check.names=FALSE)
## targets <- read.delim("targets.txt", comment="#")
## cmp <- readComp(file="targets.txt", format="matrix", delim="-")
## edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=FALSE, mdsplot="")


## ----custom_annot, eval=FALSE--------------------------------------------
## library("biomaRt")
## m <- useMart("plants_mart", dataset="athaliana_eg_gene", host="plants.ensembl.org")
## desc <- getBM(attributes=c("tair_locus", "description"), mart=m)
## desc <- desc[!duplicated(desc[,1]),]
## descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
## edgeDF <- data.frame(edgeDF, Desc=descv[rownames(edgeDF)], check.names=FALSE)
## write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)


## ----filter_degs, eval=FALSE---------------------------------------------
## edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names=1, check.names=FALSE)
## pdf("results/DEGcounts.pdf")
## DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=20))
## dev.off()
## write.table(DEG_list$Summary, "./results/DEGcounts.xls", quote=FALSE, sep="\t", row.names=FALSE)


## ----venn_diagram, eval=FALSE--------------------------------------------
## vennsetup <- overLapper(DEG_list$Up[6:9], type="vennsets")
## vennsetdown <- overLapper(DEG_list$Down[6:9], type="vennsets")
## pdf("results/vennplot.pdf")
## vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))
## dev.off()


## ----get_go_annot, eval=FALSE--------------------------------------------
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


## ----go_enrich, eval=FALSE-----------------------------------------------
## library("biomaRt")
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


## ----go_plot, eval=FALSE-------------------------------------------------
## gos <- BatchResultslim[grep("M6-V6_up_down", BatchResultslim$CLID), ]
## gos <- BatchResultslim
## pdf("GOslimbarplotMF.pdf", height=8, width=10); goBarplot(gos, gocat="MF"); dev.off()
## goBarplot(gos, gocat="BP")
## goBarplot(gos, gocat="CC")


## ----heatmap, eval=FALSE-------------------------------------------------
## library(pheatmap)
## geneids <- unique(as.character(unlist(DEG_list[[1]])))
## y <- assay(rlog(dds))[geneids, ]
## pdf("heatmap1.pdf")
## pheatmap(y, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
## dev.off()


## ----sessionInfo---------------------------------------------------------
sessionInfo()

