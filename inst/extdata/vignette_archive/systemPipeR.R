## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")))

## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE-------------------
suppressPackageStartupMessages({
    library(systemPipeR)
    library(systemPipeRdata)
    library(BiocParallel)
    library(Biostrings)
    library(Rsamtools)
    library(GenomicRanges)
    library(ggplot2)
    library(GenomicAlignments)
    library(ShortRead)
    library(ape)
})

## ----install, eval=FALSE-------------------------------------------------
## source("http://bioconductor.org/biocLite.R") # Sources the biocLite.R installation script
## biocLite("systemPipeR") # Installs systemPipeR from Bioconductor
## biocLite("tgirke/systemPipeRdata", build_vignettes=TRUE, dependencies=TRUE) # From github

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

## ----comment_lines, eval=TRUE--------------------------------------------
readLines(targetspath)[1:4]

## ----targetscomp, eval=TRUE----------------------------------------------
readComp(file=targetspath, format="vector", delim="-")

## ----param_structure, eval=TRUE------------------------------------------
parampath <- system.file("extdata", "tophat.param", package="systemPipeR")
read.delim(parampath, comment.char = "#")

## ----param_import, eval=TRUE---------------------------------------------
args <- suppressWarnings(systemArgs(sysma=parampath, mytargets=targetspath))
args

## ----sysarg_access, eval=TRUE--------------------------------------------
names(args)
modules(args)
cores(args)
outpaths(args)[1]
sysargs(args)[1]

## ----sysarg_json, eval=TRUE----------------------------------------------
systemArgs(sysma=parampath, mytargets=targetspath, type="json")

## ----load_package, eval=FALSE--------------------------------------------
## library(systemPipeR)
## library(systemPipeRdata)
## genWorkenvir(workflow="rnaseq")
## setwd("rnaseq")

## ----construct_sysargs, eval=FALSE---------------------------------------
## args <- systemArgs(sysma="param/trim.param", mytargets="targets.txt")

## ----preprocessing, eval=FALSE-------------------------------------------
## preprocessReads(args=args, Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)",
##                 batchsize=100000, overwrite=TRUE, compress=TRUE)
## writeTargetsout(x=args, file="targets_trim.txt")

## ----custom_preprocessing, eval=FALSE------------------------------------
## args <- systemArgs(sysma="param/trimPE.param", mytargets="targetsPE.txt")
## filterFct <- function(fq, cutoff=20, Nexceptions=0) {
##     qcount <- rowSums(as(quality(fq), "matrix") <= cutoff)
##     fq[qcount <= Nexceptions] # Retains reads where Phred scores are >= cutoff with N exceptions
## }
## preprocessReads(args=args, Fct="filterFct(fq, cutoff=20, Nexceptions=0)", batchsize=100000)
## writeTargetsout(x=args, file="targets_PEtrim.txt")

## ----fastq_quality, eval=FALSE-------------------------------------------
## fqlist <- seeFastq(fastq=infile1(args), batchsize=10000, klength=8)
## pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
## seeFastqPlot(fqlist)
## dev.off()

## ----fastq_quality_parallel_single, eval=FALSE---------------------------
## args <- systemArgs(sysma="param/tophat.param", mytargets="targets.txt")
## f <- function(x) seeFastq(fastq=infile1(args)[x], batchsize=100000, klength=8)
## fqlist <- bplapply(seq(along=args), f, BPPARAM = MulticoreParam(workers=8))
## seeFastqPlot(unlist(fqlist, recursive=FALSE))

## ----fastq_quality_parallel_cluster, eval=FALSE--------------------------
## library(BiocParallel); library(BatchJobs)
## f <- function(x) {
##     library(systemPipeR)
##     args <- systemArgs(sysma="param/tophat.param", mytargets="targets.txt")
##     seeFastq(fastq=infile1(args)[x], batchsize=100000, klength=8)
## }
## funs <- makeClusterFunctionsTorque("torque.tmpl")
## param <- BatchJobsParam(length(args), resources=list(walltime="20:00:00", nodes="1:ppn=1", memory="6gb"), cluster.functions=funs)
## register(param)
## fqlist <- bplapply(seq(along=args), f)
## seeFastqPlot(unlist(fqlist, recursive=FALSE))

## ----bowtie_index, eval=FALSE--------------------------------------------
## args <- systemArgs(sysma="param/tophat.param", mytargets="targets.txt")
## moduleload(modules(args)) # Skip if module system is not available
## system("bowtie2-build ./data/tair10.fasta ./data/tair10.fasta")

## ----run_bowtie_single, eval=FALSE---------------------------------------
## bampaths <- runCommandline(args=args)

## ----run_bowtie_parallel, eval=FALSE-------------------------------------
## resources <- list(walltime="20:00:00", nodes=paste0("1:ppn=", cores(args)), memory="10gb")
## reg <- clusterRun(args, conffile=".BatchJobs.R", template="torque.tmpl", Njobs=18, runid="01",
##                   resourceList=resources)
## waitForJobs(reg)

## ----process_monitoring, eval=FALSE--------------------------------------
## showStatus(reg)
## file.exists(outpaths(args))
## sapply(1:length(args), function(x) loadResult(reg, x)) # Works after job completion

## ----align_stats1, eval=FALSE--------------------------------------------
## read_statsDF <- alignStats(args)
## write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")

## ----align_stats2, eval=TRUE---------------------------------------------
read.table(system.file("extdata", "alignStats.xls", package="systemPipeR"), header=TRUE)[1:4,]

## ----align_stats_parallel, eval=FALSE------------------------------------
## f <- function(x) alignStats(args[x])
## read_statsList <- bplapply(seq(along=args), f, BPPARAM = MulticoreParam(workers=8))
## read_statsDF <- do.call("rbind", read_statsList)

## ----align_stats_parallel_cluster, eval=FALSE----------------------------
## library(BiocParallel); library(BatchJobs)
## f <- function(x) {
##     library(systemPipeR)
##     args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
##     alignStats(args[x])
## }
## funs <- makeClusterFunctionsTorque("torque.tmpl")
## param <- BatchJobsParam(length(args), resources=list(walltime="20:00:00", nodes="1:ppn=1", memory="6gb"), cluster.functions=funs)
## register(param)
## read_statsList <- bplapply(seq(along=args), f)
## read_statsDF <- do.call("rbind", read_statsList)

## ----igv, eval=FALSE-----------------------------------------------------
## symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"),
##             urlbase="http://myserver.edu/~username/",
##         urlfile="IGVurl.txt")

## ----bowtie2, eval=FALSE-------------------------------------------------
## args <- systemArgs(sysma="bowtieSE.param", mytargets="targets.txt")
## moduleload(modules(args)) # Skip if module system is not available
## bampaths <- runCommandline(args=args)

## ----bowtie2_cluster, eval=FALSE-----------------------------------------
## resources <- list(walltime="20:00:00", nodes=paste0("1:ppn=", cores(args)), memory="10gb")
## reg <- clusterRun(args, conffile=".BatchJobs.R", template="torque.tmpl", Njobs=18, runid="01",
##                   resourceList=resources)
## waitForJobs(reg)

## ----bwamem_cluster, eval=FALSE------------------------------------------
## args <- systemArgs(sysma="param/bwa.param", mytargets="targets.txt")
## moduleload(modules(args)) # Skip if module system is not available
## system("bwa index -a bwtsw ./data/tair10.fasta") # Indexes reference genome
## bampaths <- runCommandline(args=args[1:2])

## ----rsubread, eval=FALSE------------------------------------------------
## library(Rsubread)
## args <- systemArgs(sysma="param/rsubread.param", mytargets="targets.txt")
## buildindex(basename=reference(args), reference=reference(args)) # Build indexed reference genome
## align(index=reference(args), readfile1=infile1(args)[1:4], input_format="FASTQ",
##       output_file=outfile1(args)[1:4], output_format="SAM", nthreads=8, indels=1, TH1=2)
## for(i in seq(along=outfile1(args))) asBam(file=outfile1(args)[i], destination=gsub(".sam", "", outfile1(args)[i]), overwrite=TRUE, indexDestination=TRUE)

## ----gsnap, eval=FALSE---------------------------------------------------
## library(gmapR); library(BiocParallel); library(BatchJobs)
## args <- systemArgs(sysma="param/gsnap.param", mytargets="targetsPE.txt")
## gmapGenome <- GmapGenome(reference(args), directory="data", name="gmap_tair10chr/", create=TRUE)
## f <- function(x) {
##     library(gmapR); library(systemPipeR)
##     args <- systemArgs(sysma="gsnap.param", mytargets="targetsPE.txt")
##     gmapGenome <- GmapGenome(reference(args), directory="data", name="gmap_tair10chr/", create=FALSE)
##     p <- GsnapParam(genome=gmapGenome, unique_only=TRUE, molecule="DNA", max_mismatches=3)
##     o <- gsnap(input_a=infile1(args)[x], input_b=infile2(args)[x], params=p, output=outfile1(args)[x])
## }
## funs <- makeClusterFunctionsTorque("torque.tmpl")
## param <- BatchJobsParam(length(args), resources=list(walltime="20:00:00", nodes="1:ppn=1", memory="6gb"), cluster.functions=funs)
## register(param)
## d <- bplapply(seq(along=args), f)

## ----create_txdb, eval=FALSE---------------------------------------------
## library(gmapR); library(BiocParallel); library(BatchJobs)
## library(GenomicFeatures)
## txdb <- makeTxDbFromGFF(file="data/tair10.gff", format="gff", dataSource="TAIR", organism="A. thaliana")
## saveDb(txdb, file="./data/tair10.sqlite")

## ----read_counting_multicore, eval=FALSE---------------------------------
## library(BiocParallel)
## txdb <- loadDb("./data/tair10.sqlite")
## eByg <- exonsBy(txdb, by="gene")
## bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
## multicoreParam <- MulticoreParam(workers=4); register(multicoreParam); registered()
## counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", ignore.strand=TRUE, inter.feature=TRUE, singleEnd=TRUE)) # Note: for strand-specific RNA-Seq set 'ignore.strand=FALSE' and for PE data set 'singleEnd=FALSE'
## countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
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
##     args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
##     bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
##     summarizeOverlaps(eByg, bfl[x], mode="Union", ignore.strand=TRUE, inter.feature=TRUE, singleEnd=TRUE)
## }
## funs <- makeClusterFunctionsTorque("torque.tmpl")
## param <- BatchJobsParam(length(args), resources=list(walltime="20:00:00", nodes="1:ppn=1", memory="6gb"), cluster.functions=funs)
## register(param)
## counteByg <- bplapply(seq(along=args), f)
## countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
## rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(outpaths(args))

## ----read_counting_mirna, eval=FALSE-------------------------------------
## system("wget ftp://mirbase.org/pub/mirbase/19/genomes/My_species.gff3 -P ./data/")
## gff <- import.gff("./data/My_species.gff3", asRangedData=FALSE)
## gff <- split(gff, elementMetadata(gff)$ID)
## bams <- names(bampaths); names(bams) <- targets$SampleName
## bfl <- BamFileList(bams, yieldSize=50000, index=character())
## countDFmiR <- summarizeOverlaps(gff, bfl, mode="Union", ignore.strand=FALSE, inter.feature=FALSE) # Note: inter.feature=FALSE important since pre and mature miRNA ranges overlap
## rpkmDFmiR <- apply(countDFmiR, 2, function(x) returnRPKM(counts=x, gffsub=gff))
## write.table(assays(countDFmiR)$counts, "results/countDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")
## write.table(rpkmDFmiR, "results/rpkmDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")

## ----sample_tree_, eval=TRUE---------------------------------------------
library(ape,  warn.conflicts=FALSE)
rpkmDFeBygpath <- system.file("extdata", "rpkmDFeByg.xls", package="systemPipeR")
rpkmDFeByg <- read.table(rpkmDFeBygpath, check.names=FALSE)
rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
d <- cor(rpkmDFeByg, method="spearman")
hc <- hclust(as.dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)

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
## m <- useMart("ENSEMBL_MART_PLANT"); listDatasets(m)
## m <- useMart("ENSEMBL_MART_PLANT", dataset="athaliana_eg_gene")
## listAttributes(m) # Choose data types you want to download
## go <- getBM(attributes=c("go_accession", "tair_locus", "go_namespace_1003"), mart=m)
## go <- go[go[,3]!="",]; go[,3] <- as.character(go[,3])
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
## library("biomaRt"); m <- useMart("ENSEMBL_MART_PLANT", dataset="athaliana_eg_gene")
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
## y <- rpkmDFeByg[geneids, ]
## pdf("heatmap1.pdf")
## pheatmap(y, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
## dev.off()

## ----genVar_workflow_single, eval=FALSE----------------------------------
## setwd("../")
## genWorkenvir(workflow="varseq")
## setwd("varseq")

## ----sessionInfo---------------------------------------------------------
sessionInfo()

