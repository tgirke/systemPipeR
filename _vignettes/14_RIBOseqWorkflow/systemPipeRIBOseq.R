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


## ----genRibo_workflow, eval=FALSE-------------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="riboseq", bam=TRUE)
## setwd("riboseq")


## ----download_latest, eval=FALSE--------------------------
## download.file("https://raw.githubusercontent.com/tgirke/systemPipeRdata/master/vignettes/systemPipeRIBOseq.Rmd", "systemPipeRIBOseq.Rmd")


## ----load_systempiper, eval=TRUE,  messages=FALSE---------
library(systemPipeR)


## ----source_helper_fcts, eval=FALSE-----------------------
## source("systemPipeRIBOseq_Fct.R")


## ----load_targets, eval=TRUE------------------------------
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")[,1:4]
targets


## ----construct_SYSargs2_trim-se, eval=FALSE---------------
## dir_path <- system.file("extdata/cwl/preprocessReads/trim-se", package="systemPipeR")
## trim <- loadWorkflow(targets=targetspath, wf_file="trim-se.cwl", input_file="trim-se.yml", dir_path=dir_path)
## trim <- renderWF(trim, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
## trim
## output(trim)[1:2]


## ----fastq_trimming, eval=FALSE---------------------------
## fctpath <- system.file("extdata", "custom_Fct.R", package="systemPipeR")
## source(fctpath)
## iterTrim <- ".iterTrimbatch1(fq, pattern='ACACGTCT', internalmatch=FALSE, minpatternlength=6, Nnumber=1, polyhomo=50, minreadlength=16, maxreadlength=101)"
## preprocessReads(args=trim, Fct=iterTrim, batchsize=100000, overwrite=TRUE, compress=TRUE)
## writeTargetsout(x=trim, file="targets_trim.txt", step = 1,
##                 new_col = "FileName", new_col_output_index = 1, overwrite = TRUE)


## ----fastq_report, eval=FALSE-----------------------------
## fqlist <- seeFastq(fastq=infile1(trim), batchsize=10000, klength=8)
## png("./results/fastqReport.png", height=18, width=4*length(fqlist), units="in", res=72)
## seeFastqPlot(fqlist)
## dev.off()


## ----bowtie_index, eval=FALSE-----------------------------
## dir_path <- system.file("extdata/cwl/bowtie2/bowtie2-idx", package="systemPipeR")
## idx <- loadWorkflow(targets=NULL, wf_file="bowtie2-index.cwl", input_file="bowtie2-index.yml", dir_path=dir_path)
## idx <- renderWF(idx)
## idx
## cmdlist(idx)
## 
## ## Run in single machine
## runCommandline(idx, make_bam = FALSE)


## ----tophat2-se, eval=FALSE-------------------------------
## dir_path <- system.file("extdata/cwl/tophat2/tophat2-se", package="systemPipeR")
## args <- loadWorkflow(targets = targetspath, wf_file = "tophat2-mapping-se.cwl",
##     input_file = "tophat2-mapping-se.yml", dir_path = dir_path)
## args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", SampleName = "_SampleName_"))
## args
## cmdlist(args)[1:2]
## output(args)[1:2]
## 
## ## Run in single machine
## args <- runCommandline(args, make_bam = TRUE)


## ----tophat_se_cluster, eval=FALSE------------------------
## ## Run on the cluster
## moduleload(modules(args))
## resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024)
## reg <- clusterRun(args, FUN = runCommandline, more.args = list(args=args, make_bam=TRUE, dir=FALSE),
##                   conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl",
##                   Njobs=18, runid="01", resourceList=resources)
## waitForJobs(reg=reg)


## ----hisat_index, eval=FALSE------------------------------
## dir_path <- system.file("extdata/cwl/hisat2/hisat2-idx", package="systemPipeR")
## idx <- loadWorkflow(targets=NULL, wf_file="hisat2-index.cwl", input_file="hisat2-index.yml", dir_path=dir_path)
## idx <- renderWF(idx)
## idx
## cmdlist(idx)
## 
## ## Run
## runCommandline(idx, make_bam = FALSE)


## ----hisat_SYSargs2_object, eval=FALSE--------------------
## dir_path <- system.file("extdata/cwl/hisat2/hisat2-se", package="systemPipeR")
## args <- loadWorkflow(targets=targetspath, wf_file="hisat2-mapping-se.cwl",
##                      input_file="hisat2-mapping-se.yml", dir_path=dir_path)
## args <- renderWF(args, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
## args
## cmdlist(args)[1:2]
## output(args)[1:2]
## 
## ## Run
## args <- runCommandline(args)


## ----hisat2_clusterRun, eval=FALSE------------------------
## library(batchtools)
## resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024)
## reg <- clusterRun(args, FUN = runCommandline, more.args = list(args=args, make_bam=TRUE, dir=FALSE),
##                   conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl",
##                   Njobs=18, runid="01", resourceList=resources)
## getStatus(reg=reg)
## waitForJobs(reg=reg)


## ----check_files_exist, eval=FALSE------------------------
## outpaths <- subsetWF(args , slot="output", subset=1, index=1)
## file.exists(outpaths)


## ----align_stats, eval=FALSE------------------------------
## read_statsDF <- alignStats(args=args)
## write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")


## ----align_stats_view, eval=TRUE--------------------------
read.table(system.file("extdata", "alignStats.xls", package="systemPipeR"), header=TRUE)[1:4,]


## ----bam_urls, eval=FALSE---------------------------------
## symLink2bam(sysargs=args, htmldir=c("~/.html/", "projects/tests/"),
##             urlbase="http://biocluster.ucr.edu/~tgirke/",
##             urlfile="./results/IGVurl.txt")


## ----genFeatures, eval=FALSE------------------------------
## library(GenomicFeatures)
## txdb <- makeTxDbFromGFF(file="data/tair10.gff", format="gff3", organism="Arabidopsis")
## feat <- genFeatures(txdb, featuretype="all", reduce_ranges=TRUE, upstream=1000,
##                     downstream=0, verbose=TRUE)


## ----featuretypeCounts, eval=FALSE------------------------
## library(ggplot2); library(grid)
## outpaths <- subsetWF(args , slot="output", subset=1, index=1)
## fc <- featuretypeCounts(bfl=BamFileList(outpaths, yieldSize=50000), grl=feat,
##                         singleEnd=TRUE, readlength=NULL, type="data.frame")
## p <- plotfeaturetypeCounts(x=fc, graphicsfile="results/featureCounts.png",
##                            graphicsformat="png", scales="fixed", anyreadlength=TRUE,
##                            scale_length_val=NULL)


## ----featuretypeCounts_length, eval=FALSE-----------------
## fc2 <- featuretypeCounts(bfl=BamFileList(outpaths, yieldSize=50000), grl=feat,
##                          singleEnd=TRUE, readlength=c(74:76,99:102), type="data.frame")
## p2 <- plotfeaturetypeCounts(x=fc2, graphicsfile="results/featureCounts2.png",
##                             graphicsformat="png", scales="fixed", anyreadlength=FALSE,
##                             scale_length_val=NULL)


## ----pred_orf, eval=FALSE---------------------------------
## library(systemPipeRdata); library(GenomicFeatures); library(rtracklayer)
## txdb <- makeTxDbFromGFF(file="data/tair10.gff", format="gff3", organism="Arabidopsis")
## futr <- fiveUTRsByTranscript(txdb, use.names=TRUE)
## dna <- extractTranscriptSeqs(FaFile("data/tair10.fasta"), futr)
## uorf <- predORF(dna, n="all", mode="orf", longest_disjoint=TRUE, strand="sense")


## ----scale_ranges, eval=FALSE-----------------------------
## grl_scaled <- scaleRanges(subject=futr, query=uorf, type="uORF", verbose=TRUE)
## export.gff3(unlist(grl_scaled), "results/uorf.gff")


## ----translate1, eval=FALSE-------------------------------
## translate(unlist(getSeq(FaFile("data/tair10.fasta"), grl_scaled[[7]])))


## ----add_features, eval=FALSE-----------------------------
## feat <- genFeatures(txdb, featuretype="all", reduce_ranges=FALSE)
## feat <- c(feat, GRangesList("uORF"=unlist(grl_scaled)))


## ----pred_sorf, eval=FALSE--------------------------------
## feat <- genFeatures(txdb, featuretype="intergenic", reduce_ranges=TRUE)
## intergenic <- feat$intergenic
## strand(intergenic) <- "+"
## dna <- getSeq(FaFile("data/tair10.fasta"), intergenic)
## names(dna) <- mcols(intergenic)$feature_by
## sorf <- predORF(dna, n="all", mode="orf", longest_disjoint=TRUE, strand="both")
## sorf <- sorf[width(sorf) > 60] # Remove sORFs below length cutoff, here 60bp
## intergenic <- split(intergenic, mcols(intergenic)$feature_by)
## grl_scaled_intergenic <- scaleRanges(subject=intergenic, query=sorf, type="sORF", verbose=TRUE)
## export.gff3(unlist(grl_scaled_intergenic), "sorf.gff")
## translate(getSeq(FaFile("data/tair10.fasta"), unlist(grl_scaled_intergenic)))


## ----coverage_binned1, eval=FALSE-------------------------
## grl <- cdsBy(txdb, "tx", use.names=TRUE)
## fcov <- featureCoverage(bfl=BamFileList(outpaths[1:2]), grl=grl[1:4],
##                         resizereads=NULL, readlengthrange=NULL, Nbins=20, method=mean,
##                         fixedmatrix=FALSE, resizefeatures=TRUE, upstream=20,
##                         downstream=20, outfile="results/featureCoverage.xls",
##                         overwrite=TRUE)


## ----coverage_binned2, eval=FALSE-------------------------
## fcov <- featureCoverage(bfl=BamFileList(outpaths[1:4]), grl=grl[1:12], resizereads=NULL, readlengthrange=NULL, Nbins=NULL, method=mean, fixedmatrix=TRUE,
##                         resizefeatures=TRUE, upstream=20, downstream=20,
##                         outfile="results/featureCoverage.xls", overwrite=TRUE)
## plotfeatureCoverage(covMA=fcov, method=mean, scales="fixed", extendylim=2,
##                     scale_count_val=10^6)


## ----coverage_binned3, eval=FALSE-------------------------
## library(ggplot2); library(grid)
## fcov <- featureCoverage(bfl=BamFileList(outpaths[1:4]), grl=grl[1:4],
##                         resizereads=NULL, readlengthrange=NULL, Nbins=20, method=mean,
##                         fixedmatrix=TRUE, resizefeatures=TRUE, upstream=20,
##                         downstream=20,outfile="results/featureCoverage.xls",
##                         overwrite=TRUE)
## png("./results/featurePlot.png", height=12, width=24, units="in", res=72)
## plotfeatureCoverage(covMA=fcov, method=mean, scales="fixed", extendylim=2, scale_count_val=10^6)
## dev.off()


## ----coverage_nuc_level, eval=FALSE-----------------------
## fcov <- featureCoverage(bfl=BamFileList(outpaths[1:2]), grl=grl[1],
##                         resizereads=NULL, readlengthrange=NULL, Nbins=NULL, method=mean,
##                         fixedmatrix=FALSE, resizefeatures=TRUE, upstream=20,
##                         downstream=20, outfile=NULL)


## ----read_counting, eval=FALSE----------------------------
## library("GenomicFeatures"); library(BiocParallel)
## txdb <- loadDb("./data/tair10.sqlite")
## eByg <- exonsBy(txdb, by=c("gene"))
## bfl <- BamFileList(outpaths, yieldSize=50000, index=character())
## multicoreParam <- MulticoreParam(workers=8); register(multicoreParam); registered()
## counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union",
##                                                ignore.strand=TRUE,
##                                                inter.feature=FALSE,
##                                                singleEnd=TRUE))
## countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
## rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
## rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
## write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
## write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")


## ----read_counting_view, eval=FALSE-----------------------
## read.delim("results/countDFeByg.xls", row.names=1, check.names=FALSE)[1:4,1:5]


## ----read_rpkm_view, eval=FALSE---------------------------
## read.delim("results/rpkmDFeByg.xls", row.names=1, check.names=FALSE)[1:4,1:4]


## ----sample_tree, eval=FALSE------------------------------
## library(DESeq2, quietly=TRUE); library(ape,  warn.conflicts=FALSE)
## countDF <- as.matrix(read.table("./results/countDFeByg.xls"))
## colData <- data.frame(row.names=targets.as.df(targets(args))$SampleName, condition=targets.as.df(targets(args))$Factor)
## dds <- DESeq2::DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
## d <- cor(assay(DESeq2::rlog(dds)), method="spearman")
## hc <- hclust(dist(1-d))
## png("results/sample_tree.pdf")
## ape::plot.phylo(ape::as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE,
##            no.margin=TRUE)
## dev.off()


## ----deg_edger, eval=FALSE--------------------------------
## library(edgeR)
## countDF <- read.delim("results/countDFeByg.xls", row.names=1, check.names=FALSE)
## targets <- read.delim("targets.txt", comment="#")
## cmp <- readComp(file="targets.txt", format="matrix", delim="-")
## edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=FALSE, mdsplot="")


## ----add_descr, eval=FALSE--------------------------------
## library("biomaRt")
## m <- useMart("plants_mart", dataset="athaliana_eg_gene", host="plants.ensembl.org")
## desc <- getBM(attributes=c("tair_locus", "description"), mart=m)
## desc <- desc[!duplicated(desc[,1]),]
## descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
## edgeDF <- data.frame(edgeDF, Desc=descv[rownames(edgeDF)], check.names=FALSE)
## write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)


## ----filter_degs, eval=FALSE------------------------------
## edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names=1, check.names=FALSE)
## png("./results/DEGcounts.png", height=10, width=10, units="in", res=72)
## DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=20))
## dev.off()
## write.table(DEG_list$Summary, "./results/DEGcounts.xls", quote=FALSE, sep="\t", row.names=FALSE)


## ----venn_diagram, eval=FALSE-----------------------------
## vennsetup <- overLapper(DEG_list$Up[6:9], type="vennsets")
## vennsetdown <- overLapper(DEG_list$Down[6:9], type="vennsets")
## png("results/vennplot.png")
## vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))
## dev.off()


## ----get_go_annot, eval=FALSE-----------------------------
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


## ----go_enrich, eval=FALSE, warning=FALSE, message=FALSE----
## library("biomaRt")
## library(BBmisc) # Defines suppressAll()
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


## ----go_plot, eval=FALSE----------------------------------
## gos <- BatchResultslim[grep("M6-V6_up_down", BatchResultslim$CLID), ]
## gos <- BatchResultslim
## png("./results/GOslimbarplotMF.png", height=12, width=12, units="in", res=72)
## goBarplot(gos, gocat="MF")
## dev.off()
## goBarplot(gos, gocat="BP")
## goBarplot(gos, gocat="CC")


## ----diff_loading, eval=FALSE-----------------------------
## library(DESeq2)
## countDFeBygpath <- system.file("extdata", "countDFeByg.xls", package="systemPipeR")
## countDFeByg <- read.delim(countDFeBygpath, row.names=1)
## coldata <- DataFrame(assay=factor(rep(c("Ribo","mRNA"), each=4)),
##                      condition=factor(rep(as.character(targets.as.df(targets(args))$Factor[1:4]), 2)),
##                      row.names=as.character(targets.as.df(targets(args))$SampleName)[1:8])
## coldata


## ----diff_translational_eff, eval=FALSE-------------------
## dds <- DESeq2::DESeqDataSetFromMatrix(countData=as.matrix(countDFeByg[,rownames(coldata)]),
##                             colData = coldata,
##                             design = ~ assay + condition + assay:condition)
## # model.matrix(~ assay + condition + assay:condition, coldata) # Corresponding design matrix
## dds <- DESeq2::DESeq(dds, test="LRT", reduced = ~ assay + condition)
## res <- DESeq2::results(dds)
## head(res[order(res$padj),],4)
## # write.table(res, file="transleff.xls", quote=FALSE, col.names = NA, sep="\t")


## ----heatmap, eval=FALSE----------------------------------
## library(pheatmap)
## geneids <- unique(as.character(unlist(DEG_list[[1]])))
## y <- assay(DESeq2::rlog(dds))[geneids, ]
## png("heatmap1.png")
## pheatmap(y, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
## dev.off()


## ----render_report, eval=FALSE----------------------------
## rmarkdown::render("systemPipeRIBOseq.Rmd", "html_document")
## rmarkdown::render("systemPipeRIBOseq.Rmd", "pdf_document")


## ----sessionInfo------------------------------------------
sessionInfo()

