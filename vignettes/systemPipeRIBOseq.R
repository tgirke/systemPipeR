## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"---------------------------------
BiocStyle::latex(use.unsrturl=FALSE)

## ----setup, include=FALSE, cache=FALSE-------------------------------------------------------
library(knitr)
# set global chunk options for knitr
opts_chunk$set(comment=NA, warning=FALSE, message=FALSE, fig.path='figure/systemPipeR-')
options(formatR.arrow=TRUE, width=95)
unlink("test.db")

## ----eval=TRUE-------------------------------------------------------------------------------
library(systemPipeR)

## ----eval=FALSE------------------------------------------------------------------------------
#  library(systemPipeRdata)
#  genWorkenvir(workflow="chipseq")
#  setwd("chipseq")

## ----eval=FALSE------------------------------------------------------------------------------
#  source("systemPipeRIBOseq_Fct.R")

## ----eval=TRUE-------------------------------------------------------------------------------
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")[,1:4]
targets

## ----eval=FALSE------------------------------------------------------------------------------
#  args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
#  fqlist <- seeFastq(fastq=infile1(args), batchsize=100000, klength=8)
#  pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
#  seeFastqPlot(fqlist)
#  dev.off()

## ----eval=FALSE------------------------------------------------------------------------------
#  args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
#  sysargs(args)[1] # Command-line parameters for first FASTQ file

## ----eval=FALSE------------------------------------------------------------------------------
#  moduleload(modules(args))
#  system("bowtie2-build ./data/tair10.fasta ./data/tair10.fasta")
#  resources <- list(walltime="20:00:00", nodes=paste0("1:ppn=", cores(args)), memory="10gb")
#  reg <- clusterRun(args, conffile=".BatchJobs.R", template="torque.tmpl", Njobs=18, runid="01",
#                    resourceList=resources)
#  waitForJobs(reg)

## ----eval=FALSE------------------------------------------------------------------------------
#  file.exists(outpaths(args))

## ----eval=FALSE------------------------------------------------------------------------------
#  read_statsDF <- alignStats(args=args)
#  write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")

## ----eval=TRUE-------------------------------------------------------------------------------
read.table(system.file("extdata", "alignStats.xls", package="systemPipeR"), header=TRUE)[1:4,]

## ----eval=FALSE------------------------------------------------------------------------------
#  symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"),
#              urlbase="http://biocluster.ucr.edu/~tgirke/",
#  	    urlfile="./results/IGVurl.txt")

## ----eval=FALSE------------------------------------------------------------------------------
#  library(GenomicFeatures)
#  file <- system.file("extdata/annotation", "tair10.gff", package="systemPipeRdata")
#  txdb <- makeTxDbFromGFF(file=file, format="gff3", organism="Arabidopsis")
#  feat <- genFeatures(txdb, featuretype="all", reduce_ranges=TRUE, upstream=1000, downstream=0, verbose=TRUE)

## ----eval=FALSE------------------------------------------------------------------------------
#  fc <- featuretypeCounts(bfl=BamFileList(outpaths(args), yieldSize=50000), grl=feat, singleEnd=TRUE, readlength=NULL, type="data.frame")
#  p <- plotfeaturetypeCounts(x=featureCounts2, graphicsfile="featureCounts.pdf", graphicsformat="pdf", scales="fixed", anyreadlength=TRUE, scale_length_val=NULL)

## ----eval=FALSE------------------------------------------------------------------------------
#  fc2 <- featuretypeCounts(bfl=BamFileList(outpaths(args), yieldSize=50000), grl=feat, singleEnd=TRUE, readlength=c(74:76,99:102), type="data.frame")
#  p2 <- plotfeaturetypeCounts(x=fc2, graphicsfile="featureCounts2.pdf", graphicsformat="pdf", scales="fixed", anyreadlength=FALSE, scale_length_val=NULL)

## ----eval=FALSE------------------------------------------------------------------------------
#  library("GenomicFeatures"); library(BiocParallel)
#  txdb <- loadDb("./data/tair10.sqlite")
#  eByg <- exonsBy(txdb, by=c("gene"))
#  bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
#  multicoreParam <- MulticoreParam(workers=8); register(multicoreParam); registered()
#  counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union",
#                                                 ignore.strand=TRUE,
#                                                 inter.feature=FALSE,
#                                                 singleEnd=TRUE))
#  countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
#  rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
#  rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
#  write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
#  write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")

## ----eval=FALSE------------------------------------------------------------------------------
#  read.delim("results/countDFeByg.xls", row.names=1, check.names=FALSE)[1:4,1:5]

## ----eval=FALSE------------------------------------------------------------------------------
#  read.delim("results/rpkmDFeByg.xls", row.names=1, check.names=FALSE)[1:4,1:4]

## ----eval=FALSE------------------------------------------------------------------------------
#  library(DESeq2, quietly=TRUE); library(ape,  warn.conflicts=FALSE)
#  countDF <- as.matrix(read.table("./results/countDFeByg.xls"))
#  colData <- data.frame(row.names=targetsin(args)$SampleName, condition=targetsin(args)$Factor)
#  dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
#  d <- cor(assay(rlog(dds)), method="spearman")
#  hc <- hclust(dist(1-d))
#  pdf("results/sample_tree.pdf")
#  plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
#  dev.off()

## ----eval=FALSE------------------------------------------------------------------------------
#  library(edgeR)
#  countDF <- read.delim("countDFeByg.xls", row.names=1, check.names=FALSE)
#  targets <- read.delim("targets.txt", comment="#")
#  cmp <- readComp(file="targets.txt", format="matrix", delim="-")
#  edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=FALSE, mdsplot="")

## ----eval=FALSE------------------------------------------------------------------------------
#  desc <- read.delim("data/desc.xls")
#  desc <- desc[!duplicated(desc[,1]),]
#  descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
#  edgeDF <- data.frame(edgeDF, Desc=descv[rownames(edgeDF)], check.names=FALSE)
#  write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)

## ----eval=FALSE------------------------------------------------------------------------------
#  edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names=1, check.names=FALSE)
#  pdf("results/DEGcounts.pdf")
#  DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=1))
#  dev.off()
#  write.table(DEG_list$Summary, "./results/DEGcounts.xls", quote=FALSE, sep="\t", row.names=FALSE)

## ----eval=FALSE------------------------------------------------------------------------------
#  vennsetup <- overLapper(DEG_list$Up[6:9], type="vennsets")
#  vennsetdown <- overLapper(DEG_list$Down[6:9], type="vennsets")
#  pdf("results/vennplot.pdf")
#  vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))
#  dev.off()

## ----eval=FALSE------------------------------------------------------------------------------
#  library("biomaRt")
#  listMarts() # To choose BioMart database
#  m <- useMart("ENSEMBL_MART_PLANT"); listDatasets(m)
#  m <- useMart("ENSEMBL_MART_PLANT", dataset="athaliana_eg_gene")
#  listAttributes(m) # Choose data types you want to download
#  go <- getBM(attributes=c("go_accession", "tair_locus", "go_namespace_1003"), mart=m)
#  go <- go[go[,3]!="",]; go[,3] <- as.character(go[,3])
#  go[go[,3]=="molecular_function", 3] <- "F"; go[go[,3]=="biological_process", 3] <- "P"; go[go[,3]=="cellular_component", 3] <- "C"
#  go[1:4,]
#  dir.create("./data/GO")
#  write.table(go, "data/GO/GOannotationsBiomart_mod.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
#  catdb <- makeCATdb(myfile="data/GO/GOannotationsBiomart_mod.txt", lib=NULL, org="", colno=c(1,2,3), idconv=NULL)
#  save(catdb, file="data/GO/catdb.RData")

## ----eval=FALSE------------------------------------------------------------------------------
#  load("data/GO/catdb.RData")
#  DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=50), plot=FALSE)
#  up_down <- DEG_list$UporDown; names(up_down) <- paste(names(up_down), "_up_down", sep="")
#  up <- DEG_list$Up; names(up) <- paste(names(up), "_up", sep="")
#  down <- DEG_list$Down; names(down) <- paste(names(down), "_down", sep="")
#  DEGlist <- c(up_down, up, down)
#  DEGlist <- DEGlist[sapply(DEGlist, length) > 0]
#  BatchResult <- GOCluster_Report(catdb=catdb, setlist=DEGlist, method="all", id_type="gene", CLSZ=2, cutoff=0.9, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
#  library("biomaRt"); m <- useMart("ENSEMBL_MART_PLANT", dataset="athaliana_eg_gene")
#  goslimvec <- as.character(getBM(attributes=c("goslim_goa_accession"), mart=m)[,1])
#  BatchResultslim <- GOCluster_Report(catdb=catdb, setlist=DEGlist, method="slim", id_type="gene", myslimv=goslimvec, CLSZ=10, cutoff=0.01, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)

## ----eval=FALSE------------------------------------------------------------------------------
#  gos <- BatchResultslim[grep("M6-V6_up_down", BatchResultslim$CLID), ]
#  gos <- BatchResultslim
#  pdf("GOslimbarplotMF.pdf", height=8, width=10); goBarplot(gos, gocat="MF"); dev.off()
#  goBarplot(gos, gocat="BP")
#  goBarplot(gos, gocat="CC")

## ----eval=FALSE------------------------------------------------------------------------------
#  library(pheatmap)
#  geneids <- unique(as.character(unlist(DEG_list[[1]])))
#  y <- rpkmDFeByg[geneids, ]
#  pdf("heatmap1.pdf")
#  pheatmap(y, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
#  dev.off()

## ----sessionInfo, results='asis'-------------------------------------------------------------
toLatex(sessionInfo())

