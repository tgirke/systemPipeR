### R code from vignette source 'systemPipeRNAseq.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex(use.unsrturl=FALSE)


###################################################
### code chunk number 2: systemPipeRNAseq.Rnw:42-44
###################################################
options(width=95)
unlink("test.db")


###################################################
### code chunk number 3: systemPipeRNAseq.Rnw:71-72
###################################################
library(systemPipeR)


###################################################
### code chunk number 4: systemPipeRNAseq.Rnw:76-77 (eval = FALSE)
###################################################
## source("systemPipeRNAseq_Fct.R")


###################################################
### code chunk number 5: systemPipeRNAseq.Rnw:82-85
###################################################
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")[,1:4]
targets


###################################################
### code chunk number 6: systemPipeRNAseq.Rnw:97-102 (eval = FALSE)
###################################################
## args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
## fqlist <- seeFastq(fastq=infile1(args), batchsize=100000, klength=8)
## pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
## seeFastqPlot(fqlist)
## dev.off()


###################################################
### code chunk number 7: systemPipeRNAseq.Rnw:114-116 (eval = FALSE)
###################################################
## args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
## sysargs(args)[1] # Command-line parameters for first FASTQ file


###################################################
### code chunk number 8: systemPipeRNAseq.Rnw:119-124 (eval = FALSE)
###################################################
## moduleload(modules(args))
## system("bowtie2-build ./data/tair10.fasta ./data/tair10.fasta")
## resources <- list(walltime="20:00:00", nodes=paste0("1:ppn=", cores(args)), memory="10gb")
## reg <- clusterRun(args, conffile=".BatchJobs.R", template="torque.tmpl", Njobs=18, runid="01", 
##                   resourceList=resources)


###################################################
### code chunk number 9: systemPipeRNAseq.Rnw:127-128 (eval = FALSE)
###################################################
## file.exists(outpaths(args))


###################################################
### code chunk number 10: systemPipeRNAseq.Rnw:133-135 (eval = FALSE)
###################################################
## read_statsDF <- alignStats(args=args) 
## write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")


###################################################
### code chunk number 11: systemPipeRNAseq.Rnw:137-138
###################################################
read.table(system.file("extdata", "alignStats.xls", package="systemPipeR"), header=TRUE)[1:4,]


###################################################
### code chunk number 12: systemPipeRNAseq.Rnw:143-146 (eval = FALSE)
###################################################
## symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"), 
##             urlbase="http://biocluster.ucr.edu/~tgirke/", 
## 	    urlfile="./results/IGVurl.txt")


###################################################
### code chunk number 13: systemPipeRNAseq.Rnw:152-166 (eval = FALSE)
###################################################
## library("GenomicFeatures"); library(BiocParallel)
## txdb <- loadDb("./data/tair10.sqlite")
## eByg <- exonsBy(txdb, by=c("gene"))
## bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
## multicoreParam <- MulticoreParam(workers=8); register(multicoreParam); registered()
## counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", 
##                                                ignore.strand=TRUE, 
##                                                inter.feature=FALSE, 
##                                                singleEnd=TRUE)) 
## countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
## rownames(countDFeByg) <- names(rowData(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
## rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
## write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
## write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")


###################################################
### code chunk number 14: systemPipeRNAseq.Rnw:169-170 (eval = FALSE)
###################################################
## read.delim("results/countDFeByg.xls", row.names=1, check.names=FALSE)[1:4,1:5]


###################################################
### code chunk number 15: systemPipeRNAseq.Rnw:173-174 (eval = FALSE)
###################################################
## read.delim("results/rpkmDFeByg.xls", row.names=1, check.names=FALSE)[1:4,1:4]


###################################################
### code chunk number 16: systemPipeRNAseq.Rnw:179-187 (eval = FALSE)
###################################################
## library(ape)
## rpkmDFeByg <- read.delim("./results/rpkmDFeByg.xls", row.names=1, check.names=FALSE)[,-19]
## rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
## d <- cor(rpkmDFeByg, method="spearman")
## hc <- hclust(as.dist(1-d))
## pdf("results/sample_tree.pdf")
## plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
## dev.off()


###################################################
### code chunk number 17: systemPipeRNAseq.Rnw:198-203 (eval = FALSE)
###################################################
## library(edgeR)
## countDF <- read.delim("countDFeByg.xls", row.names=1, check.names=FALSE) 
## targets <- read.delim("targets.txt", comment="#")
## cmp <- readComp(file="targets.txt", format="matrix", delim="-")
## edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=FALSE, mdsplot="")


###################################################
### code chunk number 18: systemPipeRNAseq.Rnw:207-212 (eval = FALSE)
###################################################
## desc <- read.delim("data/desc.xls") 
## desc <- desc[!duplicated(desc[,1]),]
## descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
## edgeDF <- data.frame(edgeDF, Desc=descv[rownames(edgeDF)], check.names=FALSE)
## write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)


###################################################
### code chunk number 19: systemPipeRNAseq.Rnw:216-221 (eval = FALSE)
###################################################
## edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names=1, check.names=FALSE) 
## pdf("results/DEGcounts.pdf")
## DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=1))
## dev.off()
## write.table(DEG_list$Summary, "./results/DEGcounts.xls", quote=FALSE, sep="\t", row.names=FALSE)


###################################################
### code chunk number 20: systemPipeRNAseq.Rnw:231-236 (eval = FALSE)
###################################################
## vennsetup <- overLapper(DEG_list$Up[6:9], type="vennsets")
## vennsetdown <- overLapper(DEG_list$Down[6:9], type="vennsets")
## pdf("results/vennplot.pdf")
## vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))
## dev.off()


###################################################
### code chunk number 21: systemPipeRNAseq.Rnw:248-261 (eval = FALSE)
###################################################
## library("biomaRt")
## listMarts() # To choose BioMart database
## m <- useMart("ENSEMBL_MART_PLANT"); listDatasets(m) 
## m <- useMart("ENSEMBL_MART_PLANT", dataset="athaliana_eg_gene")
## listAttributes(m) # Choose data types you want to download
## go <- getBM(attributes=c("go_accession", "tair_locus", "go_namespace_1003"), mart=m)
## go <- go[go[,3]!="",]; go[,3] <- as.character(go[,3])
## go[go[,3]=="molecular_function", 3] <- "F"; go[go[,3]=="biological_process", 3] <- "P"; go[go[,3]=="cellular_component", 3] <- "C"
## go[1:4,]
## dir.create("./data/GO")
## write.table(go, "data/GO/GOannotationsBiomart_mod.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
## catdb <- makeCATdb(myfile="data/GO/GOannotationsBiomart_mod.txt", lib=NULL, org="", colno=c(1,2,3), idconv=NULL)
## save(catdb, file="data/GO/catdb.RData") 


###################################################
### code chunk number 22: systemPipeRNAseq.Rnw:266-277 (eval = FALSE)
###################################################
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


###################################################
### code chunk number 23: systemPipeRNAseq.Rnw:282-287 (eval = FALSE)
###################################################
## gos <- BatchResultslim[grep("M6-V6_up_down", BatchResultslim$CLID), ]
## gos <- BatchResultslim
## pdf("GOslimbarplotMF.pdf", height=8, width=10); goBarplot(gos, gocat="MF"); dev.off()
## goBarplot(gos, gocat="BP")
## goBarplot(gos, gocat="CC")


###################################################
### code chunk number 24: systemPipeRNAseq.Rnw:300-306 (eval = FALSE)
###################################################
## library(pheatmap)
## geneids <- unique(as.character(unlist(DEG_list[[1]])))
## y <- rpkmDFeByg[geneids, ]
## pdf("heatmap1.pdf")
## pheatmap(y, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
## dev.off()


###################################################
### code chunk number 25: sessionInfo
###################################################
toLatex(sessionInfo())


