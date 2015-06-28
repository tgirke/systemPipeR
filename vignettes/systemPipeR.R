### R code from vignette source 'systemPipeR.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex(use.unsrturl=FALSE)


###################################################
### code chunk number 2: systemPipeR.Rnw:42-44
###################################################
options(width=95)
unlink("test.db")


###################################################
### code chunk number 3: systemPipeR.Rnw:69-71 (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R") # Sources the biocLite.R installation script 
## biocLite("systemPipeR") # Installs the package


###################################################
### code chunk number 4: systemPipeR.Rnw:76-79 (eval = FALSE)
###################################################
## library("systemPipeR") # Loads the package
## library(help="systemPipeR") # Lists all functions and classes 
## vignette("systemPipeR") # Opens this PDF manual from R


###################################################
### code chunk number 5: systemPipeR.Rnw:87-90
###################################################
library(systemPipeR)
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
read.delim(targetspath, comment.char = "#")


###################################################
### code chunk number 6: systemPipeR.Rnw:94-96
###################################################
targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
read.delim(targetspath, comment.char = "#")[1:2,1:6]


###################################################
### code chunk number 7: systemPipeR.Rnw:100-101
###################################################
readComp(file=targetspath, format="vector", delim="-")


###################################################
### code chunk number 8: systemPipeR.Rnw:106-108
###################################################
parampath <- system.file("extdata", "tophat.param", package="systemPipeR")
read.delim(parampath, comment.char = "#")


###################################################
### code chunk number 9: systemPipeR.Rnw:111-113
###################################################
args <- systemArgs(sysma=parampath, mytargets=targetspath)
args


###################################################
### code chunk number 10: systemPipeR.Rnw:116-121
###################################################
names(args)
modules(args)
cores(args)
outpaths(args)[1]
sysargs(args)[1]


###################################################
### code chunk number 11: systemPipeR.Rnw:124-125
###################################################
systemArgs(sysma=parampath, mytargets=targetspath, type="json")


###################################################
### code chunk number 12: systemPipeR.Rnw:131-132 (eval = FALSE)
###################################################
## library(systemPipeR)


###################################################
### code chunk number 13: systemPipeR.Rnw:136-137 (eval = FALSE)
###################################################
## args <- systemArgs(sysma="trim.param", mytargets="targets.txt")


###################################################
### code chunk number 14: systemPipeR.Rnw:155-158 (eval = FALSE)
###################################################
## preprocessReads(args=args, Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)", 
##                 batchsize=100000, overwrite=TRUE, compress=TRUE)
## writeTargetsout(x=args, file="targets_trim.txt")


###################################################
### code chunk number 15: systemPipeR.Rnw:164-171 (eval = FALSE)
###################################################
## filterFct <- function(fq) {
## 	filter1 <- nFilter(threshold=1) # Keeps only reads without Ns
## 	filter2 <- polynFilter(threshold=20, nuc=c("A","T","G","C")) # Removes low complexity reads 
## 	filter <- compose(filter1, filter2)
## 	fq[filter(fq)]
## }
## preprocessReads(args=args, Fct="filterFct(fq)", batchsize=100000)


###################################################
### code chunk number 16: systemPipeR.Rnw:180-184 (eval = FALSE)
###################################################
## fqlist <- seeFastq(fastq=infile1(args), batchsize=10000, klength=8)
## pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
## seeFastqPlot(fqlist)
## dev.off()


###################################################
### code chunk number 17: systemPipeR.Rnw:194-198 (eval = FALSE)
###################################################
## args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
## f <- function(x) seeFastq(fastq=infile1(args)[x], batchsize=100000, klength=8)
## fqlist <- bplapply(seq(along=args), f, BPPARAM = MulticoreParam(workers=8))
## seeFastqPlot(unlist(fqlist, recursive=FALSE))


###################################################
### code chunk number 18: systemPipeR.Rnw:202-213 (eval = FALSE)
###################################################
## library(BiocParallel); library(BatchJobs)
## f <- function(x) {
##     library(systemPipeR)
##     args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
##     seeFastq(fastq=infile1(args)[x], batchsize=100000, klength=8)
## }
## funs <- makeClusterFunctionsTorque("torque.tmpl")
## param <- BatchJobsParam(length(args), resources=list(walltime="20:00:00", nodes="1:ppn=1", memory="6gb"), cluster.functions=funs)
## register(param)
## fqlist <- bplapply(seq(along=args), f)
## seeFastqPlot(unlist(fqlist, recursive=FALSE))


###################################################
### code chunk number 19: systemPipeR.Rnw:218-221 (eval = FALSE)
###################################################
## args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
## moduleload(modules(args)) # Skip if module system is not available
## system("bowtie2-build ./data/tair10.fasta ./data/tair10.fasta")


###################################################
### code chunk number 20: systemPipeR.Rnw:225-226 (eval = FALSE)
###################################################
## bampaths <- runCommandline(args=args)


###################################################
### code chunk number 21: systemPipeR.Rnw:230-236 (eval = FALSE)
###################################################
## file.copy(system.file("extdata", ".BatchJobs.R", package="systemPipeR"), ".")
## file.copy(system.file("extdata", "torque.tmpl", package="systemPipeR"), ".")
## resources <- list(walltime="20:00:00", nodes=paste0("1:ppn=", cores(args)), memory="10gb")
## reg <- clusterRun(args, conffile=".BatchJobs.R", template="torque.tmpl", Njobs=18, runid="01", 
##                   resourceList=resources)
## waitForJobs(reg)


###################################################
### code chunk number 22: systemPipeR.Rnw:240-243 (eval = FALSE)
###################################################
## showStatus(reg)
## file.exists(outpaths(args))
## sapply(1:length(args), function(x) loadResult(reg, x)) # Works after job completion


###################################################
### code chunk number 23: systemPipeR.Rnw:248-250 (eval = FALSE)
###################################################
## read_statsDF <- alignStats(args) 
## write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")


###################################################
### code chunk number 24: systemPipeR.Rnw:254-255
###################################################
read.table(system.file("extdata", "alignStats.xls", package="systemPipeR"), header=TRUE)[1:4,]


###################################################
### code chunk number 25: systemPipeR.Rnw:259-262 (eval = FALSE)
###################################################
## f <- function(x) alignStats(args[x])
## read_statsList <- bplapply(seq(along=args), f, BPPARAM = MulticoreParam(workers=8))
## read_statsDF <- do.call("rbind", read_statsList)


###################################################
### code chunk number 26: systemPipeR.Rnw:266-277 (eval = FALSE)
###################################################
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


###################################################
### code chunk number 27: systemPipeR.Rnw:282-285 (eval = FALSE)
###################################################
## symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"), 
##             urlbase="http://myserver.edu/~username/", 
## 	    urlfile="IGVurl.txt")


###################################################
### code chunk number 28: systemPipeR.Rnw:292-295 (eval = FALSE)
###################################################
## args <- systemArgs(sysma="bowtieSE.param", mytargets="targets.txt")
## moduleload(modules(args)) # Skip if module system is not available
## bampaths <- runCommandline(args=args)


###################################################
### code chunk number 29: systemPipeR.Rnw:299-301 (eval = FALSE)
###################################################
## qsubargs <- getQsubargs(queue="batch", cores=cores(args), memory="mem=10gb", time="walltime=20:00:00")
## (joblist <- qsubRun(args=args, qsubargs=qsubargs, Nqsubs=18, package="systemPipeR"))


###################################################
### code chunk number 30: systemPipeR.Rnw:306-310 (eval = FALSE)
###################################################
## args <- systemArgs(sysma="bwa.param", mytargets="targets.txt")
## moduleload(modules(args)) # Skip if module system is not available
## system("bwa index -a bwtsw ./data/tair10.fasta") # Indexes reference genome
## bampaths <- runCommandline(args=args)


###################################################
### code chunk number 31: systemPipeR.Rnw:316-322 (eval = FALSE)
###################################################
## library(Rsubread)
## args <- systemArgs(sysma="rsubread.param", mytargets="targets.txt")
## buildindex(basename=reference(args), reference=reference(args)) # Build indexed reference genome
## align(index=reference(args), readfile1=infile1(args), input_format="FASTQ", 
##       output_file=outfile1(args), output_format="SAM", nthreads=8, indels=1, TH1=2)
## for(i in seq(along=outfile1(args))) asBam(file=outfile1(args)[i], destination=gsub(".sam", "", outfile1(args)[i]), overwrite=TRUE, indexDestination=TRUE)


###################################################
### code chunk number 32: systemPipeR.Rnw:327-330 (eval = FALSE)
###################################################
## library(GenomicFeatures)
## txdb <- makeTranscriptDbFromGFF(file="data/tair10.gff", format="gff", dataSource="TAIR", species="A. thaliana")
## saveDb(txdb, file="./data/tair10.sqlite")


###################################################
### code chunk number 33: systemPipeR.Rnw:334-345 (eval = FALSE)
###################################################
## library(BiocParallel)
## txdb <- loadDb("./data/tair10.sqlite")
## eByg <- exonsBy(txdb, by="gene")
## bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
## multicoreParam <- MulticoreParam(workers=4); register(multicoreParam); registered()
## counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", ignore.strand=TRUE, inter.feature=TRUE, singleEnd=TRUE)) # Note: for strand-specific RNA-Seq set 'ignore.strand=FALSE' and for PE data set 'singleEnd=FALSE'
## countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
## rownames(countDFeByg) <- names(rowData(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
## rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
## write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
## write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")


###################################################
### code chunk number 34: systemPipeR.Rnw:349-364 (eval = FALSE)
###################################################
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


###################################################
### code chunk number 35: systemPipeR.Rnw:369-378 (eval = FALSE)
###################################################
## system("wget ftp://mirbase.org/pub/mirbase/19/genomes/My_species.gff3 -P ./data/")
## gff <- import.gff("./data/My_species.gff3", asRangedData=FALSE)
## gff <- split(gff, elementMetadata(gff)$ID)
## bams <- names(bampaths); names(bams) <- targets$SampleName
## bfl <- BamFileList(bams, yieldSize=50000, index=character())
## countDFmiR <- summarizeOverlaps(gff, bfl, mode="Union", ignore.strand=FALSE, inter.feature=FALSE) # Note: inter.feature=FALSE important since pre and mature miRNA ranges overlap
## rpkmDFmiR <- apply(countDFmiR, 2, function(x) returnRPKM(counts=x, gffsub=gff))
## write.table(assays(countDFmiR)$counts, "results/countDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")
## write.table(rpkmDFmiR, "results/rpkmDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")


###################################################
### code chunk number 36: sample_tree
###################################################
library(ape)
rpkmDFeBygpath <- system.file("extdata", "rpkmDFeByg.xls", package="systemPipeR")
rpkmDFeByg <- read.table(rpkmDFeBygpath, check.names=FALSE)
rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
d <- cor(rpkmDFeByg, method="spearman")
hc <- hclust(as.dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)


###################################################
### code chunk number 37: systemPipeR.Rnw:408-411
###################################################
targets <- read.delim(targetspath, comment="#")
cmp <- readComp(file=targetspath, format="matrix", delim="-")
cmp[[1]]


###################################################
### code chunk number 38: systemPipeR.Rnw:412-415
###################################################
countDFeBygpath <- system.file("extdata", "countDFeByg.xls", package="systemPipeR")
countDFeByg <- read.delim(countDFeBygpath, row.names=1)
edgeDF <- run_edgeR(countDF=countDFeByg, targets=targets, cmp=cmp[[1]], independent=FALSE, mdsplot="")


###################################################
### code chunk number 39: DEGcounts
###################################################
DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=10))


###################################################
### code chunk number 40: systemPipeR.Rnw:427-429
###################################################
names(DEG_list)
DEG_list$Summary[1:4,]


###################################################
### code chunk number 41: systemPipeR.Rnw:440-441
###################################################
degseqDF <- run_DESeq2(countDF=countDFeByg, targets=targets, cmp=cmp[[1]], independent=FALSE)


###################################################
### code chunk number 42: DEGcounts2
###################################################
DEG_list2 <- filterDEGs(degDF=degseqDF, filter=c(Fold=2, FDR=10))


###################################################
### code chunk number 43: vennplot
###################################################
vennsetup <- overLapper(DEG_list$Up[6:9], type="vennsets")
vennsetdown <- overLapper(DEG_list$Down[6:9], type="vennsets")
vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))


###################################################
### code chunk number 44: systemPipeR.Rnw:471-482 (eval = FALSE)
###################################################
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


###################################################
### code chunk number 45: systemPipeR.Rnw:487-498 (eval = FALSE)
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
### code chunk number 46: systemPipeR.Rnw:503-508 (eval = FALSE)
###################################################
## gos <- BatchResultslim[grep("M6-V6_up_down", BatchResultslim$CLID), ]
## gos <- BatchResultslim
## pdf("GOslimbarplotMF.pdf", height=8, width=10); goBarplot(gos, gocat="MF"); dev.off()
## goBarplot(gos, gocat="BP")
## goBarplot(gos, gocat="CC")


###################################################
### code chunk number 47: systemPipeR.Rnw:521-527 (eval = FALSE)
###################################################
## library(pheatmap)
## geneids <- unique(as.character(unlist(DEG_list[[1]])))
## y <- rpkmDFeByg[geneids, ]
## pdf("heatmap1.pdf")
## pheatmap(y, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
## dev.off()


###################################################
### code chunk number 48: sessionInfo
###################################################
toLatex(sessionInfo())


