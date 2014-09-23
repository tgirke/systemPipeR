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
### code chunk number 3: systemPipeR.Rnw:67-69 (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R") # Sources the biocLite.R installation script 
## biocLite("systemPipeR") # Installs the package


###################################################
### code chunk number 4: systemPipeR.Rnw:74-77 (eval = FALSE)
###################################################
## library("systemPipeR") # Loads the package
## library(help="systemPipeR") # Lists all functions and classes 
## vignette("systemPipeR") # Opens this PDF manual from R


###################################################
### code chunk number 5: systemPipeR.Rnw:85-88
###################################################
library(systemPipeR)
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
read.delim(targetspath, comment.char = "#")


###################################################
### code chunk number 6: systemPipeR.Rnw:92-95
###################################################
library(systemPipeR)
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
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
## args <- systemArgs(sysma="tophat.param", mytargets="targetsPE.txt")


###################################################
### code chunk number 14: systemPipeR.Rnw:146-150 (eval = FALSE)
###################################################
## fqlist <- seeFastq(fastq=infile1(args), batchsize=10000, klength=8)
## pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
## seeFastqPlot(fqlist)
## dev.off()


###################################################
### code chunk number 15: systemPipeR.Rnw:161-163 (eval = FALSE)
###################################################
## moduleload(modules(args)) # Skip if module system is not available
## system("bowtie2-build ./data/tair10.fasta ./data/tair10.fasta")


###################################################
### code chunk number 16: systemPipeR.Rnw:167-168 (eval = FALSE)
###################################################
## bampaths <- runCommandline(args=args)


###################################################
### code chunk number 17: systemPipeR.Rnw:172-174 (eval = FALSE)
###################################################
## qsubargs <- getQsubargs(queue="batch", cores=cores(args), memory="mem=10gb", time="walltime=20:00:00")
## (joblist <- qsubRun(args=args, qsubargs=qsubargs, Nqsubs=18, package="systemPipeR"))


###################################################
### code chunk number 18: systemPipeR.Rnw:178-180 (eval = FALSE)
###################################################
## read_statsDF <- alignStats(args) 
## write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")


###################################################
### code chunk number 19: systemPipeR.Rnw:185-188 (eval = FALSE)
###################################################
## symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"), 
##             urlbase="http://myserver.edu/~username/", 
## 	    urlfile="IGVurl.txt")


###################################################
### code chunk number 20: systemPipeR.Rnw:193-195 (eval = FALSE)
###################################################
## args <- systemArgs(sysma="bowtieSE.param", mytargets="targets.txt")
## bampaths <- runCommandline(args=args)


###################################################
### code chunk number 21: systemPipeR.Rnw:199-201 (eval = FALSE)
###################################################
## qsubargs <- getQsubargs(queue="batch", cores=cores(args), memory="mem=10gb", time="walltime=20:00:00")
## (joblist <- qsubRun(args=args, qsubargs=qsubargs, Nqsubs=18, package="systemPipeR"))


###################################################
### code chunk number 22: systemPipeR.Rnw:206-209 (eval = FALSE)
###################################################
## library(GenomicFeatures)
## txdb <- makeTranscriptDbFromGFF(file="data/tair10.gff", format="gff", dataSource="TAIR", species="A. thaliana")
## saveDb(txdb, file="./data/tair10.sqlite")


###################################################
### code chunk number 23: systemPipeR.Rnw:213-224 (eval = FALSE)
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
### code chunk number 24: systemPipeR.Rnw:229-238 (eval = FALSE)
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
### code chunk number 25: systemPipeR.Rnw:243-249 (eval = FALSE)
###################################################
## library(ape)
## rpkmDFeByg <- read.table("./results/rpkmDFeByg.xls", check.names=FALSE)
## rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
## d <- cor(rpkmDFeByg, method="spearman")
## hc <- hclust(as.dist(1-d))
## plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)


###################################################
### code chunk number 26: systemPipeR.Rnw:259-262
###################################################
targets <- read.delim(targetspath, comment="#")
cmp <- readComp(file=targetspath, format="matrix", delim="-")
cmp[[1]]


###################################################
### code chunk number 27: systemPipeR.Rnw:265-267 (eval = FALSE)
###################################################
## countDFeByg <- read.delim("./results/countDFeByg.xls", row.names=1)
## edgeDF <- run_edgeR(countDF=countDFeByg, targets=targets, cmp=cmp[[1]], independent=FALSE, mdsplot="")


###################################################
### code chunk number 28: systemPipeR.Rnw:270-271 (eval = FALSE)
###################################################
## DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=10))


###################################################
### code chunk number 29: systemPipeR.Rnw:279-281 (eval = FALSE)
###################################################
## names(DEG_list)
## DEG_list$Summary


###################################################
### code chunk number 30: systemPipeR.Rnw:287-298 (eval = FALSE)
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
### code chunk number 31: systemPipeR.Rnw:303-314 (eval = FALSE)
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
### code chunk number 32: systemPipeR.Rnw:319-324 (eval = FALSE)
###################################################
## gos <- BatchResultslim[grep("M6-V6_up_down", BatchResultslim$CLID), ]
## gos <- BatchResultslim
## pdf("GOslimbarplotMF.pdf", height=8, width=10); goBarplot(gos, gocat="MF"); dev.off()
## goBarplot(gos, gocat="BP")
## goBarplot(gos, gocat="CC")


###################################################
### code chunk number 33: systemPipeR.Rnw:337-343 (eval = FALSE)
###################################################
## library(pheatmap)
## geneids <- unique(as.character(unlist(DEG_list[[1]])))
## y <- rpkmDFeByg[geneids, ]
## pdf("heatmap1.pdf")
## pheatmap(y, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
## dev.off()


###################################################
### code chunk number 34: sessionInfo
###################################################
toLatex(sessionInfo())


