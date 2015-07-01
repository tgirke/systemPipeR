## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"---------------------------------
BiocStyle::latex(use.unsrturl=FALSE)

## ----setup, include=FALSE, cache=FALSE-------------------------------------------------------
library(knitr)
# set global chunk options for knitr
opts_chunk$set(comment=NA, warning=FALSE, message=FALSE, fig.path='figure/systemPipeR-')
options(formatR.arrow=TRUE, width=95)
unlink("test.db")

## ----eval=FALSE------------------------------------------------------------------------------
#  source("http://bioconductor.org/biocLite.R") # Sources the biocLite.R installation script
#  biocLite("systemPipeR") # Installs the package

## ----eval=FALSE------------------------------------------------------------------------------
#  library("systemPipeR") # Loads the package
#  library(help="systemPipeR") # Lists all functions and classes
#  vignette("systemPipeR") # Opens this PDF manual from R

## ----eval=TRUE-------------------------------------------------------------------------------
library(systemPipeR)
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
read.delim(targetspath, comment.char = "#")

## ----eval=TRUE-------------------------------------------------------------------------------
targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
read.delim(targetspath, comment.char = "#")[1:2,1:6]

## ----eval=TRUE-------------------------------------------------------------------------------
readComp(file=targetspath, format="vector", delim="-")

## ----eval=TRUE-------------------------------------------------------------------------------
parampath <- system.file("extdata", "tophat.param", package="systemPipeR")
read.delim(parampath, comment.char = "#")

## ----eval=TRUE-------------------------------------------------------------------------------
args <- systemArgs(sysma=parampath, mytargets=targetspath)
args

## ----eval=TRUE-------------------------------------------------------------------------------
names(args)
modules(args)
cores(args)
outpaths(args)[1]
sysargs(args)[1]

## ----eval=TRUE-------------------------------------------------------------------------------
systemArgs(sysma=parampath, mytargets=targetspath, type="json")

## ----eval=FALSE------------------------------------------------------------------------------
#  library(systemPipeR)

## ----eval=FALSE------------------------------------------------------------------------------
#  args <- systemArgs(sysma="trim.param", mytargets="targets.txt")

## ----eval=FALSE------------------------------------------------------------------------------
#  preprocessReads(args=args, Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)",
#                  batchsize=100000, overwrite=TRUE, compress=TRUE)
#  writeTargetsout(x=args, file="targets_trim.txt")

## ----eval=FALSE------------------------------------------------------------------------------
#  filterFct <- function(fq) {
#  	filter1 <- nFilter(threshold=1) # Keeps only reads without Ns
#  	filter2 <- polynFilter(threshold=20, nuc=c("A","T","G","C")) # Removes low complexity reads
#  	filter <- compose(filter1, filter2)
#  	fq[filter(fq)]
#  }
#  preprocessReads(args=args, Fct="filterFct(fq)", batchsize=100000)

## ----eval=FALSE------------------------------------------------------------------------------
#  fqlist <- seeFastq(fastq=infile1(args), batchsize=10000, klength=8)
#  pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
#  seeFastqPlot(fqlist)
#  dev.off()

## ----eval=FALSE------------------------------------------------------------------------------
#  args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
#  f <- function(x) seeFastq(fastq=infile1(args)[x], batchsize=100000, klength=8)
#  fqlist <- bplapply(seq(along=args), f, BPPARAM = MulticoreParam(workers=8))
#  seeFastqPlot(unlist(fqlist, recursive=FALSE))

## ----eval=FALSE------------------------------------------------------------------------------
#  library(BiocParallel); library(BatchJobs)
#  f <- function(x) {
#      library(systemPipeR)
#      args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
#      seeFastq(fastq=infile1(args)[x], batchsize=100000, klength=8)
#  }
#  funs <- makeClusterFunctionsTorque("torque.tmpl")
#  param <- BatchJobsParam(length(args), resources=list(walltime="20:00:00", nodes="1:ppn=1", memory="6gb"), cluster.functions=funs)
#  register(param)
#  fqlist <- bplapply(seq(along=args), f)
#  seeFastqPlot(unlist(fqlist, recursive=FALSE))

## ----eval=FALSE------------------------------------------------------------------------------
#  args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
#  moduleload(modules(args)) # Skip if module system is not available
#  system("bowtie2-build ./data/tair10.fasta ./data/tair10.fasta")

## ----eval=FALSE------------------------------------------------------------------------------
#  bampaths <- runCommandline(args=args)

## ----eval=FALSE------------------------------------------------------------------------------
#  file.copy(system.file("extdata", ".BatchJobs.R", package="systemPipeR"), ".")
#  file.copy(system.file("extdata", "torque.tmpl", package="systemPipeR"), ".")
#  resources <- list(walltime="20:00:00", nodes=paste0("1:ppn=", cores(args)), memory="10gb")
#  reg <- clusterRun(args, conffile=".BatchJobs.R", template="torque.tmpl", Njobs=18, runid="01",
#                    resourceList=resources)
#  waitForJobs(reg)

## ----eval=FALSE------------------------------------------------------------------------------
#  showStatus(reg)
#  file.exists(outpaths(args))
#  sapply(1:length(args), function(x) loadResult(reg, x)) # Works after job completion

## ----eval=FALSE------------------------------------------------------------------------------
#  read_statsDF <- alignStats(args)
#  write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")

## ----eval=TRUE-------------------------------------------------------------------------------
read.table(system.file("extdata", "alignStats.xls", package="systemPipeR"), header=TRUE)[1:4,]

## ----eval=FALSE------------------------------------------------------------------------------
#  f <- function(x) alignStats(args[x])
#  read_statsList <- bplapply(seq(along=args), f, BPPARAM = MulticoreParam(workers=8))
#  read_statsDF <- do.call("rbind", read_statsList)

## ----eval=FALSE------------------------------------------------------------------------------
#  library(BiocParallel); library(BatchJobs)
#  f <- function(x) {
#      library(systemPipeR)
#      args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
#      alignStats(args[x])
#  }
#  funs <- makeClusterFunctionsTorque("torque.tmpl")
#  param <- BatchJobsParam(length(args), resources=list(walltime="20:00:00", nodes="1:ppn=1", memory="6gb"), cluster.functions=funs)
#  register(param)
#  read_statsList <- bplapply(seq(along=args), f)
#  read_statsDF <- do.call("rbind", read_statsList)

## ----eval=FALSE------------------------------------------------------------------------------
#  symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"),
#              urlbase="http://myserver.edu/~username/",
#  	    urlfile="IGVurl.txt")

## ----eval=FALSE------------------------------------------------------------------------------
#  args <- systemArgs(sysma="bowtieSE.param", mytargets="targets.txt")
#  moduleload(modules(args)) # Skip if module system is not available
#  bampaths <- runCommandline(args=args)

## ----eval=FALSE------------------------------------------------------------------------------
#  qsubargs <- getQsubargs(queue="batch", cores=cores(args), memory="mem=10gb", time="walltime=20:00:00")
#  (joblist <- qsubRun(args=args, qsubargs=qsubargs, Nqsubs=18, package="systemPipeR"))

## ----eval=FALSE------------------------------------------------------------------------------
#  args <- systemArgs(sysma="bwa.param", mytargets="targets.txt")
#  moduleload(modules(args)) # Skip if module system is not available
#  system("bwa index -a bwtsw ./data/tair10.fasta") # Indexes reference genome
#  bampaths <- runCommandline(args=args)

## ----eval=FALSE------------------------------------------------------------------------------
#  library(Rsubread)
#  args <- systemArgs(sysma="rsubread.param", mytargets="targets.txt")
#  buildindex(basename=reference(args), reference=reference(args)) # Build indexed reference genome
#  align(index=reference(args), readfile1=infile1(args), input_format="FASTQ",
#        output_file=outfile1(args), output_format="SAM", nthreads=8, indels=1, TH1=2)
#  for(i in seq(along=outfile1(args))) asBam(file=outfile1(args)[i], destination=gsub(".sam", "", outfile1(args)[i]), overwrite=TRUE, indexDestination=TRUE)

## ----eval=FALSE------------------------------------------------------------------------------
#  library(GenomicFeatures)
#  txdb <- makeTranscriptDbFromGFF(file="data/tair10.gff", format="gff", dataSource="TAIR", species="A. thaliana")
#  saveDb(txdb, file="./data/tair10.sqlite")

## ----eval=FALSE------------------------------------------------------------------------------
#  library(BiocParallel)
#  txdb <- loadDb("./data/tair10.sqlite")
#  eByg <- exonsBy(txdb, by="gene")
#  bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
#  multicoreParam <- MulticoreParam(workers=4); register(multicoreParam); registered()
#  counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", ignore.strand=TRUE, inter.feature=TRUE, singleEnd=TRUE)) # Note: for strand-specific RNA-Seq set 'ignore.strand=FALSE' and for PE data set 'singleEnd=FALSE'
#  countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
#  rownames(countDFeByg) <- names(rowData(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
#  rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
#  write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
#  write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")

## ----eval=FALSE------------------------------------------------------------------------------
#  library(BiocParallel)
#  f <- function(x) {
#      library(systemPipeR); library(BiocParallel); library(GenomicFeatures)
#      txdb <- loadDb("./data/tair10.sqlite")
#      eByg <- exonsBy(txdb, by="gene")
#      args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
#      bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
#      summarizeOverlaps(eByg, bfl[x], mode="Union", ignore.strand=TRUE, inter.feature=TRUE, singleEnd=TRUE)
#  }
#  funs <- makeClusterFunctionsTorque("torque.tmpl")
#  param <- BatchJobsParam(length(args), resources=list(walltime="20:00:00", nodes="1:ppn=1", memory="6gb"), cluster.functions=funs)
#  register(param)
#  counteByg <- bplapply(seq(along=args), f)
#  countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
#  rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(outpaths(args))

## ----eval=FALSE------------------------------------------------------------------------------
#  system("wget ftp://mirbase.org/pub/mirbase/19/genomes/My_species.gff3 -P ./data/")
#  gff <- import.gff("./data/My_species.gff3", asRangedData=FALSE)
#  gff <- split(gff, elementMetadata(gff)$ID)
#  bams <- names(bampaths); names(bams) <- targets$SampleName
#  bfl <- BamFileList(bams, yieldSize=50000, index=character())
#  countDFmiR <- summarizeOverlaps(gff, bfl, mode="Union", ignore.strand=FALSE, inter.feature=FALSE) # Note: inter.feature=FALSE important since pre and mature miRNA ranges overlap
#  rpkmDFmiR <- apply(countDFmiR, 2, function(x) returnRPKM(counts=x, gffsub=gff))
#  write.table(assays(countDFmiR)$counts, "results/countDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")
#  write.table(rpkmDFmiR, "results/rpkmDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")

## ----sample_tree, eval=TRUE, include=FALSE---------------------------------------------------
library(ape)
rpkmDFeBygpath <- system.file("extdata", "rpkmDFeByg.xls", package="systemPipeR")
rpkmDFeByg <- read.table(rpkmDFeBygpath, check.names=FALSE)
rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
d <- cor(rpkmDFeByg, method="spearman")
hc <- hclust(as.dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)

## ----eval=TRUE-------------------------------------------------------------------------------
targets <- read.delim(targetspath, comment="#")
cmp <- readComp(file=targetspath, format="matrix", delim="-")
cmp[[1]]

## ----eval=TRUE-------------------------------------------------------------------------------
countDFeBygpath <- system.file("extdata", "countDFeByg.xls", package="systemPipeR")
countDFeByg <- read.delim(countDFeBygpath, row.names=1)
edgeDF <- run_edgeR(countDF=countDFeByg, targets=targets, cmp=cmp[[1]], independent=FALSE, mdsplot="")

## ----DEGcounts, eval=TRUE, include=FALSE-----------------------------------------------------
DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=10))

## ----eval=TRUE-------------------------------------------------------------------------------
names(DEG_list)
DEG_list$Summary[1:4,]

## ----eval=TRUE-------------------------------------------------------------------------------
degseqDF <- run_DESeq2(countDF=countDFeByg, targets=targets, cmp=cmp[[1]], independent=FALSE)

## ----DEGcounts2, eval=TRUE, include=FALSE----------------------------------------------------
DEG_list2 <- filterDEGs(degDF=degseqDF, filter=c(Fold=2, FDR=10))

## ----vennplot, eval=TRUE, include=FALSE------------------------------------------------------
vennsetup <- overLapper(DEG_list$Up[6:9], type="vennsets")
vennsetdown <- overLapper(DEG_list$Down[6:9], type="vennsets")
vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))

## ----eval=FALSE------------------------------------------------------------------------------
#  library("biomaRt")
#  listMarts() # To choose BioMart database
#  m <- useMart("ENSEMBL_MART_PLANT"); listDatasets(m)
#  m <- useMart("ENSEMBL_MART_PLANT", dataset="athaliana_eg_gene")
#  listAttributes(m) # Choose data types you want to download
#  go <- getBM(attributes=c("go_accession", "tair_locus", "go_namespace_1003"), mart=m)
#  go <- go[go[,3]!="",]; go[,3] <- as.character(go[,3])
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

