---
title: 3. Workflow overview
last_updated: Sun Oct 22 17:33:57 2017
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeR_3.html
---

## Define environment settings and samples

A typical workflow starts with generating the expected working environment
containing the proper directory structure, input files and parameter settings.
To simplify this task, one can load one of the existing NGS workflows templates
provided by _`systemPipeRdata`_ into the current working directory. The
following does this for the _`rnaseq`_ template. The name of the resulting
workflow directory can be specified under the _`mydirname`_ argument. The
default _`NULL`_ uses the name of the chosen workflow. An error is issued if a
directory of the same name and path exists already. On Linux and OS X systems
one can also create new workflow instances from the command-line of a terminal as shown
[here](http://bioconductor.org/packages/devel/data/experiment/vignettes/systemPipeRdata/inst/doc/systemPipeRdata.html#generate-workflow-template).
To apply workflows to custom data, the user needs to modify the _`targets`_ file and if
necessary update the corresponding _`param`_ file(s). A collection of pre-generated _`param`_
files is provided in the _`param`_ subdirectory of each workflow template. They
are also viewable in the GitHub repository of _`systemPipeRdata`_ ([see
here](https://github.com/tgirke/systemPipeRdata/tree/master/inst/extdata/param)).


```r
library(systemPipeR)
library(systemPipeRdata)
genWorkenvir(workflow="rnaseq", mydirname=NULL)
setwd("rnaseq")
```

Construct _`SYSargs`_ object from _`param`_ and _`targets`_ files.

```r
args <- systemArgs(sysma="param/trim.param", mytargets="targets.txt")
```

## Read Preprocessing
The function _`preprocessReads`_ allows to apply predefined or custom
read preprocessing functions to all FASTQ files referenced in a
_`SYSargs`_ container, such as quality filtering or adaptor trimming
routines. The paths to the resulting output FASTQ files are stored in the
_`outpaths`_ slot of the _`SYSargs`_ object. Internally,
_`preprocessReads`_ uses the _`FastqStreamer`_ function from
the _`ShortRead`_ package to stream through large FASTQ files in a
memory-efficient manner. The following example performs adaptor trimming with
the _`trimLRPatterns`_ function from the _`Biostrings`_ package.
After the trimming step a new targets file is generated (here
_`targets_trim.txt`_) containing the paths to the trimmed FASTQ files.
The new targets file can be used for the next workflow step with an updated
_`SYSargs`_ instance, _e.g._ running the NGS alignments with the
trimmed FASTQ files.
 

```r
preprocessReads(args=args, Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)", 
                batchsize=100000, overwrite=TRUE, compress=TRUE)
writeTargetsout(x=args, file="targets_trim.txt")
```

The following example shows how one can design a custom read preprocessing function 
using utilities provided by the _`ShortRead`_ package, and then run it
in batch mode with the _'preprocessReads'_ function (here on paired-end reads).

```r
args <- systemArgs(sysma="param/trimPE.param", mytargets="targetsPE.txt")
filterFct <- function(fq, cutoff=20, Nexceptions=0) {
    qcount <- rowSums(as(quality(fq), "matrix") <= cutoff)
    fq[qcount <= Nexceptions] # Retains reads where Phred scores are >= cutoff with N exceptions
}
preprocessReads(args=args, Fct="filterFct(fq, cutoff=20, Nexceptions=0)", batchsize=100000)
writeTargetsout(x=args, file="targets_PEtrim.txt")
```

## FASTQ quality report
The following _`seeFastq`_ and _`seeFastqPlot`_ functions generate and plot a series of
useful quality statistics for a set of FASTQ files including per cycle quality
box plots, base proportions, base-level quality trends, relative k-mer
diversity, length and occurrence distribution of reads, number of reads above
quality cutoffs and mean quality distribution.  
The function _`seeFastq`_ computes the quality statistics and stores the results in a
relatively small list object that can be saved to disk with _`save()`_ and
reloaded with _`load()`_ for later plotting. The argument _`klength`_ specifies the
k-mer length and _`batchsize`_ the number of reads to random sample from each
FASTQ file.


```r
fqlist <- seeFastq(fastq=infile1(args), batchsize=10000, klength=8)
pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
seeFastqPlot(fqlist)
dev.off()
```

![](./pages/mydoc/systemPipeR_files/fastqReport.png)
<div align="center"><b>Figure 2:</b> FASTQ quality report. To zoom in, rigth click image and open it in a separate browser tab. </div>


Parallelization of QC report on single machine with multiple cores


```r
args <- systemArgs(sysma="param/tophat.param", mytargets="targets.txt")
f <- function(x) seeFastq(fastq=infile1(args)[x], batchsize=100000, klength=8)
fqlist <- bplapply(seq(along=args), f, BPPARAM = MulticoreParam(workers=8))
seeFastqPlot(unlist(fqlist, recursive=FALSE))
```

Parallelization of QC report via scheduler (_e.g._ Torque) across several compute nodes


```r
library(BiocParallel); library(BatchJobs)
f <- function(x) {
    library(systemPipeR)
    args <- systemArgs(sysma="param/tophat.param", mytargets="targets.txt")
    seeFastq(fastq=infile1(args)[x], batchsize=100000, klength=8)
}
funs <- makeClusterFunctionsSLURM("slurm.tmpl")
param <- BatchJobsParam(length(args), resources=list(walltime="20:00:00", ntasks=1, ncpus=1, memory="6gb"), cluster.functions=funs)
register(param)
fqlist <- bplapply(seq(along=args), f)
seeFastqPlot(unlist(fqlist, recursive=FALSE))
```

## Alignment with _`Tophat2`_
Build _`Bowtie2`_ index.


```r
args <- systemArgs(sysma="param/tophat.param", mytargets="targets.txt")
moduleload(modules(args)) # Skip if module system is not available
system("bowtie2-build ./data/tair10.fasta ./data/tair10.fasta")
```

Execute _`SYSargs`_ on a single machine without submitting to a queuing system of a compute cluster. This way the input FASTQ files will be processed sequentially. If available, multiple CPU cores can be used for processing each file. The number of CPU cores (here 4) to use for each process is defined in the _`*.param`_ file. With _`cores(args)`_ one can return this value from the _`SYSargs`_ object. Note, if a module system is not installed or used, then the corresponding _`*.param`_ file needs to be edited accordingly by either providing an empty field in the line(s) starting with _`module`_ or by deleting these lines.


```r
bampaths <- runCommandline(args=args)
```

Alternatively, the computation can be greatly accelerated by processing many files in parallel using several compute nodes of a cluster, where a scheduling/queuing system is used for load balancing. To avoid over-subscription of CPU cores on the compute nodes, the value from _`cores(args)`_ is passed on to the submission command, here _`nodes`_ in the _`resources`_ list object. The number of independent parallel cluster processes is defined under the _`Njobs`_ argument. The following example will run 18 processes in parallel using for each 4 CPU cores. If the resources available on a cluster allow to run all 18 processes at the same time then the shown sample submission will utilize in total 72 CPU cores. Note, _`clusterRun`_ can be used with most queueing systems as it is based on utilities from the _`BatchJobs`_ package which supports the use of template files (_`*.tmpl`_) for defining the run parameters of different schedulers. To run the following code, one needs to have both a conf file (see _`.BatchJob`_ samples [here](https://code.google.com/p/batchjobs/wiki/DortmundUsage)) and a template file (see _`*.tmpl`_ samples [here](https://github.com/tudo-r/BatchJobs/tree/master/examples)) for the queueing available on a system. The following example uses the sample conf and template files for the Torque scheduler provided by this package.  


```r
resources <- list(walltime="20:00:00", ntasks=1, ncpus=cores(args), memory="10G")
reg <- clusterRun(args, conffile=".BatchJobs.R", template="slurm.tmpl", Njobs=18, runid="01", 
		  resourceList=resources)
waitForJobs(reg)
```

Useful commands for monitoring progress of submitted jobs

```r
showStatus(reg)
file.exists(outpaths(args))
sapply(1:length(args), function(x) loadResult(reg, x)) # Works after job completion
```

## Read and alignment count stats
Generate table of read and alignment counts for all samples. 

```r
read_statsDF <- alignStats(args) 
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
```

The following shows the first four lines of the sample alignment stats file provided by the _`systemPipeR`_ package. For simplicity the number of PE reads is multiplied here by 2 to approximate proper alignment frequencies where each read in a pair is counted. 

```r
read.table(system.file("extdata", "alignStats.xls", package="systemPipeR"), header=TRUE)[1:4,]
```

```
##   FileName Nreads2x Nalign Perc_Aligned Nalign_Primary Perc_Aligned_Primary
## 1      M1A   192918 177961     92.24697         177961             92.24697
## 2      M1B   197484 159378     80.70426         159378             80.70426
## 3      A1A   189870 176055     92.72397         176055             92.72397
## 4      A1B   188854 147768     78.24457         147768             78.24457
```

Parallelization of read/alignment stats on single machine with multiple cores

```r
f <- function(x) alignStats(args[x])
read_statsList <- bplapply(seq(along=args), f, BPPARAM = MulticoreParam(workers=8))
read_statsDF <- do.call("rbind", read_statsList)
```

Parallelization of read/alignment stats via scheduler (_e.g._ Torque) across several compute nodes 

```r
library(BiocParallel); library(BatchJobs)
f <- function(x) {
    library(systemPipeR)
    args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
    alignStats(args[x])
}
funs <- makeClusterFunctionsSLURM("slurm.tmpl")
param <- BatchJobsParam(length(args), resources=list(walltime="20:00:00", ntasks=1, ncpus=1, memory="6gb"), cluster.functions=funs)
register(param)
read_statsList <- bplapply(seq(along=args), f)
read_statsDF <- do.call("rbind", read_statsList)
```

## Create symbolic links for viewing BAM files in IGV
The genome browser IGV supports reading of indexed/sorted BAM files via web URLs. This way it can be avoided to create unnecessary copies of these large files. To enable this approach, an HTML directory with http access needs to be available in the user account (_e.g._ _`home/publichtml`_) of a system. If this is not the case then the BAM files need to be moved or copied to the system where IGV runs. In the following, _`htmldir`_ defines the path to the HTML directory with http access where the symbolic links to the BAM files will be stored. The corresponding URLs will be written to a text file specified under the `_urlfile`_ argument. 

```r
symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"), 
            urlbase="http://myserver.edu/~username/", 
        urlfile="IGVurl.txt")
```

## Alternative NGS Aligners

### Alignment with _`Bowtie2`_ (_e.g._ for miRNA profiling)
The following example runs _`Bowtie2`_ as a single process without submitting it to a cluster.

```r
args <- systemArgs(sysma="bowtieSE.param", mytargets="targets.txt")
moduleload(modules(args)) # Skip if module system is not available
bampaths <- runCommandline(args=args)
```

Alternatively, submit the job to compute nodes.

```r
resources <- list(walltime="20:00:00", ntasks=1, ncpus=cores(args), memory="10G")
reg <- clusterRun(args, conffile=".BatchJobs.R", template="slurm.tmpl", Njobs=18, runid="01", 
		  resourceList=resources)
waitForJobs(reg)
```

### Alignment with _`BWA-MEM`_ (_e.g._ for VAR-Seq)
The following example runs BWA-MEM as a single process without submitting it to a cluster.

```r
args <- systemArgs(sysma="param/bwa.param", mytargets="targets.txt")
moduleload(modules(args)) # Skip if module system is not available
system("bwa index -a bwtsw ./data/tair10.fasta") # Indexes reference genome
bampaths <- runCommandline(args=args[1:2])
```

### Alignment with _`Rsubread`_ (_e.g._ for RNA-Seq)
The following example shows how one can use within the `systemPipeR` environment the R-based 
aligner `Rsubread` or other R-based functions that read from input files and write to output files.

```r
library(Rsubread)
args <- systemArgs(sysma="param/rsubread.param", mytargets="targets.txt")
buildindex(basename=reference(args), reference=reference(args)) # Build indexed reference genome
align(index=reference(args), readfile1=infile1(args)[1:4], input_format="FASTQ", 
      output_file=outfile1(args)[1:4], output_format="SAM", nthreads=8, indels=1, TH1=2)
for(i in seq(along=outfile1(args))) asBam(file=outfile1(args)[i], destination=gsub(".sam", "", outfile1(args)[i]), overwrite=TRUE, indexDestination=TRUE)
```

### Alignment with _`gsnap`_ (_e.g._ for VAR-Seq and RNA-Seq)
Another R-based short read aligner is _`gsnap`_ from the _`gmapR`_ package (Wu et al., 2010). 
The code sample below introduces how to run this aligner on multiple nodes of a compute cluster. 

```r
library(gmapR); library(BiocParallel); library(BatchJobs)
args <- systemArgs(sysma="param/gsnap.param", mytargets="targetsPE.txt")
gmapGenome <- GmapGenome(reference(args), directory="data", name="gmap_tair10chr/", create=TRUE)
f <- function(x) {
    library(gmapR); library(systemPipeR)
    args <- systemArgs(sysma="gsnap.param", mytargets="targetsPE.txt")
    gmapGenome <- GmapGenome(reference(args), directory="data", name="gmap_tair10chr/", create=FALSE)
    p <- GsnapParam(genome=gmapGenome, unique_only=TRUE, molecule="DNA", max_mismatches=3)
    o <- gsnap(input_a=infile1(args)[x], input_b=infile2(args)[x], params=p, output=outfile1(args)[x])
}
funs <- makeClusterFunctionsSLURM("slurm.tmpl")
param <- BatchJobsParam(length(args), resources=list(walltime="20:00:00", ntasks=1, ncpus=1, memory="6gb"), cluster.functions=funs)
register(param)
d <- bplapply(seq(along=args), f)
```

## Read counting for mRNA profiling experiments
Create _`txdb`_ (needs to be done only once)

```r
library(GenomicFeatures)
txdb <- makeTxDbFromGFF(file="data/tair10.gff", format="gff", dataSource="TAIR", organism="A. thaliana")
saveDb(txdb, file="./data/tair10.sqlite")
```

The following performs read counting with _`summarizeOverlaps`_ in parallel mode with multiple cores. 

```r
library(BiocParallel)
txdb <- loadDb("./data/tair10.sqlite")
eByg <- exonsBy(txdb, by="gene")
bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
multicoreParam <- MulticoreParam(workers=4); register(multicoreParam); registered()
counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", ignore.strand=TRUE, inter.feature=TRUE, singleEnd=TRUE)) # Note: for strand-specific RNA-Seq set 'ignore.strand=FALSE' and for PE data set 'singleEnd=FALSE'
countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
```
Please note, in addition to read counts this step generates RPKM normalized expression values. For most statistical differential expression or abundance analysis methods, such as _`edgeR`_ or _`DESeq2`_, the raw count values should be used as input. The usage of RPKM values should be restricted to specialty applications required by some users, _e.g._ manually comparing the expression levels of different genes or features. 


Read counting with _`summarizeOverlaps`_ using multiple nodes of a cluster

```r
library(BiocParallel)
f <- function(x) {
    library(systemPipeR); library(BiocParallel); library(GenomicFeatures)
    txdb <- loadDb("./data/tair10.sqlite")
    eByg <- exonsBy(txdb, by="gene")
    args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
    bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
    summarizeOverlaps(eByg, bfl[x], mode="Union", ignore.strand=TRUE, inter.feature=TRUE, singleEnd=TRUE)
}
funs <- makeClusterFunctionsSLURM("slurm.tmpl")
param <- BatchJobsParam(length(args), resources=list(walltime="20:00:00", ntasks=1, ncpus=1, memory="6gb"), cluster.functions=funs)
register(param)
counteByg <- bplapply(seq(along=args), f) 
countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(outpaths(args))
```

## Read counting for miRNA profiling experiments
Download miRNA genes from miRBase

```r
system("wget ftp://mirbase.org/pub/mirbase/19/genomes/My_species.gff3 -P ./data/")
gff <- import.gff("./data/My_species.gff3")
gff <- split(gff, elementMetadata(gff)$ID)
bams <- names(bampaths); names(bams) <- targets$SampleName
bfl <- BamFileList(bams, yieldSize=50000, index=character())
countDFmiR <- summarizeOverlaps(gff, bfl, mode="Union", ignore.strand=FALSE, inter.feature=FALSE) # Note: inter.feature=FALSE important since pre and mature miRNA ranges overlap
rpkmDFmiR <- apply(countDFmiR, 2, function(x) returnRPKM(counts=x, gffsub=gff))
write.table(assays(countDFmiR)$counts, "results/countDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")
write.table(rpkmDFmiR, "results/rpkmDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")
```

## Correlation analysis of samples
The following computes the sample-wise Spearman correlation coefficients from the _`rlog`_ (regularized-logarithm) transformed expression values generated with the _`DESeq2`_ package. After transformation to a distance matrix, hierarchical clustering is performed with the _`hclust`_ function and the result is plotted as a dendrogram ([sample\_tree.pdf](./results/sample_tree.pdf)). 


```r
library(DESeq2, warn.conflicts=FALSE, quietly=TRUE); library(ape, warn.conflicts=FALSE)
countDFpath <- system.file("extdata", "countDFeByg.xls", package="systemPipeR")
countDF <- as.matrix(read.table(countDFpath))
colData <- data.frame(row.names=targetsin(args)$SampleName, condition=targetsin(args)$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
d <- cor(assay(rlog(dds)), method="spearman")
hc <- hclust(dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)
```

<img src="./pages/mydoc/systemPipeR_files/sample_tree_rlog-1.png" width="100%" />

<div align="center"><b>Figure 3:</b> Correlation dendrogram of samples for <i>rlog</i> values. </div>

Alternatively, the clustering can be performed with _`RPKM`_ normalized expression values. In combination with Spearman correlation the results of the two clustering methods are often relatively similar. 


```r
rpkmDFeBygpath <- system.file("extdata", "rpkmDFeByg.xls", package="systemPipeR")
rpkmDFeByg <- read.table(rpkmDFeBygpath, check.names=FALSE)
rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
d <- cor(rpkmDFeByg, method="spearman")
hc <- hclust(as.dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
```


## DEG analysis with _`edgeR`_
The following _`run_edgeR`_ function is a convenience wrapper for
identifying differentially expressed genes (DEGs) in batch mode with
_`edgeR`_'s GML method (Robinson et al., 2010) for any number of
pairwise sample comparisons specified under the _`cmp`_ argument. Users
are strongly encouraged to consult the 
[_`edgeR`_](\href{http://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) vignette 
for more detailed information on this topic and how to properly run _`edgeR`_ 
on data sets with more complex experimental designs. 

```r
targets <- read.delim(targetspath, comment="#")
cmp <- readComp(file=targetspath, format="matrix", delim="-")
cmp[[1]]
```

```
##       [,1]  [,2] 
##  [1,] "M1"  "A1" 
##  [2,] "M1"  "V1" 
##  [3,] "A1"  "V1" 
##  [4,] "M6"  "A6" 
##  [5,] "M6"  "V6" 
##  [6,] "A6"  "V6" 
##  [7,] "M12" "A12"
##  [8,] "M12" "V12"
##  [9,] "A12" "V12"
```

```r
countDFeBygpath <- system.file("extdata", "countDFeByg.xls", package="systemPipeR")
countDFeByg <- read.delim(countDFeBygpath, row.names=1)
edgeDF <- run_edgeR(countDF=countDFeByg, targets=targets, cmp=cmp[[1]], independent=FALSE, mdsplot="")
```

```
## Disp = 0.20653 , BCV = 0.4545
```

Filter and plot DEG results for up and down regulated genes. Because of the small size of the toy data set used by this vignette, the _FDR_ value has been set to a relatively high threshold (here 10%). More commonly used _FDR_ cutoffs are 1% or 5%. The definition of '_up_' and '_down_' is given in the corresponding help file. To open it, type _`?filterDEGs`_ in the R console. 

```r
DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=10))
```

<img src="./pages/mydoc/systemPipeR_files/edger_deg_counts-1.png" width="100%"  class="widefigure" />
<div align="center"><b>Figure 4:</b> Up and down regulated DEGs identified by <i>edgeR</i>. </div>


```r
names(DEG_list)
```

```
## [1] "UporDown" "Up"       "Down"     "Summary"
```

```r
DEG_list$Summary[1:4,]
```

```
##       Comparisons Counts_Up_or_Down Counts_Up Counts_Down
## M1-A1       M1-A1                 0         0           0
## M1-V1       M1-V1                 1         1           0
## A1-V1       A1-V1                 1         1           0
## M6-A6       M6-A6                 0         0           0
```

## DEG analysis with _`DESeq2`_ 
The following _`run_DESeq2`_ function is a convenience wrapper for
identifying DEGs in batch mode with _`DESeq2`_ (Love et al., 2014) for any number of
pairwise sample comparisons specified under the _`cmp`_ argument. Users
are strongly encouraged to consult the 
[_`DESeq2`_](http://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf) vignette
for more detailed information on this topic and how to properly run _`DESeq2`_ 
on data sets with more complex experimental designs. 

```r
degseqDF <- run_DESeq2(countDF=countDFeByg, targets=targets, cmp=cmp[[1]], independent=FALSE)
```

Filter and plot DEG results for up and down regulated genes. 

```r
DEG_list2 <- filterDEGs(degDF=degseqDF, filter=c(Fold=2, FDR=10))
```

<img src="./pages/mydoc/systemPipeR_files/deseq2_deg_counts-1.png" width="100%"  class="widefigure" />
<div align="center"><b>Figure 5:</b> Up and down regulated DEGs identified by <i>DESeq2</i>. </div>


## Venn Diagrams
The function _`overLapper`_ can compute Venn intersects for large numbers of sample sets (up to 20 or more) and _`vennPlot`_ can plot 2-5 way Venn diagrams. A useful feature is the possiblity to combine the counts from several Venn comparisons with the same number of sample sets in a single Venn diagram (here for 4 up and down DEG sets).

```r
vennsetup <- overLapper(DEG_list$Up[6:9], type="vennsets")
vennsetdown <- overLapper(DEG_list$Down[6:9], type="vennsets")
vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))
```

<img src="./pages/mydoc/systemPipeR_files/vennplot-1.png" width="100%" />
<div align="center"><b>Figure 6:</b> Venn Diagram for 4 Up and Down DEG Sets. </div>


## GO term enrichment analysis of DEGs
### Obtain gene-to-GO mappings
The following shows how to obtain gene-to-GO mappings from _`biomaRt`_ (here for _A. thaliana_) and how to organize them for the downstream GO term enrichment analysis. Alternatively, the gene-to-GO mappings can be obtained for many organisms from Bioconductor's  _`*.db`_ genome annotation packages or GO annotation files provided by various genome databases. For each annotation this relatively slow preprocessing step needs to be performed only once. Subsequently, the preprocessed data can be loaded with the _`load`_ function as shown in the next subsection. 

```r
library("biomaRt")
listMarts() # To choose BioMart database
m <- useMart("ENSEMBL_MART_PLANT"); listDatasets(m) 
m <- useMart("ENSEMBL_MART_PLANT", dataset="athaliana_eg_gene")
listAttributes(m) # Choose data types you want to download
go <- getBM(attributes=c("go_accession", "tair_locus", "go_namespace_1003"), mart=m)
go <- go[go[,3]!="",]; go[,3] <- as.character(go[,3])
dir.create("./data/GO")
write.table(go, "data/GO/GOannotationsBiomart_mod.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
catdb <- makeCATdb(myfile="data/GO/GOannotationsBiomart_mod.txt", lib=NULL, org="", colno=c(1,2,3), idconv=NULL)
save(catdb, file="data/GO/catdb.RData") 
```

### Batch GO term enrichment analysis
Apply the enrichment analysis to the DEG sets obtained in the above differential expression analysis. Note, in the following example the _FDR_ filter is set here to an unreasonably high value, simply because of the small size of the toy data set used in this vignette. Batch enrichment analysis of many gene sets is performed with the _`GOCluster_Report`_ function. When _`method="all"`_, it returns all GO terms passing the p-value cutoff specified under the _`cutoff`_ arguments. When _`method="slim"`_, it returns only the GO terms specified under the _`myslimv`_ argument. The given example shows how one can obtain such a GO slim vector from BioMart for a specific organism.  

```r
load("data/GO/catdb.RData")
DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=50), plot=FALSE)
up_down <- DEG_list$UporDown; names(up_down) <- paste(names(up_down), "_up_down", sep="")
up <- DEG_list$Up; names(up) <- paste(names(up), "_up", sep="")
down <- DEG_list$Down; names(down) <- paste(names(down), "_down", sep="")
DEGlist <- c(up_down, up, down)
DEGlist <- DEGlist[sapply(DEGlist, length) > 0]
BatchResult <- GOCluster_Report(catdb=catdb, setlist=DEGlist, method="all", id_type="gene", CLSZ=2, cutoff=0.9, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
library("biomaRt"); m <- useMart("ENSEMBL_MART_PLANT", dataset="athaliana_eg_gene")
goslimvec <- as.character(getBM(attributes=c("goslim_goa_accession"), mart=m)[,1])
BatchResultslim <- GOCluster_Report(catdb=catdb, setlist=DEGlist, method="slim", id_type="gene", myslimv=goslimvec, CLSZ=10, cutoff=0.01, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
```

### Plot batch GO term results
The _`data.frame`_ generated by _`GOCluster_Report`_ can be plotted with the _`goBarplot`_ function. Because of the variable size of the sample sets, it may not always be desirable to show the results from different DEG sets in the same bar plot. Plotting single sample sets is achieved by subsetting the input data frame as shown in the first line of the following example. 

```r
gos <- BatchResultslim[grep("M6-V6_up_down", BatchResultslim$CLID), ]
gos <- BatchResultslim
pdf("GOslimbarplotMF.pdf", height=8, width=10); goBarplot(gos, gocat="MF"); dev.off()
goBarplot(gos, gocat="BP")
goBarplot(gos, gocat="CC")
```
![](./pages/mydoc/systemPipeR_files/GOslimbarplotMF.png)
<div align="center"><b>Figure 7:</b> GO Slim Barplot for MF Ontology.</div>



## Clustering and heat maps
The following example performs hierarchical clustering on the _`rlog`_ transformed expression matrix subsetted by the DEGs identified in the 
above differential expression analysis. It uses a Pearson correlation-based distance measure and complete linkage for cluster joining.

```r
library(pheatmap)
geneids <- unique(as.character(unlist(DEG_list[[1]])))
y <- assay(rlog(dds))[geneids, ]
pdf("heatmap1.pdf")
pheatmap(y, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
dev.off()
```
![](./pages/mydoc/systemPipeR_files/heatmap1.png)
<div align="center"><b>Figure 8:</b> Heat map with hierarchical clustering dendrograms of DEGs.</div>


<br><br><center><a href="mydoc_systemPipeR_2.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeR_4.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
