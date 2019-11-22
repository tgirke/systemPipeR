---
title: 8. Read quantification per annotation range
last_updated: Thu Nov 21 16:49:12 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeRIBOseq_08.html
---

## Read counting with `summarizeOverlaps` in parallel mode using multiple cores

Reads overlapping with annotation ranges of interest are counted for each
sample using the `summarizeOverlaps` function (Lawrence et al., 2013). The
read counting is preformed for exonic gene regions in a non-strand-specific
manner while ignoring overlaps among different genes. Subsequently, the
expression count values are normalized by \textit{reads per kp per million
mapped reads} (RPKM). The raw read count table (`countDFeByg.xls`) and the corresponding
RPKM table (`rpkmDFeByg.xls`) are written to
separate files in the `results` directory of this project.
Parallelization is achieved with the `BiocParallel` package, here
using 8 CPU cores.


```r
library("GenomicFeatures")
library(BiocParallel)
txdb <- loadDb("./data/tair10.sqlite")
eByg <- exonsBy(txdb, by = c("gene"))
bfl <- BamFileList(outpaths, yieldSize = 50000, index = character())
multicoreParam <- MulticoreParam(workers = 8)
register(multicoreParam)
registered()
counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, 
    x, mode = "Union", ignore.strand = TRUE, inter.feature = FALSE, 
    singleEnd = TRUE))
countDFeByg <- sapply(seq(along = counteByg), function(x) assays(counteByg[[x]])$counts)
rownames(countDFeByg) <- names(rowRanges(counteByg[[1]]))
colnames(countDFeByg) <- names(bfl)
rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts = x, 
    ranges = eByg))
write.table(countDFeByg, "results/countDFeByg.xls", col.names = NA, 
    quote = FALSE, sep = "\t")
write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names = NA, 
    quote = FALSE, sep = "\t")
```

Sample of data slice of count table


```r
read.delim("results/countDFeByg.xls", row.names = 1, check.names = FALSE)[1:4, 
    1:5]
```

Sample of data slice of RPKM table


```r
read.delim("results/rpkmDFeByg.xls", row.names = 1, check.names = FALSE)[1:4, 
    1:4]
```

Note, for most statistical differential expression or abundance analysis
methods, such as `edgeR` or `DESeq2`, the raw count values
should be used as input. The usage of RPKM values should be restricted to
specialty applications required by some users, _e.g._ manually comparing
the expression levels among different genes or features. 

## Sample-wise correlation analysis

The following computes the sample-wise Spearman correlation coefficients from
the `rlog` transformed expression values generated with the
`DESeq2` package. After transformation to a distance matrix,
hierarchical clustering is performed with the `hclust` function and
the result is plotted as a dendrogram and written to a file named `sample_tree.png`
in the `results` directory. 


```r
library(DESeq2, quietly = TRUE)
library(ape, warn.conflicts = FALSE)
countDF <- as.matrix(read.table("./results/countDFeByg.xls"))
colData <- data.frame(row.names = targets.as.df(targets(args))$SampleName, 
    condition = targets.as.df(targets(args))$Factor)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = countDF, colData = colData, 
    design = ~condition)
d <- cor(assay(DESeq2::rlog(dds)), method = "spearman")
hc <- hclust(dist(1 - d))
png("results/sample_tree.pdf")
ape::plot.phylo(ape::as.phylo(hc), type = "p", edge.col = "blue", 
    edge.width = 2, show.node.label = TRUE, no.margin = TRUE)
dev.off()
```

![](./pages/mydoc/systemPipeRIBOseq_files/sample_tree.png)
<div align="center">Figure 5: Correlation dendrogram of samples.</div>


<br><br><center><a href="mydoc_systemPipeRIBOseq_07.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeRIBOseq_09.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
