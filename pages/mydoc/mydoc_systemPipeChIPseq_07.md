---
title: 7. Annotate peaks with genomic context
last_updated: Thu Nov 21 15:49:32 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeChIPseq_07.html
---

## Annotation with `ChIPpeakAnno` package

The following annotates the identified peaks with genomic context information using the `ChIPpeakAnno` and `ChIPseeker` packages, respectively (Zhu et al., 2010; Yu et al., 2015).


```r
library(ChIPpeakAnno)
library(GenomicFeatures)
dir_path <- system.file("extdata/cwl/annotate_peaks", package = "systemPipeR")
args <- loadWF(targets = "targets_macs.txt", wf_file = "annotate-peaks.cwl", 
    input_file = "annotate-peaks.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))

txdb <- makeTxDbFromGFF(file = "data/tair10.gff", format = "gff", 
    dataSource = "TAIR", organism = "Arabidopsis thaliana")
ge <- genes(txdb, columns = c("tx_name", "gene_id", "tx_type"))
for (i in seq(along = args)) {
    peaksGR <- as(read.delim(infile1(args)[i], comment = "#"), 
        "GRanges")
    annotatedPeak <- annotatePeakInBatch(peaksGR, AnnotationData = genes(txdb))
    df <- data.frame(as.data.frame(annotatedPeak), as.data.frame(values(ge[values(annotatedPeak)$feature, 
        ])))
    outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)
    write.table(df, outpaths[i], quote = FALSE, row.names = FALSE, 
        sep = "\t")
}
writeTargetsout(x = args, file = "targets_peakanno.txt", step = 1, 
    new_col = "FileName", new_col_output_index = 1, overwrite = TRUE)
```



The peak annotation results are written for each peak set to separate
files in the `results` directory. They are named after the corresponding peak
files with extensions specified in the `annotate_peaks.param` file, 
here `*.peaks.annotated.xls`.

## Annotation with `ChIPseeker` package

Same as in previous step but using the `ChIPseeker` package for annotating the peaks.


```r
library(ChIPseeker)
for (i in seq(along = args)) {
    peakAnno <- annotatePeak(infile1(args)[i], TxDb = txdb, verbose = FALSE)
    df <- as.data.frame(peakAnno)
    outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)
    write.table(df, outpaths[i], quote = FALSE, row.names = FALSE, 
        sep = "\t")
}
writeTargetsout(x = args, file = "targets_peakanno.txt", step = 1, 
    new_col = "FileName", new_col_output_index = 1, overwrite = TRUE)
```

Summary plots provided by the `ChIPseeker` package. Here applied only to one sample
for demonstration purposes.


```r
peak <- readPeakFile(infile1(args)[1])
covplot(peak, weightCol = "X.log10.pvalue.")
outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)
peakHeatmap(outpaths[1], TxDb = txdb, upstream = 1000, downstream = 1000, 
    color = "red")
plotAvgProf2(outpaths[1], TxDb = txdb, upstream = 1000, downstream = 1000, 
    xlab = "Genomic Region (5'->3')", ylab = "Read Count Frequency")
```

<br><br><center><a href="mydoc_systemPipeChIPseq_06.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeChIPseq_08.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
