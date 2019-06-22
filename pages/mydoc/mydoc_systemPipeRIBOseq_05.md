---
title: 5. Read distribution across genomic features
last_updated: Fri Jun 21 16:34:15 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeRIBOseq_05.html
---

The `genFeatures` function generates a variety of feature types from
`TxDb` objects using utilities provided by the `GenomicFeatures` package. 

## Obtain feature types

The first step is the generation of the feature type ranges based on
annotations provided by a GFF file that can be transformed into a
`TxDb` object. This includes ranges for mRNAs, exons, introns, UTRs,
CDSs, miRNAs, rRNAs, tRNAs, promoter and intergenic regions. In addition, any
number of custom annotations can be included in this routine.


```r
library(GenomicFeatures)
txdb <- makeTxDbFromGFF(file = "data/tair10.gff", format = "gff3", 
    organism = "Arabidopsis")
feat <- genFeatures(txdb, featuretype = "all", reduce_ranges = TRUE, 
    upstream = 1000, downstream = 0, verbose = TRUE)
```

## Count and plot reads of any length

The `featuretypeCounts` function counts how many reads in short read
alignment files (BAM format) overlap with entire annotation categories. This
utility is useful for analyzing the distribution of the read mappings across
feature types, _e.g._ coding versus non-coding genes. By default the
read counts are reported for the sense and antisense strand of each feature
type separately. To minimize memory consumption, the BAM files are processed in
a stream using utilities from the `Rsamtools` and
`GenomicAlignment` packages.  The counts can be reported for each read
length separately or as a single value for reads of any length.  Subsequently,
the counting results can be plotted with the associated
`plotfeaturetypeCounts` function.

The following generates and plots feature counts for any read length.


```r
library(ggplot2)
library(grid)
fc <- featuretypeCounts(bfl = BamFileList(outpaths(args), yieldSize = 50000), 
    grl = feat, singleEnd = TRUE, readlength = NULL, type = "data.frame")
p <- plotfeaturetypeCounts(x = fc, graphicsfile = "results/featureCounts.png", 
    graphicsformat = "png", scales = "fixed", anyreadlength = TRUE, 
    scale_length_val = NULL)
```

![](./pages/mydoc/systemPipeRIBOseq_files/featureCounts.png)
<div align="center">Figure 2: Read distribution plot across annotation features for any read length.</div>

## Count and plot reads of specific lengths

To determine the approximate read length of ribosome footprints in Ribo-Seq experiments, one can generate and plot the feature counts for specific read lengths separately. Typically, the most abundant read length obtained for translated features corresponds to the approximate footprint length occupied by the ribosomes.


```r
fc2 <- featuretypeCounts(bfl = BamFileList(outpaths(args), yieldSize = 50000), 
    grl = feat, singleEnd = TRUE, readlength = c(74:76, 99:102), 
    type = "data.frame")
p2 <- plotfeaturetypeCounts(x = fc2, graphicsfile = "results/featureCounts2.png", 
    graphicsformat = "png", scales = "fixed", anyreadlength = FALSE, 
    scale_length_val = NULL)
```

![](./pages/mydoc/systemPipeRIBOseq_files/featureCounts2.png)
<div align="center">Figure 3: Read distribution plot across annotation features for specific read lengths.</div>

<br><br><center><a href="mydoc_systemPipeRIBOseq_04.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeRIBOseq_06.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
