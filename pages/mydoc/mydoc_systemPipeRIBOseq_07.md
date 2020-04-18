---
title: 7. Genomic read coverage along transripts or CDSs
last_updated: Sat Apr 18 12:33:39 2020
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeRIBOseq_07.html
---

The `featureCoverage` function computes the read coverage along
single and multi component features based on genomic alignments. The coverage
segments of component features are spliced to continuous ranges, such as exons
to transcripts or CDSs to ORFs. The results can be obtained with single
nucleotide resolution (_e.g._ around start and stop codons) or as mean coverage
of relative bin sizes, such as 100 bins for each feature. The latter allows
comparisons of coverage trends among transcripts of variable length. Additionally, 
the results can be obtained for single or many features (_e.g._ any number of
transcripts) at once. Visualization of the coverage results is facilitated by
the downstream `plotfeatureCoverage` function. 

## Binned CDS coverage to compare many transcripts


```r
grl <- cdsBy(txdb, "tx", use.names = TRUE)
fcov <- featureCoverage(bfl = BamFileList(outpaths[1:2]), grl = grl[1:4], 
    resizereads = NULL, readlengthrange = NULL, Nbins = 20, method = mean, 
    fixedmatrix = FALSE, resizefeatures = TRUE, upstream = 20, 
    downstream = 20, outfile = "results/featureCoverage.xls", 
    overwrite = TRUE)
```

## Coverage upstream and downstream of start and stop codons


```r
fcov <- featureCoverage(bfl = BamFileList(outpaths[1:4]), grl = grl[1:12], 
    resizereads = NULL, readlengthrange = NULL, Nbins = NULL, 
    method = mean, fixedmatrix = TRUE, resizefeatures = TRUE, 
    upstream = 20, downstream = 20, outfile = "results/featureCoverage.xls", 
    overwrite = TRUE)
plotfeatureCoverage(covMA = fcov, method = mean, scales = "fixed", 
    extendylim = 2, scale_count_val = 10^6)
```

## Combined coverage for both binned CDS and start/stop codons


```r
library(ggplot2)
library(grid)
fcov <- featureCoverage(bfl = BamFileList(outpaths[1:4]), grl = grl[1:4], 
    resizereads = NULL, readlengthrange = NULL, Nbins = 20, method = mean, 
    fixedmatrix = TRUE, resizefeatures = TRUE, upstream = 20, 
    downstream = 20, outfile = "results/featureCoverage.xls", 
    overwrite = TRUE)
png("./results/featurePlot.png", height = 12, width = 24, units = "in", 
    res = 72)
plotfeatureCoverage(covMA = fcov, method = mean, scales = "fixed", 
    extendylim = 2, scale_count_val = 10^6)
dev.off()
```

![](./pages/mydoc/systemPipeRIBOseq_files/featurePlot.png)
<div align="center">Figure 4: Feature coverage plot with single nucleotide resolution around start and stop codons and binned coverage between them.</div>

## Nucleotide level coverage along entire transcripts/CDSs


```r
fcov <- featureCoverage(bfl = BamFileList(outpaths[1:2]), grl = grl[1], 
    resizereads = NULL, readlengthrange = NULL, Nbins = NULL, 
    method = mean, fixedmatrix = FALSE, resizefeatures = TRUE, 
    upstream = 20, downstream = 20, outfile = NULL)
```

<br><br><center><a href="mydoc_systemPipeRIBOseq_06.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeRIBOseq_08.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
