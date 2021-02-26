---
title: 12. Clustering and heat maps
last_updated: Sat Apr 18 12:33:39 2020
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeRIBOseq_12.html
---

The following example performs hierarchical clustering on the `rlog` transformed expression matrix subsetted by the DEGs identified in the 
above differential expression analysis. It uses a Pearson correlation-based distance measure and complete linkage for cluster joining.


```r
library(pheatmap)
geneids <- unique(as.character(unlist(DEG_list[[1]])))
y <- assay(DESeq2::rlog(dds))[geneids, ]
png("heatmap1.png")
pheatmap(y, scale = "row", clustering_distance_rows = "correlation", 
    clustering_distance_cols = "correlation")
dev.off()
```

![](./pages/mydoc/systemPipeRIBOseq_files/heatmap1.png)
<div align="center">Figure 9: Heat map with hierarchical clustering dendrograms of DEGs</div>

<br><br><center><a href="mydoc_systemPipeRIBOseq_11.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeRIBOseq_13.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
