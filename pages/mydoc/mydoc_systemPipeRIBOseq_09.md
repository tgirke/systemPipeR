---
title: 9. Analysis of differentially expressed genes with `edgeR`
last_updated: Thu Nov 21 16:49:12 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeRIBOseq_09.html
---

The analysis of differentially expressed genes (DEGs) is performed with the glm
method from the `edgeR` package (Robinson et al., 2010). The sample
comparisons used by this analysis are defined in the header lines of the
`targets.txt` file starting with `<CMP>`.


```r
library(edgeR)
countDF <- read.delim("results/countDFeByg.xls", row.names = 1, 
    check.names = FALSE)
targets <- read.delim("targets.txt", comment = "#")
cmp <- readComp(file = "targets.txt", format = "matrix", delim = "-")
edgeDF <- run_edgeR(countDF = countDF, targets = targets, cmp = cmp[[1]], 
    independent = FALSE, mdsplot = "")
```

Add functional gene descriptions, here from `biomaRt`. 


```r
library("biomaRt")
m <- useMart("plants_mart", dataset = "athaliana_eg_gene", host = "plants.ensembl.org")
desc <- getBM(attributes = c("tair_locus", "description"), mart = m)
desc <- desc[!duplicated(desc[, 1]), ]
descv <- as.character(desc[, 2])
names(descv) <- as.character(desc[, 1])
edgeDF <- data.frame(edgeDF, Desc = descv[rownames(edgeDF)], 
    check.names = FALSE)
write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote = FALSE, 
    sep = "\t", col.names = NA)
```

Filter and plot DEG results for up and down regulated genes. The definition of
`up` and `down` is given in the corresponding help file. To
open it, type `?filterDEGs` in the R console. 


```r
edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names = 1, 
    check.names = FALSE)
png("./results/DEGcounts.png", height = 10, width = 10, units = "in", 
    res = 72)
DEG_list <- filterDEGs(degDF = edgeDF, filter = c(Fold = 2, FDR = 20))
dev.off()
write.table(DEG_list$Summary, "./results/DEGcounts.xls", quote = FALSE, 
    sep = "\t", row.names = FALSE)
```

![](./pages/mydoc/systemPipeRIBOseq_files/DEGcounts.png)
<div align="center">Figure 6: Up and down regulated DEGs.</div>

The function `overLapper` can compute Venn intersects for large
numbers of sample sets (up to 20 or more) and `vennPlot` can plot 2-5
way Venn diagrams. A useful feature is the possiblity to combine the counts
from several Venn comparisons with the same number of sample sets in a single
Venn diagram (here for 4 up and down DEG sets).


```r
vennsetup <- overLapper(DEG_list$Up[6:9], type = "vennsets")
vennsetdown <- overLapper(DEG_list$Down[6:9], type = "vennsets")
png("results/vennplot.png")
vennPlot(list(vennsetup, vennsetdown), mymain = "", mysub = "", 
    colmode = 2, ccol = c("blue", "red"))
dev.off()
```

![](./pages/mydoc/systemPipeRIBOseq_files/vennplot.png)
<div align="center">Figure 7: Venn Diagram for 4 Up and Down DEG Sets</div>

<br><br><center><a href="mydoc_systemPipeRIBOseq_08.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeRIBOseq_10.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
