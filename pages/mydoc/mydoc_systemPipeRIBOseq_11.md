---
title: 11. Differential ribosome loading analysis (translational efficiency)
last_updated: Sat Apr 18 12:33:39 2020
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeRIBOseq_11.html
---

Combined with mRNA-Seq data, Ribo-Seq or polyRibo-Seq experiments can be used
to study changes in translational efficiencies of genes and/or transcripts for
different treatments. For test purposes the following generates a small test
data set from the sample data used in this vignette, where two types of RNA
samples (`assays`) are considered: polyribosomal mRNA (`Ribo`)
and total mRNA (`mRNA`). In addition, there are two treatments
(`conditions`): `M1` and `A1`. 


```r
library(DESeq2)
countDFeBygpath <- system.file("extdata", "countDFeByg.xls", 
    package = "systemPipeR")
countDFeByg <- read.delim(countDFeBygpath, row.names = 1)
coldata <- DataFrame(assay = factor(rep(c("Ribo", "mRNA"), each = 4)), 
    condition = factor(rep(as.character(targets.as.df(targets(args))$Factor[1:4]), 
        2)), row.names = as.character(targets.as.df(targets(args))$SampleName)[1:8])
coldata
```

Differences in translational efficiencies can be calculated by ratios of ratios
for the two conditions: 

$$(Ribo\_A1 / mRNA\_A1) / (Ribo\_M1 / mRNA\_M1)$$


The latter can be modeled with the `DESeq2` package using the design $\sim assay + condition + assay:condition$, where the interaction term $assay:condition$ represents the ratio of ratios. Using the likelihood ratio test of `DESeq2`, which removes the interaction term in the reduced model, one can test whether the translational efficiency (ribosome loading) is different in condition `A1` than in `M1`.


```r
dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(countDFeByg[, 
    rownames(coldata)]), colData = coldata, design = ~assay + 
    condition + assay:condition)
# model.matrix(~ assay + condition + assay:condition,
# coldata) # Corresponding design matrix
dds <- DESeq2::DESeq(dds, test = "LRT", reduced = ~assay + condition)
res <- DESeq2::results(dds)
head(res[order(res$padj), ], 4)
# write.table(res, file='transleff.xls', quote=FALSE,
# col.names = NA, sep='\t')
```

<br><br><center><a href="mydoc_systemPipeRIBOseq_10.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeRIBOseq_12.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
