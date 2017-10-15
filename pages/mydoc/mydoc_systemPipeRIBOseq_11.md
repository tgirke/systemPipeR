---
title: 11. Differential ribosome loading analysis (translational efficiency)
last_updated: Sun Oct 15 13:21:42 2017
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeRIBOseq_11.html
---

Combinded with mRNA-Seq data, Ribo-Seq or polyRibo-Seq experiments can be used
to study changes in translational efficiencies of genes and/or transcripts for
different treatments. For test purposes the following generates a small test
data set from the sample data used in this vignette, where two types of RNA
samples (`assays`) are considered: polyribosomal mRNA (`Ribo`)
and total mRNA (`mRNA`). In addition, there are two treatments
(`conditions`): `M1` and `A1`. 


```r
library(DESeq2)
```

```
## 
## Attaching package: 'DESeq2'
```

```
## The following object is masked from 'package:systemPipeR':
## 
##     results
```

```r
targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
parampath <- system.file("extdata", "tophat.param", package="systemPipeR")
countDFeBygpath <- system.file("extdata", "countDFeByg.xls", package="systemPipeR")
args <- suppressWarnings(systemArgs(sysma=parampath, mytargets=targetspath))
countDFeByg <- read.delim(countDFeBygpath, row.names=1)
coldata <- DataFrame(assay=factor(rep(c("Ribo","mRNA"), each=4)), 
                condition=factor(rep(as.character(targetsin(args)$Factor[1:4]), 2)), 
                row.names=as.character(targetsin(args)$SampleName)[1:8])
coldata
```

```
## DataFrame with 8 rows and 2 columns
##        assay condition
##     <factor>  <factor>
## M1A     Ribo        M1
## M1B     Ribo        M1
## A1A     Ribo        A1
## A1B     Ribo        A1
## V1A     mRNA        M1
## V1B     mRNA        M1
## M6A     mRNA        A1
## M6B     mRNA        A1
```

Differences in translational efficiencies can be calculated by ratios of ratios
for the two conditions: 

$$(Ribo\_A1 / mRNA\_A1) / (Ribo\_M1 / mRNA\_M1)$$


The latter can be modeled with the `DESeq2` package using the design $\sim assay + condition + assay:condition$, 
where the interaction term $assay:condition$
represents the ratio of ratios. Using the likelihood ratio test of
`DESeq2`, which removes the interaction term in the reduced model, one
can test whether the translational efficiency (ribosome loading) is different
in condition `A1` than in `M1`.


```r
dds <- DESeqDataSetFromMatrix(countData=as.matrix(countDFeByg[,rownames(coldata)]), 
                            colData = coldata, 
                            design = ~ assay + condition + assay:condition)
# model.matrix(~ assay + condition + assay:condition, coldata) # Corresponding design matrix
dds <- DESeq(dds, test="LRT", reduced = ~ assay + condition)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```

```r
res <- DESeq2::results(dds)
head(res[order(res$padj),],4)
```

```
## log2 fold change (MLE): assayRibo.conditionM1 
## LRT p-value: '~ assay + condition + assay:condition' vs '~ assay + condition' 
## DataFrame with 4 rows and 6 columns
##             baseMean log2FoldChange     lfcSE      stat       pvalue         padj
##            <numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
## AT5G01040   94.03820      -5.824665 1.1621008  26.97122 2.065071e-07 2.106372e-05
## AT5G01015   25.22693       5.408013 1.1955448  19.37582 1.073577e-05 5.475243e-04
## AT5G01100  540.71442      -3.610594 0.8187969  18.44376 1.749925e-05 5.949746e-04
## AT2G01021 9227.94454       4.863866 1.2581171  13.45362 2.445341e-04 6.235619e-03
```

```r
# write.table(res, file="transleff.xls", quote=FALSE, col.names = NA, sep="\t")
```

<br><br><center><a href="mydoc_systemPipeRIBOseq_10.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeRIBOseq_12.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
