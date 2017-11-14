---
title: 9. Differential binding analysis
last_updated: Mon Nov 13 16:23:52 2017
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeChIPseq_09.html
---

The `runDiff` function performs differential binding analysis in batch mode for
several count tables using `edgeR` or `DESeq2` (Robinson et al., 2010; Love et al., 2014).
Internally, it calls the functions `run_edgeR` and `run_DESeq2`. It also returns 
the filtering results and plots from the downstream `filterDEGs` function using 
the fold change and FDR cutoffs provided under the `dbrfilter` argument.


```r
args_diff <- systemArgs(sysma="param/rundiff.param", mytargets="targets_countDF.txt")
cmp <- readComp(file=args_bam, format="matrix") 
dbrlist <- runDiff(args=args_diff, diffFct=run_edgeR, targets=targetsin(args_bam), 
                    cmp=cmp[[1]], independent=TRUE, dbrfilter=c(Fold=2, FDR=1))
writeTargetsout(x=args_diff, file="targets_rundiff.txt", overwrite=TRUE)
```


<br><br><center><a href="mydoc_systemPipeChIPseq_08.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeChIPseq_10.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
