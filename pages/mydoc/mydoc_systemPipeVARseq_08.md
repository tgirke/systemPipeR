---
title: 8. Combine annotation results among samples
last_updated: Sat Feb  2 11:47:13 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeVARseq_08.html
---

To simplify comparisons among samples, the `combineVarReports`
function combines all variant annotation reports referenced in a
`SYSargs` instance (here `args`). At the same time the function
allows to consider only certain feature types of interest. For instance, the
below setting `filtercol=c(Consequence="nonsynonymous")` will include
only nonsysynonymous variances listed in the `Consequence` column of
the annotation reports. To omit filtering, one can use the setting
`filtercol="All"`.

## Combine results from `GATK`  


```r
args <- systemArgs(sysma = "param/annotate_vars.param", mytargets = "targets_gatk_filtered.txt")
combineDF <- combineVarReports(args, filtercol = c(Consequence = "nonsynonymous"))
write.table(combineDF, "./results/combineDF_nonsyn_gatk.xls", 
    quote = FALSE, row.names = FALSE, sep = "\t")
```

## Combine results from `BCFtools`  


```r
args <- systemArgs(sysma = "param/annotate_vars.param", mytargets = "targets_sambcf_filtered.txt")
combineDF <- combineVarReports(args, filtercol = c(Consequence = "nonsynonymous"))
write.table(combineDF, "./results/combineDF_nonsyn_sambcf.xls", 
    quote = FALSE, row.names = FALSE, sep = "\t")
```

## Combine results from `VariantTools`


```r
args <- systemArgs(sysma = "param/annotate_vars.param", mytargets = "targets_vartools_filtered.txt")
combineDF <- combineVarReports(args, filtercol = c(Consequence = "nonsynonymous"))
write.table(combineDF, "./results/combineDF_nonsyn_vartools.xls", 
    quote = FALSE, row.names = FALSE, sep = "\t")
combineDF[2:4, ]
```

<br><br><center><a href="mydoc_systemPipeVARseq_07.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeVARseq_09.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
