---
title: 9. Summary statistics of variants
last_updated: Fri Jun 21 16:33:06 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeVARseq_09.html
---

The `varSummary` function counts the number of variants for each feature type
included in the anntation reports.

## Summary for `GATK`


```r
args <- systemArgs(sysma = "param/annotate_vars.param", mytargets = "targets_gatk_filtered.txt")
varSummary(args)
write.table(varSummary(args), "./results/variantStats_gatk.xls", 
    quote = FALSE, col.names = NA, sep = "\t")
```

## Summary for `BCFtools`


```r
args <- systemArgs(sysma = "param/annotate_vars.param", mytargets = "targets_sambcf_filtered.txt")
varSummary(args)
write.table(varSummary(args), "./results/variantStats_sambcf.xls", 
    quote = FALSE, col.names = NA, sep = "\t")
```

## Summary for `VariantTools`  


```r
args <- systemArgs(sysma = "param/annotate_vars.param", mytargets = "targets_vartools_filtered.txt")
varSummary(args)
write.table(varSummary(args), "./results/variantStats_vartools.xls", 
    quote = FALSE, col.names = NA, sep = "\t")
```

<br><br><center><a href="mydoc_systemPipeVARseq_08.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeVARseq_10.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
