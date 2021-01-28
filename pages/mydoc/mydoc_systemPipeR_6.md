---
title: 6. Version information
last_updated: Thu Jan 28 13:38:31 2021
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeR_6.html
---

**Note:** the most recent version of this tutorial can be found <a href="http://www.bioconductor.org/packages/devel/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html">here</a>.


```r
sessionInfo()
```

```
## R version 4.0.3 (2020-10-10)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.1 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
## LAPACK: /home/dcassol/src/R-4.0.3/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] kableExtra_1.3.1            dplyr_1.0.2                
##  [3] DESeq2_1.30.0               magrittr_2.0.1             
##  [5] batchtools_0.9.14           ape_5.4-1                  
##  [7] ggplot2_3.3.2               systemPipeR_1.25.3         
##  [9] ShortRead_1.48.0            GenomicAlignments_1.26.0   
## [11] SummarizedExperiment_1.20.0 DelayedArray_0.16.0        
## [13] MatrixGenerics_1.2.0        Matrix_1.2-18              
## [15] matrixStats_0.57.0          Biobase_2.50.0             
## [17] BiocParallel_1.24.1         Rsamtools_2.6.0            
## [19] Biostrings_2.58.0           XVector_0.30.0             
## [21] GenomicRanges_1.42.0        GenomeInfoDb_1.26.2        
## [23] IRanges_2.24.1              S4Vectors_0.28.1           
## [25] BiocGenerics_0.36.0         BiocStyle_2.18.1           
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_2.0-0         rjson_0.2.20             hwriter_1.3.2           
##   [4] ellipsis_0.3.1           rstudioapi_0.13          farver_2.0.3            
##   [7] bit64_4.0.5              AnnotationDbi_1.52.0     xml2_1.3.2              
##  [10] codetools_0.2-18         splines_4.0.3            geneplotter_1.68.0      
##  [13] knitr_1.30               jsonlite_1.7.2           annotate_1.68.0         
##  [16] GO.db_3.12.1             dbplyr_2.0.0             png_0.1-7               
##  [19] pheatmap_1.0.12          graph_1.68.0             BiocManager_1.30.10     
##  [22] compiler_4.0.3           httr_1.4.2               GOstats_2.56.0          
##  [25] backports_1.2.1          assertthat_0.2.1         limma_3.46.0            
##  [28] formatR_1.7              htmltools_0.5.1.1        prettyunits_1.1.1       
##  [31] tools_4.0.3              gtable_0.3.0             glue_1.4.2              
##  [34] GenomeInfoDbData_1.2.4   Category_2.56.0          rsvg_2.1                
##  [37] rappdirs_0.3.1           V8_3.4.0                 Rcpp_1.0.5              
##  [40] vctrs_0.3.6              nlme_3.1-151             rtracklayer_1.50.0      
##  [43] xfun_0.20                stringr_1.4.0            rvest_0.3.6             
##  [46] lifecycle_0.2.0          XML_3.99-0.5             edgeR_3.32.0            
##  [49] zlibbioc_1.36.0          scales_1.1.1             BSgenome_1.58.0         
##  [52] VariantAnnotation_1.36.0 hms_0.5.3                RBGL_1.66.0             
##  [55] RColorBrewer_1.1-2       yaml_2.2.1               curl_4.3                
##  [58] memoise_1.1.0            biomaRt_2.46.0           latticeExtra_0.6-29     
##  [61] stringi_1.5.3            RSQLite_2.2.1            genefilter_1.72.0       
##  [64] checkmate_2.0.0          GenomicFeatures_1.42.1   DOT_0.1                 
##  [67] rlang_0.4.10             pkgconfig_2.0.3          bitops_1.0-6            
##  [70] evaluate_0.14            lattice_0.20-41          purrr_0.3.4             
##  [73] labeling_0.4.2           bit_4.0.4                tidyselect_1.1.0        
##  [76] GSEABase_1.52.1          AnnotationForge_1.32.0   bookdown_0.21           
##  [79] R6_2.5.0                 generics_0.1.0           base64url_1.4           
##  [82] DBI_1.1.0                pillar_1.4.7             withr_2.3.0             
##  [85] survival_3.2-7           RCurl_1.98-1.2           tibble_3.0.4            
##  [88] crayon_1.3.4             BiocFileCache_1.14.0     rmarkdown_2.6           
##  [91] jpeg_0.1-8.1             progress_1.2.2           locfit_1.5-9.4          
##  [94] grid_4.0.3               data.table_1.13.4        blob_1.2.1              
##  [97] Rgraphviz_2.34.0         webshot_0.5.2            digest_0.6.27           
## [100] xtable_1.8-4             brew_1.0-6               openssl_1.4.3           
## [103] munsell_0.5.0            viridisLite_0.3.0        askpass_1.1
```

<br><br><center><a href="mydoc_systemPipeR_5.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeR_7.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
