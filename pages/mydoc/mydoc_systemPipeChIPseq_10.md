---
title: 10. GO term enrichment analysis
last_updated: Sat Apr 18 12:27:46 2020
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeChIPseq_10.html
---

The following performs GO term enrichment analysis for each annotated peak set.


```r
dir_path <- system.file("extdata/cwl/annotate_peaks", package = "systemPipeR")
args <- loadWF(targets = "targets_bam_ref.txt", wf_file = "annotate-peaks.cwl", 
    input_file = "annotate-peaks.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName1 = "_FASTQ_PATH1_", 
    FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_"))

args_anno <- loadWF(targets = "targets_macs.txt", wf_file = "annotate-peaks.cwl", 
    input_file = "annotate-peaks.yml", dir_path = dir_path)
args_anno <- renderWF(args_anno, inputvars = c(FileName = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))
annofiles <- subsetWF(args_anno, slot = "output", subset = 1, 
    index = 1)
gene_ids <- sapply(names(annofiles), function(x) unique(as.character(read.delim(annofiles[x])[, 
    "geneId"])), simplify = FALSE)
load("data/GO/catdb.RData")
BatchResult <- GOCluster_Report(catdb = catdb, setlist = gene_ids, 
    method = "all", id_type = "gene", CLSZ = 2, cutoff = 0.9, 
    gocats = c("MF", "BP", "CC"), recordSpecGO = NULL)
```

<br><br><center><a href="mydoc_systemPipeChIPseq_09.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeChIPseq_11.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
