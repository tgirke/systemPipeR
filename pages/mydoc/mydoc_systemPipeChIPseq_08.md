---
title: 8. Count reads overlapping peaks
last_updated: Sat Apr 18 12:27:46 2020
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeChIPseq_08.html
---

The `countRangeset` function is a convenience wrapper to perform read counting
iteratively over serveral range sets, here peak range sets. Internally,
the read counting is performed with the `summarizeOverlaps` function from the 
`GenomicAlignments` package. The resulting count tables are directly saved to 
files, one for each peak set.


```r
library(GenomicRanges)
dir_path <- system.file("extdata/cwl/count_rangesets", package = "systemPipeR")
args <- loadWF(targets = "targets_macs.txt", wf_file = "count_rangesets.cwl", 
    input_file = "count_rangesets.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))

## Bam Files
targets <- system.file("extdata", "targets_chip.txt", package = "systemPipeR")
dir_path <- system.file("extdata/cwl/bowtie2/bowtie2-se", package = "systemPipeR")
args_bam <- loadWF(targets = targets, wf_file = "bowtie2-mapping-se.cwl", 
    input_file = "bowtie2-mapping-se.yml", dir_path = dir_path)
args_bam <- renderWF(args_bam, inputvars = c(FileName = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))
args_bam <- output_update(args_bam, dir = FALSE, replace = TRUE, 
    extension = c(".sam", ".bam"))
outpaths <- subsetWF(args_bam, slot = "output", subset = 1, index = 1)

bfl <- BamFileList(outpaths, yieldSize = 50000, index = character())
countDFnames <- countRangeset(bfl, args, mode = "Union", ignore.strand = TRUE)
writeTargetsout(x = args, file = "targets_countDF.txt", step = 1, 
    new_col = "FileName", new_col_output_index = 1, overwrite = TRUE)
```

<br><br><center><a href="mydoc_systemPipeChIPseq_07.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeChIPseq_09.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
