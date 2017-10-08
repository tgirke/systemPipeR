---
title: 7. Transcript Ranges
last_updated: Wed Jun  7 19:39:39 2017
sidebar: mydoc_sidebar
permalink: mydoc_Rsequences_07.html
---

Storing annotation ranges in `TranscriptDb` databases makes many operations more robust and convenient.

```r
library(GenomicFeatures)
download.file("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/gff3.gff", "data/gff3.gff")
txdb <- makeTxDbFromGFF(file="data/gff3.gff", format="gff", dataSource="TAIR", organism="Arabidopsis thaliana")
```

```
## Warning in sort(classes): closing unused connection 7 (data/gff3.gff)
```

```
## Warning in sort(classes): closing unused connection 5 (data/gff3.gff)
```

```
## Warning in .extract_exons_from_GRanges(cds_IDX, gr, ID, Name, Parent, feature = "cds", : The following orphan CDS were dropped (showing only the 6 first):
##   seqid start  end strand   ID              Parent Name
## 1  Chr1  3760 3913      + <NA> AT1G01010.1-Protein <NA>
## 2  Chr1  3996 4276      + <NA> AT1G01010.1-Protein <NA>
## 3  Chr1  4486 4605      + <NA> AT1G01010.1-Protein <NA>
## 4  Chr1  4706 5095      + <NA> AT1G01010.1-Protein <NA>
## 5  Chr1  5174 5326      + <NA> AT1G01010.1-Protein <NA>
## 6  Chr1  5439 5630      + <NA> AT1G01010.1-Protein <NA>
```

```r
saveDb(txdb, file="./data/TAIR10.sqlite")
```

```
## TxDb object:
## # Db type: TxDb
## # Supporting package: GenomicFeatures
## # Data source: TAIR
## # Organism: Arabidopsis thaliana
## # Taxonomy ID: 3702
## # miRBase build ID: NA
## # Genome: NA
## # transcript_nrow: 28
## # exon_nrow: 113
## # cds_nrow: 99
## # Db created by: GenomicFeatures package from Bioconductor
## # Creation time: 2017-04-26 10:51:39 -0700 (Wed, 26 Apr 2017)
## # GenomicFeatures version at creation time: 1.28.0
## # RSQLite version at creation time: 1.1-2
## # DBSCHEMAVERSION: 1.1
```

```r
txdb <- loadDb("./data/TAIR10.sqlite")
transcripts(txdb)
```

```
## GRanges object with 28 ranges and 2 metadata columns:
##        seqnames         ranges strand |     tx_id     tx_name
##           <Rle>      <IRanges>  <Rle> | <integer> <character>
##    [1]     Chr1 [ 3631,  5899]      + |         1 AT1G01010.1
##    [2]     Chr1 [ 5928,  8737]      - |         2 AT1G01020.1
##    [3]     Chr1 [ 6790,  8737]      - |         3 AT1G01020.2
##    [4]     Chr1 [11649, 13714]      - |         4 AT1G01030.1
##    [5]     Chr2 [ 1025,  2810]      + |         5 AT2G01008.1
##    ...      ...            ...    ... .       ...         ...
##   [24]     ChrC [  383,  1444]      - |        24 ATCG00020.1
##   [25]     ChrC [ 1717,  4347]      - |        25 ATCG00030.1
##   [26]     ChrM [11918, 12241]      + |        26 ATMG00030.1
##   [27]     ChrM [  273,   734]      - |        27 ATMG00010.1
##   [28]     ChrM [ 8848, 11415]      - |        28 ATMG00020.1
##   -------
##   seqinfo: 7 sequences (2 circular) from an unspecified genome; no seqlengths
```

```r
transcriptsBy(txdb, by = "gene")
```

```
## GRangesList object of length 22:
## $AT1G01010 
## GRanges object with 1 range and 2 metadata columns:
##       seqnames       ranges strand |     tx_id     tx_name
##          <Rle>    <IRanges>  <Rle> | <integer> <character>
##   [1]     Chr1 [3631, 5899]      + |         1 AT1G01010.1
## 
## $AT1G01020 
## GRanges object with 2 ranges and 2 metadata columns:
##       seqnames       ranges strand | tx_id     tx_name
##   [1]     Chr1 [5928, 8737]      - |     2 AT1G01020.1
##   [2]     Chr1 [6790, 8737]      - |     3 AT1G01020.2
## 
## $AT1G01030 
## GRanges object with 1 range and 2 metadata columns:
##       seqnames         ranges strand | tx_id     tx_name
##   [1]     Chr1 [11649, 13714]      - |     4 AT1G01030.1
## 
## ...
## <19 more elements>
## -------
## seqinfo: 7 sequences (2 circular) from an unspecified genome; no seqlengths
```

```r
exonsBy(txdb, by = "gene")
```

```
## GRangesList object of length 22:
## $AT1G01010 
## GRanges object with 6 ranges and 2 metadata columns:
##       seqnames       ranges strand |   exon_id   exon_name
##          <Rle>    <IRanges>  <Rle> | <integer> <character>
##   [1]     Chr1 [3631, 3913]      + |         1        <NA>
##   [2]     Chr1 [3996, 4276]      + |         2        <NA>
##   [3]     Chr1 [4486, 4605]      + |         3        <NA>
##   [4]     Chr1 [4706, 5095]      + |         4        <NA>
##   [5]     Chr1 [5174, 5326]      + |         5        <NA>
##   [6]     Chr1 [5439, 5899]      + |         6        <NA>
## 
## $AT1G01020 
## GRanges object with 12 ranges and 2 metadata columns:
##        seqnames       ranges strand | exon_id exon_name
##    [1]     Chr1 [5928, 6263]      - |       7      <NA>
##    [2]     Chr1 [6437, 7069]      - |       8      <NA>
##    [3]     Chr1 [6790, 7069]      - |       9      <NA>
##    [4]     Chr1 [7157, 7232]      - |      10      <NA>
##    [5]     Chr1 [7157, 7450]      - |      11      <NA>
##    ...      ...          ...    ... .     ...       ...
##    [8]     Chr1 [7762, 7835]      - |      14      <NA>
##    [9]     Chr1 [7942, 7987]      - |      15      <NA>
##   [10]     Chr1 [8236, 8325]      - |      16      <NA>
##   [11]     Chr1 [8417, 8464]      - |      17      <NA>
##   [12]     Chr1 [8571, 8737]      - |      18      <NA>
## 
## $AT1G01030 
## GRanges object with 2 ranges and 2 metadata columns:
##       seqnames         ranges strand | exon_id exon_name
##   [1]     Chr1 [11649, 13173]      - |      19      <NA>
##   [2]     Chr1 [13335, 13714]      - |      20      <NA>
## 
## ...
## <19 more elements>
## -------
## seqinfo: 7 sequences (2 circular) from an unspecified genome; no seqlengths
```

## `txdb` from BioMart

Alternative sources for creating `txdb` databases are BioMart, Bioc annotation packages, UCSC, etc. The following shows how to create a `txdb` from BioMart.

```r
library(GenomicFeatures); library("biomaRt")
txdb <- makeTxDbFromBiomart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host="plants.ensembl.org")
```

The following steps are useful to find out what is availble in BioMart. 

```r
listMarts() # Lists BioMart databases
listMarts(host="plants.ensembl.org")
mymart <- useMart("plants_mart", host="plants.ensembl.org") # Select one, here plants_mart_25
listDatasets(mymart) # List datasets available in the selected BioMart database
mymart <- useMart("plants_mart", dataset="athaliana_eg_gene", host="plants.ensembl.org")
listAttributes(mymart) # List available features 
getBM(attributes=c("ensembl_gene_id", "description"), mart=mymart)[1:4,]
```

## Efficient Sequence Parsing 

### `getSeq`

The following parses all annotation ranges provided by a `GRanges` object (e.g. `gff`) from a genome sequence stored in a local file.

```r
gff <- gff[values(gff)$type != "chromosome"] # Remove chromosome ranges
rand <- DNAStringSet(sapply(unique(as.character(seqnames(gff))), function(x) paste(sample(c("A","T","G","C"), 200000, replace=T), collapse="")))
writeXStringSet(DNAStringSet(rand), "./data/test")
getSeq(FaFile("./data/test"), gff)
```

```
##   A DNAStringSet instance of length 442
##       width seq                                                                 names               
##   [1]  2269 TAAGTAAGGTGAGCAAGGAAAGCGCTACACGC...GATGGAGGCGCCCACACGATTAGAATGTTACA Chr1
##   [2]  2269 TAAGTAAGGTGAGCAAGGAAAGCGCTACACGC...GATGGAGGCGCCCACACGATTAGAATGTTACA Chr1
##   [3]  1871 CTAATGATACTGGGAGTACCGTTGCGATGACC...GCAAGCCCTATAAAATTGAAGGGAGGTGTATA Chr1
##   [4]   283 TAAGTAAGGTGAGCAAGGAAAGCGCTACACGC...CAGCAGAACAAACAGATACAGCCGGCCCAAGC Chr1
##   [5]   129 TAAGTAAGGTGAGCAAGGAAAGCGCTACACGC...AGACGTCCCGTCGAGCACGCCTTCCGATTATA Chr1
##   ...   ... ...
## [438]   324 AGTGTAGCTACCAGTTTATGTGGGTGCGCTTG...GCCCGCTCGAGACTGGAATGAGCAGTTGGGCT ChrM
## [439]   324 AGTGTAGCTACCAGTTTATGTGGGTGCGCTTG...GCCCGCTCGAGACTGGAATGAGCAGTTGGGCT ChrM
## [440]   324 AGTGTAGCTACCAGTTTATGTGGGTGCGCTTG...GCCCGCTCGAGACTGGAATGAGCAGTTGGGCT ChrM
## [441]   324 AGTGTAGCTACCAGTTTATGTGGGTGCGCTTG...GCCCGCTCGAGACTGGAATGAGCAGTTGGGCT ChrM
## [442]   324 AGTGTAGCTACCAGTTTATGTGGGTGCGCTTG...GCCCGCTCGAGACTGGAATGAGCAGTTGGGCT ChrM
```

### `extractTranscriptSeqs`

Sequences composed of several ranges, such as transcripts (or CDSs) with several exons, can be parsed with `extractTranscriptSeqs`. 
Note: the following expects the genome sequence in a file called `mygenome.fasta` and a valid `txdb` defining the ranges for that
genome.

```r
library(GenomicFeatures); library(Biostrings); library(Rsamtools)
extractTranscriptSeqs(FaFile("mygenome.fasta"), exonsBy(txdb, "tx", use.names=TRUE)) 
```

<br><br><center><a href="mydoc_Rsequences_06.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_Rsequences_08.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
