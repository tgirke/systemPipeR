---
title: 4. Sequences in Bioconductor
last_updated: Wed Jun  7 19:39:39 2017
sidebar: mydoc_sidebar
permalink: mydoc_Rsequences_04.html
---

## Important Data Objects of Biostrings

### `XString` for single sequence

* `DNAString`: for DNA
* `RNAString`: for RNA
* `AAString`: for amino acid 
* `BString`: for any string

### `XStringSet` for many sequences
        
* `DNAStringSet``: for DNA
* `RNAStringSet`: for RNA
* `AAStringSet`: for amino acid 
* `BStringSet`: for any string

### `QualityScaleXStringSet` for sequences with quality data

* `QualityScaledDNAStringSet`: for DNA
* `QualityScaledRNAStringSet`: for RNA
* `QualityScaledAAStringSet`: for amino acid 
* `QualityScaledBStringSet`: for any string

## Sequence Import and Export

Download the following sequences to your current working directory and then import them into R: 
[ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.ffn](ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.ffn)


```r
dir.create("data", showWarnings = FALSE)
# system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.ffn")
download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.ffn", "data/AE004437.ffn")
```

Import FASTA file with `readDNAStringSet`

```r
myseq <- readDNAStringSet("data/AE004437.ffn")
myseq[1:3]
```

```
##   A DNAStringSet instance of length 3
##     width seq                                                                   names               
## [1]  1206 ATGACTCGGCGGTCTCGTGTCGGTGCCGGCCTC...GTCGTCGTTGTTCGACGCTGGCGGAACCCATGA gi|12057215|gb|AE...
## [2]   666 ATGAGCATCATCGAACTCGAAGGCGTGGTCAAA...GTCAACCTCGTCGATGGGGTGTTACACACGTGA gi|12057215|gb|AE...
## [3]  1110 ATGGCGTGGCGGAACCTCGGGCGGAACCGCGTG...AACGATCCGCCCGTCGAGGCGCTCGGCGAATGA gi|12057215|gb|AE...
```

Subset sequences with regular expression on sequence name field

```r
sub <- myseq[grep("99.*", names(myseq))]
length(sub)
```

```
## [1] 170
```

Export subsetted sequences to FASTA file

```r
writeXStringSet(sub, file="./data/AE004437sub.ffn", width=80)
```

Now inspect exported sequence file `AE004437sub.ffn` in a text editor
	
    
## Working with `XString` Containers

The `XString` stores the different types of biosequences in dedicated containers

```r
library(Biostrings)
d <- DNAString("GCATAT-TAC")
d
```

```
##   10-letter "DNAString" instance
## seq: GCATAT-TAC
```

```r
d[1:4]
```

```
##   4-letter "DNAString" instance
## seq: GCAT
```

RNA sequences

```r
r <- RNAString("GCAUAU-UAC") 
r <- RNAString(d) # Converts d to RNAString object
r
```

```
##   10-letter "RNAString" instance
## seq: GCAUAU-UAC
```
Protein sequences

```r
p <- AAString("HCWYHH")
p
```

```
##   6-letter "AAString" instance
## seq: HCWYHH
```

Any type of character strings

```r
b <- BString("I store any set of characters. Other XString objects store only the IUPAC characters.")
b
```

```
##   85-letter "BString" instance
## seq: I store any set of characters. Other XString objects store only the IUPAC characters.
```

## Working with `XStringSet` Containers

`XStringSet` containers allow to store many biosequences in one object

```r
dset <- DNAStringSet(c("GCATATTAC", "AATCGATCC", "GCATATTAC")) 
names(dset) <- c("seq1", "seq2", "seq3") # Assigns names
dset[1:2]
```

```
##   A DNAStringSet instance of length 2
##     width seq                                                                   names               
## [1]     9 GCATATTAC                                                             seq1
## [2]     9 AATCGATCC                                                             seq2
```

Important utilities for `XStringSet` containers

```r
width(dset) # Returns the length of each sequences
```

```
## [1] 9 9 9
```

```r
d <- dset[[1]] # The [[ subsetting operator returns a single entry as XString object
dset2 <- c(dset, dset) # Appends/concatenates two XStringSet objects
dsetchar <- as.character(dset) # Converts XStringSet to named vector 
dsetone <- unlist(dset) # Collapses many sequences to a single one stored in a DNAString container
```

Sequence subsetting by positions:

```r
DNAStringSet(dset, start=c(1,2,3), end=c(4,8,5)) 
```

```
##   A DNAStringSet instance of length 3
##     width seq                                                                   names               
## [1]     4 GCAT                                                                  seq1
## [2]     7 ATCGATC                                                               seq2
## [3]     3 ATA                                                                   seq3
```

## Multiple Alignment Class

The `XMultipleAlignment` class stores the different types of multiple sequence alignments:


```r
origMAlign <- readDNAMultipleAlignment(filepath = system.file("extdata",
              "msx2_mRNA.aln", package = "Biostrings"), format = "clustal")
origMAlign
```

```
## DNAMultipleAlignment with 8 rows and 2343 columns
##      aln                                                                        names               
## [1] -----TCCCGTCTCCGCAGCAAAAAAGTTTGAGTCG...TTGTCCAAACTCACAATTAAAAAAAAAAAAAAAAA gi|84452153|ref|N...
## [2] ------------------------------------...----------------------------------- gi|208431713|ref|...
## [3] ------------------------------------...----------------------------------- gi|118601823|ref|...
## [4] ----------------------AAAAGTTGGAGTCT...----------------------------------- gi|114326503|ref|...
## [5] ------------------------------------...----------------------------------- gi|119220589|ref|...
## [6] ------------------------------------...----------------------------------- gi|148540149|ref|...
## [7] --------------CGGCTCCGCAGCGCCTCACTCG...----------------------------------- gi|45383056|ref|N...
## [8] GGGGGAGACTTCAGAAGTTGTTGTCCTCTCCGCTGA...----------------------------------- gi|213515133|ref|...
```

## Basic Sequence Manipulations

### Reverse and Complement


```r
randset <- DNAStringSet(rand)
complement(randset[1:2])
```

```
##   A DNAStringSet instance of length 2
##     width seq
## [1]    10 CCCATGCTGT
## [2]    14 TAGAAATACTGGAG
```

```r
reverse(randset[1:2])
```

```
##   A DNAStringSet instance of length 2
##     width seq
## [1]    10 ACAGCATGGG
## [2]    14 CTCCAGTATTTCTA
```

```r
reverseComplement(randset[1:2])
```

```
##   A DNAStringSet instance of length 2
##     width seq
## [1]    10 TGTCGTACCC
## [2]    14 GAGGTCATAAAGAT
```

## Translate DNA into Protein

```r
translate(randset[1:2])
```

```
## Warning in .Call2("DNAStringSet_translate", x, skip_code, dna_codes[codon_alphabet], : in 'x[[1]]':
## last base was ignored
```

```
## Warning in .Call2("DNAStringSet_translate", x, skip_code, dna_codes[codon_alphabet], : in 'x[[2]]':
## last 2 bases were ignored
```

```
##   A AAStringSet instance of length 2
##     width seq
## [1]     3 GYD
## [2]     4 IFMT
```

## Pattern Matching

### Pattern matching with mismatches

Find pattern matches in reference 

```r
myseq1 <- readDNAStringSet("./data/AE004437.ffn") 
mypos <- matchPattern("ATGGTG", myseq1[[1]], max.mismatch=1) 
```

Count only the corresponding matches

```r
countPattern("ATGGCT", myseq1[[1]], max.mismatch=1) 
```

```
## [1] 3
```

Count matches in many sequences

```r
vcountPattern("ATGGCT", myseq1, max.mismatch=1)[1:20]
```

```
##  [1] 3 0 5 4 1 2 2 1 4 3 0 0 1 2 0 1 4 0 0 1
```

Results shown in DNAStringSet object

```r
tmp <- c(DNAStringSet("ATGGTG"), DNAStringSet(mypos)) 
```

Return a consensus  matrix for query and hits

```r
consensusMatrix(tmp)[1:4,] 
```

```
##   [,1] [,2] [,3] [,4] [,5] [,6]
## A    3    0    0    0    0    0
## C    1    1    0    0    0    0
## G    0    0    4    4    1    4
## T    0    3    0    0    3    0
```

Find all pattern matches in reference

```r
myvpos <- vmatchPattern("ATGGCT", myseq1, max.mismatch=1) 
myvpos # The results are stored as MIndex object.
```

```
## MIndex object of length 2058
## $`gi|12057215|gb|AE004437.1|:248-1453 Halobacterium sp. NRC-1, complete genome`
## IRanges object with 3 ranges and 0 metadata columns:
##           start       end     width
##       <integer> <integer> <integer>
##   [1]         1         6         6
##   [2]       383       388         6
##   [3]       928       933         6
## 
## $`gi|12057215|gb|AE004437.1|:1450-2115 Halobacterium sp. NRC-1, complete genome`
## IRanges object with 0 ranges and 0 metadata columns:
##        start       end     width
##    <integer> <integer> <integer>
## 
## $`gi|12057215|gb|AE004437.1|:2145-3254 Halobacterium sp. NRC-1, complete genome`
## IRanges object with 5 ranges and 0 metadata columns:
##           start       end     width
##       <integer> <integer> <integer>
##   [1]         1         6         6
##   [2]        94        99         6
##   [3]       221       226         6
##   [4]       535       540         6
##   [5]       601       606         6
## 
## ...
## <2055 more elements>
```

```r
Views(myseq1[[1]], start(myvpos[[1]]), end(myvpos[[1]])) # Retrieves the result for single entry
```

```
##   Views on a 1206-letter DNAString subject
## subject: ATGACTCGGCGGTCTCGTGTCGGTGCCGGCCTCGCAGCCATTGT...TTGCGATCGTCGTCGTCGTTGTTCGACGCTGGCGGAACCCATGA
## views:
##     start end width
## [1]     1   6     6 [ATGACT]
## [2]   383 388     6 [ATGGCA]
## [3]   928 933     6 [ATGACT]
```

Return all matches

```r
sapply(seq(along=myseq1), function(x) 
       as.character(Views(myseq1[[x]], start(myvpos[[x]]), end(myvpos[[x]]))))[1:4] 
```

### Pattern matching with regular expression support


```r
myseq <- DNAStringSet(c("ATGCAGACATAGTG", "ATGAACATAGATCC", "GTACAGATCAC"))
myseq[grep("^ATG", myseq, perl=TRUE)] # String searching with regular expression support
```

```
##   A DNAStringSet instance of length 2
##     width seq
## [1]    14 ATGCAGACATAGTG
## [2]    14 ATGAACATAGATCC
```

```r
pos1 <- regexpr("AT", myseq) # Searches 'myseq' for first match of pattern "AT"
as.numeric(pos1); attributes(pos1)$match.length # Returns position information of matches
```

```
## [1] 1 1 7
```

```
## [1] 2 2 2
```

```r
pos2 <- gregexpr("AT", myseq) # Searches 'myseq' for all matches of pattern "AT"
as.numeric(pos2[[1]]); attributes(pos2[[1]])$match.length # Match positions in first sequence
```

```
## [1] 1 9
```

```
## [1] 2 2
```

```r
DNAStringSet(gsub("^ATG", "NNN", myseq)) # String substitution with regular expression support
```

```
##   A DNAStringSet instance of length 3
##     width seq
## [1]    14 NNNCAGACATAGTG
## [2]    14 NNNAACATAGATCC
## [3]    11 GTACAGATCAC
```

## PWM Viewing and Searching

### Plot with `seqLogo`


```r
library(seqLogo) 
pwm <- PWM(DNAStringSet(c("GCT", "GGT", "GCA"))) 
pwm
```

```
##        [,1]      [,2]      [,3]
## A 0.0000000 0.0000000 0.2312611
## C 0.0000000 0.3157205 0.0000000
## G 0.3685591 0.2312611 0.0000000
## T 0.0000000 0.0000000 0.3157205
```

```r
seqLogo(t(t(pwm) * 1/colSums(pwm)))
```

<img src="./pages/mydoc/Rsequences_files/pwm_logo-1.png" width="672" />

### Plot with `ggseqlogo`

The `ggseqlogo` package ([manual](https://omarwagih.github.io/ggseqlogo/))
provides many customization options for plotting sequence logos. It also supports
various alphabets including sequence logos for amino acid sequences.



```r
library(ggplot2); library(ggseqlogo)
pwm <- PWM(DNAStringSet(c("GCT", "GGT", "GCA"))) 
ggseqlogo(pwm)
```

<img src="./pages/mydoc/Rsequences_files/pwm_logo2-1.png" width="672" />

Search sequence for PWM matches with score better than `min.score`

```r
chr <- DNAString("AAAGCTAAAGGTAAAGCAAAA") 
matchPWM(pwm, chr, min.score=0.9) 
```

```
##   Views on a 21-letter DNAString subject
## subject: AAAGCTAAAGGTAAAGCAAAA
## views:
##     start end width
## [1]     4   6     3 [GCT]
## [2]    10  12     3 [GGT]
## [3]    16  18     3 [GCA]
```

<br><br><center><a href="mydoc_Rsequences_03.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_Rsequences_05.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
