---
title: 3. Strings in R Base
last_updated: Wed Jun  7 19:39:39 2017
sidebar: mydoc_sidebar
permalink: mydoc_Rsequences_03.html
---

## Basic String Matching and Parsing

### String matching

Generate sample sequence data set


```r
myseq <- c("ATGCAGACATAGTG", "ATGAACATAGATCC", "GTACAGATCAC")
```

String searching with regular expression support

```r
myseq[grep("ATG", myseq)] 
```

```
## [1] "ATGCAGACATAGTG" "ATGAACATAGATCC"
```

Searches `myseq` for first match of pattern "AT"

```r
pos1 <- regexpr("AT", myseq) 
as.numeric(pos1); attributes(pos1)$match.length # Returns position information of matches
```

```
## [1] 1 1 7
```

```
## [1] 2 2 2
```

Searches `myseq` for all matches of pattern "AT"

```r
pos2 <- gregexpr("AT", myseq) 
as.numeric(pos2[[1]]); attributes(pos2[[1]])$match.length # Returns positions of matches in first sequence
```

```
## [1] 1 9
```

```
## [1] 2 2
```

String substitution with regular expression support

```r
gsub("^ATG", "atg", myseq) 
```

```
## [1] "atgCAGACATAGTG" "atgAACATAGATCC" "GTACAGATCAC"
```

### Positional parsing

```r
nchar(myseq) # Computes length of strings
```

```
## [1] 14 14 11
```

```r
substring(myseq[1], c(1,3), c(2,5)) # Positional parsing of several fragments from one string
```

```
## [1] "AT"  "GCA"
```

```r
substring(myseq, c(1,4,7), c(2,6,10)) # Positional parsing of many strings
```

```
## [1] "AT"   "AAC"  "ATCA"
```

## Random Sequence Generation

### Random DNA sequences of any length


```r
rand <- sapply(1:100, function(x) paste(sample(c("A","T","G","C"), sample(10:20), replace=T), collapse=""))
rand[1:3]
```

```
## [1] "GGGTACGACA"       "ATCTTTATGACCTC"   "ATATCCTTGTGGGACA"
```

### Count identical sequences


```r
table(c(rand[1:4], rand[1]))
```

```
## 
## ATATCCTTGTGGGACA   ATCTTTATGACCTC       GGGTACGACA     TTTTCATACAAG 
##                1                1                2                1
```

### Extract reads from reference

Note: this requires `Biostrings` package.


```r
library(Biostrings)
ref <- DNAString(paste(sample(c("A","T","G","C"), 100000, replace=T), collapse=""))
randstart <- sample(1:(length(ref)-15), 1000)
randreads <- Views(ref, randstart, width=15)
rand_set <- DNAStringSet(randreads)
unlist(rand_set)
```

```
##   15000-letter "DNAString" instance
## seq: GGATTGTCGATTGCCCAGCGCCAAGAACTGTTATTTCCCTAACCCA...CATGATGCTAACCCGGTCGGTGTAATATGTATTCCCACGCCAGGCC
```

<br><br><center><a href="mydoc_Rsequences_02.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_Rsequences_04.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
