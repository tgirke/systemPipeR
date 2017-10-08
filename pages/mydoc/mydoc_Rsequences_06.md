---
title: 6. Range Operations  
last_updated: Wed Jun  7 19:39:39 2017
sidebar: mydoc_sidebar
permalink: mydoc_Rsequences_06.html
---

## Important Data Objects for Range Operations

* `IRanges`: stores range data only (IRanges library)
* `GRanges`: stores ranges and annotations (GenomicRanges library)
* `GRangesList`: list version of GRanges container (GenomicRanges library)

## Range Data Are Stored in `IRanges` and `GRanges` Containers

### Construct `GRanges` Object 


```r
library(GenomicRanges); library(rtracklayer)
gr <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)), ranges = IRanges(1:10, end = 7:16, names = head(letters, 10)), strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)), score = 1:10, GC = seq(1, 0, length = 10)) # Example of creating a GRanges object with its constructor function.
```

### Import GFF into `GRanges` Object

```r
gff <- import.gff("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/gff3.gff") # Imports a simplified GFF3 genome annotation file.
seqlengths(gff) <- end(ranges(gff[which(values(gff)[,"type"]=="chromosome"),])) 
names(gff) <- 1:length(gff) # Assigns names to corresponding slot
gff[1:4,]
```

```
## GRanges object with 4 ranges and 10 metadata columns:
##     seqnames           ranges strand |   source       type     score     phase                  ID
##        <Rle>        <IRanges>  <Rle> | <factor>   <factor> <numeric> <integer>         <character>
##   1     Chr1 [   1, 30427671]      + |   TAIR10 chromosome      <NA>      <NA>                Chr1
##   2     Chr1 [3631,     5899]      + |   TAIR10       gene      <NA>      <NA>           AT1G01010
##   3     Chr1 [3631,     5899]      + |   TAIR10       mRNA      <NA>      <NA>         AT1G01010.1
##   4     Chr1 [3760,     5630]      + |   TAIR10    protein      <NA>      <NA> AT1G01010.1-Protein
##            Name                Note          Parent       Index Derives_from
##     <character>     <CharacterList> <CharacterList> <character>  <character>
##   1        Chr1                                            <NA>         <NA>
##   2   AT1G01010 protein_coding_gene                        <NA>         <NA>
##   3 AT1G01010.1                           AT1G01010           1         <NA>
##   4 AT1G01010.1                                            <NA>  AT1G01010.1
##   -------
##   seqinfo: 7 sequences from an unspecified genome
```

### Coerce `GRanges` object to `data.frame`

```r
as.data.frame(gff)[1:4, 1:7]
```

```
##   seqnames start      end    width strand source       type
## 1     Chr1     1 30427671 30427671      + TAIR10 chromosome
## 2     Chr1  3631     5899     2269      + TAIR10       gene
## 3     Chr1  3631     5899     2269      + TAIR10       mRNA
## 4     Chr1  3760     5630     1871      + TAIR10    protein
```

### Coerce `GRanges` to `RangedData` object and vice versa

```r
gff_rd <- as(gff, "RangedData") 
```

```
## Warning in cntxt$tailcall: closing unused connection 7 (/tmp/Rtmp7jsa8n/file78fc7c51c885)
```

```
## Warning in cntxt$tailcall: closing unused connection 5 (/tmp/Rtmp7jsa8n/file78fc7c51c885)
```

```r
gff_gr <- as(gff_rd, "GRanges") 
```

## Utilities for Range Containers

### Accessor and subsetting methods for GRanges objects

Subsetting and replacement

```r
gff[1:4]
```

```
## GRanges object with 4 ranges and 10 metadata columns:
##     seqnames           ranges strand |   source       type     score     phase                  ID
##        <Rle>        <IRanges>  <Rle> | <factor>   <factor> <numeric> <integer>         <character>
##   1     Chr1 [   1, 30427671]      + |   TAIR10 chromosome      <NA>      <NA>                Chr1
##   2     Chr1 [3631,     5899]      + |   TAIR10       gene      <NA>      <NA>           AT1G01010
##   3     Chr1 [3631,     5899]      + |   TAIR10       mRNA      <NA>      <NA>         AT1G01010.1
##   4     Chr1 [3760,     5630]      + |   TAIR10    protein      <NA>      <NA> AT1G01010.1-Protein
##            Name                Note          Parent       Index Derives_from
##     <character>     <CharacterList> <CharacterList> <character>  <character>
##   1        Chr1                                            <NA>         <NA>
##   2   AT1G01010 protein_coding_gene                        <NA>         <NA>
##   3 AT1G01010.1                           AT1G01010           1         <NA>
##   4 AT1G01010.1                                            <NA>  AT1G01010.1
##   -------
##   seqinfo: 7 sequences from an unspecified genome
```

```r
gff[1:4, c("type", "ID")] 
```

```
## GRanges object with 4 ranges and 2 metadata columns:
##     seqnames           ranges strand |       type                  ID
##        <Rle>        <IRanges>  <Rle> |   <factor>         <character>
##   1     Chr1 [   1, 30427671]      + | chromosome                Chr1
##   2     Chr1 [3631,     5899]      + |       gene           AT1G01010
##   3     Chr1 [3631,     5899]      + |       mRNA         AT1G01010.1
##   4     Chr1 [3760,     5630]      + |    protein AT1G01010.1-Protein
##   -------
##   seqinfo: 7 sequences from an unspecified genome
```

```r
gff[2] <- gff[3] 
```

GRanges objects can be concatenated with the `c` function

```r
c(gff[1:2], gff[401:402]) 
```

```
## GRanges object with 4 ranges and 10 metadata columns:
##       seqnames           ranges strand |   source           type     score     phase
##          <Rle>        <IRanges>  <Rle> | <factor>       <factor> <numeric> <integer>
##     1     Chr1 [   1, 30427671]      + |   TAIR10     chromosome      <NA>      <NA>
##     2     Chr1 [3631,     5899]      + |   TAIR10           mRNA      <NA>      <NA>
##   401     Chr5 [5516,     5769]      - |   TAIR10        protein      <NA>      <NA>
##   402     Chr5 [5770,     5801]      - |   TAIR10 five_prime_UTR      <NA>      <NA>
##                        ID        Name            Note          Parent       Index Derives_from
##               <character> <character> <CharacterList> <CharacterList> <character>  <character>
##     1                Chr1        Chr1                                        <NA>         <NA>
##     2         AT1G01010.1 AT1G01010.1                       AT1G01010           1         <NA>
##   401 AT5G01015.2-Protein AT5G01015.2                                        <NA>  AT5G01015.2
##   402                <NA>        <NA>                     AT5G01015.2        <NA>         <NA>
##   -------
##   seqinfo: 7 sequences from an unspecified genome
```

Acessor functions

```r
seqnames(gff)
```

```
## factor-Rle of length 449 with 7 runs
##   Lengths:   72   22   38  118  172   13   14
##   Values : Chr1 Chr2 Chr3 Chr4 Chr5 ChrC ChrM
## Levels(7): Chr1 Chr2 Chr3 Chr4 Chr5 ChrC ChrM
```

```r
ranges(gff)
```

```
## IRanges object with 449 ranges and 0 metadata columns:
##           start       end     width
##       <integer> <integer> <integer>
##     1         1  30427671  30427671
##     2      3631      5899      2269
##     3      3631      5899      2269
##     4      3760      5630      1871
##     5      3631      3913       283
##   ...       ...       ...       ...
##   445     11918     12241       324
##   446     11918     12241       324
##   447     11918     12241       324
##   448     11918     12241       324
##   449     11918     12241       324
```

```r
strand(gff)
```

```
## factor-Rle of length 449 with 13 runs
##   Lengths:  18  54  28  21  12 117   1 171   1  12   1   8   5
##   Values :   +   -   +   -   +   -   +   -   +   -   +   -   +
## Levels(3): + - *
```

```r
seqlengths(gff) 
```

```
##     Chr1     Chr2     Chr3     Chr4     Chr5     ChrC     ChrM 
## 30427671 19698289 23459830 18585056 26975502   154478   366924
```

```r
start(gff[1:4])
```

```
## [1]    1 3631 3631 3760
```

```r
end(gff[1:4])
```

```
## [1] 30427671     5899     5899     5630
```

```r
width(gff[1:4]) 
```

```
## [1] 30427671     2269     2269     1871
```

Accessing metadata component

```r
values(gff) # or elementMetadata(gff)
```

```
## DataFrame with 449 rows and 10 columns
##       source       type     score     phase                  ID        Name                Note
##     <factor>   <factor> <numeric> <integer>         <character> <character>     <CharacterList>
## 1     TAIR10 chromosome        NA        NA                Chr1        Chr1                    
## 2     TAIR10       mRNA        NA        NA         AT1G01010.1 AT1G01010.1                    
## 3     TAIR10       mRNA        NA        NA         AT1G01010.1 AT1G01010.1                    
## 4     TAIR10    protein        NA        NA AT1G01010.1-Protein AT1G01010.1                    
## 5     TAIR10       exon        NA        NA                  NA          NA                    
## ...      ...        ...       ...       ...                 ...         ...                 ...
## 445   TAIR10       gene        NA        NA           ATMG00030   ATMG00030 protein_coding_gene
## 446   TAIR10       mRNA        NA        NA         ATMG00030.1 ATMG00030.1                    
## 447   TAIR10    protein        NA        NA ATMG00030.1-Protein ATMG00030.1                    
## 448   TAIR10       exon        NA        NA                  NA          NA                    
## 449   TAIR10        CDS        NA         0                  NA          NA                    
##                              Parent       Index Derives_from
##                     <CharacterList> <character>  <character>
## 1                                            NA           NA
## 2                         AT1G01010           1           NA
## 3                         AT1G01010           1           NA
## 4                                            NA  AT1G01010.1
## 5                       AT1G01010.1          NA           NA
## ...                             ...         ...          ...
## 445                                          NA           NA
## 446                       ATMG00030           1           NA
## 447                                          NA  ATMG00030.1
## 448                     ATMG00030.1          NA           NA
## 449 ATMG00030.1,ATMG00030.1-Protein          NA           NA
```

```r
values(gff)[, "type"][1:20] 
```

```
##  [1] chromosome      mRNA            mRNA            protein         exon            five_prime_UTR 
##  [7] CDS             exon            CDS             exon            CDS             exon           
## [13] CDS             exon            CDS             exon            CDS             three_prime_UTR
## [19] gene            mRNA           
## Levels: chromosome gene mRNA protein exon five_prime_UTR CDS three_prime_UTR rRNA tRNA
```

```r
gff[values(gff)[ ,"type"] == "gene"] 
```

```
## GRanges object with 21 ranges and 10 metadata columns:
##       seqnames         ranges strand |   source     type     score     phase          ID
##          <Rle>      <IRanges>  <Rle> | <factor> <factor> <numeric> <integer> <character>
##    19     Chr1 [ 5928,  8737]      - |   TAIR10     gene      <NA>      <NA>   AT1G01020
##    64     Chr1 [11649, 13714]      - |   TAIR10     gene      <NA>      <NA>   AT1G01030
##    74     Chr2 [ 1025,  2810]      + |   TAIR10     gene      <NA>      <NA>   AT2G01008
##    84     Chr2 [ 3706,  5513]      + |   TAIR10     gene      <NA>      <NA>   AT2G01010
##    87     Chr2 [ 5782,  5945]      + |   TAIR10     gene      <NA>      <NA>   AT2G01020
##   ...      ...            ...    ... .      ...      ...       ...       ...         ...
##   427     ChrC [  383,  1444]      - |   TAIR10     gene      <NA>      <NA>   ATCG00020
##   432     ChrC [ 1717,  4347]      - |   TAIR10     gene      <NA>      <NA>   ATCG00030
##   437     ChrM [  273,   734]      - |   TAIR10     gene      <NA>      <NA>   ATMG00010
##   442     ChrM [ 8848, 11415]      - |   TAIR10     gene      <NA>      <NA>   ATMG00020
##   445     ChrM [11918, 12241]      + |   TAIR10     gene      <NA>      <NA>   ATMG00030
##              Name                Note          Parent       Index Derives_from
##       <character>     <CharacterList> <CharacterList> <character>  <character>
##    19   AT1G01020 protein_coding_gene                        <NA>         <NA>
##    64   AT1G01030 protein_coding_gene                        <NA>         <NA>
##    74   AT2G01008 protein_coding_gene                        <NA>         <NA>
##    84   AT2G01010                rRNA                        <NA>         <NA>
##    87   AT2G01020                rRNA                        <NA>         <NA>
##   ...         ...                 ...             ...         ...          ...
##   427   ATCG00020 protein_coding_gene                        <NA>         <NA>
##   432   ATCG00030                tRNA                        <NA>         <NA>
##   437   ATMG00010 protein_coding_gene                        <NA>         <NA>
##   442   ATMG00020                rRNA                        <NA>         <NA>
##   445   ATMG00030 protein_coding_gene                        <NA>         <NA>
##   -------
##   seqinfo: 7 sequences from an unspecified genome
```

### Useful utilities for GRanges objects

Remove chromosome ranges

```r
gff <- gff[values(gff)$type != "chromosome"] 
```

Erase the strand information

```r
strand(gff) <- "*" 
```

Collapses overlapping ranges to continuous ranges.

```r
reduce(gff) 
```

```
## GRanges object with 22 ranges and 0 metadata columns:
##        seqnames         ranges strand
##           <Rle>      <IRanges>  <Rle>
##    [1]     Chr1 [ 3631,  5899]      *
##    [2]     Chr1 [ 5928,  8737]      *
##    [3]     Chr1 [11649, 13714]      *
##    [4]     Chr2 [ 1025,  2810]      *
##    [5]     Chr2 [ 3706,  5513]      *
##    ...      ...            ...    ...
##   [18]     ChrC [  383,  1444]      *
##   [19]     ChrC [ 1717,  4347]      *
##   [20]     ChrM [  273,   734]      *
##   [21]     ChrM [ 8848, 11415]      *
##   [22]     ChrM [11918, 12241]      *
##   -------
##   seqinfo: 7 sequences from an unspecified genome
```

Return uncovered regions

```r
gaps(gff) 
```

```
## GRanges object with 43 ranges and 0 metadata columns:
##        seqnames           ranges strand
##           <Rle>        <IRanges>  <Rle>
##    [1]     Chr1 [   1, 30427671]      +
##    [2]     Chr1 [   1, 30427671]      -
##    [3]     Chr1 [   1,     3630]      *
##    [4]     Chr1 [5900,     5927]      *
##    [5]     Chr1 [8738,    11648]      *
##    ...      ...              ...    ...
##   [39]     ChrM  [    1, 366924]      -
##   [40]     ChrM  [    1,    272]      *
##   [41]     ChrM  [  735,   8847]      *
##   [42]     ChrM  [11416,  11917]      *
##   [43]     ChrM  [12242, 366924]      *
##   -------
##   seqinfo: 7 sequences from an unspecified genome
```

More intuitive way to get uncovered regions

```r
setdiff(as(seqinfo(gff), "GRanges"), gff) 
```

```
## GRanges object with 29 ranges and 0 metadata columns:
##        seqnames            ranges strand
##           <Rle>         <IRanges>  <Rle>
##    [1]     Chr1 [    1,     3630]      *
##    [2]     Chr1 [ 5900,     5927]      *
##    [3]     Chr1 [ 8738,    11648]      *
##    [4]     Chr1 [13715, 30427671]      *
##    [5]     Chr2 [    1,     1024]      *
##    ...      ...               ...    ...
##   [25]     ChrC   [ 4348, 154478]      *
##   [26]     ChrM   [    1,    272]      *
##   [27]     ChrM   [  735,   8847]      *
##   [28]     ChrM   [11416,  11917]      *
##   [29]     ChrM   [12242, 366924]      *
##   -------
##   seqinfo: 7 sequences from an unspecified genome
```

Return disjoint ranges

```r
disjoin(gff)
```

```
## GRanges object with 211 ranges and 0 metadata columns:
##         seqnames         ranges strand
##            <Rle>      <IRanges>  <Rle>
##     [1]     Chr1   [3631, 3759]      *
##     [2]     Chr1   [3760, 3913]      *
##     [3]     Chr1   [3914, 3995]      *
##     [4]     Chr1   [3996, 4276]      *
##     [5]     Chr1   [4277, 4485]      *
##     ...      ...            ...    ...
##   [207]     ChrC [ 1752,  4310]      *
##   [208]     ChrC [ 4311,  4347]      *
##   [209]     ChrM [  273,   734]      *
##   [210]     ChrM [ 8848, 11415]      *
##   [211]     ChrM [11918, 12241]      *
##   -------
##   seqinfo: 7 sequences from an unspecified genome
```

Returns coverage of ranges

```r
coverage(gff)
```

```
## RleList of length 7
## $Chr1
## integer-Rle of length 30427671 with 45 runs
##   Lengths:     3630      129      154       82      281 ...      233      161      380 30413957
##   Values :        0        4        5        3        5 ...        4        2        4        0
## 
## $Chr2
## integer-Rle of length 19698289 with 14 runs
##   Lengths:     1024      248      185       53      362 ...      164      625      102 19691617
##   Values :        0        5        3        5        3 ...        3        0        5        0
## 
## $Chr3
## integer-Rle of length 23459830 with 29 runs
##   Lengths:     1652      145      139      111       95 ...      155      148      156 23453781
##   Values :        0        4        5        3        5 ...        3        5        4        0
## 
## $Chr4
## integer-Rle of length 18585056 with 72 runs
##   Lengths:     1179      357     1358      128      872 ...      212      114       74 18571697
##   Values :        0        5        0        5        3 ...        3        5        4        0
## 
## $Chr5
## integer-Rle of length 26975502 with 64 runs
##   Lengths:     1222       28       28      109       72 ...       76       55      174 26967058
##   Values :        0        4        7       13       16 ...        3        5        4        0
## 
## ...
## <2 more elements>
```

Return the index pairings for overlapping ranges

```r
findOverlaps(gff, gff[1:4]) 
```

```
## Hits object with 55 hits and 0 metadata columns:
##        queryHits subjectHits
##        <integer>   <integer>
##    [1]         1           1
##    [2]         1           2
##    [3]         1           4
##    [4]         1           3
##    [5]         2           1
##    ...       ...         ...
##   [51]        16           1
##   [52]        16           2
##   [53]        16           3
##   [54]        17           1
##   [55]        17           2
##   -------
##   queryLength: 442 / subjectLength: 4
```

Counts overlapping ranges 

```r
countOverlaps(gff, gff[1:4])[1:40]
```

```
##  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 
##  4  4  4  4  3  4  3  3  3  3  3  3  3  3  3  3  2  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 35 36 37 38 39 40 41 
##  0  0  0  0  0  0  0
```

Return only overlapping ranges 

```r
subsetByOverlaps(gff, gff[1:4]) 
```

```
## GRanges object with 17 ranges and 10 metadata columns:
##      seqnames       ranges strand |   source            type     score     phase
##         <Rle>    <IRanges>  <Rle> | <factor>        <factor> <numeric> <integer>
##    2     Chr1 [3631, 5899]      * |   TAIR10            mRNA      <NA>      <NA>
##    3     Chr1 [3631, 5899]      * |   TAIR10            mRNA      <NA>      <NA>
##    4     Chr1 [3760, 5630]      * |   TAIR10         protein      <NA>      <NA>
##    5     Chr1 [3631, 3913]      * |   TAIR10            exon      <NA>      <NA>
##    6     Chr1 [3631, 3759]      * |   TAIR10  five_prime_UTR      <NA>      <NA>
##   ..      ...          ...    ... .      ...             ...       ...       ...
##   14     Chr1 [5174, 5326]      * |   TAIR10            exon      <NA>      <NA>
##   15     Chr1 [5174, 5326]      * |   TAIR10             CDS      <NA>         0
##   16     Chr1 [5439, 5899]      * |   TAIR10            exon      <NA>      <NA>
##   17     Chr1 [5439, 5630]      * |   TAIR10             CDS      <NA>         0
##   18     Chr1 [5631, 5899]      * |   TAIR10 three_prime_UTR      <NA>      <NA>
##                       ID        Name            Note                          Parent       Index
##              <character> <character> <CharacterList>                 <CharacterList> <character>
##    2         AT1G01010.1 AT1G01010.1                                       AT1G01010           1
##    3         AT1G01010.1 AT1G01010.1                                       AT1G01010           1
##    4 AT1G01010.1-Protein AT1G01010.1                                                        <NA>
##    5                <NA>        <NA>                                     AT1G01010.1        <NA>
##    6                <NA>        <NA>                                     AT1G01010.1        <NA>
##   ..                 ...         ...             ...                             ...         ...
##   14                <NA>        <NA>                                     AT1G01010.1        <NA>
##   15                <NA>        <NA>                 AT1G01010.1,AT1G01010.1-Protein        <NA>
##   16                <NA>        <NA>                                     AT1G01010.1        <NA>
##   17                <NA>        <NA>                 AT1G01010.1,AT1G01010.1-Protein        <NA>
##   18                <NA>        <NA>                                     AT1G01010.1        <NA>
##      Derives_from
##       <character>
##    2         <NA>
##    3         <NA>
##    4  AT1G01010.1
##    5         <NA>
##    6         <NA>
##   ..          ...
##   14         <NA>
##   15         <NA>
##   16         <NA>
##   17         <NA>
##   18         <NA>
##   -------
##   seqinfo: 7 sequences from an unspecified genome
```

## GRangesList Objects


```r
sp <- split(gff, seq(along=gff)) # Stores every range in separate component of a GRangesList object
split(gff, seqnames(gff)) # Stores ranges of each chromosome in separate component.
```

```
## GRangesList object of length 7:
## $Chr1 
## GRanges object with 71 ranges and 10 metadata columns:
##      seqnames         ranges strand |   source            type     score     phase
##         <Rle>      <IRanges>  <Rle> | <factor>        <factor> <numeric> <integer>
##    2     Chr1   [3631, 5899]      * |   TAIR10            mRNA      <NA>      <NA>
##    3     Chr1   [3631, 5899]      * |   TAIR10            mRNA      <NA>      <NA>
##    4     Chr1   [3760, 5630]      * |   TAIR10         protein      <NA>      <NA>
##    5     Chr1   [3631, 3913]      * |   TAIR10            exon      <NA>      <NA>
##    6     Chr1   [3631, 3759]      * |   TAIR10  five_prime_UTR      <NA>      <NA>
##   ..      ...            ...    ... .      ...             ...       ...       ...
##   68     Chr1 [13335, 13714]      * |   TAIR10            exon      <NA>      <NA>
##   69     Chr1 [12941, 13173]      * |   TAIR10  five_prime_UTR      <NA>      <NA>
##   70     Chr1 [11864, 12940]      * |   TAIR10             CDS      <NA>         0
##   71     Chr1 [11649, 11863]      * |   TAIR10 three_prime_UTR      <NA>      <NA>
##   72     Chr1 [11649, 13173]      * |   TAIR10            exon      <NA>      <NA>
##                       ID        Name            Note                          Parent       Index
##              <character> <character> <CharacterList>                 <CharacterList> <character>
##    2         AT1G01010.1 AT1G01010.1                                       AT1G01010           1
##    3         AT1G01010.1 AT1G01010.1                                       AT1G01010           1
##    4 AT1G01010.1-Protein AT1G01010.1                                                        <NA>
##    5                <NA>        <NA>                                     AT1G01010.1        <NA>
##    6                <NA>        <NA>                                     AT1G01010.1        <NA>
##   ..                 ...         ...             ...                             ...         ...
##   68                <NA>        <NA>                                     AT1G01030.1        <NA>
##   69                <NA>        <NA>                                     AT1G01030.1        <NA>
##   70                <NA>        <NA>                 AT1G01030.1,AT1G01030.1-Protein        <NA>
##   71                <NA>        <NA>                                     AT1G01030.1        <NA>
##   72                <NA>        <NA>                                     AT1G01030.1        <NA>
##      Derives_from
##       <character>
##    2         <NA>
##    3         <NA>
##    4  AT1G01010.1
##    5         <NA>
##    6         <NA>
##   ..          ...
##   68         <NA>
##   69         <NA>
##   70         <NA>
##   71         <NA>
##   72         <NA>
## 
## ...
## <6 more elements>
## -------
## seqinfo: 7 sequences from an unspecified genome
```

```r
unlist(sp) # Returns data as GRanges object
```

```
## GRanges object with 442 ranges and 10 metadata columns:
##           seqnames         ranges strand |   source           type     score     phase
##              <Rle>      <IRanges>  <Rle> | <factor>       <factor> <numeric> <integer>
##       1.2     Chr1   [3631, 5899]      * |   TAIR10           mRNA      <NA>      <NA>
##       2.3     Chr1   [3631, 5899]      * |   TAIR10           mRNA      <NA>      <NA>
##       3.4     Chr1   [3760, 5630]      * |   TAIR10        protein      <NA>      <NA>
##       4.5     Chr1   [3631, 3913]      * |   TAIR10           exon      <NA>      <NA>
##       5.6     Chr1   [3631, 3759]      * |   TAIR10 five_prime_UTR      <NA>      <NA>
##       ...      ...            ...    ... .      ...            ...       ...       ...
##   438.445     ChrM [11918, 12241]      * |   TAIR10           gene      <NA>      <NA>
##   439.446     ChrM [11918, 12241]      * |   TAIR10           mRNA      <NA>      <NA>
##   440.447     ChrM [11918, 12241]      * |   TAIR10        protein      <NA>      <NA>
##   441.448     ChrM [11918, 12241]      * |   TAIR10           exon      <NA>      <NA>
##   442.449     ChrM [11918, 12241]      * |   TAIR10            CDS      <NA>         0
##                            ID        Name                Note                          Parent
##                   <character> <character>     <CharacterList>                 <CharacterList>
##       1.2         AT1G01010.1 AT1G01010.1                                           AT1G01010
##       2.3         AT1G01010.1 AT1G01010.1                                           AT1G01010
##       3.4 AT1G01010.1-Protein AT1G01010.1                                                    
##       4.5                <NA>        <NA>                                         AT1G01010.1
##       5.6                <NA>        <NA>                                         AT1G01010.1
##       ...                 ...         ...                 ...                             ...
##   438.445           ATMG00030   ATMG00030 protein_coding_gene                                
##   439.446         ATMG00030.1 ATMG00030.1                                           ATMG00030
##   440.447 ATMG00030.1-Protein ATMG00030.1                                                    
##   441.448                <NA>        <NA>                                         ATMG00030.1
##   442.449                <NA>        <NA>                     ATMG00030.1,ATMG00030.1-Protein
##                 Index Derives_from
##           <character>  <character>
##       1.2           1         <NA>
##       2.3           1         <NA>
##       3.4        <NA>  AT1G01010.1
##       4.5        <NA>         <NA>
##       5.6        <NA>         <NA>
##       ...         ...          ...
##   438.445        <NA>         <NA>
##   439.446           1         <NA>
##   440.447        <NA>  ATMG00030.1
##   441.448        <NA>         <NA>
##   442.449        <NA>         <NA>
##   -------
##   seqinfo: 7 sequences from an unspecified genome
```

```r
sp[1:4, "type"] # Subsetting of GRangesList objects is similar to GRanges objects.
```

```
## GRangesList object of length 4:
## $1 
## GRanges object with 1 range and 1 metadata column:
##     seqnames       ranges strand |     type
##        <Rle>    <IRanges>  <Rle> | <factor>
##   2     Chr1 [3631, 5899]      * |     mRNA
## 
## $2 
## GRanges object with 1 range and 1 metadata column:
##     seqnames       ranges strand | type
##   3     Chr1 [3631, 5899]      * | mRNA
## 
## $3 
## GRanges object with 1 range and 1 metadata column:
##     seqnames       ranges strand |    type
##   4     Chr1 [3760, 5630]      * | protein
## 
## ...
## <1 more element>
## -------
## seqinfo: 7 sequences from an unspecified genome
```

```r
lapply(sp[1:4], length) # Looping over GRangesList objects similar to lists
```

```
## $`1`
## [1] 1
## 
## $`2`
## [1] 1
## 
## $`3`
## [1] 1
## 
## $`4`
## [1] 1
```

<br><br><center><a href="mydoc_Rsequences_05.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_Rsequences_07.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
