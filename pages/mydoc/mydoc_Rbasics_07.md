---
title: 7. Data Objects
last_updated: Sun Jun 25 17:41:27 2017
sidebar: mydoc_sidebar
permalink: mydoc_Rbasics_07.html
---

## Object types

## Vectors (1D)

Definition: `numeric` or `character`


```r
myVec <- 1:10; names(myVec) <- letters[1:10]  
myVec[1:5]
```

```
## a b c d e 
## 1 2 3 4 5
```

```r
myVec[c(2,4,6,8)]
```

```
## b d f h 
## 2 4 6 8
```

```r
myVec[c("b", "d", "f")]
```

```
## b d f 
## 2 4 6
```

## Factors (1D)

Definition: vectors with grouping information


```r
factor(c("dog", "cat", "mouse", "dog", "dog", "cat"))
```

```
## [1] dog   cat   mouse dog   dog   cat  
## Levels: cat dog mouse
```

## Matrices (2D)

Definition: two dimensional structures with data of same type


```r
myMA <- matrix(1:30, 3, 10, byrow = TRUE) 
class(myMA)
```

```
## [1] "matrix"
```

```r
myMA[1:2,]
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
## [1,]    1    2    3    4    5    6    7    8    9    10
## [2,]   11   12   13   14   15   16   17   18   19    20
```

```r
myMA[1, , drop=FALSE]
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
## [1,]    1    2    3    4    5    6    7    8    9    10
```

## Data Frames (2D)

Definition: two dimensional objects with data of variable types


```r
myDF <- data.frame(Col1=1:10, Col2=10:1) 
myDF[1:2, ]
```

```
##   Col1 Col2
## 1    1   10
## 2    2    9
```

## Arrays

Definition: data structure with one, two or more dimensions


## Lists

Definition: containers for any object type


```r
myL <- list(name="Fred", wife="Mary", no.children=3, child.ages=c(4,7,9)) 
myL
```

```
## $name
## [1] "Fred"
## 
## $wife
## [1] "Mary"
## 
## $no.children
## [1] 3
## 
## $child.ages
## [1] 4 7 9
```

```r
myL[[4]][1:2] 
```

```
## [1] 4 7
```

## Functions

Definition: piece of code


```r
myfct <- function(arg1, arg2, ...) { 
	function_body 
}
```

## Subsetting of data objects

__(1.) Subsetting by positive or negative index/position numbers__


```r
myVec <- 1:26; names(myVec) <- LETTERS 
myVec[1:4]
```

```
## A B C D 
## 1 2 3 4
```

__(2.) Subsetting by same length logical vectors__


```r
myLog <- myVec > 10
myVec[myLog] 
```

```
##  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z 
## 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
```

__(3.) Subsetting by field names__


```r
myVec[c("B", "K", "M")]
```

```
##  B  K  M 
##  2 11 13
```

__(4.) Subset with `$` sign__: references a single column or list component by its name 


```r
iris$Species[1:8]
```

```
## [1] setosa setosa setosa setosa setosa setosa setosa setosa
## Levels: setosa versicolor virginica
```

<br><br><center><a href="mydoc_Rbasics_06.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_Rbasics_08.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
