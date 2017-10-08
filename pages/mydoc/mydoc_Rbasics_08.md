---
title: 8. Important Utilities
last_updated: Sun Jun 25 17:41:27 2017
sidebar: mydoc_sidebar
permalink: mydoc_Rbasics_08.html
---
	
## Combining Objects

The `c` function combines vectors and lists


```r
c(1, 2, 3)
```

```
## [1] 1 2 3
```

```r
x <- 1:3; y <- 101:103
c(x, y)
```

```
## [1]   1   2   3 101 102 103
```

```r
iris$Species[1:8]
```

```
## [1] setosa setosa setosa setosa setosa setosa setosa setosa
## Levels: setosa versicolor virginica
```

The `cbind` and `rbind` functions can be used to append columns and rows, respecively.

```r
ma <- cbind(x, y)
ma
```

```
##      x   y
## [1,] 1 101
## [2,] 2 102
## [3,] 3 103
```

```r
rbind(ma, ma)
```

```
##      x   y
## [1,] 1 101
## [2,] 2 102
## [3,] 3 103
## [4,] 1 101
## [5,] 2 102
## [6,] 3 103
```

## Accessing Dimensions of Objects

Length and dimension information of objects


```r
length(iris$Species)
```

```
## [1] 150
```

```r
dim(iris)
```

```
## [1] 150   5
```

## Accessing Name Slots of Objects

Accessing row and column names of 2D objects

```r
rownames(iris)[1:8]
```

```
## [1] "1" "2" "3" "4" "5" "6" "7" "8"
```

```r
colnames(iris)
```

```
## [1] "Sepal.Length" "Sepal.Width"  "Petal.Length" "Petal.Width"  "Species"
```

Return name field of vectors and lists

```r
names(myVec)
```

```
##  [1] "A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" "N" "O" "P" "Q" "R" "S" "T" "U" "V" "W" "X"
## [25] "Y" "Z"
```

```r
names(myL)
```

```
## [1] "name"        "wife"        "no.children" "child.ages"
```

## Sorting Objects

The function `sort` returns a vector in ascending or descending order

```r
sort(10:1)
```

```
##  [1]  1  2  3  4  5  6  7  8  9 10
```

The function `order` returns a sorting index for sorting an object

```r
sortindex <- order(iris[,1], decreasing = FALSE)
sortindex[1:12]
```

```
##  [1] 14  9 39 43 42  4  7 23 48  3 30 12
```

```r
iris[sortindex,][1:2,]
```

```
##    Sepal.Length Sepal.Width Petal.Length Petal.Width Species
## 14          4.3         3.0          1.1         0.1  setosa
## 9           4.4         2.9          1.4         0.2  setosa
```

```r
sortindex <- order(-iris[,1]) # Same as decreasing=TRUE
```
Sorting multiple columns

```r
iris[order(iris$Sepal.Length, iris$Sepal.Width),][1:2,]
```

```
##    Sepal.Length Sepal.Width Petal.Length Petal.Width Species
## 14          4.3         3.0          1.1         0.1  setosa
## 9           4.4         2.9          1.4         0.2  setosa
```

<br><br><center><a href="mydoc_Rbasics_07.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_Rbasics_09.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
