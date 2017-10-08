---
title: 11. Useful R Functions
last_updated: Sun Jun 25 17:41:27 2017
sidebar: mydoc_sidebar
permalink: mydoc_Rbasics_11.html
---

## Unique entries

Make vector entries unique with `unique`


```r
length(iris$Sepal.Length)
```

```
## [1] 150
```

```r
length(unique(iris$Sepal.Length))
```

```
## [1] 35
```

## Count occurrences

Count occurrences of entries with `table`

```r
table(iris$Species)
```

```
## 
##     setosa versicolor  virginica 
##         50         50         50
```

## Aggregate data

Compute aggregate statistics with `aggregate`

```r
aggregate(iris[,1:4], by=list(iris$Species), FUN=mean, na.rm=TRUE)
```

```
##      Group.1 Sepal.Length Sepal.Width Petal.Length Petal.Width
## 1     setosa        5.006       3.428        1.462       0.246
## 2 versicolor        5.936       2.770        4.260       1.326
## 3  virginica        6.588       2.974        5.552       2.026
```

## Intersect data

Compute intersect between two vectors with `%in%`

```r
month.name %in% c("May", "July")
```

```
##  [1] FALSE FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE
```

## Merge data frames

Join two data frames by common field entries with `merge` (here row names `by.x=0`). To obtain only the common rows, change `all=TRUE` to `all=FALSE`. To merge on specific columns, refer to them by their position numbers or their column names.

```r
frame1 <- iris[sample(1:length(iris[,1]), 30), ]
frame1[1:2,]
```

```
##     Sepal.Length Sepal.Width Petal.Length Petal.Width    Species
## 79           6.0         2.9          4.5         1.5 versicolor
## 124          6.3         2.7          4.9         1.8  virginica
```

```r
dim(frame1)
```

```
## [1] 30  5
```

```r
my_result <- merge(frame1, iris, by.x = 0, by.y = 0, all = TRUE)
dim(my_result)
```

```
## [1] 150  11
```

<br><br><center><a href="mydoc_Rbasics_10.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_Rbasics_12.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
