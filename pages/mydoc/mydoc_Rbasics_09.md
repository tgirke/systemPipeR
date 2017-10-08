---
title: 9. Operators and Calculations
last_updated: Sun Jun 25 17:41:27 2017
sidebar: mydoc_sidebar
permalink: mydoc_Rbasics_09.html
---

## Comparison Operators

Comparison operators: `==`, `!=`, `<`, `>`, `<=`, `>=`

```r
1==1
```

```
## [1] TRUE
```
Logical operators: AND: `&`, OR: `|`, NOT: `!`

```r
x <- 1:10; y <- 10:1
x > y & x > 5
```

```
##  [1] FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE
```

## Basic Calculations

To look up math functions, see Function Index [here](http://cran.at.r-project.org/doc/manuals/R-intro.html#Function-and-variable-index)

```r
x + y
```

```
##  [1] 11 11 11 11 11 11 11 11 11 11
```

```r
sum(x)
```

```
## [1] 55
```

```r
mean(x)
```

```
## [1] 5.5
```

```r
apply(iris[1:6,1:3], 1, mean) 
```

```
##        1        2        3        4        5        6 
## 3.333333 3.100000 3.066667 3.066667 3.333333 3.666667
```

<br><br><center><a href="mydoc_Rbasics_08.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_Rbasics_10.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
