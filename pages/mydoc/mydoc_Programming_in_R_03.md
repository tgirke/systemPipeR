---
title: 3. Loops
last_updated: Mon May  1 17:40:07 2017
sidebar: mydoc_sidebar
permalink: mydoc_Programming_in_R_03.html
---

## `for` loop

`for` loops iterate over elements of a looping vector.

__Syntax__

```r
for(variable in sequence) { 
	statements 
}
```
__Example__

```r
mydf <- iris
myve <- NULL
for(i in seq(along=mydf[,1])) {
	myve <- c(myve, mean(as.numeric(mydf[i,1:3])))
}
myve[1:8]
```

```
## [1] 3.333333 3.100000 3.066667 3.066667 3.333333 3.666667 3.133333 3.300000
```

__Note:__ Inject into objecs is much faster than append approach with `c`, `cbind`, etc.

__Example__

```r
myve <- numeric(length(mydf[,1]))
for(i in seq(along=myve)) {
	myve[i] <- mean(as.numeric(mydf[i,1:3]))
}
myve[1:8]
```

```
## [1] 3.333333 3.100000 3.066667 3.066667 3.333333 3.666667 3.133333 3.300000
```

### Conditional Stop of Loops

The `stop` function can be used to break out of a loop (or a function) when a condition becomes `TRUE`. In addition, an error message will be printed.

__Example__

```r
x <- 1:10
z <- NULL
for(i in seq(along=x)) { 
	if(x[i] < 5) { 
		z <- c(z, x[i]-1)  
	} else { 
		stop("values need to be < 5") 
	}
}
```

## `while` loop

Iterates as long as a condition is true.

__Syntax__

```r
while(condition) {
	statements
}
```

__Example__

```r
z <- 0
while(z<5) { 
	z <- z + 2
	print(z)  
}
```

```
## [1] 2
## [1] 4
## [1] 6
```

## The `apply` Function Family

### `apply`

__Syntax__

```r
apply(X, MARGIN, FUN, ARGs)
```

__Arguments__

* `X`: `array`, `matrix` or `data.frame`
* `MARGIN`: `1` for rows, `2` for columns
* `FUN`: one or more functions
* `ARGs`: possible arguments for functions

__Example__

```r
apply(iris[1:8,1:3], 1, mean)
```

```
##        1        2        3        4        5        6        7        8 
## 3.333333 3.100000 3.066667 3.066667 3.333333 3.666667 3.133333 3.300000
```

### `tapply`

Applies a function to vector components that are defined by a factor.

__Syntax__

```r
tapply(vector, factor, FUN)
```

__Example__

```r
iris[1:2,]
```

```
##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
## 1          5.1         3.5          1.4         0.2  setosa
## 2          4.9         3.0          1.4         0.2  setosa
```

```r
tapply(iris$Sepal.Length, iris$Species, mean)
```

```
##     setosa versicolor  virginica 
##      5.006      5.936      6.588
```

### `sapply` and `lapply`

Both apply a function to vector or list objects. The `lapply` function always returns a list object, while `sapply` returns `vector` or `matrix` objects when it is possible. 

__Examples__

```r
x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))
lapply(x, mean)
```

```
## $a
## [1] 5.5
## 
## $beta
## [1] 4.535125
## 
## $logic
## [1] 0.5
```

```r
sapply(x, mean)
```

```
##        a     beta    logic 
## 5.500000 4.535125 0.500000
```

Often used in combination with a function definition:

```r
lapply(names(x), function(x) mean(x))
sapply(names(x), function(x) mean(x))
```

<br><br><center><a href="mydoc_Programming_in_R_02.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_Programming_in_R_04.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
