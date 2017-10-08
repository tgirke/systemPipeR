---
title: 2. Control Structures
last_updated: Mon May  1 17:40:07 2017
sidebar: mydoc_sidebar
permalink: mydoc_Programming_in_R_02.html
---

## Important Operators

### Comparison operators

* `==` (equal)
* `!=` (not equal)
* `>` (greater than)
* `>=` (greater than or equal)
* `<` (less than)
* `<=` (less than or equal)

### Logical operators
		
* `&` (and)
* `|` (or) 
* `!` (not)

## Conditional Executions: `if` Statements

An `if` statement operates on length-one logical vectors.

__Syntax__

```r
if(TRUE) { 
	statements_1 
} else { 
	statements_2 
}
```

__Example__

```r
if(1==0) { 
	print(1) 
} else { 
	print(2) 
}
```

```
## [1] 2
```

## Conditional Executions: `ifelse` Statements

The `ifelse` statement operates on vectors.

__Syntax__

```r
ifelse(test, true_value, false_value)
```
__Example__

```r
x <- 1:10 
ifelse(x<5, x, 0)
```

```
##  [1] 1 2 3 4 0 0 0 0 0 0
```

<br><br><center><a href="mydoc_Programming_in_R_01.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_Programming_in_R_03.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
