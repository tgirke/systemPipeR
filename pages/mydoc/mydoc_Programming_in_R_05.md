---
title: 5. Useful Utilities
last_updated: Mon May  1 17:40:07 2017
sidebar: mydoc_sidebar
permalink: mydoc_Programming_in_R_05.html
---

## Debugging Utilities

Several debugging utilities are available for R. They include:

* `traceback`
* `browser`
* `options(error=recover)`
* `options(error=NULL)`
* `debug`

The [Debugging in R page](http://www.stats.uwo.ca/faculty/murdoch/software/debuggingR/) provides an overview of the available resources.

## Regular Expressions

R's regular expression utilities work similar as in other languages. To learn how to use them in R, one can consult the main help page on this topic with `?regexp`.

### String matching with `grep`

The grep function can be used for finding patterns in strings, here letter `A` in vector `month.name`.

```r
month.name[grep("A", month.name)] 
```

```
## [1] "April"  "August"
```

### String substitution with `gsub`

Example for using regular expressions to substitute a pattern by another one using a back reference. Remember: single escapes `\` need to be double escaped `\\` in R.


```r
gsub('(i.*a)', 'xxx_\\1', "virginica", perl = TRUE) 
```

```
## [1] "vxxx_irginica"
```

## Interpreting a Character String as Expression

Some useful examples

Generates vector of object names in session

```r
mylist <- ls()
mylist[1] 
```

```
## [1] "i"
```

Executes 1st entry as expression


```r
get(mylist[1])
```

```
## [1] 150
```

Alternative approach 

```r
eval(parse(text=mylist[1])) 
```

```
## [1] 150
```

## Replacement, Split and Paste Functions for Strings

__Selected examples__

Substitution with back reference which inserts in this example `_` character

```r
x <- gsub("(a)","\\1_", month.name[1], perl=T) 
x
```

```
## [1] "Ja_nua_ry"
```

Split string on inserted character from above

```r
strsplit(x,"_")
```

```
## [[1]]
## [1] "Ja"  "nua" "ry"
```

Reverse a character string by splitting first all characters into vector fields


```r
paste(rev(unlist(strsplit(x, NULL))), collapse="") 
```

```
## [1] "yr_aun_aJ"
```

## Time, Date and Sleep

__Selected examples__

Return CPU (and other) times that an expression used (here ls)

```r
system.time(ls()) 
```

```
##    user  system elapsed 
##       0       0       0
```

Return the current system date and time

```r
date() 
```

```
## [1] "Sun Apr 16 08:11:08 2017"
```

Pause execution of R expressions for a given number of seconds (e.g. in loop)

```r
Sys.sleep(1) 
```

### Example

#### Import of Specific File Lines with Regular Expression

The following example demonstrates the retrieval of specific lines from an external file with a regular expression. First, an external file is created with the `cat` function, all lines of this file are imported into a vector with `readLines`, the specific elements (lines) are then retieved with the `grep` function, and the resulting lines are split into vector fields with `strsplit`.


```r
cat(month.name, file="zzz.txt", sep="\n")
x <- readLines("zzz.txt")
x[1:6] 
```

```
## [1] "January"  "February" "March"    "April"    "May"      "June"
```

```r
x <- x[c(grep("^J", as.character(x), perl = TRUE))]
t(as.data.frame(strsplit(x, "u")))
```

```
##                 [,1]  [,2] 
## c..Jan....ary.. "Jan" "ary"
## c..J....ne..    "J"   "ne" 
## c..J....ly..    "J"   "ly"
```
## Calling External Software

External command-line software can be called with `system`. The following example calls `blastall` from R

```r
system("blastall -p blastp -i seq.fasta -d uniprot -o seq.blastp")
```

<br><br><center><a href="mydoc_Programming_in_R_04.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_Programming_in_R_06.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
