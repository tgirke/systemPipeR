---
title: 13. SQLite Databases
last_updated: Sun Jun 25 17:41:27 2017
sidebar: mydoc_sidebar
permalink: mydoc_Rbasics_13.html
---

`SQLite` is a lightweight relational database solution. The `RSQLite` package provides an easy to use interface to create, manage and query `SQLite` databases directly from R. Basic instructions
for using `SQLite` from the command-line are available [here](https://www.sqlite.org/cli.html). A short introduction to `RSQLite` is available [here](https://github.com/rstats-db/RSQLite/blob/master/vignettes/RSQLite.Rmd).

## Loading data into SQLite databases

The following loads two `data.frames` derived from the `iris` data set (here `mydf1` and `mydf2`) 
into an SQLite database (here `test.db`).


```r
library(RSQLite)
mydb <- dbConnect(SQLite(), "test.db") # Creates database file test.db
mydf1 <- data.frame(ids=paste0("id", seq_along(iris[,1])), iris)
mydf2 <- mydf1[sample(seq_along(mydf1[,1]), 10),]
dbWriteTable(mydb, "mydf1", mydf1)
```

```
## [1] TRUE
```

```r
dbWriteTable(mydb, "mydf2", mydf2)
```

```
## [1] TRUE
```

## List names of tables in database


```r
dbListTables(mydb)
```

```
## [1] "mydf1" "mydf2"
```

## Import table into `data.frame`


```r
dbGetQuery(mydb, 'SELECT * FROM mydf2')
```

```
##      ids Sepal.Length Sepal.Width Petal.Length Petal.Width    Species
## 1   id99          5.1         2.5          3.0         1.1 versicolor
## 2   id16          5.7         4.4          1.5         0.4     setosa
## 3   id51          7.0         3.2          4.7         1.4 versicolor
## 4  id141          6.7         3.1          5.6         2.4  virginica
## 5  id135          6.1         2.6          5.6         1.4  virginica
## 6  id127          6.2         2.8          4.8         1.8  virginica
## 7   id73          6.3         2.5          4.9         1.5 versicolor
## 8  id148          6.5         3.0          5.2         2.0  virginica
## 9   id13          4.8         3.0          1.4         0.1     setosa
## 10  id40          5.1         3.4          1.5         0.2     setosa
```

## Query database


```r
dbGetQuery(mydb, 'SELECT * FROM mydf1 WHERE "Sepal.Length" < 4.6')
```

```
##    ids Sepal.Length Sepal.Width Petal.Length Petal.Width Species
## 1  id9          4.4         2.9          1.4         0.2  setosa
## 2 id14          4.3         3.0          1.1         0.1  setosa
## 3 id39          4.4         3.0          1.3         0.2  setosa
## 4 id42          4.5         2.3          1.3         0.3  setosa
## 5 id43          4.4         3.2          1.3         0.2  setosa
```

## Join tables

The two tables can be joined on the shared `ids` column as follows. 


```r
dbGetQuery(mydb, 'SELECT * FROM mydf1, mydf2 WHERE mydf1.ids = mydf2.ids')
```

```
##      ids Sepal.Length Sepal.Width Petal.Length Petal.Width    Species   ids Sepal.Length
## 1   id13          4.8         3.0          1.4         0.1     setosa  id13          4.8
## 2   id16          5.7         4.4          1.5         0.4     setosa  id16          5.7
## 3   id40          5.1         3.4          1.5         0.2     setosa  id40          5.1
## 4   id51          7.0         3.2          4.7         1.4 versicolor  id51          7.0
## 5   id73          6.3         2.5          4.9         1.5 versicolor  id73          6.3
## 6   id99          5.1         2.5          3.0         1.1 versicolor  id99          5.1
## 7  id127          6.2         2.8          4.8         1.8  virginica id127          6.2
## 8  id135          6.1         2.6          5.6         1.4  virginica id135          6.1
## 9  id141          6.7         3.1          5.6         2.4  virginica id141          6.7
## 10 id148          6.5         3.0          5.2         2.0  virginica id148          6.5
##    Sepal.Width Petal.Length Petal.Width    Species
## 1          3.0          1.4         0.1     setosa
## 2          4.4          1.5         0.4     setosa
## 3          3.4          1.5         0.2     setosa
## 4          3.2          4.7         1.4 versicolor
## 5          2.5          4.9         1.5 versicolor
## 6          2.5          3.0         1.1 versicolor
## 7          2.8          4.8         1.8  virginica
## 8          2.6          5.6         1.4  virginica
## 9          3.1          5.6         2.4  virginica
## 10         3.0          5.2         2.0  virginica
```


<br><br><center><a href="mydoc_Rbasics_12.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_Rbasics_14.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
