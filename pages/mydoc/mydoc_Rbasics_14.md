---
title: 14. Graphics in R
last_updated: Sun Jun 25 17:41:27 2017
sidebar: mydoc_sidebar
permalink: mydoc_Rbasics_14.html
---

## Advantages

- Powerful environment for visualizing scientific data
- Integrated graphics and statistics infrastructure
- Publication quality graphics
- Fully programmable 
- Highly reproducible
- Full [LaTeX](http://www.latex-project.org/) and Markdown support via `knitr` and `R markdown`
- Vast number of R packages with graphics utilities

## Documentation for R Graphics

__General__

- Graphics Task Page - [URL](http://cran.r-project.org/web/views/Graphics.html)
- R Graph Gallery - [URL](http://addictedtor.free.fr/graphiques/allgraph.php)
- R Graphical Manual - [URL](http://cged.genes.nig.ac.jp/RGM2/index.php)
- Paul Murrell's book R (Grid) Graphics - [URL](http://www.stat.auckland.ac.nz/~paul/RGraphics/rgraphics.html)

__Interactive graphics__

- rggobi` (GGobi) - [URL](http://www.ggobi.org/)
- `iplots` - [URL](http://www.rosuda.org/iplots/)
- Open GL (`rgl`) - [URL](http://rgl.neoscientists.org/gallery.shtml)

## Graphics Environments

__Viewing and saving graphics in R__

- On-screen graphics
- postscript, pdf, svg
- jpeg, png, wmf, tiff, ...

__Four major graphic environments__

(a) Low-level infrastructure

- R Base Graphics (low- and high-level)
- `grid`: [Manual](http://www.stat.auckland.ac.nz/~paul/grid/grid.html)
        
(b) High-level infrastructure
        \begin{itemize}
- `lattice`: [Manual](http://lmdvr.r-forge.r-project.org), [Intro](http://www.his.sunderland.ac.uk/~cs0her/Statistics/UsingLatticeGraphicsInR.htm), [Book](http://www.amazon.com/Lattice-Multivariate-Data-Visualization-Use/dp/0387759689)
- `ggplot2`: [Manual](http://had.co.nz/ggplot2/), [Intro](http://www.ling.upenn.edu/~joseff/rstudy/summer2010_ggplot2_intro.html), [Book](http://had.co.nz/ggplot2/book/)

## Base Graphics: Overview

__Important high-level plotting functions__

- `plot`: generic x-y plotting
- `barplot`: bar plots
- `boxplot`: box-and-whisker plot
- `hist`: histograms
- `pie`: pie charts
- `dotchart`: cleveland dot plots
- `image, heatmap, contour, persp`: functions to generate image-like plots
- `qqnorm, qqline, qqplot`: distribution comparison plots
- `pairs, coplot`: display of multivariant data

__Help on graphics functions__

- `?myfct`
- `?plot`
- `?par`

### Preferred Object Types

- Matrices and data frames
- Vectors
- Named vectors

## Scatter Plots

### Basic Scatter Plot

Sample data set for subsequent plots


```r
set.seed(1410)
y <- matrix(runif(30), ncol=3, dimnames=list(letters[1:10], LETTERS[1:3]))
```

Plot data

```r
plot(y[,1], y[,2]) 
```

<img src="./pages/mydoc/Rbasics_files/basic_scatter_plot-1.png" width="672" />

### All pairs


```r
pairs(y) 
```

<img src="./pages/mydoc/Rbasics_files/pairs_scatter_plot-1.png" width="672" />

### With labels


```r
plot(y[,1], y[,2], pch=20, col="red", main="Symbols and Labels")
text(y[,1]+0.03, y[,2], rownames(y))
```

<img src="./pages/mydoc/Rbasics_files/labels_scatter_plot-1.png" width="672" />

## More examples

__Print instead of symbols the row names__


```r
plot(y[,1], y[,2], type="n", main="Plot of Labels")
text(y[,1], y[,2], rownames(y)) 
```

<img src="./pages/mydoc/Rbasics_files/row_scatter_plot-1.png" width="672" />

__Usage of important plotting parameters__


```r
grid(5, 5, lwd = 2) 
op <- par(mar=c(8,8,8,8), bg="lightblue")
plot(y[,1], y[,2], type="p", col="red", cex.lab=1.2, cex.axis=1.2, 
     cex.main=1.2, cex.sub=1, lwd=4, pch=20, xlab="x label", 
     ylab="y label", main="My Main", sub="My Sub")
par(op)
```
__Important arguments_

- `mar`: specifies the margin sizes around the plotting area in order: `c(bottom, left, top, right)` 
- `col`: color of symbols
- `pch`: type of symbols, samples: `example(points)`
- `lwd`: size of symbols
- `cex.*`: control font sizes
- For details see `?par`


### Add regression line 


```r
plot(y[,1], y[,2])
myline <- lm(y[,2]~y[,1]); abline(myline, lwd=2) 
```

<img src="./pages/mydoc/Rbasics_files/plot_regression-1.png" width="672" />

```r
summary(myline) 
```

```
## 
## Call:
## lm(formula = y[, 2] ~ y[, 1])
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.40357 -0.17912 -0.04299  0.22147  0.46623 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)  
## (Intercept)   0.5764     0.2110   2.732   0.0258 *
## y[, 1]       -0.3647     0.3959  -0.921   0.3839  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.3095 on 8 degrees of freedom
## Multiple R-squared:  0.09589,	Adjusted R-squared:  -0.01712 
## F-statistic: 0.8485 on 1 and 8 DF,  p-value: 0.3839
```

### Log scale

Same plot as above, but on log scale


```r
plot(y[,1], y[,2], log="xy") 
```

<img src="./pages/mydoc/Rbasics_files/plot_regression_log-1.png" width="672" />

### Add a mathematical expression


```r
plot(y[,1], y[,2]); text(y[1,1], y[1,2], expression(sum(frac(1,sqrt(x^2*pi)))), cex=1.3) 
```

<img src="./pages/mydoc/Rbasics_files/plot_regression_math-1.png" width="672" />

## Homework 3B 

Homework 3B: [Scatter Plots](http://girke.bioinformatics.ucr.edu/GEN242/mydoc_homework_03.html)


## Line Plots

### Single data set


```r
plot(y[,1], type="l", lwd=2, col="blue") 
```

<img src="./pages/mydoc/Rbasics_files/plot_line_single-1.png" width="672" />

### Many Data Sets

Plots line graph for all columns in data frame `y`. The `split.screen` function is used in this example in a for loop to overlay several line graphs in the same plot. 


```r
split.screen(c(1,1)) 
```

```
## [1] 1
```

```r
plot(y[,1], ylim=c(0,1), xlab="Measurement", ylab="Intensity", type="l", lwd=2, col=1)
for(i in 2:length(y[1,])) { 
	screen(1, new=FALSE)
	plot(y[,i], ylim=c(0,1), type="l", lwd=2, col=i, xaxt="n", yaxt="n", ylab="", xlab="", main="", bty="n") 
}
```

<img src="./pages/mydoc/Rbasics_files/plot_line_many-1.png" width="672" />

```r
close.screen(all=TRUE) 
```

## Bar Plots 

### Basics


```r
barplot(y[1:4,], ylim=c(0, max(y[1:4,])+0.3), beside=TRUE, legend=letters[1:4]) 
text(labels=round(as.vector(as.matrix(y[1:4,])),2), x=seq(1.5, 13, by=1) + sort(rep(c(0,1,2), 4)), y=as.vector(as.matrix(y[1:4,]))+0.04) 
```

<img src="./pages/mydoc/Rbasics_files/plot_bar_simple-1.png" width="672" />
    
### Error Bars


```r
bar <- barplot(m <- rowMeans(y) * 10, ylim=c(0, 10))
stdev <- sd(t(y))
arrows(bar, m, bar, m + stdev, length=0.15, angle = 90)
```

<img src="./pages/mydoc/Rbasics_files/plot_bar_error-1.png" width="672" />

## Histograms


```r
hist(y, freq=TRUE, breaks=10)
```

<img src="./pages/mydoc/Rbasics_files/plot_hist-1.png" width="672" />

## Density Plots


```r
plot(density(y), col="red")
```

<img src="./pages/mydoc/Rbasics_files/plot_dens-1.png" width="672" />

## Pie Charts


```r
pie(y[,1], col=rainbow(length(y[,1]), start=0.1, end=0.8), clockwise=TRUE)
legend("topright", legend=row.names(y), cex=1.3, bty="n", pch=15, pt.cex=1.8, 
col=rainbow(length(y[,1]), start=0.1, end=0.8), ncol=1) 
```

<img src="./pages/mydoc/Rbasics_files/plot_pie-1.png" width="672" />

## Color Selection Utilities

Default color palette and how to change it


```r
palette()
```

```
## [1] "black"   "red"     "green3"  "blue"    "cyan"    "magenta" "yellow"  "gray"
```

```r
palette(rainbow(5, start=0.1, end=0.2))
palette()
```

```
## [1] "#FF9900" "#FFBF00" "#FFE600" "#F2FF00" "#CCFF00"
```

```r
palette("default")
```

The `gray` function allows to select any type of gray shades by providing values from 0 to 1

```r
gray(seq(0.1, 1, by= 0.2))
```

```
## [1] "#1A1A1A" "#4D4D4D" "#808080" "#B3B3B3" "#E6E6E6"
```

Color gradients with `colorpanel` function from `gplots` library`

```r
library(gplots)
colorpanel(5, "darkblue", "yellow", "white")
```

```
## [1] "#00008B" "#808046" "#FFFF00" "#FFFF80" "#FFFFFF"
```
Much more on colors in R see Earl Glynn's color chart [here](http://research.stowers-institute.org/efg/R/Color/Chart/)


## Saving Graphics to File

After the `pdf()` command all graphs are redirected to file `test.pdf`. Works for all common formats similarly: jpeg, png, ps, tiff, ...

```r
pdf("test.pdf")
plot(1:10, 1:10)
dev.off() 
```

Generates Scalable Vector Graphics (SVG) files that can be edited in vector graphics programs, such as InkScape.


```r
library("RSvgDevice")
devSVG("test.svg")
plot(1:10, 1:10)
dev.off() 
```



## Homework 3C

Homework 3C: [Bar Plots](http://girke.bioinformatics.ucr.edu/GEN242/mydoc_homework_03.html)

<br><br><center><a href="mydoc_Rbasics_13.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_Rbasics_15.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
