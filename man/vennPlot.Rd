\name{vennPlot}
\alias{vennPlot}
\title{
Plot 2-5 way Venn diagrams
}
\description{
Ploting function of 2-5 way Venn diagrams from 'VENNset' objects or count set
vectors. A useful feature is the possiblity to combine the counts from several
Venn comparisons with the same number of label sets in a single Venn diagram.
}
\usage{
vennPlot(x, mymain = "Venn Diagram", mysub = "default", setlabels = "default", yoffset = seq(0, 10, by = 0.34), ccol = rep(1, 31), colmode = 1, lcol = c("#FF0000", "#008B00", "#0000FF", "#FF00FF", "#CD8500"), lines = c("#FF0000", "#008B00", "#0000FF", "#FF00FF", "#CD8500"), mylwd = 3, diacol = 1, type = "ellipse", ccex = 1, lcex = 1, sepsplit = "_", ...)
}
\arguments{
  \item{x}{
\code{VENNset} or list of \code{VENNset} objects. Alternatively, a
\code{vector} of Venn counts or a \code{list} of vectors of Venn counts can be
provided as input. If several Venn comparisons are provided in a \code{list} 
then their results are combined in a single Venn diagram, where the count sets
are organized above each other.}
  \item{mymain}{
Main title of plot.
}
  \item{mysub}{
Subtitle of plot. Default \code{mysub="default"} reports the number of unique
items in all sets, as well as the number of unique items in each individual set,
respectively.
}
  \item{setlabels}{
The argument \code{setlabels} allows to provide a vector of custom sample labels.
However, assigning the proper names in the name slots of the initial \code{setlist} 
is preferred for tracking purposes.
}
  \item{yoffset}{
The results from several Venn comparisons can be combined in a single Venn diagram 
by assigning to \code{x} a \code{list} with several \code{VENNsets} or count vectors. 
The positonal offset of the count sets in the plot can be controlled with the 
\code{yoffset} argument. The argument setting \code{colmode} allows to assign 
different colors to each count set. For instance, with \code{colmode=2} one can 
assign to \code{ccol} a color vector or a list, such as \code{ccol=c("blue", "red")}
 or \code{ccol=list(1:8, 8:1)}.
}
  \item{ccol}{
Character or numeric \code{vector} to define colors of count values, e.g. 
\code{ccol=c("black","black","red")}.
}
  \item{colmode}{
See argument \code{yoffset}.
}
  \item{lcol}{
Character or numeric \code{vector} to define colors of set labels, e.g. 
\code{lcol=c("red", "green")}
}
  \item{lines}{
Character or numeric \code{vector} to define colors of lines in plot. 
}
  \item{mylwd}{
Defines line width of shapes used in plot. 
}
  \item{diacol}{
See argument \code{type}.
}
  \item{type}{
Defines shapes used to plot 4-way Venn diagram. Default \code{type="ellipse"} uses
ellipses. The setting \code{type="circle"} returns an incomplete 4-way Venn diagram as 
circles. This representation misses two overlap sectors, but is sometimes easier to 
navigate than the default ellipse version. The missing Venn intersects are reported
below the Venn diagram. Their font color can be controled with the argument \code{diacol}.
}
  \item{ccex}{
Controls font size for count values.
}
  \item{lcex}{
Controls font size for set labels.
}
  \item{sepsplit}{
Character used to separate sample labels in Venn counts.
}
  \item{\dots}{
Additional arguments to pass on.
}
}
\value{
Venn diagram plot.
}
\references{
See examples in 'The Electronic Journal of Combinatorics':
http://www.combinatorics.org/files/Surveys/ds5/VennSymmExamples.html
}
\author{
Thomas Girke
}
\note{
The functions provided here are an extension of the Venn diagram resources on this 
site: http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual#TOC-Venn-Diagrams
}

\seealso{
\code{overLapper}, \code{olBarplot}
}
\examples{
## Sample data
setlist <- list(A=sample(letters, 18), B=sample(letters, 16),
                C=sample(letters, 20), D=sample(letters, 22),
                E=sample(letters, 18), F=sample(letters, 22))

## 2-way Venn diagram
vennset <- overLapper(setlist[1:2], type="vennsets")
vennPlot(vennset)

## 3-way Venn diagram
vennset <- overLapper(setlist[1:3], type="vennsets")
vennPlot(vennset)

## 4-way Venn diagram
vennset <- overLapper(setlist[1:4], type="vennsets")
vennPlot(list(vennset, vennset))

## Pseudo 4-way Venn diagram with circles
vennPlot(vennset, type="circle")

## 5-way Venn diagram
vennset <- overLapper(setlist[1:5], type="vennsets")
vennPlot(vennset)

## Alternative Venn count input to vennPlot (not recommended!)
counts <- sapply(vennlist(vennset), length)
vennPlot(counts)

## 6-way Venn comparison as bar plot
vennset <- overLapper(setlist[1:6], type="vennsets")
olBarplot(vennset, mincount=1)

## Bar plot of standard intersect counts
interset <- overLapper(setlist, type="intersects")
olBarplot(interset, mincount=1)

## Accessor methods for VENNset/INTERSECTset objects
names(vennset)
names(interset)
setlist(vennset)
intersectmatrix(vennset)
complexitylevels(vennset)
vennlist(vennset)
intersectlist(interset)

## Coerce VENNset/INTERSECTset object to list
as.list(vennset)
as.list(interset)

## Pairwise intersect matrix and heatmap
olMA <- sapply(names(setlist), 
		function(x) sapply(names(setlist), 
		function(y) sum(setlist[[x]] \%in\% setlist[[y]])))
olMA
heatmap(olMA, Rowv=NA, Colv=NA)

## Presence-absence matrices for large numbers of sample sets
interset <- overLapper(setlist=setlist, type="intersects", complexity=2)
(paMA <- intersectmatrix(interset))
heatmap(paMA, Rowv=NA, Colv=NA, col=c("white", "gray")) 
}
\keyword{ utilities }
