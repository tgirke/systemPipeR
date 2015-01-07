\name{overLapper}
\alias{overLapper}
\alias{setlist}
\alias{setlist-methods}
\alias{intersectmatrix}
\alias{intersectmatrix-methods}
\alias{complexitylevels}
\alias{complexitylevels-methods}
\alias{vennlist}
\alias{vennlist-methods}
\alias{intersectlist}
\alias{intersectlist-methods}
\title{
Set Intersect and Venn Diagram Functions
}
\description{
Function for computing Venn intersects or standard intersects among large
numbers of label sets provided as \code{list} of \code{vectors}. The resulting
intersect objects can be used for plotting 2-5 way Venn diagrams or intersect
bar plots using the functions \code{vennPlot} or \code{olBarplot}, respetively.
The \code{overLapper} function scales to 2-20 or more label vectors for Venn
intersect calculations and to much larger sample numbers for standard
intersects. The different intersect types are explained below under the
definition of the \code{type} argument. The upper Venn limit around 20 label
sets is unavoidable because the complexity of Venn intersects increases
exponentially with the label set number \code{n} according to this
relationship: \code{2^n - 1}. The overall workflow of the
\code{overLapper} environment is to first compute for a list of label sets
their Venn intersects using the \code{overLapper} function, which organizes the
result sets in an object of class \code{VENNset}. Subsequently, the Venn counts
are computed and plotted as bar or Venn diagrams. The current implementation of
the plotting function \code{vennPlot} supports Venn diagrams for 2-5 label
sets. To visually analyze larger numbers of label sets, a variety of intersect
methods are introduced in the \code{olBarplot} help file. These methods are
much more scalable than Venn diagrams, but lack their restrictive intersect
logic. 
}

\usage{
overLapper(setlist, complexity = "default", sep = "_", cleanup = FALSE, keepdups = FALSE, type)
}

\arguments{
  \item{setlist}{
Object of class \code{list} where each list component stores a label set as
\code{vector} and the name of each label set is stored in the name slot of 
each list component. The names are used for naming the label sets in all 
downstream analysis steps and plots.
}
  \item{complexity}{
Complexity level of intersects specified as integer vector. For Venn intersects 
it needs to be assigned \code{1:length(setlist)} (default). If \code{complexity=2}
the function returns all pairwise intersects.
}
  \item{sep}{
Character used to separate set labels. 
}
  \item{cleanup}{
If set to \code{TRUE} then all characters of the label sets are set to upper case,
and leading and trailing spaces are removed. The default \code{cleanup=FALSE} omits 
this step.
}
  \item{keepdups}{
By default all duplicates are removed from the label sets. The setting
\code{keepdups=TRUE} will retain duplicates by appending a counter to each entry.
}
  \item{type}{
With the default setting \code{type="vennsets"} the overLapper function
computes the typical Venn intersects for the label sets provided under
\code{setlist}. With the setting \code{type="intersects"} the function will
compute pairwise intersects (not compatible with Venn diagrams). Venn
intersects follow the typical 'only in' intersect logic of Venn comparisons, such
as: labels present only in set A, labels present only in the intersect of A
& B, etc. Due to this restrictive intersect logic, the combined Venn sets
contain no duplicates. In contrast to this, regular intersects follow this
logic: labels present in the intersect of A & B, labels present in the
intersect of A & B & C, etc. This approach results usually in many duplications
of labels among the intersect sets.  
}
}
\details{
Additional Venn diagram resources are provided by the packages \code{limma}, 
\code{gplots}, \code{vennerable}, \code{eVenn} and \code{VennDiagram}, or
online resources such as shapes, Venn Diagram Generator and Venny.
}
\value{
\code{overLapper} returns standard intersect and Venn intersect results as
\code{INTERSECTset} or \code{VENNset} objects, respectively. These S4 objects
contain the following components:
	\item{setlist}{Original label sets accessible with \code{setlist()}.}
	\item{intersectmatrix}{Present-absent matrix accessible with \code{intersectmatrix()}, 
where each overlap set in the \code{vennlist} data component is labeled according to 
the label set names provided under \code{setlist}. For instance, the composite name 
'ABC' indicates that the entries are restricted to A, B and C. The seperator used 
for naming the intersect sets can be specified under the \code{sep} argument.} 
	\item{complexitylevels}{Complexity levels accessible with \code{complexitylevels()}.}
	\item{vennlist}{Venn intersects for \code{VENNset} objects accessible with \code{vennlist()}.}
	\item{intersectlist}{Standard intersects for \code{INTERSECTset} objects accessible with \code{intersectlist()}.}
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
\code{vennPlot}, \code{olBarplot}
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
