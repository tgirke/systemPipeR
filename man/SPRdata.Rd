\name{SPRdata}
\alias{SPRdata}
\docType{data}
\title{
SPRdata
}
\description{
Integrates targets files and count table from *systemPipeR* to a SummarizedExperiment object
}
\usage{
  SPRdata(counts = SimpleList(), rowData = NULL, rowRanges = GRangesList(),
                    cmp = FALSE, targetspath = NULL, SEobj = NULL, SEobjName = "default")}
\arguments{
  \item{counts}{
  counts A list or S4Vector::SimpleList() of matrix-like elements, or a matrix-like 
  object (e.g. an ordinary matrix, a data frame, a DataFrame object from the S4Vectors package, etc...).
  All elements of the list must have the same dimensions, and dimension names (if present)
  must be consistent across elements and with the row names of rowRanges and colData.
}
  \item{rowData}{
  A DataFrame object describing the rows. Row names, if present become
  the row names of the SummarizedExperiment object. The number of rows of the DataFrame
  must equal the number of rows of the matrices in assays.
 }
  \item{rowRanges}{
  A GRanges or GRangesList object describing the ranges of interest.
  Names, if present, become the row names of the SummarizedExperiment object. 
  The length of the GRanges or GRangesList must equal the number of rows of the
  matrices in assays. If rowRanges is missing, a SummarizedExperiment instance is returned.
 }
  \item{cmp}{
  character matrix where comparisons are defined in two columns.
  This matrix should be generated with the readComp() function from the targets file.
  Values used for comparisons need to match those in the Factor column of the targets file.
 }
   \item{targetspath}{
  targets data.frame
 }
  \item{SEobj}{
  A SummarizedExperiment object.
 }
  \item{SEobjName}{
  character vector name for the added assays
 }
}
\value{
A SummarizedExperiment object.
}
\author{
Daniela Cassol
}
\examples{
\dontrun{
library(systemPipeR)
# Targets file
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment.char="#")
cmp <- readComp(file=targetspath, format="matrix", delim="-")

## Count table file
countMatrixPath <- system.file("extdata", "countDFeByg.xls", package="systemPipeR")
countMatrix <- read.delim(countMatrixPath, row.names=1)

## Create empty SummarizedExperiment
sprSE <- SPRdata()
class(sprSE)
sprSE
metadata(sprSE)

## Create an object with targets file and comparison
sprSE <- SPRdata(cmp=TRUE, targetspath = targetspath)
sprSE
metadata(sprSE)
colData(sprSE)
assays(sprSE)

## Create an object with targets file and comparison and count table
sprSE <- SPRdata(counts = countMatrix, cmp=TRUE, targetspath = targetspath)
sprSE
metadata(sprSE)
colData(sprSE)
assays(sprSE)
}}
\keyword{ utilities }