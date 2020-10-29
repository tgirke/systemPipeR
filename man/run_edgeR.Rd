\name{run_edgeR}
\alias{run_edgeR}
\title{
Runs edgeR 
}
\description{
Convenience wrapper function to identify differentially expressed genes
(DEGs) in batch mode with the \code{edgeR} GML method for any number of pairwise
sample comparisons specified under the \code{cmp} argument.  Users are strongly
encouraged to consult the \code{edgeR} vignette for more detailed information
on this topic and how to properly run \code{edgeR} on data sets with more
complex experimental designs.
}
\usage{
run_edgeR(countDF, targets, cmp, independent = TRUE, paired = NULL, mdsplot = "")
}
\arguments{
  \item{countDF}{
\code{date.frame} containing raw read counts
}
  \item{targets}{
	targets \code{data.frame}
}
  \item{cmp}{
	\code{character matrix} where comparisons are defined in two columns. This matrix should be generated with \code{readComp()} from the targets file. Values used for comparisons need to match those in the \code{Factor} column of the targets file.
}
  \item{independent}{
	If \code{independent=TRUE} then the \code{countDF} will be subsetted for each comparison. This behavior can be useful when working with samples from unrelated studies. For samples from the same or comparable studies, the setting \code{independent=FALSE} is usually preferred.
}
  \item{paired}{
	Defines pairs (\code{character vector}) for paired analysis. Default is unpaired (\code{paired=NULL}).
}
  \item{mdsplot}{
	Directory where \code{plotMDS} should be written to. Default setting \code{mdsplot=""} will omit the plotting step.
}
}
\value{
\code{data.frame} containing \code{edgeR} results from all comparisons. Comparison labels are appended to column titles for tracking.
}
\references{
Please properly cite the \code{edgeR} papers when using this function:
http://www.bioconductor.org/packages/devel/bioc/html/edgeR.html
}
\author{
Thomas Girke
}
\seealso{
\code{run_DESeq2}, \code{readComp} and \code{edgeR} vignette
}
\examples{
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")
cmp <- readComp(file=targetspath, format="matrix", delim="-")
countfile <- system.file("extdata", "countDFeByg.xls", package="systemPipeR")
countDF <- read.delim(countfile, row.names=1)
edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=FALSE, mdsplot="")
pval <- edgeDF[, grep("_FDR$", colnames(edgeDF)), drop=FALSE]
fold <- edgeDF[, grep("_logFC$", colnames(edgeDF)), drop=FALSE]
DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=10))
names(DEG_list)
DEG_list$Summary
}
\keyword{ utilities }
