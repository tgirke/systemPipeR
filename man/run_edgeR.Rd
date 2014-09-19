\name{run_edgeR}
\alias{run_edgeR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Runs edgeR 
}
\description{
Convenience function to rund \code{edgeR GLM} on many comparisons. Analysis 
can be performed on entire count matrix or only on those columns required for each
comparison. 
}
\usage{
run_edgeR(countDF, targets, cmp, independent = TRUE, paired = NULL, mdsplot = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{countDF}{
\code{date.frame} containing raw read counts
}
  \item{targets}{
	targets \code{data.frame}
}
  \item{cmp}{
	\code{character matrix} where comparisons are defined in two columns. This matrix should be generated with \code{readComp()} from targets file. Values used for comparisons need to match those in \code{Factor} column of targets file.
}
  \item{independent}{
	If \code{independent=TRUE} then \code{countDF} will be subsetted for each comparison
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
\author{
Thomas Girke
}
\seealso{
\code{edgeR} vignette and \code{readComp}
}
\examples{
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment="#")
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
