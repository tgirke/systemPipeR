\name{run_DESeq2}
\alias{run_DESeq2}
\title{
Runs DESeq2 
}
\description{
Convenience function to rund \code{DESeq2} on many comparisons. Analysis 
can be performed on entire count matrix or only on those columns required for each
comparison.  
}
\usage{
run_DESeq2(countDF, targets, cmp, independent = FALSE)
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
}
\value{
\code{data.frame} containing \code{DESeq2} results from all comparisons. Comparison labels are appended to column titles for tracking.
}
\author{
Thomas Girke
}
\seealso{
\code{run_edgeR}, \code{readComp} and \code{DESeq2} vignette
}
\examples{
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment="#")
cmp <- readComp(file=targetspath, format="matrix", delim="-")
countfile <- system.file("extdata", "countDFeByg.xls", package="systemPipeR")
countDF <- read.delim(countfile, row.names=1)
degseqDF <- run_DESeq2(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=FALSE)
pval <- degseqDF[, grep("_FDR$", colnames(degseqDF)), drop=FALSE]
fold <- degseqDF[, grep("_logFC$", colnames(degseqDF)), drop=FALSE]
DEG_list <- filterDEGs(degDF=degseqDF, filter=c(Fold=2, FDR=10))
names(DEG_list)
DEG_list$Summary
}
\keyword{ utilities }
