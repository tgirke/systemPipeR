\name{trimbatch}
\alias{trimbatch}
\title{
Genome read coverage by transcript models
}
\description{
## Trims adaptors hierarchically from longest to shortest match from right end of read. 
## If 'internalmatch=TRUE' then internal matches will trigger the same behavior. The 
## argument minpatternlength defines shortest adaptor match to consider for reads 
## containing only partial adaptors at the right end.

}
\usage{
trimbatch(fq, pattern, internalmatch=FALSE, minpatternlength=8, 
                            Nnumber=1, polyhomo=100, minreadlength=18, 
                            maxreadlength)
}
\arguments{
  \item{fq}{
    \code{character} path to fastq file that contains the target sequences.
}
  \item{pattern}{
    \code{character} pattern used to trim the sequence.
}
  \item{internalmatch}{
    The default is \code{FALSE}. Trims adaptors hierarchically from longest to 
    shortest match from right end of read. If 'internalmatch=TRUE' then internal 
    matches will trigger the same behavior.
}
  \item{minpatternlength}{
    It defines shortest adaptor match to consider for reads containing only 
    partial adaptors at the right end.
}
  \item{Nnumber}{
    A numeric value representing a minimum criterion for the filter. It selects
    elements with fewer than \code{Nnumber} 'N' symbols in each element.
}
  \item{polyhomo}{
    A numeric value representing a maximum criterion for the filter. It selects 
    elements with fewer than threshold copies of any nucleotide.
}
  \item{minreadlength}{
\code{numeric} value representing minimum read length.
}
  \item{maxreadlength}{
  \code{numeric} value representing maximun read length.
}
}
\author{
Thomas Girke
}
\examples{
## Preprocessing of paired-end reads
dir_path <- system.file("extdata/cwl/preprocessReads/trim-pe", package="systemPipeR")
targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
trim <- loadWorkflow(targets=targetspath, wf_file="trim-pe.cwl", input_file="trim-pe.yml", dir_path=dir_path)
trim <- renderWF(trim, inputvars=c(FileName1="_FASTQ_PATH1_", FileName2="_FASTQ_PATH2_", SampleName="_SampleName_"))
trim
\dontrun{
iterTrim <- "trimbatch(fq, pattern='ACACGTCT', internalmatch=FALSE, minpatternlength=6, Nnumber=1, polyhomo=50, minreadlength=16, maxreadlength=101)" 
preprocessReads(args=trim[1], Fct=iterTrim, batchsize=100000, overwrite=TRUE, compress=TRUE)
}
}
\keyword{ utilities }
