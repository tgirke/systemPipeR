\name{preprocessReads}
\alias{preprocessReads}
\title{
Run custom read preprocessing functions
}
\description{
Applies custom read preprocessing functions to single-end or paired-end FASTQ 
files. The function uses the \code{FastqStreamer} function from the \code{ShortRead} 
package to stream through large files in a memory-efficient manner. 
}
\usage{
preprocessReads(args = NULL,
                FileName1 = NULL, FileName2 = NULL, 
                outfile1 = NULL, outfile2 = NULL, 
                Fct, batchsize = 100000, overwrite = TRUE, ...)
}
\arguments{
  \item{args}{
Object of class \code{SYSargs} or \code{SYSargs2}.
}
  \item{FileName1}{
Path to input forward fastq file. Defaul is \code{NULL}.
}
  \item{FileName2}{
Path to input reverse fastq file. Defaul is \code{NULL}
}
  \item{outfile1}{
Path to output forward fastq file. Defaul is \code{NULL}
}
  \item{outfile2}{
Path to output reverse fastq file. Defaul is \code{NULL}
}
  \item{Fct}{
\code{character} string of custom read preprocessing function call where both the 
input and output needs to be an object of class \code{ShortReadQ}. The name of the 
input \code{ShortReadQ} object needs to be \code{fq}.
}
  \item{batchsize}{
Number of reads to process in each iteration by the internally used \code{FastqStreamer}
function.
}
  \item{overwrite}{
If \code{TRUE} existing file will be overwritten.
}
  \item{\dots}{
To pass on additional arguments to the internally used \code{writeFastq} function.
}
}
\value{
Writes to files in FASTQ format. Their names are specified by \code{outpaths(args)}.
}
\author{
Thomas Girke
}
\seealso{
FastqStreamer 
}
\examples{
## Preprocessing of single-end reads
dir_path <- system.file("extdata/cwl/preprocessReads/trim-se", package="systemPipeR")
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
trim <- loadWorkflow(targets=targetspath, wf_file="trim-se.cwl", input_file="trim-se.yml", dir_path=dir_path)
trim <- renderWF(trim, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
\dontrun{
preprocessReads(args=trim[1], Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)", batchsize=100000, overwrite=TRUE, compress=TRUE)
}

## Preprocessing of paired-end reads
dir_path <- system.file("extdata/cwl/preprocessReads/trim-pe", package="systemPipeR")
targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
trim <- loadWorkflow(targets=targetspath, wf_file="trim-pe.cwl", input_file="trim-pe.yml", dir_path=dir_path)
trim <- renderWF(trim, inputvars=c(FileName1="_FASTQ_PATH1_", FileName2="_FASTQ_PATH2_", SampleName="_SampleName_"))
trim
\dontrun{
preprocessReads(args=trim[1], Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)", batchsize=100000, overwrite=TRUE, compress=TRUE)
}
}
\keyword{ utilities }
