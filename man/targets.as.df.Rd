\name{targets.as.df}
\alias{targets.as.df}
\title{
Convert targets list to data.frame
}
\description{
Convert \code{list}, which stores data from each target \code{input} file to \code{data.frame} object. 
}
\usage{
targets.as.df(x)
}
\arguments{
  \item{x}{
  An object of the class \code{"list"} that stores data from each target \code{input} file, as \code{targets(WF)}.
  }
}
\value{
\code{data.frame} containing all the \code{input} file information. 
}
\author{
Daniela Cassol and Thomas Girke
}
\examples{
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl", package="systemPipeR")
WF <- loadWorkflow(targets=targets, wf_file="hisat2/hisat2-mapping-se.cwl", 
                  input_file="hisat2/hisat2-mapping-se.yml", dir_path=dir_path)
WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
WF
targets.as.df(targets(WF))
}
\keyword{ utilities }
