\name{printParam}
\alias{printParam}
\alias{subsetParam}
\alias{replaceParam}
\alias{renameParam}
\alias{appendParam}

\title{
Accessories function to modify the Command-line 
}
\description{
Accessories function to modify the Command-line 
}
\usage{
printParam(sysargs, position, index = NULL)
subsetParam(sysargs, position, index = NULL, trim = TRUE, mute = FALSE)
replaceParam(sysargs, position, index = NULL, replace, mute = FALSE)
renameParam(sysargs, position, index = FALSE, rename, mute = FALSE)
appendParam(sysargs, position, index = NULL, append, after, mute = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sysargs}{
  Object of class \code{SYSargs2}. Output from the \code{createParamFiles} function.
}
  \item{position}{
  string, one of \code{baseCommand}, \code{inputs}, \code{outputs} to view or apply a modification.
}
  \item{index}{
  numeric or character vector, index to view or change a single item in 
  \code{baseCommand}, \code{inputs}, \code{outputs}.
}
  \item{trim}{
  logical, only keep arguments specified by \code{index}. Default is "TRUE".
}
  \item{replace}{
  named list, replace arguments in different positions. Replace list length
  must be the same as \code{index}. Different positions will have different requirements.
}
  \item{rename}{
  character vector, rename items in different positions. \code{rename} vector length
  must be the same as \code{index}.
}
  \item{append}{
  named list, same requirements as \code{replace}, however it cannot append \code{baseCommand}.
}
  \item{after}{
  a subscript, after which the values are to be appended. If \code{NULL} will be 
  after the last argument or specify a numeric integer.
}
  \item{mute}{
  logical, print the raw command-line string and output after replacing or rename.
}
}

\details{
- \code{printParam}: prints its arguments defined by \code{position} and \code{index}.

- \code{subsetParam}: returns subsets of command-line, keeping the arguments defined by \code{position} and \code{index}.

- \code{replaceParam}: replaces the values in command-line with indices given in list by those given in values

- \code{renameParam}: rename the names of the arguments.

- \code{appendParam}: Add arguments to the original command line.
}
\value{
\code{SYSargs2} object
}
\references{
For more details on \code{CWL}, please consult the following page: https://www.commonwl.org/
}
\author{
Le Zhang and Daniela Cassol
}
\seealso{
\code{writeParamFiles}
\code{createParamFiles}
\code{loadWorkflow}
\code{renderWF}
\code{showClass("SYSargs2")}
}

\examples{
command <- "
hisat2 \
    -S <F, out: ./results/M1A.sam> \
    -x <F: ./data/tair10.fasta> \
    -k <int: 1> \
    -min-intronlen <int: 30> \
    -max-intronlen <int: 3000> \
    -threads <int: 4> \
    -U <F: ./data/SRR446027_1.fastq.gz> \
    --verbose
"
cmd <- createParamFiles(command)
cmdlist(cmd)
}

\keyword{ utilities }

