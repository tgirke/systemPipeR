\name{importWF}
\alias{importWF}
\title{
Import R Markdown file as workflow
}
\description{
Import R Markdown file as workflow. Each R code chunk will be set as a step in the
workflow. This operation requires a few extra settings on the R Markdown chunk
options, to include a particular code chunk in the workflow analysis. Please check
\code{Details}.
}
\usage{
importWF(
    sysargs, file_path, ignore_eval = TRUE,
    update = FALSE, confirm = FALSE, verbose = TRUE
)
}
\arguments{
  \item{sysargs}{
  \code{\link{SYSargsList}} empty object. More information on \code{Details}.
}
  \item{file_path}{
  string, file path of the workflow file.
}
  \item{ignore_eval}{
  logical, treat all R chunks' \code{eval} option as \code{TRUE} in workflow
  \code{Rmd} file even if some chunks have \code{eval=FALSE}.
}
  \item{update}{
  logical, If you have ever changed the template and want to sync new changes,
  turn \code{TRUE} to update the workflow. This function will find the
  difference between old template and this new template, update line number
  records, update preprocessing code, and try to import new steps. The
  updated template is useful in \code{renderReport()} function.
}
  \item{confirm}{
  logical, Only useful when you combine \code{update = TRUE}, some questions
  will be asked during update. Changing this to \code{TRUE} would
  skip these questions and directly say "Yes" to all answers. Useful
  in non-interactive mode. Default is to say "No" in  non-interactive mode.
}
  \item{verbose}{
  logical, print out verbose message while function running.
}
}
\value{
\code{importWF} will return an \code{SYSargsList} object.
}
\details{
To include a particular code chunk from the R Markdown file in the workflow
analysis, please use the following code chunk options:

  - \code{spr = 'r'}: for code chunks with R code lines;
  - \code{spr = 'sysargs'}: for code chunks with an `SYSargsList` object;
  - \code{spr.dep = <StepName>}: for specify the previous dependency. If this options
  is not found, it will automaticly add the previous step.

For \code{spr = 'sysargs'}, the last object assigned needs to be the
\code{\link{SYSargsList}}.
If the \code{spr} flag is not found, the R chunk will not be included in the workflow.

It is required to start a project using \code{SPRproject()} function, and use the
object to populate the steps from R Markdown file.
}
\author{
Le Zhang and Daniela Cassol
}
\seealso{
\code{\link{SYSargsList}}
\code{\link{renderReport}}
}
\examples{
file_path <- system.file("extdata/spr_simple_lw.Rmd", package="systemPipeR")
sal <- SPRproject(overwrite = TRUE)
sal <- importWF(sal, file_path)
}
\keyword{ utilities }
