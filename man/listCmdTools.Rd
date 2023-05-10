
\name{listCmdTools}
\alias{listCmdModules}
\title{
List/check the existence of command-line tools of a workflow
}
\description{
These functions list/check whether required command-line tools/modules are installed in the PATH and are callable.
}

\usage{
listCmdTools(sal, check_path = FALSE, check_module = FALSE)

listCmdModules(sal, check_module = FALSE)
}

\arguments{
  \item{sal}{
	SPR workflow object in \code{SYSargsList} class.
  }
  \item{check_path}{
    logical, whether to check if the required tools are in PATH.
  }
  \item{check_module}{
    logical, whether to check if the required modules are installed.
  }
}

\details{
Both functions by default will not check the existence of tools or modules. The default is to list the requirement.
}

\value{
Both functions print out the list/check results as dataframe. The first column is workflow step names that
require certain tools/modules. The second column is the tool/module names. The third column is logical, \code{TRUE} for
the existence of the tool in PATH/modular system, if \code{check_path = TRUE} or \code{check_module = TRUE}. Otherwise,
the third column will be \code{NA}.

In the case of both \code{check_path = TRUE, check_module = TRUE} for \code{listCmdTools}, the returned dataframe is
still results for tool PATH checking but not module checking results. If one wish to obtain the module checking results,
please use \code{listCmdModules}.

When the current workflow has no command-line (SYSargs) step, or there is no module required, or there is no modular
system installed, the return will be \code{NULL}.

These two functions are automatically performed when \code{importWF} is called.
}

\seealso{
\code{importWF}
\code{module}
}

\author{
Le Zhang
}

\examples{
# See examples of `importWF`
}
\keyword{ utilities }
