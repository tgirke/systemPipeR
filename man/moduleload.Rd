
\name{moduleload}
\alias{module}
\alias{moduleload}
\alias{modulelist}
\alias{moduleAvail}
\alias{moduleClear}
\alias{moduleInit}
\alias{moduleUnload}
\title{
Interface to allow full use of the Environment Modules system for Unix
}
\description{
The function \code{module} enables use of the Environment Modules system
(http://modules.sourceforge.net/) from within the R environment. 
The user's login shell environment (i.e. \code{bash -l}) will be used to 
initialize the current session. The module function can also; load or unload 
specific software, list all the loaded software within the current session, 
and list all the applications available for loading from the module system.
Lastly, the module function can remove all loaded software from the current
session.
}

\usage{
module(action_type, module_name = NULL)
moduleload(module_name)
moduleUnload(module_name)
modulelist()
moduleAvail()
moduleClear()
moduleInit()
}

\arguments{
  \item{action_type}{
	Name of the action to be executed as character vector. The following switches are accepted: \code{avail}, \code{list}, \code{init}, \code{load}, \code{unload}, and \code{clear}.
  }
  \item{module_name}{
	Name of software to load as character vector. Examples: \code{"hisat2"}, 
	\code{"hisat2/2.1.0"}, \code{c("hisat2", "samtools")}.
  }
}

\details{
Partial failure would also result `FALSE`, e.g. "load" two modules, one 
successful and the other failed, then the return is `FALSE`.
For "unload" action will always return `TRUE` even if the module is not 
loaded at all or not found.
}

\author{
Tyler Backman, Jordan Hayes and Thomas Girke
}

\examples{
\dontrun{
## List all available software from the module system
avail <- moduleAvail()

## List loaded software in the current session
modulelist()

## Example for loading a software into the shell environment
moduleload("hisat2") 
moduleload("hisat2/2.2.1")

## Example for removing software from the shell environment
moduleUnload("hisat2")

## Clear all of the software from the shell's initialization files
moduleClear()

## List and load all the software loaded in users default login shell into the 
current session (default)
moduleInit()
}
}
\keyword{ utilities }
