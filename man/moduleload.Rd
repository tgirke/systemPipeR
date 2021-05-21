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
(http://modules.sourceforge.net/) from within the R environment.  By default
the user's login shell environment (ie. \code{bash -l}) will be used to
initialize the current session.  The module function can also; load or unload
specific software, list all the loaded software within the current session, and
list all the applications available for loading from the module system.
Lastly, the module function can remove all loaded software from the current
session.
}

\usage{
module(action_type, module_name="")
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
	Name of software to load as character vector.
  }
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
moduleload(modules(args)) 

## Example for removing software from the shell environment
moduleUnload("hisat2")

## Clear all of the software from the shell's initialization files
moduleClear()

## List and load all the software loaded in users default login shell into the current session (default)
moduleInit()
}
}
\keyword{ utilities }
