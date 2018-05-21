\name{module}
\alias{module}
\title{
Interface to allow full use of the Environment Modules system for Unix
}
\description{
The function \code{module} enables use of the Environment Modules system (http://modules.sourceforge.net/) from within the R environment.
By default the user's login shell environment (ie. \dQuote{bash -l}) will be used to initialize the current session.
The module function can also; load or unload specific software, list all the loaded software within the current session, and list all the applications available for loading from the module system.
Lastly, the module function can remove all loaded software from the current session.
}

\usage{
module(action_type, module_name="")
}

\arguments{
  \item{action_type}{
	Name of the action to be executed as character vector. The following switches are accepted: \dQuote{avail}, \dQuote{list}, \dQuote{init}, \dQuote{load}, \dQuote{unload}, and \dQuote{clear}.
  }
  \item{module_name}{
	Name of software to load as character vector.
  }
}

\author{
Jordan Hayes and Thomas Girke
}

\examples{
\dontrun{
## List all available software from the module system
module("avail")

## List loaded software in the current session
module("list")

## Example for loading a software into the shell environment
module("load","tophat") 

## Example for removing software from the shell environment
module("unload", "tophat")

## Clear all of the software from the shell’s initialization files
module("clear")

## List and load all the software loaded in users default login shell into the current session (default)
module("init")
}
}
\keyword{ utilities }
