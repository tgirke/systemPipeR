\name{module}
\alias{module}
\title{
New interface to use the Environment Modules package from Unix 
}
\description{
Function to be able to use the Environment Modules package (http://modules.sourceforge.net/) in R environment.
The function \code{module} is a tool that allows an easy shell initialization with the pre-defined settings,modifies the environment during the session, loading or unloading a specific application, list all the package initialized in the current session, and list all the application software available in the system. Also, the user can remove all the application for the current session.
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

## List and load to the current session all the software pre-defined in the shell's initialization file 
module("init")

## Example for loading software into the shell environment 
module("load","tophat") 

## Example for remove software from the shell environment
module("unload", "tophat")

## Clear all of the software from the shell’s initialization files
module("clear")
}
}
\keyword{ utilities }
