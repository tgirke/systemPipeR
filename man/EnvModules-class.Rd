\name{EnvModules-class}
\docType{class}
\alias{EnvModules-class}
\alias{[[,EnvModules,ANY,missing-method}
\alias{[[<-,EnvModules-method}
\alias{[,EnvModules,ANY,ANY,ANY-method}
\alias{[[,EnvModules,ANY,ANY-method}
\alias{[[<-,EnvModules,ANY,ANY,ANY-method}
\alias{$,EnvModules-method}
\alias{available_modules,EnvModules-method}
\alias{coerce,EnvModules,list-method}
\alias{coerce,list,EnvModules-method}
\alias{default_modules,EnvModules-method}
\alias{EnvModules,EnvModules-method}
\alias{loaded_modules,EnvModules-method}
\alias{modulecmd,EnvModules-method}
\alias{names,EnvModules-method}
\alias{show,EnvModules-method}

\title{Class \code{"EnvModules"}}
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
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("EnvModules", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{available_modules}:}{Object of class \code{"list"} ~~ }
    \item{\code{loaded_modules}:}{Object of class \code{"list"} ~~ }
    \item{\code{default_modules}:}{Object of class \code{"list"} ~~ }
    \item{\code{modulecmd}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "EnvModules")}: ... }
    \item{[[}{\code{signature(x = "EnvModules", i = "ANY", j = "missing")}: ... }
    \item{[[<-}{\code{signature(x = "EnvModules")}: ... }
    \item{$}{\code{signature(x = "EnvModules")}: ... }
    \item{available_modules}{\code{signature(x = "EnvModules")}: ... }
    \item{coerce}{\code{signature(from = "EnvModules", to = "list")}: ... }
    \item{coerce}{\code{signature(from = "list", to = "EnvModules")}: ... }
    \item{default_modules}{\code{signature(x = "EnvModules")}: ... }
    \item{EnvModules}{\code{signature(x = "EnvModules")}: ... }
    \item{loaded_modules}{\code{signature(x = "EnvModules")}: ... }
    \item{modulecmd}{\code{signature(x = "EnvModules")}: ... }
    \item{names}{\code{signature(x = "EnvModules")}: ... }
    \item{show}{\code{signature(object = "EnvModules")}: ... }
	 }
}
\author{
Jordan Hayes and Daniela Cassol
}

\examples{
showClass("EnvModules")
}
\keyword{classes}
