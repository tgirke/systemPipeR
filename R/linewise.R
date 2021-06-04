###############################################
## Class and Method Definitions for LINEWISE ##
###############################################
## Define LineWise class
setClass("LineWise", slots = c(
  codeLine = "expression",
  codeChunkStart = "integer",
  rmdPath = "character", 
  stepName="character", 
  dependency = "character"
))
## Methods to return LineWise components
setGeneric(name = "codeLine", def = function(x) standardGeneric("codeLine"))
setMethod(f = "codeLine", signature = "LineWise", definition = function(x) {
  return(cat(as.character(x@codeLine), sep = "\n"))
})

setGeneric(name = "rmdPath", def = function(x) standardGeneric("rmdPath"))
setMethod(f = "rmdPath", signature = "LineWise", definition = function(x) {
  return(x@rmdPath)
})

setGeneric(name = "codeChunkStart", def = function(x) standardGeneric("codeChunkStart"))
setMethod(f = "codeChunkStart", signature = "LineWise", definition = function(x) {
  return(x@codeChunkStart)
})

setGeneric(name = "stepName", def = function(x) standardGeneric("stepName"))
setMethod(f = "stepName", signature = "LineWise", definition = function(x) {
  return(x@stepName)
})

# setGeneric(name = "dependency", def = function(x) standardGeneric("dependency"))
setMethod(f = "dependency", signature = "LineWise", definition = function(x) {
  return(x@dependency)
})


## Constructor methods
## List to LineWise with: as(mylist, "LineWise")
setAs(from = "list", to = "LineWise",
  def = function(from) {
    new("LineWise",
        codeLine = from$codeLine,
        codeChunkStart = from$codeChunkStart,
        rmdPath = from$rmdPath, 
        stepName = from$stepName,
        dependency = from$dependency
    )
})

## Coerce back to list: as(LineWise, "list")
setGeneric(name = "linewise", def = function(x) standardGeneric("linewise"))
setMethod(f = "linewise", signature = "LineWise", definition = function(x) {
  linewise <- list(codeLine = x@codeLine, codeChunkStart= x@codeChunkStart, 
                   rmdPath = x@rmdPath, stepName = x@stepName, dependency = x@dependency)
  return(linewise)
})

## LineWise to list with: as(LineWise, "list")
setAs(from = "LineWise", to = "list",  def = function(from) {
  linewise(from)
})

## Define print behavior for LineWise
setMethod(f = "show", signature = "LineWise",
  definition = function(object) {
    cat(crayon::yellow$bold(paste0("Instance of '", class(object), "'")),
        paste0("    R chunk length: ", length(object)),
        sep = "\n")
})

## Extend names() method
setMethod(f = "names", signature = "LineWise",
  definition = function(x) {
    return(slotNames(x))
})

## Extend length() method
setMethod(
  f = "length", signature = "LineWise", definition = function(x) {
    return(length(x@codeLine))
})


# Behavior of "[" operator for LineWise
setMethod(f = "[", signature = "LineWise", definition = function(x, i, ..., drop) {
  if (is.logical(i)) {
    i <- which(i)
  }
  x@codeLine <- x@codeLine[i]
  x@codeChunkStart <- x@codeChunkStart[i]
  x@stepName <- x@stepName[i]
  x@dependency <- x@dependency[i]
  return(x)
})

## Behavior of "[[" operator for LineWise
setMethod(f = "[[", signature = c("LineWise", "ANY", "missing"),
  definition = function(x, i, ..., drop) {
    return(as(x, "list")[[i]])
})

## Behavior of "$" operator for LineWise
setMethod("$", signature = "LineWise",
          definition = function(x, name) {
            slot(x, name)
})

## Replacement method for LineWise using "[" operator
setReplaceMethod(f = "[[", signature = "LineWise", definition = function(x, i, j, value) {
  if (i == 1) x@codeLine <- value
  if (i == 2) x@codeChunkStart <- value
  if (i == 3) x@rmdPath <- value
  if (i == 4) x@stepName <- value
  if (i == 5) x@dependency <- value
  if (i == "codeLine") x@codeLine <- value
  if (i == "codeChunkStart") x@codeChunkStart <- value
  if (i == "rmdPath") x@rmdPath <- value
  if (i == "stepName") x@stepName <- value
  if (i == "dependency") x@stepName <- value
  return(x)
})

## Replacement methods for codeLine 

# setGeneric(name="codeLine<-", def=function(x, line, value) standardGeneric("codeLine<-"))
# setReplaceMethod("codeLine", c("LineWise"), function(x, value) {
#   x@codeLine <- value
#   x
# })

setGeneric(name="replaceCodeLine<-", def=function(x, line, ..., value) standardGeneric("replaceCodeLine<-"))
setReplaceMethod("replaceCodeLine", signature = c("LineWise"), function(x, line, ..., value) {
  if(!inherits(x, "LineWise")) stop("Provide 'LineWise' class object")
      y <- as(x, "list")
      y$codeLine <- as.character(y$codeLine)
      y$codeLine[line] <- value
      y$codeLine <- parse(text = y$codeLine)
      x <- as(y, "LineWise")
      x
})

setReplaceMethod("replaceCodeLine", signature = c("SYSargsList"), function(x, step, line, value) {
  y <- x$stepsWF[step][[1]]
  if(!inherits(y, "LineWise")) stop("Provide 'LineWise' class object")
  y <- as(y, "list")
  y$codeLine <- as.character(y$codeLine)
  y$codeLine[line] <- value
  y$codeLine <- parse(text = y$codeLine)
  y <- as(y, "LineWise")
  x <- as(x, "list")
  x$stepsWF[step][[1]] <- y
  x <- as(x, "SYSargsList")
  x
})

setGeneric(name="appendCodeLine<-", def=function(x, after=length(x),..., value) standardGeneric("appendCodeLine<-"))
setReplaceMethod("appendCodeLine", signature = c("LineWise"), function(x, after=length(x), value) {
  lengx <- length(x)
  x <- linewise(x)
  value <- parse(text=value)
  if (!after) {
    x$codeLine <- c(value, x$codeLine)
  } else if (after >= lengx) {
    x$codeLine <- c(x$codeLine, value)
  } else {
    x$codeLine <- c(x$codeLine[1L:after], value, x$codeLine[(after + 1L):lengx])
  }
  x <- as(x, "LineWise")
  x
})

setReplaceMethod("appendCodeLine", signature = c("SYSargsList"), function(x, step, after=NULL, value) {
  if(is.null(after)) after <- length(stepsWF(x[step])[[1]])
  y <- x$stepsWF[step][[1]]
  if(!inherits(y, "LineWise")) stop("Provide 'LineWise' class object")
  lengx <- length(y)
  y <- linewise(y)
  value <- parse(text=value)
  if (!after) {
    y$codeLine <- c(value, y$codeLine)
  } else if (after >= lengx) {
    y$codeLine <- c(y$codeLine, value)
  } else {
    y$codeLine <- c(y$codeLine[1L:after], value, y$codeLine[(after + 1L):lengx])
  }
  y <- as(y, "LineWise")
  x <- as(x, "list")
  x$stepsWF[step][[1]] <- y
  x <- as(x, "SYSargsList")
  x
})


setMethod(f = "codeLine", signature = "SYSargsList", definition = function(x) {
  rcode <- sapply(names(x$stepsWF), function(x) list(NULL))
  for(i in seq_along(x)){
    if(!inherits(stepsWF(x)[[i]], "LineWise")) stop("This step is 'SYSargs2'. Please provide an 'LineWise' class step")
    code_list <- x$stepsWF[[i]]$codeLine
    rcode[[i]] <- code_list
  }
  for(i in seq_along(rcode)){
    cat(crayon::blue(names(rcode[i])), paste0("    ", as.character(rcode[[i]])), sep = "\n")
  }
  #return(rcode)
})



