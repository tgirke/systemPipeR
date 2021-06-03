###############################################
## Class and Method Definitions for LINEWISE ##
###############################################
## Define LineWise class
setClass("LineWise", slots = c(
  codeLine = "expression",
  rmdStart = "integer",
  script = "character"
))
## Methods to return LineWise components
setGeneric(name = "codeLine", def = function(x) standardGeneric("codeLine"))
setMethod(f = "codeLine", signature = "LineWise", definition = function(x) {
  return(cat(as.character(x@codeLine), sep = "\n"))
})

setGeneric(name = "script", def = function(x) standardGeneric("script"))
setMethod(f = "script", signature = "LineWise", definition = function(x) {
  return(x@script)
})

setGeneric(name = "rmdStart", def = function(x) standardGeneric("rmdStart"))
setMethod(f = "rmdStart", signature = "LineWise", definition = function(x) {
  return(x@rmdStart)
})

## Constructor methods
## List to LineWise with: as(mylist, "LineWise")
setAs(from = "list", to = "LineWise",
  def = function(from) {
    new("LineWise",
        codeLine = from$codeLine,
        rmdStart = from$rmdStart,
        script = from$script
    )
})

## Coerce back to list: as(LineWise, "list")
setGeneric(name = "linewise", def = function(x) standardGeneric("linewise"))
setMethod(f = "linewise", signature = "LineWise", definition = function(x) {
  linewise <- list(codeLine = x@codeLine, rmdStart= x@rmdStart, script = x@script)
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
  x@rmdStart <- x@rmdStart[i]
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
  if (i == 2) x@rmdStart <- value
  if (i == 3) x@script <- value
  if (i == "codeLine") x@codeLine <- value
  if (i == "rmdStart") x@rmdStart <- value
  if (i == "script") x@script <- value
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
      y <- as(x, "list")
      y$codeLine <- as.character(y$codeLine)
      y$codeLine[line] <- value
      y$codeLine <- parse(text = y$codeLine)
      x <- as(y, "LineWise")
      x
})

setReplaceMethod("replaceCodeLine", signature = c("SYSargsList"), function(x, step, line, value) {
  y <- stepsWF(x)[step][[1]]
  y <- as(y, "list")
  y$codeLine <- as.character(y$codeLine)
  y$codeLine[line] <- value
  y$codeLine <- parse(text = y$codeLine)
  y <- as(y, "LineWise")
  stepsWF(x, step) <- y
  x
})

setGeneric(name="appendCodeLine<-", def=function(x, ..., after=length(x), value) standardGeneric("appendCodeLine<-"))
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
  y <- stepsWF(x)[step][[1]]
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
  stepsWF(x, step) <- y
  x
})


setMethod(f = "codeLine", signature = "SYSargsList", definition = function(x) {
  rcode <- sapply(names(x$stepsWF), function(x) list(NULL))
  for(i in seq_along(x)){
    if(!inherits(stepsWF(x)[[i]], "LineWise")) stop("Provide a stepWF with a 'LineWise' class")
    code_list <- x$stepsWF[[i]]$codeLine
    rcode[[i]] <- code_list
  }
  for(i in seq_along(rcode)){
    cat(names(rcode[i]), paste0("    ", as.character(rcode[[i]])), sep = "\n")
  }
  #return(rcode)
})



