#####################################
## Method Definitions for LineWise ##
#####################################
## Methods to return LineWise components
setMethod(f = "codeLine", signature = "LineWise", definition = function(x) {
    return(cat(as.character(x@codeLine), sep = "\n"))
})
setMethod(f = "rmdPath", signature = "LineWise", definition = function(x) {
    return(x@files$rmdPath)
})
setMethod(f = "codeChunkStart", signature = "LineWise", definition = function(x) {
    return(x@codeChunkStart)
})
setMethod(f = "stepName", signature = "LineWise", definition = function(x) {
    return(x@stepName)
})
setMethod(f = "dependency", signature = "LineWise", definition = function(x) {
    return(x@dependency)
})
setMethod(f = "status", signature = "LineWise", definition = function(x) {
    return(x@status)
})
setMethod(f = "files", signature = "LineWise", definition = function(x) {
    return(x@files)
})
setMethod(f = "runInfo", signature = "LineWise", definition = function(x) {
    return(x@runInfo)
})

#######################################
## Constructor Methods for LineWise ##
#######################################
## Constructor methods
## List to LineWise with: as(mylist, "LineWise")
setAs(
    from = "list", to = "LineWise",
    def = function(from) {
        new("LineWise",
            codeLine = from$codeLine,
            codeChunkStart = from$codeChunkStart,
            stepName = from$stepName,
            dependency = from$dependency,
            status = from$status,
            files = from$files,
            runInfo = from$runInfo
        )
    }
)

## Coerce back to list: as(LineWise, "list")
setMethod(f = "linewise", signature = "LineWise", definition = function(x) {
    linewise <- list(
        codeLine = x@codeLine, codeChunkStart = x@codeChunkStart,
        # rmdPath = x@rmdPath,
        stepName = x@stepName, dependency = x@dependency,
        status = x@status, files = x@files, runInfo = x@runInfo
    )
    return(linewise)
})

## LineWise to list with: as(LineWise, "list")
setAs(from = "LineWise", to = "list", def = function(from) {
    linewise(from)
})

## Define print behavior for LineWise
setMethod(
    f = "show", signature = "LineWise",
    definition = function(object) {
        cat(crayon::yellow$bold(paste0("Instance of '", class(object), "'")),
            paste0("    Code Chunk length: ", length(object)),
            sep = "\n"
        )
    }
)

## Extend names() method
setMethod(
    f = "names", signature = "LineWise",
    definition = function(x) {
        return(slotNames(x))
    }
)

## Extend length() method
setMethod(f = "length", signature = "LineWise", definition = function(x) {
    return(length(x@codeLine))
})


# Behavior of "[" operator for LineWise
setMethod(f = "[", signature = "LineWise", definition = function(x, i, ..., drop) {
    if (is.logical(i)) {
        i <- which(i)
    }
    x@codeLine <- x@codeLine[i]
    x@stepName <- x@stepName
    x@dependency <- x@dependency
    x@status <- x@status
    x@files <- x@files
    x@runInfo <- x@runInfo
    return(x)
})

## Behavior of "[[" operator for LineWise
setMethod(
    f = "[[", signature = c("LineWise", "ANY", "missing"),
    definition = function(x, i, ..., drop) {
        return(as(x, "list")[[i]])
    }
)

## Behavior of "$" operator for LineWise
setMethod("$",
    signature = "LineWise",
    definition = function(x, name) {
        slot(x, name)
    }
)

###########################################
## Replacement Definitions for LineWise ##
##########################################

## Replacement method for LineWise using "[" operator
setReplaceMethod(f = "[[", signature = "LineWise", definition = function(x, i, j, value) {
    if (i == 1) x@codeLine <- value
    if (i == 2) x@codeChunkStart <- value
    if (i == 4) x@stepName <- value
    if (i == 5) x@dependency <- value
    if (i == 6) x@status <- value
    if (i == 7) x@files <- value
    if (i == 8) x@runInfo <- value
    if (i == "codeLine") x@codeLine <- value
    if (i == "codeChunkStart") x@codeChunkStart <- value
    if (i == "stepName") x@stepName <- value
    if (i == "dependency") x@dependency <- value
    if (i == "status") x@status <- value
    if (i == "files") x@files <- value
    if (i == "runInfo") x@runInfo <- value
    return(x)
})

setReplaceMethod("replaceCodeLine", signature = c("LineWise"), function(x, line, ..., value) {
    if (!inherits(x, "LineWise")) stop("Provide 'LineWise' class object")
    y <- as(x, "list")
    y$codeLine <- as.character(y$codeLine)
    y$codeLine[line] <- value
    y$codeLine <- parse(text = y$codeLine)
    x <- as(y, "LineWise")
    x
})

setReplaceMethod("appendCodeLine", signature = c("LineWise"), function(x, after = length(x), value) {
    lengx <- length(x)
    x <- linewise(x)
    value <- parse(text = value)
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
