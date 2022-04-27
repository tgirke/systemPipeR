#####################################
## Method Definitions for SYSargs2 ##
#####################################
## Methods to return SYSargs2 components
setMethod(f = "targets", signature = "SYSargs2", definition = function(x) {
    return(x@targets)
})
setMethod(f = "targetsheader", signature = "SYSargs2", definition = function(x) {
    return(x@targetsheader)
})
setMethod(f = "modules", signature = "SYSargs2", definition = function(x) {
    return(setNames(as.character(x@modules), names(x@modules)))
})
setMethod(f = "wf", signature = "SYSargs2", definition = function(x) {
    return(x@wf)
})
setMethod(f = "clt", signature = "SYSargs2", definition = function(x) {
    return(x@clt)
})
setMethod(f = "yamlinput", signature = "SYSargs2", definition = function(x) {
    return(x@yamlinput)
})
setMethod(f = "cmdlist", signature = "SYSargs2", definition = function(x) {
    return(x@cmdlist)
})
setMethod(f = "input", signature = "SYSargs2", definition = function(x) {
    return(x@input)
})
setMethod(f = "output", signature = "SYSargs2", definition = function(x) {
    return(x@output)
})
setMethod(f = "files", signature = "SYSargs2", definition = function(x) {
    return(x@files)
})
setMethod(f = "inputvars", signature = "SYSargs2", definition = function(x) {
    return(x@inputvars)
})
setMethod(f = "cmdToCwl", signature = "SYSargs2", definition = function(x) {
    return(x@cmdToCwl)
})
setMethod(f = "status", signature = "SYSargs2", definition = function(x) {
    return(x@status)
})

## Constructor methods
## List to SYSargs2 with: as(mylist, "SYSargs2")
setAs(
    from = "list", to = "SYSargs2",
    def = function(from) {
        new("SYSargs2",
            targets = from$targets,
            targetsheader = from$targetsheader,
            modules = from$modules,
            wf = from$wf,
            clt = from$clt,
            yamlinput = from$yamlinput,
            cmdlist = from$cmdlist,
            input = from$input,
            output = from$output,
            files = from$files,
            inputvars = from$inputvars,
            cmdToCwl = from$cmdToCwl,
            status = from$status,
            internal_outfiles = from$internal_outfiles
        )
    }
)

setMethod(f = "sysargs2", signature = "SYSargs2", definition = function(x) {
    sysargs2 <- list(
        targets = x@targets, targetsheader = x@targetsheader, modules = x@modules, wf = x@wf,
        clt = x@clt, yamlinput = x@yamlinput, cmdlist = x@cmdlist, input = x@input, output = x@output,
        files = x@files, inputvars = x@inputvars, cmdToCwl = x@cmdToCwl,
        status = x@status, internal_outfiles = x@internal_outfiles
    )
    return(sysargs2)
})

## SYSargs2 to list with: as(SYSargs2, "list")
setAs(from = "SYSargs2", to = "list", def = function(from) {
    sysargs2(from)
})

## Define print behavior for SYSargs2
setMethod(
    f = "show", signature = "SYSargs2",
    definition = function(object) {
        cat(crayon::green$bold(paste0("Instance of '", class(object), "':")),
            paste0("   Slot names/accessors: "),
            paste0(
                "      targets: ", length(object@targets),
                " (", head(names(object@targets), 1), "...",
                tail(names(object@targets), 1), ")",
                ", targetsheader: ", length(unlist(object@targetsheader)), " (lines)"
            ),
            paste0("      modules: ", length(object@modules)),
            paste0(
                "      wf: ", length(object@wf$steps),
                ", clt: ", length(object@clt),
                ", yamlinput: ", length(object@yamlinput), " (inputs)"
            ),
            paste0(
                "      input: ", length(object@input),
                ", output: ", length(object@output)
            ),
            paste0("      cmdlist: ", length(object@cmdlist)),
            "   Sub Steps:",
            paste0(
                "      ", seq_along(object@clt), ". ", object@files$steps,
                " (rendered: ", length(object@cmdlist[[1]]) != 0, ")"
            ),
            "\n",
            sep = "\n"
        )
    }
)

## Extend names() method
setMethod(
    f = "names", signature = "SYSargs2",
    definition = function(x) {
        return(slotNames(x))
    }
)

## Extend infile1() method
setMethod(f = "infile1", signature = "SYSargs2", definition = function(x, input=c("FileName", "FileName1")) {
    subset_input <- input(x)
    input_sub <- input[input %in% names(input(x)[[1]])]
    subset_sample <- sapply(names(subset_input), function(y) subset_input[[y]][[input_sub]])
    subset_sample <- sapply(names(subset_sample), function(y) ifelse(is.null(subset_sample[[y]]), "", subset_sample[y]))
    return(subset_sample)
})

## Extend infile2() method
setMethod(f = "infile2", signature = "SYSargs2", definition = function(x, input="FileName2") {
    subset_input <- input(x)
    input_sub <- input[input %in% names(input(x)[[1]])]
    if(length(input_sub)==0) input_sub <- ""
    subset_sample <- sapply(names(subset_input), function(y) subset_input[[y]][[input_sub]])
    subset_sample <- sapply(names(subset_sample), function(y) ifelse(is.null(subset_sample[[y]]), "", subset_sample[y]))
    return(subset_sample)
})

## Extend length() method
setMethod(
    f = "length", signature = "SYSargs2", definition = function(x) {
        return(length(x@cmdlist))
    }
)

# Behavior of "[" operator for SYSargs2
setMethod(f = "[", signature = "SYSargs2", definition = function(x, i, ..., drop) {
    if (is.logical(i)) {
        i <- which(i)
    }
    x@targets <- x@targets[i]
    x@input <- x@input[i]
    x@output <- x@output[i]
    x@internal_outfiles <- x@internal_outfiles[i]
    x@cmdlist <- x@cmdlist[i]
    return(x)
})

## Behavior of "[[" operator for SYSargs2
setMethod(
    f = "[[", signature = c("SYSargs2", "ANY", "missing"),
    definition = function(x, i, ..., drop) {
        return(as(x, "list")[[i]])
    }
)

## Behavior of "$" operator for SYSargs2
setMethod("$",
          signature = "SYSargs2",
          definition = function(x, name) {
              slot(x, name)
          }
)

## Convert targets data.frame to list
targets.as.list <- function(x, id="SampleName") {
    targetslist <- yaml::yaml.load(yaml::as.yaml(x, column.major = FALSE))
    names(targetslist) <- x[[id]]
    return(targetslist)
}

## Usage:
# targets <- read.delim("targets.txt", comment.char = "#")
# targetslist <- targets.as.list(x=targets)

## Convert targets list to data.frame
targets.as.df <- function(x) {
    targetstmp <- sapply(x, as.character, simplify = FALSE)
    targetsDF <- as.data.frame(do.call("rbind", targetstmp))
    rownames(targetsDF) <- NULL
    colnames(targetsDF) <- names(x[[1]])
    return(targetsDF)
}

## Usage:
# targets.as.df(x=targetslist)

## targets slot from a SYSargs2 obj to df with: as(SYSargs2, "data.frame")
setAs(from = "SYSargs2", to = "DataFrame", def = function(from) {
    S4Vectors::DataFrame(targets.as.df(targets(from)))
})

setMethod("baseCommand", signature = "SYSargs2", definition = function(x) {
    return(x@clt[[1]]$baseCommand[[1]])
})

setMethod("SampleName", signature = "SYSargs2", definition = function(x) {
    targets_x <- targets(x)
    if (length(targets_x) > 0) {
        sample_name_x <- as(x, "DataFrame")
        id <- x[['files']][['id']]
        return(sample_name_x[[id]])
    } else if (length(targets_x) == 0) {
        message("This step doesn't contain multiple samples.")
    }
})

## Replacement method for SYSargs2 using "[" operator
setReplaceMethod(f = "[[", signature = "SYSargs2", definition = function(x, i, j, value) {
    if (i == 1) x@targets <- value
    if (i == 2) x@targetsheader <- value
    if (i == 3) x@modules <- value
    if (i == 4) x@wf <- value
    if (i == 5) x@clt <- value
    if (i == 6) x@yamlinput <- value
    if (i == 7) x@cmdlist <- value
    if (i == 8) x@input <- value
    if (i == 9) x@output <- value
    if (i == 10) x@files <- value
    if (i == 11) x@status <- value
    if (i == 12) x@internal_outfiles <- value
    if (i == "targets") x@targets <- value
    if (i == "targetsheader") x@targetsheader <- value
    if (i == "modules") x@modules <- value
    if (i == "wf") x@wf <- value
    if (i == "clt") x@clt <- value
    if (i == "yamlinput") x@yamlinput <- value
    if (i == "cmdlist") x@cmdlist <- value
    if (i == "input") x@input <- value
    if (i == "output") x@output <- value
    if (i == "files") x@files <- value
    if (i == "cmdToCwl") x@cmdToCwl <- value
    if (i == "status") x@status <- value
    if (i == "internal_outfiles") x@internal_outfiles <- value
    return(x)
})

## Replacement method
setReplaceMethod("yamlinput", c("SYSargs2"), function(x, paramName, value) {
    x <- as(x, "list")
    ## Check paramName
    if (!paramName %in% names(x$yamlinput))
        stop("'paramName' argument must be one of the following: ", "\n",
             paste0(names(x$yamlinput), collapse = ", "))
    ## Check class of value
    if (!identical(class(x$yamlinput[[paramName]]), class(value))) 
        stop("'value' argument must be the same class of the 'paramName': ", "\n",
                    class(x$yamlinput[[paramName]]))
    x$yamlinput[paramName] <- value
    x <- as(x, "SYSargs2")
    x <- updateWF(x)
    x
})

setReplaceMethod("cmdToCwl", c("SYSargs2"), function(x, ..., value) {
    x@cmdToCwl <- value
    x
})
