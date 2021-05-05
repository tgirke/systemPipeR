###############################################
## Class and Method Definitions for SYSargs2 ##
###############################################
## Define SYSargs2 class
setClass("SYSargs2", representation(
  targets = "list",
  targetsheader = "list",
  modules = "list",
  wf = "list",
  clt = "list",
  yamlinput = "list",
  cmdlist = "list",
  input = "list",
  output = "list",
  cwlfiles = "list",
  inputvars = "list"
))
## Methods to return SYSargs2 components
setGeneric(name = "targets", def = function(x) standardGeneric("targets"))
setMethod(f = "targets", signature = "SYSargs2", definition = function(x) {
  return(x@targets)
})
setMethod(f = "targetsheader", signature = "SYSargs2", definition = function(x) {
  return(x@targetsheader)
})
setMethod(f = "modules", signature = "SYSargs2", definition = function(x) {
  return(setNames(as.character(x@modules), names(x@modules)))
})
setGeneric(name = "wf", def = function(x) standardGeneric("wf"))
setMethod(f = "wf", signature = "SYSargs2", definition = function(x) {
  return(x@wf)
})
setGeneric(name = "clt", def = function(x) standardGeneric("clt"))
setMethod(f = "clt", signature = "SYSargs2", definition = function(x) {
  return(x@clt)
})
setGeneric(name = "yamlinput", def = function(x) standardGeneric("yamlinput"))
setMethod(f = "yamlinput", signature = "SYSargs2", definition = function(x) {
  return(x@yamlinput)
})
setGeneric(name = "cmdlist", def = function(x) standardGeneric("cmdlist"))
setMethod(f = "cmdlist", signature = "SYSargs2", definition = function(x) {
  return(x@cmdlist)
})
setGeneric(name = "input", def = function(x) standardGeneric("input"))
setMethod(f = "input", signature = "SYSargs2", definition = function(x) {
  return(x@input)
})
setGeneric(name = "output", def = function(x) standardGeneric("output"))
setMethod(f = "output", signature = "SYSargs2", definition = function(x) {
  return(x@output)
})
setGeneric(name = "cwlfiles", def = function(x) standardGeneric("cwlfiles"))
setMethod(f = "cwlfiles", signature = "SYSargs2", definition = function(x) {
  return(x@cwlfiles)
})
setGeneric(name = "inputvars", def = function(x) standardGeneric("inputvars"))
setMethod(f = "inputvars", signature = "SYSargs2", definition = function(x) {
  return(x@inputvars)
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
      cwlfiles = from$cwlfiles,
      inputvars = from$inputvars
    )
  }
)

## Coerce back to list: as(SYSargs2, "list")
setGeneric(name = "sysargs2", def = function(x) standardGeneric("sysargs2"))
setMethod(f = "sysargs2", signature = "SYSargs2", definition = function(x) {
  sysargslist <- list(targets = x@targets, targetsheader = x@targetsheader, modules = x@modules, wf = x@wf, clt = x@clt, yamlinput = x@yamlinput, cmdlist = x@cmdlist, input = x@input, output = x@output, cwlfiles = x@cwlfiles, inputvars = x@inputvars)
  return(sysargslist)
})

## SYSargs2 to list with: as("SYSargs2", list)
setAs(
  from = "SYSargs2", to = "list",
  def = function(from) {
    sysargs2(from)
  }
)

## Define print behavior for SYSargs2
setMethod(
  f = "show", signature = "SYSargs2",
  definition = function(object) {
    cat(paste0("Instance of '", class(object), "':"),
      paste0("   Slot names/accessors: "),
      paste0(
        "      targets: ", length(object@targets),
        " (", head(names(object@targets), 1), "...",
        tail(names(object@targets), 1), ")",
        ", targetsheader: ", length(unlist(object@targetsheader)), " (lines)"
      ),
      paste0("      modules: ", length(object@modules)),
      paste0(
        "      wf: ", length(object@wf),
        ", clt: ", length(object@clt),
        ", yamlinput: ", length(object@yamlinput), " (components)"
      ),
      paste0(
        "      input: ", length(object@input),
        ", output: ", length(object@output)
      ),
      paste0("      cmdlist: ", length(object@cmdlist)),
      "   WF Steps:",
      paste0(
        "      ", seq_along(object@clt), ". ", object@cwlfiles$steps,
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
setMethod(f = "infile1", signature = "SYSargs2", definition = function(x) {
  subset_input <- input(x)
  subset_sample <- sapply(names(subset_input), function(x) list(NULL))
  for (i in seq_along(names(subset_input))) {
    if ("FileName" %in% names(subset_input[[i]])) {
      subset_sample[[i]] <- normalizePath(subset_input[[i]][["FileName"]])
    } else {
      subset_sample[[i]] <- normalizePath(subset_input[[i]][["FileName1"]])
    }
    subset_sample <- as.character(subset_sample)
    names(subset_sample) <- names(subset_input)
  }
  return(subset_sample)
})

## Extend infile2() method
setMethod(f = "infile2", signature = "SYSargs2", definition = function(x) {
  subset_input <- input(x)
  subset_sample <- sapply(names(subset_input), function(x) list(NULL))
  for (i in seq_along(names(subset_input))) {
    if ("FileName2" %in% names(subset_input[[i]])) {
      subset_sample[[i]] <- normalizePath(subset_input[[i]][["FileName2"]])
    }
    subset_sample <- as.character(subset_sample)
    names(subset_sample) <- names(subset_input)
  }
  return(subset_sample)
})

## Extend length() method
setMethod(
  f = "length", signature = "SYSargs2",
  definition = function(x) {
    return(length(x@targets))
  }
)

## Convert targets data.frame to list
targets.as.list <- function(x) {
  targetslist <- yaml::yaml.load(yaml::as.yaml(x, column.major = FALSE))
  names(targetslist) <- x$SampleName
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

# Behavior of "[" operator for SYSargs2
setMethod(f = "[", signature = "SYSargs2", definition = function(x, i, ..., drop) {
  if (is.logical(i)) {
    i <- which(i)
  }
  x@targets <- x@targets[i]
  x@input <- x@input[i]
  x@output <- x@output[i]
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
  if (i == 10) x@cwlfiles <- value
  if (i == "targets") x@targets <- value
  if (i == "targetsheader") x@targetsheader <- value
  if (i == "modules") x@modules <- value
  if (i == "wf") x@wf <- value
  if (i == "clt") x@clt <- value
  if (i == "yamlinput") x@yamlinput <- value
  if (i == "cmdlist") x@cmdlist <- value
  if (i == "input") x@input <- value
  if (i == "output") x@output <- value
  if (i == "cwlfiles") x@cwlfiles <- value
  return(x)
})

####################################################
## Class and Method Definitions for SYSargsList ##
###################################################
## Define SYSargsList class
setClass("SYSargsList", slots = c(
  sysconfig = "list",
  codeSteps = "list",
  stepsWF = "character",
  dataWF = "data.frame",
  SYSargs2_steps = "list",
  statusWF = "list",
  projectWF = "list"
))

## Methods to return SYSargsList components
setGeneric(name = "sysconfig", def = function(x) standardGeneric("sysconfig"))
setMethod(f = "sysconfig", signature = "SYSargsList", definition = function(x) {
  return(x@sysconfig)
})
setGeneric(name = "codeSteps", def = function(x) standardGeneric("codeSteps"))
setMethod(f = "codeSteps", signature = "SYSargsList", definition = function(x) {
  return(x@codeSteps)
})
setGeneric(name = "stepsWF", def = function(x) standardGeneric("stepsWF"))
setMethod(f = "stepsWF", signature = "SYSargsList", definition = function(x) {
  return(names(x@stepsWF))
})
setGeneric(name = "dataWF", def = function(x) standardGeneric("dataWF"))
setMethod(f = "dataWF", signature = "SYSargsList", definition = function(x) {
  return(x@dataWF)
})
setGeneric(name = "SYSargs2_steps", def = function(x) standardGeneric("SYSargs2_steps"))
setMethod(f = "SYSargs2_steps", signature = "SYSargsList", definition = function(x) {
  return(x@SYSargs2_steps)
})
setGeneric(name = "statusWF", def = function(x) standardGeneric("statusWF"))
setMethod(f = "statusWF", signature = "SYSargsList", definition = function(x) {
  return(x@statusWF)
})
setGeneric(name = "projectWF", def = function(x) standardGeneric("projectWF"))
setMethod(f = "projectWF", signature = "SYSargsList", definition = function(x) {
  return(x@projectWF)
})

## Constructor methods
## List to SYSargsList
setAs(
  from = "list", to = "SYSargsList",
  def = function(from) {
    new("SYSargsList",
      sysconfig = from$sysconfig, codeSteps = from$codeSteps, stepsWF = from$stepsWF, dataWF = from$dataWF, SYSargs2_steps = from$SYSargs2_steps,
      statusWF = from$statusWF, projectWF = from$projectWF
    )
  }
)

## Coerce back to list: as(SYSargsList, "list")
setGeneric(name = "sysargslist", def = function(x) standardGeneric("sysargslist"))
setMethod(f = "sysargslist", signature = "SYSargsList", definition = function(x) {
  sysargsset <- list(sysconfig = x@sysconfig, codeSteps = x@codeSteps, stepsWF = x@stepsWF, dataWF = x@dataWF, SYSargs2_steps = x@SYSargs2_steps, statusWF = x@statusWF, projectWF = x@projectWF)
  return(sysargsset)
})

## SYSargsList to list with: as("SYSargsList", list)
setAs(
  from = "SYSargsList", to = "list",
  def = function(from) {
    sysargslist(from)
  }
)

## Define print behavior for SYSargsList
setMethod(
  f = "show", signature = "SYSargsList",
  definition = function(object) {
    cat(paste0("Instance of '", class(object), "':"),
      "   WF Steps:",
      names(object@stepsWF),
      "\n",
      sep = "\n"
    )
  }
)

## Extend names() method
setMethod(
  f = "names", signature = "SYSargsList",
  definition = function(x) {
    return(slotNames(x))
  }
)

## Extend length() method
setMethod(
  f = "length", signature = "SYSargsList",
  definition = function(x) {
    return(length(x@stepsWF))
  }
)

# Behavior of "[" operator for SYSargsList
setMethod(f = "[", signature = "SYSargsList", definition = function(x, i, ..., drop) {
  if (is.logical(i)) {
    i <- which(i)
  }
  x@codeSteps <- x@codeSteps[i]
  return(x)
})

## Behavior of "[[" operator for SYSargsList
setMethod(
  f = "[[", signature = c("SYSargsList","ANY", "missing"),
  definition = function(x, i, ..., drop) {
    return(as(x, "list")[[i]])
  }
)

## Behavior of "$" operator for SYSargsList
setMethod("$",
  signature = "SYSargsList",
  definition = function(x, name) {
    slot(x, name)
  }
)

## Replacement method for SYSargsList using "[" operator
setReplaceMethod(f = "[[", signature = "SYSargsList", definition = function(x, i, j, value) {
  if (i == 1) x@sysconfig <- value
  if (i == 2) x@codeSteps <- value
  if (i == 3) x@stepsWF <- value
  if (i == 4) x@dataWF <- value
  if (i == 5) x@SYSargs2_steps <- value
  if (i == 6) x@statusWF <- value
  if (i == 7) x@projectWF <- value
  if (i == "sysconfig") x@sysconfig <- value
  if (i == "codeSteps") x@codeSteps <- value
  if (i == "stepsWF") x@stepsWF <- value
  if (i == "dataWF") x@dataWF <- value
  if (i == "SYSargs2_steps") x@SYSargs2_steps <- value
  if (i == "statusWF") x@statusWF <- value
  if (i == "projectWF") x@projectWF <- value
  return(x)
})
