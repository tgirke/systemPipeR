####################################
## Method Definitions for SYSargs ##
####################################
## Methods to return SYSargs components
setMethod(f = "targetsin", signature = "SYSargs", definition = function(x) {
    return(x@targetsin)
})
setMethod(f = "targetsout", signature = "SYSargs", definition = function(x) {
    return(x@targetsout)
})
setMethod(f = "targetsheader", signature = "SYSargs", definition = function(x) {
    return(x@targetsheader)
})
setMethod(f = "modules", signature = "SYSargs", definition = function(x) {
    return(as.character(x@modules))
})
setMethod(f = "software", signature = "SYSargs", definition = function(x) {
    return(as.character(x@software))
})
setMethod(f = "cores", signature = "SYSargs", definition = function(x) {
    return(x@cores)
})
setMethod(f = "other", signature = "SYSargs", definition = function(x) {
    return(x@other)
})
setMethod(f = "reference", signature = "SYSargs", definition = function(x) {
    return(x@reference)
})
setMethod(f = "results", signature = "SYSargs", definition = function(x) {
    return(as.character(x@results))
})
setMethod(f = "infile1", signature = "SYSargs", definition = function(x) {
    return(x@infile1)
})
setMethod(f = "infile2", signature = "SYSargs", definition = function(x) {
    return(x@infile2)
})
setMethod(f = "outfile1", signature = "SYSargs", definition = function(x) {
    return(x@outfile1)
})
setMethod(f = "SampleName", signature = "SYSargs", definition = function(x) {
    return(names(x@sysargs))
})
setMethod(f = "sysargs", signature = "SYSargs", definition = function(x) {
    return(x@sysargs)
})
setMethod(f = "outpaths", signature = "SYSargs", definition = function(x) {
    return(x@outpaths)
})

## Constructor methods
## List to SYSargs with: as(mylist, "SYSargs")
setAs(
    from = "list", to = "SYSargs",
    def = function(from) {
        new("SYSargs",
            targetsin = from$targetsin,
            targetsout = from$targetsout,
            targetsheader = from$targetsheader,
            modules = from$modules,
            software = from$software,
            cores = from$cores,
            other = from$other,
            reference = from$reference,
            results = from$results,
            infile1 = from$infile1,
            infile2 = from$infile2,
            outfile1 = from$outfile1,
            sysargs = from$sysargs,
            outpaths = from$outpaths
        )
    }
)

## Define print behavior for SYSargs
setMethod(
    f = "show", signature = "SYSargs",
    definition = function(object) {
        cat("An instance of '", class(object), "' for running '", object@software, "' on ", length(object@sysargs), " samples ", "\n", sep = "")
    }
)

## Extend names() method
setMethod(
    f = "names", signature = "SYSargs",
    definition = function(x) {
        return(slotNames(x))
    }
)

## Extend length() method
setMethod(
    f = "length", signature = "SYSargs",
    definition = function(x) {
        return(length(x@infile1))
    }
)

## Behavior of "[" operator for SYSargs
setMethod(f = "[", signature = "SYSargs", definition = function(x, i, ..., drop) {
    if (is.logical(i)) {
        i <- which(i)
    }
    x@targetsin <- x@targetsin[i, ]
    x@targetsout <- x@targetsout[i, ]
    x@infile1 <- x@infile1[i]
    x@infile2 <- x@infile2[i]
    x@outfile1 <- x@outfile1[i]
    x@sysargs <- x@sysargs[i]
    x@outpaths <- x@outpaths[i]
    return(x)
})

## Construct SYSargs object from param and targets files
systemArgs <- function(sysma, mytargets, type = "SYSargs") {
    ## Read sysma and convert to arglist; if NULL is assigned to sysma then a dummy version is generated
    sysmapath <- sysma
    if (length(sysmapath) != 0) {
        sysma <- as.matrix(read.delim(sysma, comment.char = "#"))
        sysma[is.na(sysma)] <- ""
    } else {
        sysma <- cbind(
            PairSet = c("software", "cores", "other", "outfile1", "reference", "infile1", "infile1", "infile2", "infile2"),
            Name = c("", "", "", "", "", "", "path", "", "path"),
            Value = c("", "1", "", "<FileName1>", "", "<FileName1>", "", "<FileName2>", "")
        )
    }
    if (any(sysma[, 1] %in% "type")) { # Detects software type: commandline or R
        iscommandline <- sysma[sysma[, 1] %in% "type", , drop = FALSE]
        iscommandline <- as.logical(iscommandline[1, "Value"])
        sysma <- sysma[!sysma[, 1] %in% "type", ] # removes type row(s)
    } else {
        iscommandline <- TRUE # If type line not present then 'commandline' will be assumed
    }
    arglist <- sapply(as.character(unique(sysma[, "PairSet"])), function(x) as.vector(t(as.matrix(sysma[sysma[, "PairSet"] == x, 2:3]))), simplify = FALSE)
    for (i in seq(along = arglist)) names(arglist[[i]]) <- paste(rep(c("n", "v"), length(arglist[[i]]) / 2), rep(1:(length(arglist[[i]]) / 2), 2), sep = "")
    if (type == "json") {
        checkPkg("rjson", quietly = FALSE)
        return(rjson::toJSON(arglist))
    }
    ## Read comment/header lines from targets file
    targetsheader <- readLines(mytargets)
    targetsheader <- targetsheader[grepl("^#", targetsheader)]
    ## Validity checks
    mytargets <- read.delim(mytargets, comment.char = "#")
    mytargetsorig <- mytargets
    if (any(duplicated(mytargets$SampleName))) stop("SampleName column of mytargets cannot contain duplicated entries!")
    ## Preprocessing of targets input
    colnames(mytargets)[1] <- "FileName1" # To support FileName column for SE data
    ## Insert empty FileName2 column if not present
    if (length(mytargets$FileName2) == 0) mytargets <- data.frame(FileName1 = mytargets$FileName1, FileName2 = "", mytargets[, !colnames(mytargets) %in% "FileName1"])
    ## Check name:value violations in arglist
    check <- sapply(names(arglist), function(x) sum(grepl("^n", names(arglist[[x]]))) == sum(grepl("^n", names(arglist[[x]]))))
    if (any(!check)) stop(paste("Name:Value violation in arglist component(s):", paste(names(check[check]), collapse = ", ")))
    ## Modify arglist object as specified in arglist and mytargets
    ## Remove module component and store values in separate container
    modules <- as.character(arglist$modules[grepl("v", names(arglist$modules))])
    arglist <- arglist[!names(arglist) %in% "modules"]
    ## Extract single value components
    software <- as.character(arglist$software[grepl("v", names(arglist$software))])
    other <- as.character(arglist$other[grepl("v", names(arglist$other))])
    if (!(is.na(arglist[["reference"]][["v1"]]) | nchar(arglist[["reference"]][["v1"]]) == 0)) {
        reference <- as.character(arglist$reference[grepl("v", names(arglist$reference))])
        if (!grepl("^/", reference)) reference <- paste0(getwd(), gsub("^\\.", "", reference)) # Turn relative into absolute path.
        arglist[["reference"]]["v1"] <- reference
    } else {
        reference <- ""
    }
    cores <- as.numeric(arglist$cores[grepl("v", names(arglist$cores))])
    ## Populate arglist$infile1
    if (any(grepl("^<.*>$", arglist$infile1))) {
        infile1 <- gsub("<|>", "", arglist$infile1[grepl("^<.*>$", arglist$infile1)][[1]])
        infile1 <- as.character(mytargets[, infile1])
        infile1 <- normalizePath(infile1)
        argname <- arglist$infile1[grep("<.*>", arglist$infile1)[1] - 1]
        path <- arglist$infile1[grep("path", arglist$infile1)[1] + 1]
        infile1back <- paste(path, infile1, sep = "")
        names(infile1back) <- as.character(mytargets$SampleName)
        infile1 <- paste(argname, " ", path, infile1, sep = "")
        arglist[["infile1"]] <- gsub("(^ {1,})|( ${1,})", "", infile1)
    } else {
        infile1back <- rep("", length(mytargets[, 1]))
        infile1 <- infile1back
        names(infile1back) <- as.character(mytargets$SampleName)
        arglist[["infile1"]] <- infile1back
    }
    ## Populate arglist$infile2
    if (any(grepl("^<.*>$", arglist$infile2))) {
        infile2 <- gsub("<|>", "", arglist$infile2[grepl("^<.*>$", arglist$infile2)][[1]])
        infile2 <- as.character(mytargets[, infile2])
        if (nchar(infile2[1]) > 0) infile2 <- normalizePath(infile2)
        argname <- arglist$infile2[grep("<.*>", arglist$infile2)[1] - 1]
        path <- arglist$infile2[grep("path", arglist$infile2)[1] + 1]
        infile2back <- paste(path, infile2, sep = "")
        names(infile2back) <- as.character(mytargets$SampleName)
        infile2 <- paste(argname, " ", path, infile2, sep = "")
        arglist[["infile2"]] <- gsub("(^ {1,})|( ${1,})", "", infile2)
    } else {
        infile2back <- rep("", length(mytargets[, 1]))
        infile2 <- infile2back
        names(infile2back) <- as.character(mytargets$SampleName)
        arglist[["infile2"]] <- infile2back
    }
    ## Populate arglist$outfile1
    outfile1 <- gsub("<|>", "", arglist$outfile1[grepl("^<.*>$", arglist$outfile1)][[1]])
    outfile1 <- as.character(mytargets[, outfile1])
    outfile1 <- gsub("^.*/", "", outfile1)
    remove <- arglist$outfile1[grep("remove", arglist$outfile1)[1] + 1]
    outfile1 <- gsub(as.character(remove), "", outfile1)
    outfile1back <- outfile1
    outpaths <- outfile1
    outextension <- as.character(arglist$outfile1[grep("outextension", arglist$outfile1) + 1])
    append <- arglist$outfile1[grep("append", arglist$outfile1)[1] + 1]
    outfile1 <- paste(outfile1, append, sep = "")
    argname <- arglist$outfile1[grep("<.*>", arglist$outfile1)[1] - 1]
    path <- arglist$outfile1[grep("path", arglist$outfile1)[1] + 1]
    path <- gsub("^\\./|^/|/$", "", path)
    resultpath <- paste(getwd(), "/", path, "/", sep = "")
    outfile1back <- paste(getwd(), "/", path, "/", outfile1, sep = "")
    names(outfile1back) <- as.character(mytargets$SampleName)
    outfile1 <- paste(argname, " ", getwd(), "/", path, "/", outfile1, sep = "")
    arglist[["outfile1"]] <- gsub("(^ {1,})|( ${1,})", "", outfile1)
    ## Populate arglist$outfile2 if it exists (usually only required for PE trimming)
    if ("outfile2" %in% names(arglist)) {
        outfile2 <- gsub("<|>", "", arglist$outfile2[grepl("^<.*>$", arglist$outfile2)][[1]])
        outfile2 <- as.character(mytargets[, outfile2])
        outfile2 <- gsub("^.*/", "", outfile2)
        remove2 <- arglist$outfile2[grep("remove", arglist$outfile2)[1] + 1]
        outfile2 <- gsub(as.character(remove2), "", outfile2)
        outfile2back <- outfile2
        outpaths2 <- outfile2
        outextension2 <- as.character(arglist$outfile2[grep("outextension", arglist$outfile2) + 1])
        append2 <- arglist$outfile2[grep("append", arglist$outfile2)[1] + 1]
        outfile2 <- paste(outfile2, append2, sep = "")
        argname2 <- arglist$outfile2[grep("<.*>", arglist$outfile2)[1] - 1]
        path2 <- arglist$outfile2[grep("path", arglist$outfile2)[1] + 1]
        path2 <- gsub("^\\./|^/|/$", "", path2)
        resultpath2 <- paste(getwd(), "/", path2, "/", sep = "")
        outfile2back <- paste(getwd(), "/", path2, "/", outfile2, sep = "")
        names(outfile2back) <- as.character(mytargets$SampleName)
        outfile2 <- paste(argname2, " ", getwd(), "/", path2, "/", outfile2, sep = "")
        arglist[["outfile2"]] <- gsub("(^ {1,})|( ${1,})", "", outfile2)
    }
    ## Generate arglist$outpaths
    outpaths <- paste(getwd(), "/", path, "/", outpaths, outextension, sep = "")
    names(outpaths) <- as.character(mytargets$SampleName)
    ## Generate targetsout
    targetsout <- mytargetsorig
    targetsout[, 1] <- outpaths[as.character(targetsout$SampleName)]
    if ("outfile2" %in% names(arglist)) {
        outpaths2 <- paste(getwd(), "/", path2, "/", outpaths2, outextension2, sep = "")
        names(outpaths2) <- as.character(mytargets$SampleName)
        targetsout[, 2] <- outpaths2[as.character(targetsout$SampleName)]
        arglist <- arglist[!names(arglist) %in% "outfile2"]
    } else {
        colnames(targetsout)[1] <- "FileName"
        targetsout <- targetsout[, !colnames(targetsout) %in% "FileName2"]
    }
    ## Collapse remaining components to single string vectors
    remaining <- names(arglist)[!names(arglist) %in% c("outfile1", "infile1", "infile2", "outpaths")]
    for (i in remaining) arglist[[i]] <- rep(gsub("(^ {1,})|( ${1,})", "", paste(arglist[[i]], collapse = " ")), length(arglist$infile1))
    args <- do.call("cbind", arglist)
    rownames(args) <- as.character(mytargets$SampleName)
    args <- apply(args, 1, paste, collapse = " ")
    if (software == "bash_commands") { # If command-line is series of bash commands
        args <- gsub("' {1,}| {1,}'", "'", args)
        args <- gsub("bash_commands {1,}", "", args)
    }
    ## When software is R-based then system commands make no sense and NA is used instead
    if (iscommandline == FALSE) args[] <- ""
    ## If sysma=NULL then make adjustments that are most reasonable
    if (length(sysmapath) == 0) {
        targetsout <- mytargetsorig
        modules <- ""
        software <- "R functions"
        other <- ""
        reference <- ""
        resultpath <- ""
        outfile1back <- infile1back
        args[] <- ""
        outpaths[] <- infile1back
    }
    ## Construct SYSargs object from components
    syslist <- list(
        targetsin = mytargetsorig,
        targetsout = targetsout,
        targetsheader = targetsheader,
        modules = modules,
        software = software,
        cores = cores,
        other = other,
        reference = reference,
        results = resultpath,
        infile1 = infile1back,
        infile2 = infile2back,
        outfile1 = outfile1back,
        sysargs = args,
        outpaths = outpaths
    )
    sysargs <- as(syslist, "SYSargs")
    if (type == "SYSargs") {
        return(sysargs)
    }
}

## Usage:
# args <- systemArgs(sysma="../inst/extdata/tophat.param", mytargets="../inst/extdata/targets.txt")
# names(args); modules(args); cores(args); outpaths(args); sysargs(args)
