########################################
## Method Definitions for SYSargsList ##
########################################
## Methods to return SYSargsList components
setMethod(f = "stepsWF", signature = "SYSargsList", definition = function(x) {
    return(x@stepsWF)
})
setMethod(f = "statusWF", signature = "SYSargsList", definition = function(x) {
    return(lapply(lapply(x@statusWF, "[[", 2), function(y) S4Vectors::DataFrame(y, check.names = FALSE)))
})
setMethod(f = "targetsWF", signature = "SYSargsList", definition = function(x) {
    return(x@targetsWF)
})
setMethod(f = "outfiles", signature = "SYSargsList", definition = function(x) {
    return(x@outfiles)
})
setMethod(f = "SEobj", signature = "SYSargsList", definition = function(x) {
    return(x@SEobj)
})
setMethod(f = "dependency", signature = "SYSargsList", definition = function(x) {
    return(x@dependency)
})
setMethod(f = "projectInfo", signature = "SYSargsList", definition = function(x) {
    return(x@projectInfo)
})
setMethod(f = "runInfo", signature = "SYSargsList", definition = function(x) {
    return(x@runInfo)
})

## Coerce back to list: as(SYSargsList, "list")
setMethod(f = "sysargslist", signature = "SYSargsList", definition = function(x) {
    sysargslist <- list(
        stepsWF = x@stepsWF,
        statusWF = x@statusWF,
        targetsWF = x@targetsWF,
        outfiles = x@outfiles,
        SEobj = x@SEobj,
        dependency = x@dependency,
        targets_connection = x@targets_connection,
        projectInfo = x@projectInfo,
        runInfo = x@runInfo
    )
    return(sysargslist)
})

## Constructor methods
## List to SYSargsList
setAs(
    from = "list", to = "SYSargsList",
    def = function(from) {
        new("SYSargsList",
            stepsWF = from$stepsWF,
            statusWF = from$statusWF,
            targetsWF = from$targetsWF,
            outfiles = from$outfiles,
            SEobj = from$SEobj,
            dependency = from$dependency,
            targets_connection = from$targets_connection,
            projectInfo = from$projectInfo,
            runInfo = from$runInfo
        )
    }
)

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
        if (length(object) > 0) {
            status_1 <- .showHelper(object)
            status <- append("   WF Steps:\n", status_1)
        } else {
            status <- "No workflow steps added"
        }
        cat(
            crayon::blue$bold(paste0("Instance of '", class(object), "':")), "\n",
            status, "\n"
        )
    }
)

.showHelper <- function(object) {
    status_1 <- as.character()
    for (i in seq_along(object@stepsWF)) {
        status_color <- switch(tolower(object@statusWF[[i]]$status.summary),
            "pending" = crayon::blue$bold,
            "warning" = crayon::make_style("orange")$bold,
            "error" = crayon::red$bold,
            "success" = crayon::green$bold
        )
        if (inherits(object@stepsWF[[i]], "SYSargs2")) {
            status_1 <- c(
                status_1,
                c(
                    paste0("      ", i, ". ", crayon::blue(names(object@stepsWF)[i], "--> Status: "), status_color(object@statusWF[[i]]$status.summary), sep = ""), "\n",
                    paste0(c("          Total Files: ", "Existing: ", "Missing: "), colSums(object@statusWF[[i]][[2]][2:4]), collapse = " | "), "\n",
                    paste0(
                        # "         Sub Steps:", "\n",
                        paste0("        ", i, ".", seq_along(object@stepsWF[[i]]@clt), ". ", crayon::green(object@stepsWF[[i]]$files$steps)), "\n",
                        paste0("             cmdlist: ", length(object@stepsWF[[i]]), " | "),
                        sapply(as.list(object@stepsWF[[i]]$files$steps), function(x) {
                            paste0(
                                paste0(names(table(unlist(object@statusWF[[i]][[2]][object@stepsWF[[i]]$files$steps][x]))), ": ",
                                    table(unlist(object@statusWF[[i]][[2]][object@stepsWF[[i]]$files$steps][x])),
                                    collapse = " | "
                                )
                            )
                        }),
                        "\n"
                    )
                )
            )
        } else if (inherits(object@stepsWF[[i]], "LineWise")) {
            status_1 <- c(
                status_1,
                paste0("      ", i, ". ", crayon::blue(names(object@stepsWF)[i], "--> Status: "), status_color(object@statusWF[[i]]$status.summary), "\n")
            )
        }
    }
    return(status_1)
}

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
    if (missing(i)) {
        i <- 1:length(x)
    }
    if (inherits(i, "character")) {
        if (!all(i %in% stepName(x))) {
            stop(paste0(
                "\n",
                "Step name doesn't exist. Please subset accordingly with the 'stepName(x)'",
                "\n",
                paste0(stepName(x), collapse = ", ")
            ))
        }
        i <- which(stepName(x) %in% i)
    }
    if (is.logical(i)) {
        i <- which(i)
    }
    for (s in seq_along(i)) {
        tryCatch(
            {
                if (i[s] < 0) {
                    ii <- i[s] * -1
                } else {
                    ii <- i[s]
                }
                x@stepsWF[[ii]]
            }, error = function(e) {
                e$message <- paste0(
                    "\n",
                    "Step number is out of range. Please subset accordingly with the 'length(x)'",
                    "\n",
                    paste0(1:length(x), collapse = ", ")
                )
                stop(e)
            }
        )
    }
    x@stepsWF <- x@stepsWF[i]
    names_tc <- names(x@stepsWF)
    x@statusWF <- x@statusWF[i]
    x@targetsWF <- x@targetsWF[i]
    x@outfiles <- x@outfiles[i]
    x@SEobj <- x@SEobj[i]
    x@dependency <- x@dependency[i]
    x@targets_connection <- x@targets_connection[names(x@targets_connection) %in% names_tc]
    x@projectInfo <- x@projectInfo
    x@runInfo$directory <- x@runInfo$directory[i]
    # x <- .check_write_SYSargsList(x)
    return(x)
})

## Behavior of "[[" operator for SYSargsList
setMethod(f = "[[", signature = c("SYSargsList", "ANY", "missing"), definition = function(x, i, ..., drop) {
    return(as(x, "list")[[i]])
})

## Behavior of "$" operator for SYSargsList
setMethod("$",
          signature = "SYSargsList",
          definition = function(x, name) {
              slot(x, name)
          }
)

setMethod("SampleName", signature = "SYSargsList", definition = function(x, step) {
    ## Check steps
    if (inherits(step, "numeric")) {
        if (!step %in% 1:length(x)) stop("We can not find this step in the Workflow")
    } else if (inherits(step, "character")) {
        if (!step %in% stepName(x)) stop("We can not find this step in the Workflow")
    }
    if (!is.null(targetsWF(x)[[step]]$SampleName)) {
        return(targetsWF(x)[[step]]$SampleName)
    } else if (is.null(targetsWF(x)[[step]]$SampleName)) {
        message("This step doesn't contain multiple samples.")
    }
})

setMethod("baseCommand", signature = "SYSargsList", definition = function(x, step) {
    return(baseCommand(x[step]$stepsWF[[1]]))
})

setMethod("stepName", signature = "SYSargsList", definition = function(x) {
    return(names(stepsWF(x)))
})

setMethod("targetsheader", signature = "SYSargsList", definition = function(x, step) {
    return(stepsWF(x)[[step]]$targetsheader)
})

## cmdlist method for SYSargslist
setMethod(f = "cmdlist", signature = "SYSargsList", definition = function(x, step, targets = NULL) {
    if (missing(step)) {
        step <- 1:length(x)
    }
    step <- .StepClass(x, class="SYSargs2", step)
    if(length(step)==0) stop("No selected step is a 'SYSargs2' object.")
    x <- x[step]
    cmd <- sapply(names(x$stepsWF), function(x) list(NULL))
    for (i in seq_along(x)) {
        if (nchar(cmdlist(x$stepsWF[[i]])[[1]][[1]]) > 0) {
            cmd_list <- cmdlist(x$stepsWF[[i]])
            if (!is.null(targets)) {
                cmd_list <- cmd_list[targets]
            }
            cmd[[i]] <- cmd_list
        }
    }
    return(cmd)
})

setMethod(f = "yamlinput", signature = "SYSargsList", definition = function(x, step) {
    if (inherits(stepsWF(x)[[step]], "LineWise")) stop("Provide a stepWF with a 'SYSargs2' class")
    stepsWF(x)[[step]]$yamlinput
})

## Behavior of "subset" method for SYSargsList
setMethod(f = "subset", signature = "SYSargsList", definition = function(x, subset_steps, input_targets, keep_steps = TRUE) {
    x_sub <- x[subset_steps]
    ## check subset_steps length
    if (length(unique(sapply(stepsWF(x_sub), function(x) length(x)))) > 1) stop("All 'subset_steps' should contain the same length.")
    if (missing(input_targets)) {
        input_targets <- 1:max(sapply(stepsWF(x_sub), function(x) length(x)))
    }
    # Check targets index, names
    if (inherits(input_targets, "numeric")) {
        if (!all(input_targets %in% sapply(stepsWF(x_sub), function(x) 1:length(x)))) {
            stop(
                "Please select the number of 'input_targets' accordingly, options are: ",
                paste0(1:length(x_sub@stepsWF[[1]]), collapse = ", ")
            )
        }
    } else if (inherits(input_targets, "character")) {
        if (!all(input_targets %in% sapply(stepsWF(x_sub), function(x) SampleName(x)))) {
            stop(
                "Please select the number of 'input_targets' accordingly, options are: ",
                paste0(SampleName(x_sub@stepsWF[[1]]), collapse = ", ")
            )
        }
        input_targets <- which(SampleName(x_sub$stepsWF[[1]]) %in% input_targets)
    }

    if (keep_steps == FALSE) {
        x <- x_sub
        subset_steps <- 1:length(x_sub)
    }
    for (s in subset_steps) {
        x@stepsWF[[s]] <- x@stepsWF[[s]][input_targets]
        x@statusWF[[s]]$status.completed <- x@statusWF[[s]]$status.completed[input_targets, ]
        x@statusWF[[s]]$status.time <- x@statusWF[[s]]$status.time[input_targets, ]
        x@targetsWF[[s]] <- x@targetsWF[[s]][input_targets, ]
        out <- DataFrame(x@outfiles[[s]][input_targets, ])
        colnames(out) <- colnames(x@outfiles[[s]])
        x@outfiles[[s]] <- out
        # x_sub@SEobj[[s]] <- x_sub@SEobj[[s]][i]
        x@dependency <- x@dependency
        x@targets_connection <- x@targets_connection
        x@projectInfo <- x@projectInfo
        x@runInfo$directory <- x@runInfo$directory
    }
    x <- .check_write_SYSargsList(x)
    x
})

setMethod("getColumn", signature = "SYSargsList", definition = function(x, step, position = c("outfiles", "targetsWF"), column = 1, names = SampleName(x, step)) {
    ## assertions
    stopifnot(inherits(x, "SYSargsList"))
    stopifnot(length(step) == 1)
    stopifnot(length(column) == 1)
    position <- match.arg(position, c("outfiles", "targetsWF"))
    ## Check steps
    if (inherits(step, "numeric")) {
        if (!step %in% 1:length(x)) stop("We can not find this step in the Workflow")
    } else if (inherits(step, "character")) {
        if (!step %in% stepName(x)) stop("We can not find this step in the Workflow")
    }
    ## Check column
    if (inherits(column, "numeric")) {
        if (!column %in% 1:ncol(x[[position]][[step]])) stop("We can not find this column in the Workflow")
    } else if (inherits(column, "character")) {
        if (!column %in% colnames(x[[position]][[step]])) stop("We can not find this column in the Workflow")
    }
    ## Check names
    if (!length(names) == length(x[[position]][[step]][[column]])) stop("'names' argument needs to have the same length of desired output")
    ##
    if (!is.null(x[[position]][[step]][[column]])) {
        subset <- x[[position]][[step]][[column]]
        names(subset) <- if (is.null(names)) names <- rep("", length(subset)) else names
    } else {
        message("This step doesn't contain expected outfiles.")
    }
    return(subset)
})


#' modify values of outfiles or input targets df columns
#' @param x SAL object
#' @param step
#' @param df a dataframe that must have the same rows as the modifying dataframe.
#' However, if there is no column in the original dataframe. This new df will replace the empty one.
#' If there is a non-empty dataframe, any existing columns with the same name as this new df
#' will be replaced. Any columns that do not exist in the original column but this new df will
#' be added to the original df.
#' @param position which slot in SAL to modify, one of "outfiles", "targetsWF"
#'
#' @return returns a SAL object
#' @export
#'
#' @examples
setMethod("updateColumn", signature = "SYSargsList", definition = function(x, step, df, position = c("outfiles", "targetsWF")) {
    ## assertions
    stopifnot(inherits(x, "SYSargsList"))
    stopifnot(length(step) == 1)
    stopifnot(inherits(df, c("DFrame", "data.frame")))
    position <- match.arg(position, c("outfiles", "targetsWF"))
    ## Check steps
    if (inherits(step, "numeric")) {
        if (!step %in% 1:length(x)) stop("We can not find this step in the Workflow")
    } else if (inherits(step, "character")) {
        if (!step %in% stepName(x)) stop("We can not find this step in the Workflow")
    }
    ## get some info
    df_names <- names(df)
    df_rows <- nrow(df)
    sal_name <- as.character(match.call()$x)
    ## if empty original df
    # if(nrow(x[[position]][[step]]) == 0) {x[[position]][[step]] <- as(df, "DataFrame"); return(.updateSAL(x, sal_name))}
    if (nrow(x[[position]][[step]]) == 0) {
        x[[position]][[step]] <- as(df, "DataFrame")
        return(x)
    }
    ## if not empty
    if (nrow(x[[position]][[step]]) != df_rows) stop("updateColumn: Original dataframe has different rows than the new dataframe.")
    x[[position]][[step]][, df_names] <- df
    x
    # .updateSAL(x, sal_name)
})

# .updateSAL <- function(x, sal_name){
#   assign(sal_name, x, globalenv())
# }

## Print accessor for codeLine slot
setMethod(f = "codeLine", signature = "SYSargsList", definition = function(x, step) {
    if (missing(step)) {
        step <- 1:length(x)
    }
    step <- .StepClass(x, class="LineWise", step)
    if(length(step)==0) stop("No selected step is a 'LineWise' object.")
    x <- x[step]
    rcode <- sapply(names(x$stepsWF), function(x) list(NULL))
    for (i in seq_along(x)) {
        if (!inherits(stepsWF(x)[[i]], "LineWise")) stop("This step is 'SYSargs2'. Please provide an 'LineWise' class step")
        code_list <- x$stepsWF[[i]]$codeLine
        rcode[[i]] <- code_list
    }
    for (i in seq_along(rcode)) {
        cat(crayon::blue(names(rcode[i])), paste0("    ", as.character(rcode[[i]])), sep = "\n")
    }
    # return(rcode)
})

.StepClass <- function(x, class=c("SYSargs2", "LineWise"), step){
    x_class <- sapply(stepsWF(x), function(y) class(y))
    if(any(class(step)==c("numeric", "integer"))){
        step <- stepName(x)[step]
    }
    select <- x_class[names(x_class) %in% step]
    if(all(select %in% class)){
       return(names(select)) 
    } else if(!all(select %in% class)){
        message(paste0(paste0(names(select[!select %in% class]), collapse = " AND "), " step have been dropped because it is not a ", class, " object."), "\n")
        return(names(select[select %in% class]))
        }
}

## viewEnvir() methods for SYSargslist
setMethod(f = "viewEnvir", signature = "SYSargsList", definition = function(x) {
    print(x@runInfo$env)
    print(ls(x@runInfo$env, all.names = TRUE))
})

## copyEnvir() methods for SYSargslist
setMethod(f = "copyEnvir", signature = "SYSargsList", definition = function(x, list = character(), new.env = globalenv()) {
    envir <- x@runInfo$env
    print(envir)
    if (length(list) == 0) {
        list <- ls(envir, all.names = TRUE)
    } else {
        list <- list
    }
    for (l in list) {
        assign(l, get(l, envir), new.env)
    }
    cat(paste0("Copying to 'new.env': ", "\n", paste0(list, collapse = ", ")))
})


## Replacement method for SYSargsList using "[[" operator
setReplaceMethod(f = "[[", signature = "SYSargsList", definition = function(x, i, j, value) {
    if (i == 1) x@stepsWF <- value
    if (i == 2) x@statusWF <- value
    if (i == 3) x@targetsWF <- value
    if (i == 4) x@outfiles <- value
    if (i == 5) x@SEobj <- value
    if (i == 6) x@dependency <- value
    if (i == 7) x@targets_connection <- value
    if (i == 8) x@projectInfo <- value
    if (i == 9) x@runInfo <- value
    if (i == "stepsWF") x@stepsWF <- value
    if (i == "statusWF") x@statusWF <- value
    if (i == "targetsWF") x@targetsWF <- value
    if (i == "outfiles") x@outfiles <- value
    if (i == "SEobj") x@SEobj <- value
    if (i == "dependency") x@dependency <- value
    if (i == "targets_connection") x@targets_connection <- value
    if (i == "projectInfo") x@projectInfo <- value
    if (i == "runInfo") x@runInfo <- value
    return(x)
})

## Replacement method
setReplaceMethod("appendStep", c("SYSargsList"), function(x, after = length(x), ..., value) {
    ## used in `importWF`
    on.exit(options(spr_importing = FALSE))
    ## used in `+.SYSargsList`
    on.exit(options(appendPlus = FALSE))
    ## append position
    lengx <- length(x)
    after <- after
    if (stepName(value) %in% stepName(x)) stop("Steps Names need to be unique.")
    ## Dependency
    if (all(is.na(dependency(value)) && length(x) > 0) && !getOption("spr_importing") && !getOption("appendPlus")) stop("'dependency' argument is required to append a step in the workflow.")
    if (dependency(value) == "") value[["dependency"]][[1]] <- NA
    ## Append
    if (inherits(value, "SYSargsList")) {
        value <- .validationStepConn(x, value)
        x <- sysargslist(x)
        if (names(value$stepsWF) == "Step_x") {
            step_name <- paste0("Step_", after + 1L)
            renameStep(value, 1) <- step_name
        }
        if (!after) {
            x$stepsWF <- c(value$stepsWF, x$stepsWF)
            x$targetsWF <- c(targetsWF(value), x$targetsWF)
            x$statusWF <- c(value$statusWF, x$statusWF)
            x$dependency <- c(dependency(value), x$dependency)
            x$outfiles <- c(outfiles(value), x$outfiles)
            x$targets_connection <- c(value$targets_connection, x$targets_connection)
            x$runInfo$directory <- c(value$runInfo$directory, x$runInfo$directory)
        } else if (after >= lengx) {
            x$stepsWF <- c(x$stepsWF, value$stepsWF)
            x$targetsWF <- c(x$targetsWF, targetsWF(value))
            x$statusWF <- c(x$statusWF, value$statusWF)
            x$dependency <- c(x$dependency, dependency(value))
            x$outfiles <- c(x$outfiles, outfiles(value))
            x$targets_connection <- c(x$targets_connection, value$targets_connection)
            x$runInfo$directory <- c(x$runInfo$directory, value$runInfo$directory)
        } else {
            after_tc <- names(x$stepsWF)[1L:after]
            before_tc <- names(x$stepsWF)[(after + 1L):lengx]
            x$targets_connection <- c(x$targets_connection[names(x$targets_connection) %in% after_tc], value$targets_connection, x$targets_connection[names(x$targets_connection) %in% before_tc])
            x$stepsWF <- c(x$stepsWF[1L:after], value$stepsWF, x$stepsWF[(after + 1L):lengx])
            x$targetsWF <- c(x$targetsWF[1L:after], targetsWF(value), x$targetsWF[(after + 1L):lengx])
            x$statusWF <- c(x$statusWF[1L:after], value$statusWF, x$statusWF[(after + 1L):lengx])
            x$dependency <- c(x$dependency[1L:after], dependency(value), x$dependency[(after + 1L):lengx])
            x$outfiles <- c(x$outfiles[1L:after], outfiles(value), x$outfiles[(after + 1L):lengx])
            x$runInfo$directory <- c(x$runInfo$directory[1L:after], value$runInfo$directory, x$runInfo$directory[(after + 1L):lengx])
        }
        x <- as(x, "SYSargsList")
    } else if (inherits(value, "LineWise")) {
        if (value$stepName == "Step_x") {
            step_name <- paste0("Step_", after + 1L)
        } else {
            step_name <- value$stepName
        }
        x <- sysargslist(x)
        if (!after) {
            x$stepsWF <- c(value, x$stepsWF)
            x$targetsWF <- c(list(DataFrame()), x$targetsWF)
            x$statusWF <- c(list(value$status), x$statusWF)
            x$dependency <- c(value$dependency, x$dependency)
            x$outfiles <- c(list(DataFrame()), x$outfiles)
            x$targets_connection <- c(list(NULL), x$targets_connection)
            x$runInfo$directory <- c(list(FALSE), x$runInfo$directory)
        } else if (after >= lengx) {
            x$stepsWF <- c(x$stepsWF, value)
            x$targetsWF <- c(x$targetsWF, list(DataFrame()))
            x$statusWF <- c(x$statusWF, list(value$status))
            x$dependency <- c(x$dependency, value$dependency)
            x$outfiles <- c(x$outfiles, list(DataFrame()))
            x$targets_connection <- c(x$targets_connection, list(NULL))
            x$runInfo$directory <- c(x$runInfo$directory, list(FALSE))
        } else {
            after_tc <- names(x$stepsWF)[1L:after]
            before_tc <- names(x$stepsWF)[(after + 1L):lengx]
            x$targets_connection <- c(x$targets_connection[names(x$targets_connection) %in% after_tc], list(NULL), x$targets_connection[names(x$targets_connection) %in% before_tc])
            x$stepsWF <- c(x$stepsWF[1L:after], value, x$stepsWF[(after + 1L):lengx])
            x$targetsWF <- c(x$targetsWF[1L:after], list(DataFrame()), x$targetsWF[(after + 1L):lengx])
            x$statusWF <- c(x$statusWF[1L:after], list(value$status), x$statusWF[(after + 1L):lengx])
            x$dependency <- c(x$dependency[1L:after], value$dependency, x$dependency[(after + 1L):lengx])
            x$outfiles <- c(x$outfiles[1L:after], list(DataFrame()), x$outfiles[(after + 1L):lengx])
            x$runInfo$directory <- c(x$runInfo$directory[1L:after], list(FALSE), x$runInfo$directory[(after + 1L):lengx])
        }
        names(x$stepsWF)[after + 1L] <- step_name
        names(x$statusWF)[after + 1L] <- step_name
        names(x$dependency)[after + 1L] <- step_name
        names(x$targetsWF)[after + 1L] <- step_name
        names(x$outfiles)[after + 1L] <- step_name
        names(x$targets_connection)[after + 1L] <- step_name
        names(x$runInfo$directory)[after + 1L] <- step_name
        x <- as(x, "SYSargsList")
    } else {
        stop("Argument 'value' needs to be assigned an object of class 'SYSargsList' OR 'LineWise'.")
    }
    # if(any(duplicated(names(stepsWF(x))))) warning("Duplication is found in names(stepsWF(x)). Consider renaming the steps.")
    x <- .check_write_SYSargsList(x)
    x
})

## Internal functions ## 
.check_write_SYSargsList <- function(x) {
    if (!inherits(x, "SYSargsList")) stop("Argument 'x' needs to be assigned an object of class 'SYSargsList'.")
    sys.file <- projectInfo(x)$sysargslist
    if (is.null(sys.file)) {
        if (interactive()) {
            init <- readline(cat(
                cat(crayon::bgMagenta("** Object 'x' was NOT initialized with SPRproject function! **\n")),
                "Do you like to initialize the project now using the default values?", "Type a number: \n",
                "\n 1. Please create the project \n 2. I don't need the project information \n 3. Quit"
            ))
        } else {
            ## For an non-interactive session
            init <- "1"
        }
        if (init == "1") {
            init <- SPRproject()
            x[["projectInfo"]] <- init$projectInfo
            sys.file <- projectInfo(x)$sysargslist
            write_SYSargsList(x, sys.file, silent = TRUE)
            return(x)
        } else if (init == "2") {
            print("For more details, check help(SRRproject)")
            return(x)
        } else if (init == "3") {
            stop("Quiting...")
        }
    } else if (!is.null(sys.file)) {
        sys.file <- projectInfo(x)$sysargslist
        write_SYSargsList(x, sys.file, silent = TRUE)
        return(x)
    }
}

.validationStepConn <- function(x, value) {
    ## Check outfiles names
    if (any(duplicated(unlist(append(
            lapply(outfiles(x), function(y) names(y)), lapply(outfiles(value), function(y) names(y))))
            ))) stop("'outfiles' columns names need to be unique")
    ## Check value length
    if (length(value) > 1) stop("One step can be appended in each operation.")
    targesCon <- value$targets_connection[[1]]
    if (!is.null(targesCon[[1]])) {
        step <- targesCon[[1]][[1]]
        if (any(!step %in% names(stepsWF(x)))) {
            stop(paste0(
                "'targets' argument needs to be assigned as valid targets file OR the names of a previous step, for example: ", "\n",
                paste0(names(stepsWF(x)), collapse = " OR ")
            ))
        }
        if (length(step) == 1) {
            targets_name <- paste(colnames(targetsWF(x)[step][[1]]), collapse = "|")
            new_targets_col <- targesCon[[2]][[1]][-c(which(grepl(targets_name, targesCon[[2]][[1]])))]
            ## add skip
            if (all(!new_targets_col %in% colnames(x$outfiles[[step]]))) {
                stop(paste0(
                    "'targets_column' argument needs to be assigned as valid column names of a previous step, for example: ", "\n",
                    paste0(colnames(x$outfiles[[step]]), collapse = " OR \n")
                ))
            }
            ## Check outfiles names X targets names
            if (any(new_targets_col %in% targets_name)) warning("We found duplication on 'outfiles' colnames and targetsWF colnames... Would you please make sure you are connecting the right steps? ")
            if (is.null(targesCon[[3]][[1]])) {
                old_targets <- x$targetsWF[[step]]
            } else {
                old_targets <- x$targetsWF[[step]][-c(which(grepl(paste(targesCon[[3]][[1]], collapse = "|"), colnames(x$targetsWF[[step]]))))]
            }
            new_targets <- cbind(x$outfiles[[step]][new_targets_col], old_targets)
            new_targetsheader <- targetsheader(x, step)
            ## DOUBLE CONNECTION
        } else if (length(step) > 1) {
            targets_list <- sapply(step, function(y) targetsWF(x)[[y]])
            targets_list_name <- unique(unlist(lapply(targets_list, function(y) names(y))))
            old_targets <- Reduce(function(x, y) merge(x, y, by = targets_list_name, all = TRUE), targets_list)
            targets_name <- paste(targets_list_name, collapse = "|")
            new_targets_col <- targesCon[[2]][[1]][-c(which(grepl(targets_name, targesCon[[2]][[1]])))]
            ## Check outfiles names X targets names
            if (any(new_targets_col %in% targets_name)) warning("We found duplication on 'outfiles' colnames and targetsWF colnames... Would you please make sure you are connecting the right steps? ")
            colnames_outfiles <- sapply(outfiles(x), function(y) names(y))
            if (!all(new_targets_col %in% colnames_outfiles)) {
                stop(paste0(
                    "'targets_column' argument needs to be assigned as valid column names of a previous step, for example: ", "\n",
                    paste0(colnames_outfiles, collapse = " OR \n")
                ))
            }
            if (is.null(targesCon[[3]][[1]])) {
                old_targets <- old_targets
            } else {
                old_targets <- old_targets[-c(which(grepl(paste(targesCon[[3]][[1]], collapse = "|"), colnames(old_targets))))]
            }
            new_col_list <- lapply(step, function(y) outfiles(x)[[y]])
            new_targets <- cbind(new_col_list, old_targets)
            new_targetsheader <- sapply(step, function(y) targetsheader(x, y))[1]
            names(new_targetsheader) <- "targetsheader"
        }
        WF <- value$stepsWF[[1]]
        # inputvars_v <- unlist(WF$inputvars)
        ## TODO: check inputvars...
        WF2 <- updateWF(WF, new_targets = targets.as.list(data.frame(new_targets)), new_targetsheader = new_targetsheader, inputvars = WF$inputvars, write.yaml = FALSE)
        value <- sysargslist(value)
        value$stepsWF[[1]] <- WF2
        value$targetsWF[[1]] <- as(WF2, "DataFrame")
        value$outfiles[[1]] <- output.as.df(WF2)
        value$statusWF[[1]] <- WF2$status
        value <- as(value, "SYSargsList")
    } # else if (is.null(targesCon[[1]])){
    #   value <- sysargslist(value)
    #   value$targets_connection <- list(NULL)
    #   names(value$targets_connection) <- names(value$stepsWF)
    #   value <- as(value, "SYSargsList")
    # }
    ## for SYSargs2...
    ## outfiles... if dir=TRUE
    if (inherits(value, "SYSargs2")) {
        value[["statusWF"]][[1]]$status.completed <- cbind(check.output(value)[[1]], value$statusWF[[1]]$status.completed[5:ncol(value$statusWF[[1]]$status.completed)])
    }
    if (all(!is.na(dependency(value)) && !getOption("spr_importing"))) {
        print("a")
        dep <- dependency(value)[[1]]
        if (inherits(dep, "character")) {
            if (all(!dep %in% names(stepsWF(x)))) {
                stop(
                    "'dependency' argument needs to be assigned as valid previous Step Name, for example: ", "\n",
                    paste0(names(stepsWF(x)), collapse = " OR ")
                )
            }
        } else {
            if (inherits(dep, "numeric")) {
                if (all(!dep %in% 1:length(stepsWF(x)))) {
                    stop(
                        "'dependency' argument needs to be assigned as valid previous Step Index, for example: ", "\n",
                        paste0(1:length(stepsWF(x)), collapse = " OR ")
                    )
                }
            }
        }
    }
    ## runInfo
    if ("env" %in% names(value$runInfo)) {
        value[["runInfo"]] <- value$runInfo$directory
    }
    return(value)
}

## Usage:
# appendStep(sal) <- SYSargsList(WF)
# appendStep(sal, after=0) <- SYSargsList(WF)
# appendStep(sal, after=0, step_index="test_11") <- SYSargsList(WF)

setReplaceMethod("yamlinput", c("SYSargsList"), function(x, step, paramName, value) {
    x_sub <- x[step]
    args <- x_sub@stepsWF[[1]]
    yamlinput(args, paramName) <- value
    x <- sysargslist(x)
    x$stepsWF[[1]] <- args
    x <- as(x, "SYSargsList")
    x <- .check_write_SYSargsList(x)
    x
})

setReplaceMethod("replaceStep", c("SYSargsList"), function(x, step, step_name = "default", value) {
    if (any(!inherits(value, "SYSargsList") && !inherits(value, "LineWise"))) 
        stop("Argument 'value' needs to be assigned an object of class 'SYSargsList' or 'LineWise'.")
    if (all(inherits(value, "SYSargsList") && length(value) > 1)) stop("Argument 'value' cannot have 'length(value) > 1")
    if (inherits(step, "numeric")) {
        if (step > length(x)) stop(paste0("Argument 'step' cannot be greater than ", length(x)))
    } else if (inherits(step, "character")) {
        if (!step %in% names(stepsWF(x))) {
            stop(paste0(
                "Argument 'step' needs to be assigned one of the following: ",
                paste(names(stepsWF(x)), collapse = " OR ")
            ))
        }
    }
    x <- sysargslist(x)
    if (inherits(value, "SYSargsList")){
        x$stepsWF[step] <- value$stepsWF
        x$statusWF[step] <- value$statusWF
    } else if(inherits(value, "LineWise")){
        x$stepsWF[[step]] <- value
        x$statusWF[[step]] <- value$status
    }
    if (step_name == "default") {
        name <- stepName(value)
        if (name %in% names(x$stepsWF)) {
            names(x$stepsWF)[step] <- paste0("Step_", step)
            names(x$statusWF)[step] <- paste0("Step_", step)
            names(x$dependency)[step] <- paste0("Step_", step)
            names(x$targetsWF)[step] <- paste0("Step_", step)
            names(x$outfiles)[step] <- paste0("Step_", step)
            cat(paste0("Index name of x", "[", step, "]", " was rename to ", paste0("Step_", step), " to avoid duplications."))
        } else {
            names(x$stepsWF)[step] <- stepName(value)
            names(x$statusWF)[step] <- stepName(value)
            names(x$dependency)[step] <- stepName(value)
            names(x$targetsWF)[step] <- stepName(value)
            names(x$outfiles)[step] <- stepName(value)
        }
    } else {
        names(x$stepsWF)[step] <- step_name
        names(x$statusWF)[step] <- step_name
        names(x$dependency)[step] <- step_name
        names(x$targetsWF)[step] <- step_name
        names(x$outfiles)[step] <- step_name
    }
    x <- as(x, "SYSargsList")
    if (any(duplicated(names(stepsWF(x))))) warning("Duplication is found in names(stepsWF(x)). Consider renaming the steps.")
    x <- .check_write_SYSargsList(x)
    x
})

setReplaceMethod("renameStep", c("SYSargsList"), function(x, step, ..., value) {
    if (length(step) != length(value)) stop("value argument needs to be the same length of the step for rename")
    if (inherits(value, "character")) {
        if (any(grepl("[[:space:]]", value))) message("Spaces found in the Step Name has been replaced by `_`")
        value <- gsub("[[:space:]]", "_", value)
        names(x@stepsWF)[step] <- value
        names(x@statusWF)[step] <- value
        names(x@dependency)[step] <- value
        names(x@targetsWF)[step] <- value
        names(x@outfiles)[step] <- value
        # names(x@SEobj)[step] <- value
        names(x@targets_connection)[step] <- value
        if (!is.null(x$runInfo$directory)) {
            names(x@runInfo$directory)[step] <- value
        }
    } else {
        stop("Replace value needs to be assigned an 'character' name for the workflow step.")
    }
    if (!is.null(x$projectInfo$sysargslist)) {
        x <- .check_write_SYSargsList(x)
    }
    x
})

setReplaceMethod("dependency", c("SYSargsList"), function(x, step, ..., value) {
    x@dependency[[step]] <- value
    x <- .check_write_SYSargsList(x)
    x
})

setReplaceMethod("replaceCodeLine", signature = c("SYSargsList"), function(x, step, line, value) {
    y <- x$stepsWF[step][[1]]
    if (!inherits(y, "LineWise")) stop("The step argument needs to be assigned a 'LineWise' object")
    y <- as(y, "list")
    y$codeLine <- as.character(y$codeLine)
    y$codeLine[line] <- value
    y$codeLine <- parse(text = y$codeLine)
    y <- as(y, "LineWise")
    x <- as(x, "list")
    x$stepsWF[step][[1]] <- y
    x <- as(x, "SYSargsList")
    sys.file <- projectInfo(x)$sysargslist
    write_SYSargsList(x, sys.file, silent = TRUE)
    x
})

setReplaceMethod("appendCodeLine", signature = c("SYSargsList"), function(x, step, after = NULL, value) {
    if (is.null(after)) after <- length(stepsWF(x[step])[[1]])
    y <- x$stepsWF[step][[1]]
    if (!inherits(y, "LineWise")) stop("Provide 'LineWise' class object")
    lengx <- length(y)
    y <- linewise(y)
    value <- parse(text = value)
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
    sys.file <- projectInfo(x)$sysargslist
    write_SYSargsList(x, sys.file, silent = TRUE)
    x
})

setReplaceMethod("stepsWF", c("SYSargsList"), function(x, step, ..., value) {
    x@stepsWF[[step]] <- value
    x <- .check_write_SYSargsList(x)
    x
})

setReplaceMethod("statusWF", c("SYSargsList"), function(x, step, ..., value) {
    x@statusWF[[step]] <- value
    x <- .check_write_SYSargsList(x)
    x
})