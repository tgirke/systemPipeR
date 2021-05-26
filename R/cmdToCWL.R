## createParamFiles ##
createParamFiles <- function(commandline, cwlVersion = "v1.1", class = "CommandLineTool", 
                             results_path = "./results", module_load = "baseCommand", 
                             file = "default", writeParamFiles = FALSE, 
                             overwrite = FALSE, silent = FALSE) {
    commandline <- .modifyCwlParse(commandline = commandline)
    if (interactive()) {
        correct <- readline(cat(
            cat(crayon::bgMagenta("*****Checking Output*****\n")),
            "Is it all correct? ", "\n",
            "Would you like to proceed now? Type a number: \n 1. Yes \n 2. No \n"
        ))
    } else {
        ## For an non-interactive session
        correct <- "1"
    }
    if (correct == "1") {
        WF <- createWF(
            targets = NULL, commandline, results_path = results_path,
            module_load = module_load, overwrite = TRUE, cwlVersion = cwlVersion, class = class
        )
        WF <- renderWF(WF)
        if (writeParamFiles) {
            ## Action to save files
            writeParamFiles(sysargs = WF, file = file, overwrite = overwrite, silent = silent)
        } 
        ## Return
        return(WF)
    } else if (correct == "2") {
        WF <- createWF(
            targets = NULL, commandline, results_path = results_path,
            module_load = module_load, overwrite = TRUE, cwlVersion = cwlVersion, class = class
        )
        return(WF)
    }
}

## Usage: 
# command <- "
# hisat2 \
#     -S <F, out: ./results/M1A.sam> \
#     -x <F: ./data/tair10.fasta> \
#     -k <int: 1> \
#     -min-intronlen <int: 30> \
#     -max-intronlen <int: 3000> \
#     -threads <int: 4> \
#     -U <F: ./data/SRR446027_1.fastq.gz> \
#     --verbose
# "
# cmd <- createParamFiles(command)
# cmdlist(cmd)


## writeParamFiles ##
writeParamFiles <- function(sysargs, file = "default", overwrite = TRUE, silent = FALSE) {
    if (any(length(cmdlist(sysargs)[[1]]) == 0)) stop("Argument 'sysargs' needs to be assigned an object of class 'SYSargs2' fully rendered.")
    baseCommand <- clt(sysargs)[[1]]$baseCommand
    if ("default" %in% file) {
        if (dir.exists(file.path("param/cwl", baseCommand)) == FALSE) dir.create(path = file.path("param/cwl", baseCommand), recursive = TRUE)
        file.cwl <- file.path("param/cwl", baseCommand, paste0(baseCommand, ".cwl"))
        file.yml <- file.path("param/cwl", baseCommand, paste0(baseCommand, ".yml"))
    } else {
        for (i in seq_along(file)) {
            extension <- sub(".*\\.", "", file[[i]])
            if (!c("cwl") %in% extension & !c("yml") %in% extension) stop("Argument 'file' needs to be assigned as a character vector with the names of the two param file. For example, 'test.cwl' and 'test.yml'.")
            if (c("yml") %in% extension) {
                file.yml <- file[[i]]
            } else if (c("cwl") %in% extension) {
                file.cwl <- file[[i]]
            }
        }
    }
    if (file.exists(file.cwl) & overwrite == FALSE) {
          stop(paste(
              "I am not allowed to overwrite files; please delete existing file:",
              file.cwl, "or set 'overwrite=TRUE', or provide a different name in the 'file' argument"
          ))
      }
    if (file.exists(file.yml) & overwrite == FALSE) {
          stop(paste(
              "I am not allowed to overwrite files; please delete existing file:",
              file.yml, "or set 'overwrite=TRUE', or provide a different name in the 'file' argument"
          ))
      }
    args <- sysargs2(sysargs) ## TODO
    cwlVersion <- args$clt[[1]]$cwlVersion
    class <- args$clt[[1]]$class
    module_load <- args$yamlinput$ModulesToLoad[[1]]
    results_path <- args$yamlinput$results_path$path
    clt <- write.clt(args$cmdToCwl, cwlVersion = cwlVersion, class = class, file.cwl = file.cwl, writeout = TRUE, silent = silent)
    yamlinput <- write.yml(args$cmdToCwl, file.yml = file.yml, results_path = results_path, module_load = module_load, writeout = TRUE, silent = silent)
}

## Usage: 
# cmd <- createParamFiles(command)
# writeParamFiles(cmd)

## printParam ##
printParam <- function(sysargs, position, index = NULL) {
    if (inherits(sysargs, "SYSargs2")) {
        param <- cmdToCwl(sysargs)
        newParam <- .modifyCwlParse(
            commandline = param, position = position,
            index = index
        )
    } else {
        if (inherits(sysargs, "cwlParse")) {
            commandline <- .modifyCwlParse(
                commandline = sysargs, position = position,
                index = index
            )
        }
    }
}

## Usage
# printParam(cmd, position = "baseCommand")
# printParam(cmd, position = "inputs")
# printParam(cmd, position = "outputs")
# printParam(cmd, position = "inputs", index = 1:2)
# # or use character to index
# printParam(cmd, position = "inputs", index = c("S", "U"))
# printParam(cmd, position = "inputs", index = -1:-2)


## subsetParam  ##
subsetParam <- function(sysargs, position, index = NULL, trim = TRUE, mute = FALSE) {
    if (inherits(sysargs, "SYSargs2")) {
        param <- cmdToCwl(sysargs)
        newParam <- .modifyCwlParse(
            commandline = param, position = position,
            index = index, trim = trim, mute = mute
        )
        cmdToCwl(sysargs) <- newParam
        sysargs <- updateWF(sysargs)
        return(invisible(sysargs))
    } else {
        if (inherits(sysargs, "cwlParse")) {
          sysargs <- .modifyCwlParse(
                commandline = commandline, position = position,
                index = index, trim = trim, mute = mute
            )
            return(invisible(sysargs))
        }
    }
}

## Usage:
# cmd2 <- subsetParam(cmd, position = "inputs", index = 1:2, trim = TRUE)
# cmdlist(cmd2)
# cmd3 <- subsetParam(cmd, position = "inputs", index = -1, trim = TRUE)
# cmdlist(cmd3)


## replaceParam  ##
replaceParam <- function(sysargs, position, index = NULL, replace, mute = FALSE) {
    if (inherits(sysargs, "SYSargs2")) {
        param <- cmdToCwl(sysargs)
        newParam <- .modifyCwlParse(
            commandline = param, position = position,
            index = index, replacing = replace, mute = mute
        )
        cmdToCwl(sysargs) <- newParam
        sysargs <- updateWF(sysargs)
        return(invisible(sysargs))
    } else {
        if (inherits(sysargs, "cwlParse")) {
          sysargs <- .modifyCwlParse(
                commandline = commandline, position = position,
                index = index, replacing = replacing, mute = mute
            )
            return(invisible(sysargs))
        }
    }
}

## Usage:
# cmd4 <- replaceParam(cmd, "base", index = 1, replace = list(baseCommand = "bwa"))
# cmdlist(cmd4)
# 
# newIn <- new_inputs <- list(
#     "new_input1" = list(type = "File", preF="-b", yml ="myfile"),
#     "new_input2" = "-L <int: 4>"
# )
# cmd5 <- replaceParam(cmd, "inputs", index = 1:2, replace = new_inputs)
# cmdlist(cmd5)
# 
# # string format
# new_outs <- list(
#     "sam_out" = "<F: $(inputs.results_path)/test.sam>"
# )
# cmd6 <- replaceParam(cmd, "outputs", index = 1, replace = new_outs) 
# output(cmd6)
# 
# # list format
# new_outs <- list(
#     "sam_out" = list(type = "File", value = "$(inputs.results_path)/test2.sam")
# )
# cmd7 <- replaceParam(cmd, "outputs", index = 1, replace = new_outs)
# output(cmd7)

## renameParam  ##
renameParam <- function(sysargs, position, index = FALSE, rename, mute = FALSE) {
    if (inherits(sysargs, "SYSargs2")) {
        param <- cmdToCwl(sysargs)
        newParam <- .modifyCwlParse(
            commandline = param, position = position,
            index = index, rename = rename, mute = mute
        )
        cmdToCwl(sysargs) <- newParam
        sysargs <- updateWF(sysargs)
        return(invisible(sysargs))
    } else {
        if (inherits(sysargs, "cwlParse")) {
          sysargs <- .modifyCwlParse(
                commandline = commandline, position = position,
                index = index, rename = rename, mute = mute
            )
            return(invisible(sysargs))
        }
    }
}

# ## Usage:
# cmd8 <- renameParam(cmd, "inputs", index = c("S"), rename = c("Sample"))

## appendParam  ##
appendParam <- function(sysargs, position, index = NULL, append, after=NULL, mute = FALSE) {
    if (inherits(sysargs, "SYSargs2")) {
        param <- cmdToCwl(sysargs)
        newParam <- .modifyCwlParse(
            commandline = param, position = position,
            index = index, appending = append, after = after, mute = mute
        )
        cmdToCwl(sysargs) <- newParam
        sysargs <- updateWF(sysargs)
        return(invisible(sysargs))
    } else {
        if (inherits(sysargs, "cwlParse")) {
          sysargs <- .modifyCwlParse(
                commandline = commandline, position = position,
                index = index, appending = append, after = after, mute = mute
            )
            return(invisible(sysargs))
        }
    }
}


# newIn <- new_inputs <- list(
#     "new_input1" = list(type = "File", preF="-b1", yml ="myfile1"),
#     "new_input2" = list(type = "File", preF="-b2", yml ="myfile2"),
#     "new_input3" = "-b3 <F: myfile3>"
# )
# cmd9 <- appendParam(cmd, "inputs", append = new_inputs)
# cmdlist(cmd9)
# 
# cmd10 <- appendParam(cmd, "inputs", append = new_inputs, after=0)
# cmdlist(cmd10)

## Internal Functions
.cmdToCwl <- function(cmd, mute = FALSE) {
    stopifnot(is.character(cmd) && length(cmd) == 1)

    cmd_split <- stringr::str_split(cmd, "\n", simplify = TRUE) %>% # split args by line
        stringr::str_remove_all("^[ ]+") %>% # remove all leading, ending spaces
        stringr::str_remove_all("[ ]+$")
    # remove empty lines
    cmd_split <- cmd_split[cmd_split != ""]
    if (any(start_check <- stringr::str_detect(cmd_split[-1], "^[^<-]"))) {
        stop("Args must start with '-' or '<' for no prefix arg: ", cmd_split[-1][start_check])
    }
    if (any(slash_check <- stringr::str_detect(cmd_split, "\\\\$"))) {
        stop("Invalid line ends with '\\\\': ", stringr::str_replace(cmd_split[slash_line], "\\\\", "\\\\\\\\"))
    }
    # base command
    cmd_base <- cmd_split[1]
    cmd_list <- cmd_split[-1] %>%
        stringr::str_replace("^<", "no_prefix <") # find no prefix args
    # find no default args and replace with keyword no_default
    no_defaults <- cmd_list[stringr::str_detect(cmd_list, "<.*>", negate = TRUE)] %>%
        stringr::str_replace("$", " <no_default: no_default>")
    # replace them
    cmd_list[stringr::str_detect(cmd_list, "<.*>", negate = TRUE)] <- no_defaults
    # grep groups
    cmd_list <- stringr::str_match(cmd_list, "(^[-]{0,2}[A-Za-z0-9-_]+) (<.*>$)")

    # get prefix
    cmd_arg_prefix <- cmd_list[, 2]
    # give unique names of no prefix args
    # use prefix without - as name
    cmd_names <- cmd_arg_prefix %>% stringr::str_remove("^[-]{0,2}")
    # give unique names of no prefix args
    cmd_names[cmd_names == "no_prefix"] <- paste0("no_prefix", seq(sum(cmd_names == "no_prefix")))
    # now safe to replace no_prefix to ""
    cmd_arg_prefix[cmd_arg_prefix == "no_prefix"] <- ""

    # get arg content
    cmd_args <- cmd_list[, 3] %>%
        stringr::str_remove("^<") %>%
        stringr::str_remove(">$") %>%
        stringr::str_split(": ")
    # get defaults
    cmd_arg_defaults <- unlist(lapply(cmd_args, `[[`, 2)) %>%
        stringr::str_remove("^[ ]+") %>%
        stringr::str_remove("[ ]+$") %>%
        stringr::str_replace("^no_default$", "")
    # get input type
    cmd_arg_info <- unlist(lapply(cmd_args, `[[`, 1)) %>% stringr::str_split(",")
    cmd_arg_type <- unlist(lapply(cmd_arg_info, `[[`, 1)) %>%
        stringr::str_remove("^[ ]+") %>%
        stringr::str_remove("[ ]+$") %>%
        stringr::str_replace("^F$", "File") %>%
        # change F to File
        stringr::str_replace("^no_default$", "") # change no_defaults to empty
    # get to know if any input will be used for output
    cmd_out_steps <- unlist(lapply(cmd_arg_info, `[`, 2)) %>%
        stringr::str_remove("^[ ]+") %>%
        stringr::str_remove("[ ]+$") %>%
        stringr::str_which("out")
    # assemble inputs
    inputs <- mapply(
        function(type, prefix, default) {
            list(type = type, preF = prefix, yml = default)
        },
        type = cmd_arg_type, prefix = cmd_arg_prefix, default = cmd_arg_defaults,
        SIMPLIFY = FALSE
    )
    names(inputs) <- cmd_names
    # assemble outputs
    outputs <- mapply(
        function(type, default) {
            list(type = type, value = default)
        },
        type = cmd_arg_type[cmd_out_steps],
        default = cmd_arg_defaults[cmd_out_steps],
        SIMPLIFY = FALSE
    )
    names(outputs) <- paste0("output", seq(length(cmd_out_steps)))
    # assemble all
    commandline <- structure(list(baseCommand = cmd_base, inputs = inputs, outputs = outputs), class = c("list", "cwlParse"))
    if (!mute) .catCmdRaw(commandline)
    invisible(commandline)
}

########## .modifyCwlParse ##########
.modifyCwlParse <- function(commandline,
                           position = NULL,
                           index = NULL,
                           trim = FALSE,
                           replacing = NULL,
                           appending = NULL,
                           after = NULL,
                           rename = NULL,
                           mute = FALSE) {
    # assertions
    if (inherits(commandline, "character")) {
        commandline <- .cmdToCwl(commandline, mute = TRUE)
    } else if (inherits(commandline, "cwlParse")) {
        commandline <- commandline
    } else {
        stop("Expecting a 'character' or 'cwlParse' object")
    }
    # If position is NULL print it and return
    if (is.null(position)) {
        .catBase(commandline$baseCommand)
        .catInputs(commandline$inputs)
        .catOutputs(commandline$outputs)
        .catCmdRaw(commandline)
        return(invisible(commandline)) # always return itself
    }

    position <- match.arg(position, c("baseCommand", "inputs", "outputs"))
    if (position == "baseCommand" && !is.null(rename)) stop("'rename' does not work on 'baseCommand' position")
    stopifnot(is.character(index) || is.numeric(index) || is.null(index))
    stopifnot(is.logical(trim) && length(trim) == 1)
    stopifnot(is.logical(mute) && length(mute) == 1)

    # selected printing and/or trimming
    if (is.null(replacing) && is.null(rename) && is.null(appending) || trim) {
        commandline <- .cwlIndexTrim(commandline, position, index, trim, mute)
        return(invisible(commandline))
    }
    # replace and rename
    if (!is.null(replacing)) commandline <- .cwlReplace(commandline, position, index, replacing, mute)
    if (!is.null(rename)) commandline <- .cwlRename(commandline, position, index, rename, mute)
    if (!is.null(appending)) commandline <- .cwlAppend(commandline, position, appending, after, mute)

    return(invisible(commandline))
}

############# helper funcs ###########
## printing methods----
.catBase <- function(baseCommand) {
    cat(crayon::bgBlue("*****BaseCommand*****\n"))
    cat(baseCommand, "\n")
}

.catInputs <- function(inputs) {
    input_names <- names(inputs)
    cat(crayon::bgBlue("*****Inputs*****\n"))
    for (i in seq_along(inputs)) {
        cat(
            input_names[i], ":\n",
            "    type: ", inputs[[i]][["type"]], "\n",
            "    preF: ", inputs[[i]][["preF"]], "\n",
            "    yml: ", inputs[[i]][["yml"]], "\n",
            sep = ""
        )
    }
}

.catOutputs <- function(outputs) {
    out_names <- names(outputs)
    cat(crayon::bgBlue("*****Outputs*****\n"))
    for (i in seq_along(outputs)) {
        cat(
            out_names[i], ":\n",
            "    type: ", outputs[[i]][["type"]], "\n",
            "    value: ", outputs[[i]][["value"]], "\n",
            sep = ""
        )
    }
}

.catCmdRaw <- function(commandline) {
    cmd_string <- lapply(commandline$inputs, `[`, c("preF", "yml")) %>%
        unlist() %>%
        paste0(collapse = " ")
    cat(crayon::bgBlue("*****Parsed raw command line*****\n"))
    cat(commandline$baseCommand, cmd_string, "\n")
}

## print certain position or trimming----
.cwlIndexTrim <- function(commandline, position, index, trim, mute) {
    switch(position,
        "baseCommand" = {
            .catBase(commandline$baseCommand)
        },
        "inputs" = {
            inputs <- commandline$inputs
            if (!is.null(index)) {
                if (is.character(index)) {
                    index <- index[index %in% names(inputs)]
                } else if (is.numeric(index)) {
                    if (sum(index > 0) > 0 && sum(index < 0) > 0) {
                        stop("Cannot mix negative and positive indexing")
                    } else if (length(index) == 1 && index == 0) {
                        stop("Cannot use a single 0 to index")
                    }
                    index <- index[index <= length(inputs)]
                } else {
                    stop("Incorrect class of index")
                }
                inputs <- commandline$inputs[index]
            } else {
                inputs <- commandline$inputs
            }

            if (trim) {
                commandline$inputs <- inputs
                if (!mute) {
                    .catInputs(inputs)
                    .catCmdRaw(commandline)
                }
            } else {
                .catInputs(inputs)
            }
        },
        "outputs" = {
            outputs <- commandline$outputs
            if (!is.null(index)) {
                if (is.character(index)) {
                    index <- index[index %in% names(outputs)]
                } else if (is.numeric(index)) {
                    if (sum(index > 0) > 0 && sum(index < 0) > 0) {
                        stop("Cannot mix negative and positive indexing")
                    } else if (length(index) == 1 && index == 0) {
                        stop("Cannot use a single 0 to index")
                    }
                    index <- index[index <= length(outputs)]
                } else {
                    stop("Incorrect class of index")
                }
                outputs <- commandline$outputs[index]
            } else {
                outputs <- commandline$outputs
            }
            if (trim) {
                commandline$outputs <- outputs
                if (!mute) {
                    .catOutputs(outputs)
                    .catCmdRaw(commandline)
                }
            } else {
                .catOutputs(outputs)
            }
        }
    )
    commandline
}

## replacing method ----------
.cwlReplace <- function(commandline, position, index, replacing, mute) {
    # many assertions
    stopifnot(is.list(replacing))
    ## assert names
    re_names <- names(replacing)
    if ("" %in% re_names) stop("All items in 'replacing' must be named")
    if (length(unique(re_names)) != length(re_names)) stop("All items in 'replacing' must be unique")

    if (!is.null(replacing) && is.null(index)) stop("'replace' specified but 'index' is empty")
    if (length(replacing) != length(index)) stop("'replace' and 'index' must have the same length")
    if (length(unique(index)) != length(replacing)) stop("You have duplicated index")
    if (length(unique(names(replacing))) != length(replacing)) stop("All replacing items must be uniquely named")
    ## start to replace
    switch(position,
        "baseCommand" = {
            if (length(replacing) != 1) stop("The list to replace baseCommand needs to be length of 1")
            if (names(replacing[1]) != "baseCommand") stop("Replace list for baseCommand must have an item named 'baseCommand'")
            commandline$baseCommand <- replacing[["baseCommand"]][[1]]
            if (!mute) {
                cat("Replacing baseCommand\n")
                .catBase(commandline$baseCommand)
                .catCmdRaw(commandline)
            }
            return(invisible(commandline))
        },
        "inputs" = {
            replacing <- .parseNewIn(replacing)
            ### verify index
            if (is.character(index)) {
                # change char to numeric
                index <- which(names(commandline$inputs) %in% index)
            } else if (is.numeric(index)) {
                if (sum(index < 0) > 0) {
                    stop("Cannot use negative indexing for replacing")
                } else if (length(index) == 1 && index == 0) {
                    stop("Cannot use a single 0 to index")
                }
                index <- index[index <= length(commandline$inputs)]
            }
            if (length(index) < 1) stop("No valid index left")
            # after some invalid index removed, validate index and replacing again
            if (length(index) != length(replacing)) {
                stop(
                    "length of index is ", length(index),
                    " but length of replacing is ", length(replacing),
                    " after removing invalid indices"
                )
            }
            ## start to replace
            commandline$inputs[index] <- replacing
            ## overwrite names
            names(commandline$inputs)[index] <- names(replacing)
            if (!mute) {
                cat("Replacing inputs\n")
                .catInputs(commandline$inputs)
                .catCmdRaw(commandline)
            }
        },
        "outputs" = {
            ## assert input type
            replacing <- .parseNewOut(replacing)
            ### verify index
            if (is.character(index)) {
                # change char to numeric
                index <- which(names(commandline$outputs) %in% index)
            } else if (is.numeric(index)) {
                if (sum(index < 0) > 0) {
                    stop("Cannot use negative indexing for replacing")
                } else if (length(index) == 1 && index == 0) {
                    stop("Cannot use a single 0 to index")
                }
                index <- index[index <= length(commandline$outputs)]
            }
            if (length(index) < 1) stop("No valid index left")
            # after some invalid index removed, validate index and replacing again
            if (length(index) != length(replacing)) {
                stop(
                    "length of index is ", length(index),
                    " but length of replacing is ", length(replacing),
                    " after removing invalid indices"
                )
            }
            ## start to replace
            commandline$outputs[index] <- replacing
            ## overwrite names
            names(commandline$outputs)[index] <- names(replacing)
            if (!mute) {
                cat("Replacing outputs\n")
                .catOutputs(commandline$outputs)
                .catCmdRaw(commandline)
            }
        }
    )
    commandline
}
## rename --------
.cwlRename <- function(commandline, position, index, rename, mute) {
    if (is.null(index)) stop("'rename' specified but 'index' is empty")
    if (length(rename) != length(index)) stop("'rename' and 'index' must have the same length")
    stopifnot(is.character(rename))
    if (length(unique(rename)) != length(rename)) stop("Duplicated reanme items detected")
    if (position == "baseCommand") stop("Cannot rename baseCommand position")

    switch(position,
        "inputs" = {
            if (is.character(index)) {
                # change char to numeric
                index <- which(names(commandline$inputs) %in% index)
            } else if (is.numeric(index)) {
                if (sum(index < 0) > 0) {
                    stop("Cannot use negative indexing for rename")
                } else if (length(index) == 1 && index == 0) {
                    stop("Cannot use a single 0 to index")
                }
                index <- index[index <= length(commandline$inputs)]
            }
            if (length(index) < 1) stop("No valid index left")
            if (length(index) != length(rename)) {
                stop(
                    "length of index is ", length(index),
                    " but length of rename is ", length(rename),
                    " after removing invalid indices"
                )
            }
            ## overwrite names
            names(commandline$inputs)[index] <- rename
            if (!mute) {
                cat("Renaming inputs\n")
                .catInputs(commandline$inputs)
                .catCmdRaw(commandline)
            }
        },
        "outputs" = {
            if (is.character(index)) {
                # change char to numeric
                index <- which(names(commandline$outputs) %in% index)
            } else if (is.numeric(index)) {
                if (sum(index < 0) > 0) {
                    stop("Cannot use negative indexing for rename")
                } else if (length(index) == 1 && index == 0) {
                    stop("Cannot use a single 0 to index")
                }
                index <- index[index <= length(commandline$outputs)]
            }
            if (length(index) < 1) stop("No valid index left")
            if (length(index) != length(rename)) {
                stop(
                    "length of index is ", length(index),
                    " but length of rename is ", length(rename),
                    " after removing invalid indices"
                )
            }
            names(commandline$outputs)[index] <- rename
            if (!mute) {
                cat("Renaming outputs\n")
                .catOutputs(commandline$outputs)
                .catCmdRaw(commandline)
            }
        }
    )
    commandline
}

## appending method ----
.cwlAppend <- function(commandline, position, appending, after, mute) {
    # many assertions
    if (position == "baseCommand") stop("Cannot appending to baseCommand")
    stopifnot(is.list(appending))
    ## assert names
    ap_names <- names(appending)
    if ("" %in% ap_names) stop("All items in 'appending' must be named")
    if (length(unique(ap_names)) != length(ap_names)) stop("All items in 'replacing' must be unique")
    after <- if (is.null(after)) {
        length(commandline[[position]])
    } else {
        stopifnot(is.numeric(after) && length(after) == 1)
        as.integer(after)
    }
    ## start to replace
    switch(position,
        "inputs" = {
            ## assert input type
            appending <- .parseNewIn(appending, print_word = "appending")
            ## start to append
            commandline$inputs <- append(commandline$inputs, appending, after = after)
            if (!mute) {
                cat("Replacing inputs\n")
                .catInputs(commandline$inputs)
                .catCmdRaw(commandline)
            }
        },
        "outputs" = {
            ## assert input type
            appending <- .parseNewOut(appending, print_word = "appending")
            ## start to append
            commandline$outputs <- append(commandline$outputs, appending, after = after)
            if (!mute) {
                cat("Replacing outputs\n")
                .catOutputs(commandline$outputs)
                .catCmdRaw(commandline)
            }
        }
    )
    commandline
}



## in/out parsing methods
.parseNewIn <- function(newIn, print_word = "replacing") {
    in_names <- names(newIn)
    newIn <- lapply(seq_along(newIn), function(x) {
        if (is.list(newIn[[x]])) {
            ## assert items
            if (length(newIn[[x]]) != 3 || !all(names(newIn[[x]]) %in% c("type", "preF", "yml"))) {
                cat("You are ", print_word, " inputs:\n")
                stop(names(newIn[x]), ' item in "', print_word, '" must be length 3 and have "type", "preF" and "yml"')
            }
            ## assert item type
            lapply(newIn[[x]], function(i) {
                if (length(i) != 1 || !is.character(i)) {
                    cat("You are ", print_word, " inputs:\n")
                    stop('In "', names(newIn[x]), '": "type", "preF", "yml" in each item must all be length 1 character string')
                }
            })
            newIn[[x]]
        } else if (is.character(newIn[[x]]) && length(newIn[[x]]) == 1) {
            .parseInStr(newIn[[x]])
        } else {
            cat("You are ", print_word, " inputs:\n")
            stop("Item ", names(newIn[x]), " in ", print_word, " must be a list with 3 sub-item or a chracter string")
        }
    })
    names(newIn) <- in_names
    return(newIn)
}

.parseInStr <- function(inStr) {
    cmd_list <- stringr::str_replace(inStr, "^<", "no_prefix <") # find no prefix args
    no_defaults <- cmd_list[stringr::str_detect(cmd_list, "<.*>", negate = TRUE)] %>%
        stringr::str_replace("$", " <no_default: no_default>")
    cmd_list[stringr::str_detect(cmd_list, "<.*>", negate = TRUE)] <- no_defaults

    cmd_list <- stringr::str_match(cmd_list, "(^[-]{0,2}[A-Za-z0-9-_]+) (<.*>$)")
    if (is.na(cmd_list[1, 1])) stop("Invalid string: ", inStr)
    cmd_arg_prefix <- cmd_list[, 2]
    cmd_arg_prefix[cmd_arg_prefix == "no_prefix"] <- ""

    cmd_args <- cmd_list[, 3] %>%
        stringr::str_remove("^<") %>%
        stringr::str_remove(">$") %>%
        stringr::str_split(": ")

    cmd_arg_defaults <- unlist(lapply(cmd_args, `[[`, 2)) %>%
        stringr::str_remove("^[ ]+") %>%
        stringr::str_remove("[ ]+$") %>%
        stringr::str_replace("^no_default$", "")

    cmd_arg_info <- unlist(lapply(cmd_args, `[[`, 1)) %>% stringr::str_split(",")
    cmd_arg_type <- unlist(lapply(cmd_arg_info, `[[`, 1)) %>%
        stringr::str_remove("^[ ]+") %>%
        stringr::str_remove("[ ]+$") %>%
        stringr::str_replace("^F$", "File") %>%
        # change F to File
        stringr::str_replace("^no_default$", "")
    list(type = cmd_arg_type, preF = cmd_arg_prefix, yml = cmd_arg_defaults)
}


.parseNewOut <- function(newOut, print_word = "replacing") {
    out_names <- names(newOut)
    newOut <- lapply(seq_along(newOut), function(x) {
        if (is.list(newOut[[x]])) {
            ## assert items
            if (length(newOut[[x]]) != 2 || !all(names(newOut[[x]]) %in% c("type", "value"))) {
                cat("You are ", print_word, " outputs:\n")
                stop(names(newOut[x]), ' item in "', print_word, '" must be length 2 and have "type", "value"')
            }
            ## assert item type
            lapply(newOut[[x]], function(i) {
                if (length(i) != 1 || !is.character(i)) {
                    cat("You are ", print_word, " outputs:\n")
                    stop('In "', names(newOut[x]), '": "type", "value" in each item must all be length 1 character string')
                }
            })
            newOut[[x]]
        } else if (is.character(newOut[[x]]) && length(newOut[[x]]) == 1) {
            .parseOutStr(newOut[[x]])
        } else {
            cat("You are ", print_word, " outputs:\n")
            stop("Item ", names(newOut[x]), " in ", print_word, " must be a list with 3 sub-item or a chracter string")
        }
    })
    names(newOut) <- out_names
    return(newOut)
}

.parseOutStr <- function(outStr) {
    cmd_list <- outStr
    if (!stringr::str_detect(cmd_list, "^<.*>$")) stop("Invalid string: '", outStr, " ', needs to be inside '<...>'")
    if (length(stringr::str_extract_all(cmd_list, ":", simplify = TRUE)) != 1) {
          stop("Invalid string: '", outStr, " ', needs to have exact 1 ':'")
      }

    cmd_args <- cmd_list %>%
        stringr::str_remove("^<") %>%
        stringr::str_remove(">$") %>%
        stringr::str_split(": ")
    if (stringr::str_detect(cmd_args[[1]][1], "\\W")) {
          stop("Invalid string: '", outStr, " ', no special character for output type")
      }

    cmd_arg_defaults <- unlist(lapply(cmd_args, `[[`, 2)) %>%
        stringr::str_remove("^[ ]+") %>%
        stringr::str_remove("[ ]+$") %>%
        stringr::str_replace("^no_default$", "")

    cmd_arg_info <- unlist(lapply(cmd_args, `[[`, 1)) %>% stringr::str_split(",")
    cmd_arg_type <- unlist(lapply(cmd_arg_info, `[[`, 1)) %>%
        stringr::str_remove("^[ ]+") %>%
        stringr::str_remove("[ ]+$") %>%
        stringr::str_replace("^F$", "File") %>%
        # change F to File
        stringr::str_replace("^no_default$", "")
    list(type = cmd_arg_type, value = cmd_arg_defaults)
}
