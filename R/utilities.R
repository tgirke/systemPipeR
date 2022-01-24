###########################################################
## Additional utilities for SYSargs and SYSargs2 objects ##
###########################################################

######################################################
## Convenience write function for targetsout(args) ##
######################################################
writeTargetsout <- function(x, file = "default", silent = FALSE, overwrite = FALSE, step = NULL, new_col = NULL,
                            new_col_output_index = NULL, remove = FALSE, ...) {
    if (all(class(x) != "SYSargs" & class(x) != "SYSargs2")) stop("Argument 'x' needs to be assigned an object of class 'SYSargs' OR 'SYSargs2")
    ## SYSargs class
    if (class(x) == "SYSargs") {
        targets <- targetsout(x)
        software <- software(x)
        if (file == "default") {
            file <- paste("targets_", software, ".txt", sep = "")
            file <- gsub(" {1,}", "_", file)
        } else {
            file <- file
        }
        headerlines <- targetsheader(x)
        ## SYSargs2 class
    } else if (class(x) == "SYSargs2") {
        if (is.null(step)) {
              stop(paste(
                  "Argument 'step' needs to be assigned one of the following values:",
                  paste(names(x$clt), collapse = ", "), "OR the corresponding position"
              ))
          }
        if (all(!is.null(step) & is.character(step) & !any(names(x$clt) %in% step))) {
              stop(paste(
                  "Argument 'step' can only be assigned one of the following values:",
                  paste(names(x$clt), collapse = ", "), "OR the corresponding position"
              ))
          }
        if (all(!is.null(step) & is.numeric(step) & !any(seq_along(names(x$clt)) %in% step))) {
              stop(paste(
                  "Argument 'step' can only be assigned one of the following position:",
                  paste(seq_along(names(x$clt)), collapse = ", "), "OR the corresponding names"
              ))
          }
        targets <- targets.as.df(targets(x))
        if (remove == TRUE) {
            targets <- targets[, -c(which(grepl("FileName", colnames(targets))))]
        }
        ## Adding the columns
        if ((!is.null(new_col) & is.null(new_col_output_index)) |
            (is.null(new_col) & !is.null(new_col_output_index)) |
            (is.null(new_col) & is.null(new_col_output_index))) {
            cat("One of 'new_col' and 'new_col_output_index' is null. It is using default column naming and adding all the output files expected, and each one will be written in a different column. \n")
            for (i in seq_len(length(output(x)[[1]][[step]]))) {
                pout <- sapply(names(output(x)), function(y) output(x)[[y]][[step]][[i]], simplify = FALSE)
                targets[[paste0(files(x)$steps, "_", i)]] <- as.character(pout)
            }
        } else if (!is.null(new_col) & !is.null(new_col_output_index)) {
            if (any(length(output(x)[[1]][[step]]) < new_col_output_index) | any(new_col_output_index < 1)) {
                stop(paste0("'new_col_output_index' argument needs to be equal or bigger than 1 and smaller than ", length(output(x)[[1]][[1]]), ", the maximum number of outputs files."))
            }
            if (length(new_col) != length(new_col_output_index)) {
                stop("'new_col' should have the same length as 'new_col_output_index'")
            }
            for (i in seq_along(new_col)) {
                pout <- sapply(names(output(x)), function(y) output(x)[[y]][[step]][[new_col_output_index[i]]], simplify = FALSE)
                targets[[as.character(new_col[i])]] <- as.character(pout)
            }
        }
        ## Workflow and Step Name
        software <- strsplit(basename(files(x)$cwl), split = "\\.")[[1]][1]
        if (is.character(step)) {
            step <- strsplit(step, split = "\\.")[[1]][1]
        } else {
            step <- strsplit(names(x$clt)[step], split = "\\.")[[1]][1]
        }
        if (file == "default") {
            file <- paste("targets_", step, ".txt", sep = "")
            file <- gsub(" {1,}", "_", file)
        } else {
            file <- file
        }
        headerlines <- targetsheader(x)[[1]]
    }
    if (file.exists(file) & overwrite == FALSE) stop(paste("I am not allowed to overwrite files; please delete existing file:", file, "or set 'overwrite=TRUE'"))
    names <- c(new_col, colnames(targets[, -c(which(grepl(paste(new_col, collapse = "|"), colnames(targets))))]))
    targets <- cbind(targets[, new_col], targets[, -c(which(grepl(paste(new_col, collapse = "|"), colnames(targets))))])
    colnames(targets) <- names
    targetslines <- c(paste(colnames(targets), collapse = "\t"), apply(targets, 1, paste, collapse = "\t"))
    writeLines(c(headerlines, targetslines), file, ...)
    if (silent != TRUE) cat("\t", "Written content of 'targetsout(x)' to file:", file, "\n")
}
## Usage:
# writeTargetsout(x=args, file="default") ## SYSargs class
# writeTargetsout(x=WF, file="default", step=1, new_col = "FileName1", new_col_output_index = 1) ## SYSargs2 class

##############################################################################
## Function to run NGS aligners including sorting and indexing of BAM files ##
##############################################################################
runCommandline <- function(args, runid = "01",
                           make_bam = FALSE, del_sam = TRUE,
                           dir = TRUE, dir.name = NULL,
                           force = FALSE, input_targets = NULL, ...) {
    ## Validation for 'args'
    if (any(!inherits(args, "SYSargs") & !inherits(args, "SYSargs2"))) stop("Argument 'args' needs to be assigned an object of class 'SYSargs' OR 'SYSargs2'")
    ## Environment Modules section
    .moduleload(args)
    ## SYSargs class ##
    if (class(args) == "SYSargs") {
        .sysargsrunCommandline(args = args, runid = runid, make_bam = make_bam, del_sam = del_sam)
    } else if (class(args) == "SYSargs2") {
        ## SYSargs2 class ##
        ## Workflow Name (Workflow OR CommandLineTool class)
        if (length(wf(args)$steps) == 0) {
            cwl.wf <- gsub("[[:space:]]", "_", paste(baseCommand(args), collapse = "_"), perl = TRUE)
        } else {
            cwl.wf <- Biostrings::strsplit(basename(files(args)$cwl), split = "\\.")[[1]][1]
        }
        ## Check if "results" folders exists...
        if (!dir.exists(file.path(yamlinput(args)$results_path$path))) {
              stop("Please check our current directory ('getwd()').
           The PATH defined at `yamlinput(args)` was not found and it is required.")
        }
        logdir <- file.path(yamlinput(args)$results_path$path)
        if (is.null(dir.name)) {
            dir.name <- cwl.wf
        } else {
            dir.name <- dir.name
        }
        ## Check if expected files exists or not
        ## Important: this happen in the new location in the case of dir=TRUE
        return <- .checkOutArgs2(args,
                                 make_bam = make_bam, dir = dir, dir.name = dir.name,
                                 force = force, del_sam = del_sam)
        args.return <- return$args_complete
        completed <- return$completed
        args <- return$args
        ## Check if the command line is exported
        if (!baseCommand(args) == c("bash")) {
            cmd_test <- tryCMD(command = baseCommand(args), silent = TRUE)
            if (cmd_test == "error") {
                stop(paste0(
                    "\n", baseCommand(args), ": command not found. ",
                    "\n", "Please make sure to configure your PATH environment variable according to the software in use."
                ), call. = FALSE)
            }
        }
        ## Create log files
        file_log <- file.path(logdir, paste0("submitargs", runid, "_", 
                                             dir.name, "_log_",
                                             format(Sys.time(), "%b%d%Y_%H%Ms%S"),
                                             paste(sample(0:9, 4), collapse = "")))
        ## Check input
        if (length(args$inputvars) >= 1) {
            inpVar <- args$inputvars
            check.inp <- colSums(sapply(inpVar, function(y) sapply(yamlinput(args), function(x) x["path"] %in% y)))
            df.targets <- targets.as.df(args$targets)[names(check.inp[check.inp > 0])]
            inp_targets2 <- FALSE
        } else {
            inp_targets2 <- TRUE
        }
        ## Targets input
        if (!is.null(input_targets)) {
            if (inherits(input_targets, c("numeric", "integer"))) {
                if (!all(input_targets %in% seq_along(cmdlist(args)))) {
                      stop("Please select the 'input_targets' number accordingly, options are: ", "\n",
                          "        ", paste0(seq_along(cmdlist(args)), collapse = ", "),
                          call. = FALSE)
                  }
            } else if (inherits(input_targets, "character")) {
                if (!all(input_targets %in% SampleName(args))) {
                      stop("Please select the 'input_targets' name accordingly, options are: ", "\n",
                          "        ", paste0(SampleName(args), collapse = ", "),
                          call. = FALSE
                      )
                  }
                input_targets <- which(SampleName(args) %in% input_targets)
            }
            input_targets <- input_targets
        } else {
            input_targets <- seq_along(cmdlist(args))
        }
        ## Sample status and Time
        sample_status <- sapply(names(cmdlist(args)), function(x) {
          list(sapply(args$files$steps, function(x) list(NULL)))
        })
        sample_status <- sample_status[input_targets]
        time_status <- args$status$status.time
        ## Progress bar
        pb <- utils::txtProgressBar(min = 0, max = length(input_targets), style = 3)
        ## LOOP
        for (i in input_targets) {
            ii <- names(cmdlist(args)[i])
            utils::setTxtProgressBar(pb, i)
            cat("## ", names(cmdlist(args)[i]), "\n", file = file_log, fill = TRUE, append = TRUE)
            ## Time
            time_status[ii, ]$time_start <- as.POSIXct(Sys.time(), origin = "1970-01-01")
            for (j in seq_along(cmdlist(args)[[i]])) {
                ## Add_command file at log file
                cat(c(
                  paste0("Time: ", paste0(format(Sys.time(), "%b%d%Y_%H%Ms%S"))), "\n",
                  "### Code: ",
                  "```{r, eval=FALSE} ",
                  cmdlist(args)[[i]][[j]][[1]],
                  "```", "\n",
                  "### Stdout: ",
                  "```{r, eval=FALSE}"
                  ), file = file_log, sep = "\n", append = TRUE)
                ## Run the commandline only for samples for which no output file is available OR force == TRUE
                if (all(force == FALSE && all(as.logical(completed[[i]][[j]])))) {
                    cat("The expected output file(s) already exist", "\n", file = file_log, fill = TRUE, append = TRUE)
                    ## close R chunk
                    cat("```", file = file_log, sep = "\n", append = TRUE)
                    sample_status[[ii]][[args$files$steps[j]]] <- "Success"
                    next()
                } else {
                    ## Create an object for executable
                    command <- gsub(" .*", "", as.character(cmdlist(args)[[i]][[j]]))
                    commandargs <- gsub("^.*? ", "", as.character(cmdlist(args)[[i]][[j]]))
                    ## Return an error if input file doesn't exist
                    if (all(inp_targets2) || all(inp_targets <- file.exists(as.character(df.targets[i, ])))) {
                        stdout <- .tryRunC(command, commandargs)
                    } else {
                        stdout <- list(
                            stdout = paste(paste0(as.character(df.targets[i, ])[!inp_targets], collapse = ", "), "\n are missing"),
                            warning = "", error = "We have an error because target(s) file doesn't exit")
                    }
                    cat(stdout$stdout, file = file_log, sep = "\n", append = TRUE)
                    ## report status
                    if (!is.null(stdout$error)) {
                        cat("## Error", file = file_log, sep = "\n", append = TRUE)
                        cat(stdout$error, file = file_log, sep = "\n", append = TRUE)
                        sample_status[[ii]][[args$files$steps[j]]] <- "Error"
                    } else if (!is.null(stdout$warning)) {
                        cat("## Warning", file = file_log, sep = "\n", append = TRUE)
                        cat(stdout$warning, file = file_log, sep = "\n", append = TRUE)
                        sample_status[[ii]][[args$files$steps[j]]] <- "Warning"
                    } else if (all(is.null(stdout$error) && is.null(stdout$warning))) {
                        sample_status[[ii]][[args$files$steps[j]]] <- "Success"
                    }
                    ## close R chunk
                    cat("```", file = file_log, sep = "\n", append = TRUE)
                }
                ## converting sam to bam using Rsamtools package
                .makeBam(output(args)[[i]][[j]], make_bam = make_bam, del_sam = del_sam)
            }
            time_status[ii, ]$time_end <- as.POSIXct(Sys.time(), origin = "1970-01-01")
        }
        args.return[["files"]][["log"]] <- file_log
        ## In the case of make_bam=TRUE
        if (make_bam == TRUE) {
            args <- .checkOutArgs2(args, make_bam = make_bam, dir = FALSE,
                                   force = force, dir.name = dir.name,
                                   del_sam = del_sam)$args_complete
        }
        ## Create recursive the subfolders
        if (dir == TRUE) {
            for (i in seq_along(names(cmdlist(args)))) {
                full_path <- file.path(logdir, dir.name)
                if (dir.exists(full_path) == FALSE) {
                    dir.create(full_path, recursive = TRUE)
                }
            }
            ## copy log files
            if (dir.exists(file.path(logdir, dir.name, "_logs")) == FALSE) {
                dir.create(file.path(logdir, dir.name, "_logs"), recursive = TRUE)
            }
            file.rename(from = file_log, to = file.path(logdir, dir.name, "_logs", basename(file_log)))
            args.return[["files"]][["log"]] <- file.path(logdir, dir.name, "_logs", basename(file_log))
            ## output FILES
            for (i in seq_along(output(args))) {
                if (length(output(args)[[i]]) > 1) {
                    for (j in seq_along(output(args)[[i]])) {
                        for (k in seq_along(output(args)[[i]][[j]])) {
                            if (file.exists(output(args)[[i]][[j]])) {
                                name <- Biostrings::strsplit(output(args)[[i]][[j]][[k]], split = "\\/")[[1]]
                                name <- name[length(name)]
                                file.rename(from = output(args)[[i]][[j]][[k]], to = file.path(logdir, dir.name, name))
                            } else if (!file.exists(output(args.return)[[i]][[j]][[k]])) {
                                dump <- "No such file or directory"
                            }
                        }
                    }
                } else if (length(output(args)[[i]]) == 1) {
                    for (j in seq_along(output(args)[[i]][[1]])) {
                        if (file.exists(output(args)[[i]][[1]][[j]])) {
                            name <- Biostrings::strsplit(output(args)[[i]][[1]][[j]], split = "\\/")[[1]]
                            name <- name[length(name)]
                            file.rename(from = output(args)[[i]][[1]][[j]], to = file.path(logdir, dir.name, name))
                        } else if (!file.exists(output(args)[[i]][[1]][[j]])) {
                            dump <- "No such file or directory"
                        }
                    }
                }
            }
        }
        ## Status and log.files
        sample_status <- lapply(sample_status, function(x) lapply(x, function(y) if (is.null(y)) y <- "Pending" else y))
        df.status <- data.frame(matrix(do.call("c", sample_status), nrow = length(sample_status), byrow = TRUE))
        colnames(df.status) <- files(args.return)$steps
        check <- check.output(args.return)[input_targets, ]
        df.status.f <- cbind(check, df.status)
        df.status.f[c(2:4)] <- sapply(df.status.f[c(2:4)], as.numeric)
        df.status.f[c(5:ncol(df.status.f))] <- sapply(df.status.f[c(5:ncol(df.status.f))], as.character)
        rownames(df.status.f) <- df.status.f$Targets
        ## updating object
        args.return[["status"]]$status.summary <- .statusSummary(df.status.f)
        args.return[["status"]]$status.completed[match(df.status.f$Targets, args.return[["status"]]$status.completed$Targets), ] <- df.status.f
        args.return[["status"]]$status.time <- time_status
        ## time
        args.return[["status"]]$status.time$time_start <- as.POSIXct(args.return[["status"]]$status.time$time_start, origin = "1970-01-01")
        args.return[["status"]]$status.time$time_end <- as.POSIXct(args.return[["status"]]$status.time$time_end, origin = "1970-01-01")
        # registering total time
        newTimeDF <- args.return[["status"]]$status.time[input_targets, ]
        args.return[["status"]]$total.time <- list(time_start = newTimeDF$time_start[1], time_end = newTimeDF$time_end[length(newTimeDF$time_end)])
        ## double check output file
        cat("\n")
        cat(crayon::blue("---- Summary ----"), "\n")
        print(df.status.f)
        close(pb)
        return(args.return)
    }
}

## Usage:
# WF <- runCommandline(WF, make_bam=TRUE) # creates the files in the ./results folder
# WF <- runCommandline(WF, dir=TRUE) # creates the files in the ./results/workflowName/Samplename folder
# WF <- runCommandline(WF, make_bam = FALSE, dir=TRUE) ## For hisat2-mapping.cwl template
# runid="01"; make_bam=FALSE; del_sam=FALSE; dir=TRUE; dir.name=NULL; force=FALSE
# runid="01"; make_bam=FALSE; del_sam=FALSE; dir=TRUE; dir.name=NULL; force=TRUE
# runid="01"; make_bam=FALSE; del_sam=FALSE; dir=TRUE; dir.name=NULL; force=TRUE

###############
## .tryRunC ##
###############
.tryRunC <- function(command, commandargs) {
    warning <- error <- NULL
    value <- withCallingHandlers(
        tryCatch(
            if (command %in% "bwa") {
                bwa_err <- tempfile()
                system2(command, args = commandargs, stdout = T, stderr = bwa_err)
                readLines(bwa_err)
            } else if (command %in% c("bash")) {
                system(paste(command, commandargs))
            } else if (isTRUE(grep("\\$", command) == 1)) {
                system(paste(command, commandargs))
            } else {
                system2(command, args = commandargs, stdout = TRUE, stderr = TRUE)
            },
            error = function(e) {
                error <<- conditionMessage(e)
                NULL
            }
        ),
        warning = function(w) {
            warning <<- append(warning, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    list(stdout = value, warning = warning, error = error)
}

## Usage:
# .tryRunC("ls", "la")


###########################################################################
## .moduleload function: internal function to module load <modules> ##
###########################################################################
.moduleload <- function(args) {
    ## Environment Modules section
    if (any(nchar(gsub(" {1,}", "", modules(args))) > 0)) {
        ## Check if "Environment Modules" is in the PATH
        # try(suppressWarnings(modulecmd_path <- system("which modulecmd", intern = TRUE, ignore.stderr = TRUE)),
        #     silent = TRUE
        # )
        modulecmd_path <- suppressMessages(is.modules.avail())
        ## "Environment Modules" is not available
        if (length(modulecmd_path) == 0) {
            dump <- "do nothing"
            ## "Environment Modules" is available and proceed the module load
        } else if (length(modulecmd_path) > 0) {
            for (j in modules(args)) moduleload(j) # loads specified software from module system
        }
    }
}


###########################################################################
## .makeBam function: internal function to convert *.sam to *.bam file ##
###########################################################################
.makeBam <- function(output_args, make_bam = TRUE, del_sam = TRUE) {
    if (all(!is.logical(c(make_bam, del_sam)))) stop("Arguments needs to be assigned 'TRUE' and 'FALSE'")
    if (make_bam == TRUE) {
        sam_files <- grepl(".sam$", output_args)
        others_files <- grepl("vcf$|bcf$|xls$|bed$", output_args)
        completed.bam <- grepl(".bam$", output_args)
        if (any(sam_files)) {
            for (k in which(sam_files)) {
                if (file.exists(output_args[k])) {
                    Rsamtools::asBam(file = output_args[k], destination = gsub("\\.sam$", "", output_args[k]), overwrite = TRUE, indexDestination = TRUE)
                    if (del_sam == TRUE) {
                        unlink(output_args[k])
                    } else if (del_sam == FALSE) {
                        dump <- "do nothing"
                    }
                } else {
                    next()
                }
            }
        } else if (any(others_files)) {
            dump <- "do nothing"
        }
        if (any(completed.bam)) { # If output is unindexed *.bam file (e.g. Tophat2)
            for (k in which(completed.bam)) {
                Rsamtools::sortBam(file = output_args[k], destination = gsub("\\.bam$", "", output_args[k]))
                Rsamtools::indexBam(output_args[k])
            }
        }
    } else if (make_bam == FALSE) {
        dump <- "do nothing"
    }
}
## Usage:
# .makeBam(output(args)[[1]][[1]], make_bam=TRUE, del_sam=FALSE)

########################################################
## .checkOutArgs2 function: internal function to check
## if the expected output has been created  ##
########################################################
.checkOutArgs2 <- function(args, make_bam, dir = dir, dir.name, force, del_sam = TRUE) {
    args_complete <- args
    if (make_bam == TRUE) {
        ## Validation for Hisat2
        if (any(grepl("samtools", names(clt(args_complete))))) {
            stop("argument 'make_bam' should be 'FALSE' when using the workflow with 'SAMtools'", call. = FALSE)
        }
        if(any(grepl(".sam", args_complete$output[[1]][[1]]))){
            args_complete <- output_update(args_complete, dir = dir, dir.name = dir.name, replace = TRUE, extension = c(".sam", ".bam"), make_bam = make_bam, del_sam = del_sam)
        }
    }
    if (dir == TRUE) {
        args_complete <- output_update(args_complete, dir = dir, dir.name = dir.name)
    }
    completed <- output(args_complete)
    for (i in seq_along(output(args_complete))) {
        for (j in seq_along(output(args_complete)[[i]])) {
            completed[[i]][[j]] <- file.exists(output(args_complete)[[i]][[j]])
            names(completed[[i]][[j]]) <- output(args_complete)[[i]][[j]]
        }
    }
    check <- check.output(args_complete, "list")
    if (force) check <- replace(check, check == TRUE, FALSE)
    for (i in seq_along(check)) {
        if (check[i] == FALSE) {
            args[["output"]][[i]] <- args$internal_outfiles[[i]]
        }
    }
    return <- list(args = args, args_complete = args_complete, completed = completed)
    return(return)
}
## Usage:
# return <- .checkOutArgs2(args, make_bam=TRUE)
# args <- return$args
# completed <- return$completed

##########################################################################################
## .sysargsrunCommandline function: Old version of runCommandline accepts SYSargs class ##
##########################################################################################
.sysargsrunCommandline <- function(args, runid = "01", make_bam = TRUE, del_sam = TRUE, ...) {
    message("This method is using the previus version of SYSargs workflow control modules.
          Please check the new version 'SYSargs2'.")
    commands <- sysargs(args)
    completed <- file.exists(outpaths(args))
    names(completed) <- outpaths(args)
    logdir <- results(args)
    for (i in seq(along = commands)) {
        ## Run alignmets only for samples for which no BAM file is available.
        if (as.logical(completed)[i]) {
            next()
        } else {
            ## Create soubmitargsID_command file
            cat(commands[i], file = paste(logdir, "submitargs", runid, sep = ""), sep = "\n", append = TRUE)
            ## Run executable
            command <- gsub(" .*", "", as.character(commands[i]))
            commandargs <- gsub("^.*? ", "", as.character(commands[i]))
            ## Execute system command; note: BWA needs special treatment in stderr handling since it writes
            ## some stderr messages to sam file if used with system2()
            if (software(args) %in% c("bwa aln", "bwa mem")) {
                stdout <- system2(command, args = commandargs, stdout = TRUE, stderr = FALSE)
            } else if (software(args) %in% c("bash_commands")) {
                stdout <- system(paste(command, commandargs))
            } else {
                stdout <- system2(command, args = commandargs, stdout = TRUE, stderr = TRUE)
            }
            ## Create submitargsID_stdout file
            cat(commands[i], file = paste(logdir, "submitargs", runid, "_log", sep = ""), sep = "\n", append = TRUE)
            cat(unlist(stdout), file = paste(logdir, "submitargs", runid, "_log", sep = ""), sep = "\n", append = TRUE)
            ## Conditional postprocessing of results
            if (make_bam == TRUE) {
                if (grepl(".sam$", outfile1(args)[i])) { # If output is *.sam file (e.g. Bowtie2)
                    Rsamtools::asBam(file = outfile1(args)[i], destination = gsub("\\.sam$", "", outfile1(args)[i]), overwrite = TRUE, indexDestination = TRUE)
                    if (del_sam == TRUE) {
                        unlink(outfile1(args)[i])
                    } else if (del_sam == FALSE) {
                        dump <- "do nothing"
                    }
                } else if (grepl("vcf$|bcf$|xls$|bed$", outpaths(args)[i])) {
                    dump <- "do nothing"
                } else { # If output is unindexed *.bam file (e.g. Tophat2)
                    Rsamtools::sortBam(file = names(completed[i]), destination = gsub("\\.bam$", "", names(completed[i])))
                    Rsamtools::indexBam(names(completed[i]))
                }
            }
        }
    }
    bamcompleted <- gsub("sam$", "bam$", file.exists(outpaths(args)))
    names(bamcompleted) <- SampleName(args)
    cat("Missing alignment results (bam files):", sum(!as.logical(bamcompleted)), "\n")
    cat("Existing alignment results (bam files):", sum(as.logical(bamcompleted)), "\n")
    return(bamcompleted)
}
## Usage:
# args <- systemArgs(sysma="param/hisat2.param", mytargets="targets.txt")
# sysargs(args)[1] # Command-line parameters for first FASTQ file
# system("hisat2-build ./data/tair10.fasta ./data/tair10.fasta")
# .sysargsrunCommandline(args=args)

############################################################################################
## batchtools-based function to submit runCommandline jobs to queuing system of a cluster ##
############################################################################################
## The advantage of this function is that it should work with most queuing/scheduling systems such as SLURM, Troque, SGE, ...
clusterRun <- function(args,
                       FUN = runCommandline,
                       more.args = list(args = args, make_bam = TRUE),
                       conffile = ".batchtools.conf.R",
                       template = "batchtools.slurm.tmpl",
                       Njobs, runid = "01", resourceList) {
    checkPkg("batchtools", quietly = FALSE)
    ## Validity checks of inputs
    if (!any(inherits(args, c("SYSargs", "SYSargs2", "SYSargsList")))) {
          stop("Argument 'args' needs to be assigned an object of class 'SYSargs' OR 'SYSargs2'")
      }
    if (!inherits(FUN, "function")) {
          stop("Value assigned to 'FUN' argument is not an object of class function.")
      }
    if (!file.exists(conffile)) {
          stop(
              "Need to point under 'conffile' argument to proper config file. ",
              "See more information here: https://mllg.github.io/batchtools/reference/makeRegistry.html. ",
              "Note: in this file *.tmpl needs to point to a valid template file."
          )
      }
    if (!file.exists(template)) {
          stop(
              "Need to point under 'template' argument to proper template file. ",
              "Sample template files for different schedulers are available here: https://github.com/mllg/batchtools/blob/master/inst/templates/"
          )
      }
    if (!class(more.args) == "list") {
          stop("'more.args' needs to be object of class 'list'.")
      }
    if (any(!names(more.args) %in% names(as.list(formals(FUN))))) {
          stop(paste(
              "The list of arguments assigned to 'more.args' can only be the ",
              "following arguments defined in the function 'FUN':",
              paste(names(as.list(formals(FUN))), collapse = ", ")
          ))
      }
    if (grepl("slurm", template)) {
        if (!grepl("slurm", Sys.getenv("PATH"))) {
            try(suppressWarnings(modulecmd_path <- system("which modulecmd", intern = TRUE, ignore.stderr = TRUE)), silent = TRUE)
            ## "Environment Modules" is not available
            if (length(modulecmd_path) == 0) {
                stop("Please export Slurm to your PATH")
                ## "Environment Modules" is available and proceed the module load
            } else if (length(modulecmd_path) > 0) {
                module_out <- moduleload("slurm") # loads specified software from module system
            }
        }
    }
    if (inherits(args, c("SYSargs", "SYSargs2"))) {
        if (!identical(FUN, runCommandline)) stop("For 'SYargs' OR 'SYSargs2' class, please use `FUN=runCommandline` function")
        ## SYSargs class
        if (inherits(args, "SYSargs")) {
            path <- normalizePath(results(args))
            args.f <- seq(along = args)
            ## SYSargs2 class
        } else if (inherits(args, "SYSargs2")) {
            path <- normalizePath(args$yamlinput$results_path$path)
            args.f <- seq(along = cmdlist(args))
        }
        ## batchtools routines
        fct <- function(i, args, ...) FUN(args = args[i], ...)
        logdir1 <- paste0(path, "/submitargs", runid, "_btdb_", paste(sample(0:9, 4), collapse = ""))
        reg <- batchtools::makeRegistry(file.dir = logdir1, conf.file = conffile, packages = "systemPipeR")
        ids <- batchtools::batchMap(fun = fct, args.f, more.args = more.args, reg = reg)
        chunk <- batchtools::chunk(ids$job.id, n.chunks = Njobs, shuffle = FALSE)
        ids$chunk <- chunk
        done <- batchtools::submitJobs(ids = ids, reg = reg, resources = resourceList)
        return(reg)
    } else if (inherits(args, c("SYSargsList"))) {
        if (!identical(FUN, runWF)) stop("For 'SYargsList' class, please use `FUN=runWF` function")
        sal <- args
        fct_sal <- function(i, sal, ...) {
            sal <- runWF(sysargs = sal, ...)
            return(sal)
        }
        path <- normalizePath(sal$projectInfo$results)
        for (j in seq_along(sal)) {
            logdir1 <- paste0(path, "/submitargs", 01, "_btdb_", paste(sample(0:9, 4), collapse = ""))
            reg <- batchtools::makeRegistry(file.dir = logdir1, conf.file = conffile, packages = "systemPipeR")
            ids <- batchtools::batchMap(fun = fct_sal, 1, more.args = list(sal = sal, steps = j), reg = reg)
            chunk <- batchtools::chunk(ids$job.id, n.chunks = 1, shuffle = FALSE)
            ids$chunk <- chunk
            done <- batchtools::submitJobs(ids = ids, reg = reg, resources = resourceList)
            batchtools::waitForJobs()
            sal <- batchtools::loadResult(reg = reg, id = ids)
        }
        return(sal)
    }
}

## Usage:
# targets <- system.file("extdata", "targets.txt", package="systemPipeR")
# dir_path <- system.file("extdata/cwl", package="systemPipeR")
# WF <- loadWorkflow(targets=targets, wf_file="hisat2/hisat2-mapping-se.cwl",
#                    input_file="hisat2/hisat2-mapping-se.yml", dir_path=dir_path)
# WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
# WF
# runCommandline(WF[1])
#
# resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024)
# reg <- clusterRun(WF[1], FUN = runCommandline,
#                   more.args = list(args = WF, make_bam = TRUE),
#                   conffile=".batchtools.conf.R",
#                   template="batchtools.slurm.tmpl",
#                   Njobs=9, runid="01", resourceList=resources)
# getStatus(reg=reg)
# waitForJobs(reg=reg)
#
# sal <- SPRproject(overwrite=TRUE)
# targetspath <- system.file("extdata/cwl/example/targets_example.txt", package="systemPipeR")
# appendStep(sal) <- SYSargsList(step_name = "echo", targets=targetspath, dir=TRUE, wf_file="example/workflow_example.cwl", input_file="example/example.yml", dir_path = system.file("extdata/cwl", package="systemPipeR"),inputvars = c(Message = "_STRING_", SampleName = "_SAMPLE_"))
#
# appendStep(sal) <- LineWise(code = {
#   hello <- lapply(getColumn(sal, step=1, 'outfiles'), function(x) yaml::read_yaml(x))
# },
# step_name = "R_getCol",
# dependency = "echo")
#
# args <- sal
# sal <- clusterRun(sal, FUN = runWF,
#                   more.args = list(),
#                   conffile=".batchtools.conf.R",
#                   template="batchtools.slurm.tmpl",
#                   Njobs=1, runid="01", resourceList=resources)

########################
## Read preprocessing ##
########################
preprocessReads <- function(args = NULL,
                            FileName1 = NULL, FileName2 = NULL,
                            outfile1 = NULL, outfile2 = NULL,
                            Fct, batchsize = 100000, overwrite = TRUE, ...) {
    ## Checking input either args or individual files
    if (all(is.null(args) && is.null(FileName1) && is.null(outfile1))) {
          stop("At least one argument must be provided, either 'args' or 'FileName1' and 'outfile1'.")
      }
    ## Checking the function argument
    if (class(Fct) != "character") stop("Argument 'Fct' needs to be of class character")
    if (!is.null(args)) {
        if (all(!inherits(args, "SYSargs") & !inherits(args, "SYSargs2"))) stop("Argument 'args' needs to be assigned an object of class 'SYSargs' OR 'SYSargs2")
        if (class(args) == "SYSargs") {
            colnames_args <- colnames(targetsout(args)) # SYSargs
            targets_in <- targetsin(args)
            if (!all(c("FileName1", "FileName2") %in% colnames_args)) {
                outpaths <- outpaths(args) # SYSargs
            } else {
                outpaths1 <- as.character(targetsout(args)$FileName1[i])
                outpaths2 <- as.character(targetsout(args)$FileName2[i])
            }
        } else if (class(args) == "SYSargs2") {
            colnames_args <- colnames(targets.as.df(args$targets)) # SYSargs2
            targets_in <- targets.as.df(args$targets)
            if (!all(c("FileName1", "FileName2") %in% colnames_args)) {
                outpaths <- subsetWF(args = args, slot = "output", subset = 1, index = 1)
            } else {
                outpaths1 <- subsetWF(args = args, slot = "output", subset = 1, index = 1)
                outpaths2 <- subsetWF(args = args, slot = "output", subset = 1, index = 2)
            }
        }
    } else if (!is.null(FileName1)) {
        if (!is.null(FileName2)) {
            colnames_args <- c("FileName1", "FileName2")
            targets_in <- data.frame(FileName1, FileName2)
            if (is.null(outfile1) | is.null(outfile2)) stop("'outfile1' and 'outfile2' are required.")
            outpaths1 <- outfile1
            outpaths2 <- outfile2
        } else {
            colnames_args <- c("FileName")
            targets_in <- data.frame(FileName1)
            if (is.null(outfile1)) stop("'outfile1' is required.")
            outpaths <- outfile1
        }
    }
    ## Run function in loop over all fastq files
    ## Single end fastq files
    if (!all(c("FileName1", "FileName2") %in% colnames_args)) {
        for (i in seq(along = outpaths)) {
            outfile <- outpaths[i]
            ## Delete existing fastq files with same names, since writeFastq will append to them
            if (overwrite == TRUE) {
                if (any(file.exists(outfile))) unlink(outfile)
            } else {
                if (any(file.exists(outfile))) stop(paste("File", outfile, "exists. Please delete file first or set overwrite=TRUE."))
            }
            ## Run preprocessor function with FastqStreamer
            counter <- 0
            f <- ShortRead::FastqStreamer(normalizePath(targets_in$FileName[i]), batchsize)
            while (length(fq <- ShortRead::yield(f))) {
                fqtrim <- eval(parse(text = Fct))
                ShortRead::writeFastq(fqtrim, outfile, mode = "a", ...)
                counter <- counter + length(fqtrim)
                cat(counter, "processed reads written to file:", outfile, "\n")
            }
            close(f)
        }
    }
    ## Paired end fastq files
    if (all(c("FileName1", "FileName2") %in% colnames_args)) {
        for (i in seq(along = outpaths1)) {
            p1 <- as.character(targets_in$FileName1[i])
            p2 <- as.character(targets_in$FileName2[i])
            p1out <- as.character(outpaths1[i])
            p2out <- as.character(outpaths2[i])
            ## Delete existing fastq files with same names, since writeFastq will append to them
            if (overwrite == TRUE) {
                if (any(file.exists(p1out))) unlink(p1out)
                if (any(file.exists(p2out))) unlink(p2out)
            } else {
                if (any(file.exists(p1out))) stop(paste("File", p1out, "exists. Please delete file first or set overwrite=TRUE."))
                if (any(file.exists(p2out))) stop(paste("File", p2out, "exists. Please delete file first or set overwrite=TRUE."))
            }
            ## Run preprocessor function with FastqStreamer
            counter1 <- 0
            counter2 <- 0
            f1 <- ShortRead::FastqStreamer(p1, batchsize)
            f2 <- ShortRead::FastqStreamer(p2, batchsize)
            while (length(fq1 <- ShortRead::yield(f1))) {
                fq2 <- ShortRead::yield(f2)
                if (length(fq1) != length(fq2)) stop("Paired end files cannot have different read numbers.")
                ## Process p1
                fq <- fq1 # for simplicity in eval
                fq1trim <- eval(parse(text = Fct))
                ## Index for p1
                index1 <- as.character(ShortRead::id(fq1)) %in% as.character(ShortRead::id(fq1trim))
                names(index1) <- seq(along = index1)
                index1 <- names(index1[index1])
                ## Process p2
                fq <- fq2 # for simplicity in eval
                fq2trim <- eval(parse(text = Fct))
                ## Index for p1
                index2 <- as.character(ShortRead::id(fq2)) %in% as.character(ShortRead::id(fq2trim))
                names(index2) <- seq(along = index2)
                index2 <- names(index2[index2])
                ## Export to processed paired files
                indexpair1 <- index1 %in% index2
                ShortRead::writeFastq(fq1trim[indexpair1], p1out, mode = "a", ...)
                indexpair2 <- index2 %in% index1
                ShortRead::writeFastq(fq2trim[indexpair2], p2out, mode = "a", ...)
                counter1 <- counter1 + sum(indexpair1)
                cat(counter1, "processed reads written to file:", p1out, "\n")
                counter2 <- counter2 + sum(indexpair2)
                cat(counter2, "processed reads written to file:", p2out, "\n")
            }
            close(f1)
            close(f2)
        }
    }
}
## Usage:
# preprocessReads(args=args, Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)", batchsize=100000, overwrite=TRUE, compress=TRUE)
# preprocessReads(FileName1 = "data/SRR446027_1.fastq.gz", FileName2 = NULL,
#                  outfile1 = "results/SRR446027_1.fastq_trim.gz", outfile2 = NULL,
#                  Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)", batchsize=100000, overwrite=TRUE, compress=TRUE)
#
# preprocessReads(FileName1 = "data/SRR446027_1.fastq.gz", FileName2 = "data/SRR446027_2.fastq.gz",
#                  outfile1 = "results/SRR446027_1.fastq_trim.gz", outfile2 = "results/SRR446027_2.fastq_trim.gz",
#                  Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)", batchsize=100000, overwrite=TRUE, compress=TRUE)

##################################################################
## Function to create sym links to bam files for viewing in IGV ##
##################################################################
symLink2bam <- function(sysargs, command = "ln -s", htmldir, ext = c(".bam", ".bai"), urlbase, urlfile) {
    ## Create URL file
    if (all(inherits(sysargs, "SYSargs") & inherits(sysargs, "SYSargs2"))) stop("Argument 'sysargs' needs to be assigned an object of class 'SYSargs' OR 'SYSargs2")
    ## SYSargs class
    if (inherits(sysargs, "SYSargs")) {
        bampaths <- outpaths(sysargs)
        symname <- SampleName(sysargs)
        ## SYSargs2 class ##
    } else if (inherits(sysargs, "SYSargs2")) {
        bampaths <- normalizePath(subsetWF(args = sysargs, slot = "output", subset = 1, index = 1))
        symname <- names(targets(sysargs))
    }
    urls <- paste(urlbase, htmldir[2], symname, ext[1], "\t", symname, sep = "")
    writeLines(urls, urlfile)
    ## Creat correspoding sym links
    dir.create(paste(htmldir, collapse = ""))
    symname <- rep(symname, each = 2)
    symname <- paste(symname, c(ext[1], paste(ext, collapse = "")), sep = "")
    bampaths2 <- as.vector(t(cbind(bampaths, paste(bampaths, ext[2], sep = ""))))
    symcommands <- paste(command, " ", bampaths2, " ", paste(htmldir, collapse = ""), symname, sep = "")
    for (i in symcommands) system(i)
}
## Usage:
# symLink2bam(sysargs=args, command="ln -s", htmldir=c("~/.html/", "somedir/"), ext=c(".bam", ".bai"), urlbase="http://cluster.hpcc.ucr.edu/~tgirke/", urlfile="IGVurl.txt")

#####################
## Alignment Stats ##
#####################
alignStats <- function(args, fqpaths, pairEnd = TRUE,
                       output_index = 1,
                       subset = "FileName1") {
    # if (class(args) %in% c("SYSargs", "SYSargs2")) {
    pairEnd <- pairEnd
    #   ## SYSargs class
    if (class(args) == "SYSargs") {
        fqpaths <- infile1(args)
        bampaths <- outpaths(args)
        if (!nchar(infile2(args))[1] > 0) pairEnd <- FALSE
        # SYSargs2 class
    } else if (class(args) == "SYSargs2") {
        fqpaths <- infile1(args)
        output.all <- subsetWF(args, slot = "output", subset = 1, index = output_index)
        bampaths <- as.character()
        for (i in seq_along(output.all)) {
            for (j in seq_along(output.all[[i]])) {
                if (grepl(".sam$", output.all[[i]][[j]]) == TRUE & grepl(".bam$", output.all[[i]][[j]]) == FALSE) {
                    stop("Please provide files in BAM format; it can be checked as 'output(args)'. Also, check 'output_update' function, if the BAM files were previously generated.")
                } else if (grepl(".bam$", output.all[[i]][[j]]) == TRUE & grepl("sorted.bam$", output.all[[i]][[j]]) == FALSE) {
                    bampaths <- c(bampaths, output.all[[i]][[j]])
                }
            }
        }
        names(bampaths) <- names(output.all)
        if (!nchar(infile2(args))[1] > 0) pairEnd <- FALSE
    } else if (class(args) == "character") {
        if (missing(fqpaths)) stop("Please provide fqpaths! It must be a named vector, names will be used as sample names")
        bampaths <- args
        fqpaths <- fqpaths
    }
    bamexists <- file.exists(bampaths)
    fqpaths <- fqpaths[bamexists]
    bampaths <- bampaths[bamexists]
    ## Obtain total read number from FASTQ files
    Nreads <- ShortRead::countLines(fqpaths) / 4
    names(Nreads) <- names(fqpaths)
    ## If reads are PE multiply by 2 as a rough approximation
    if (pairEnd) Nreads <- Nreads * 2
    ## Obtain total number of alignments from BAM files
    bfl <- BamFileList(bampaths, yieldSize = 50000, index = character())
    param <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE))
    Nalign <- Rsamtools::countBam(bfl, param = param)
    ## Obtain number of primary alignments from BAM files
    param <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isSecondaryAlignment = FALSE, isUnmappedQuery = FALSE))
    Nalignprim <- Rsamtools::countBam(bfl, param = param)
    statsDF <- data.frame(
        FileName = names(Nreads),
        Nreads = Nreads,
        Nalign = Nalign$records,
        Perc_Aligned = Nalign$records / Nreads * 100,
        Nalign_Primary = Nalignprim$records,
        Perc_Aligned_Primary = Nalignprim$records / Nreads * 100
    )
    if (pairEnd) colnames(statsDF)[which(colnames(statsDF) == "Nreads")] <- "Nreads2x"
    return(statsDF)
    # } else {
    #     stopifnot(is.character(args))
    #     stopifnot(is.character(out_dir) && length(out_dir) == 1)
    #     if(is.null(names(args))) stop("args must be a named vector, names will be used as sample names")
    #     if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    #     if(!dir.exists(out_dir)) stop("Cannot create output directory", out_dir)
    #     if(!all(check_args <- file.exists(args))) stop("Some args are missing:\n", paste0(args[!check_args], collapse = ",\n"))
    #     if(!all(check_ext <- stringr::str_detect(args, "\\.bam$"))) stop("alignStats: All args need to end with .bam\n", paste0(args[!check_ext], collapse = ",\n"))
    #
    #     file_base <- basename(args)
    #     out_args <- file.path(out_dir, gsub("\\.bam$", "_bam_stats.csv", file_base))
    #     samplenames <- names(args)
    #     bam_df <- data.frame(Sample= samplenames, File= file_base, Mapped = NA, Unmapped = NA)
    #     for (i in seq_along(args)) {
    #         bam_stats <- Rsamtools::idxstatsBam(args[i])
    #         message("Write bam stats for ", samplenames[i], " to ", basename(out_args[i]))
    #         write.csv(bam_stats, out_args[i], row.names = FALSE)
    #         map_info <- colSums(bam_stats[, c(3, 4)])
    #         bam_df[i, c(3, 4)] <- map_info
    #     }
    #     bam_df$Perc_Mapped <- round(bam_df$Mapped/(bam_df$Mapped + bam_df$Unmapped), 3) * 100
    #     return(bam_df)
    # }
}
## Usage:
# read_statsDF <- alignStats(args=args)

########################
## RPKM Normalization ##
########################
returnRPKM <- function(counts, ranges) {
    geneLengthsInKB <- sum(width(reduce(ranges))) / 1000 # Length of exon union per gene in kbp
    millionsMapped <- sum(counts) / 1e+06 # Factor for converting to million of mapped reads.
    rpm <- counts / millionsMapped # RPK: reads per kilobase of exon model.
    rpkm <- rpm / geneLengthsInKB # RPKM: reads per kilobase of exon model per million mapped reads.
    return(rpkm)
}
## Usage:
# countDFrpkm <- apply(countDF, 2, function(x) returnRPKM(counts=x, ranges=eByg))

###############################################
## Read Sample Comparisons from Targets File ##
###############################################
## Parses sample comparisons from <CMP> line(s) in targets.txt file or SYSars object.
## All possible comparisons can be specified with 'CMPset: ALL'.
readComp <- function(file, format = "vector", delim = "-") {
    if (!format %in% c("vector", "matrix")) stop("Argument format can only be assigned: vector or matrix!")
    if (class(file) == "SYSargs") {
        if (length(targetsheader(file)) == 0) stop("Input has no targets header lines.")
        comp <- targetsheader(file)
        ## SYSargs2 class
    } else if (class(file) == "SYSargs2") {
        if (length(targetsheader(file))[[1]] == 0) stop("Input has no targets header lines.")
        comp <- targetsheader(file)[[1]]
    } else {
        comp <- readLines(file)
    }
    comp <- comp[grepl("<CMP>", comp)]
    comp <- gsub("#.*<CMP>| {1,}", "", comp)
    comp <- gsub("\t", "", comp)
    comp <- gsub("^\"|\"$", "", comp) # Often required if Excel is used for editing targets file
    comp <- strsplit(comp, ":|,")
    names(comp) <- lapply(seq(along = comp), function(x) comp[[x]][1])
    comp <- sapply(names(comp), function(x) comp[[x]][-1], simplify = FALSE)
    ## Check whether all samples are present in Factor column of targets file
    checkvalues <- unique(unlist(strsplit(unlist(comp), "-")))
    checkvalues <- checkvalues[checkvalues != "ALL"]
    if (class(file) == "SYSargs") {
        all <- unique(as.character(targetsin(file)$Factor))
    } else if (class(file) == "SYSargs2") {
        all <- unique(as.character(targets.as.df(targets(file))$Factor))
    } else {
        all <- unique(as.character(read.delim(file, comment.char = "#")$Factor))
    }
    if (any(!checkvalues %in% all)) stop(paste("The following samples are not present in Factor column of targets file:", paste(checkvalues[!checkvalues %in% all], collapse = ", ")))
    ## Generate outputs
    allindex <- sapply(names(comp), function(x) any(grepl("ALL", comp[[x]])))
    if (any(allindex)) for (i in which(allindex)) comp[[i]] <- combn(all, m = 2, FUN = paste, collapse = delim)
    if (format == "vector" & delim != "-") comp <- sapply(names(comp), function(x) gsub("-", delim, comp[[x]]), simplify = FALSE)
    if (format == "vector") {
        return(comp)
    }
    if (format == "matrix") {
        return(sapply(names(comp), function(x) do.call("rbind", Biostrings::strsplit(comp[[x]], "-")), simplify = FALSE))
    }
}
## Usage:
# targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
# cmp <- readComp(targetspath, format="vector", delim="-")
# cmp <- readComp(WF, format="vector", delim="-")

#######################################################################
## Run edgeR GLM with entire count matrix or subsetted by comparison ##
#######################################################################
## If independent=TRUE then countDF will be subsetted for each comparison
run_edgeR <- function(countDF, targets, cmp, independent = TRUE, paired = NULL, mdsplot = "") {
    pkg <- c("edgeR", "limma", "grDevices")
    checkPkg(pkg, quietly = FALSE)
    ## if(class(cmp) != "matrix" & length(cmp)==2) cmp <- t(as.matrix(cmp)) # If cmp is vector of length 2, convert it to matrix.
    ## fix for _R_CHECK_LENGTH_1_LOGIC2_ error: " --- failure: the condition has length > 1 ---"
    if (all(class(cmp) != "matrix" & length(cmp) == 2)) cmp <- t(as.matrix(cmp))
    samples <- as.character(targets$Factor)
    names(samples) <- paste(as.character(targets$SampleName), "", sep = "")
    countDF <- countDF[, names(samples)]
    countDF[is.na(countDF)] <- 0
    edgeDF <- data.frame(row.names = rownames(countDF))
    group <- as.character(samples)
    if (independent == TRUE) {
        loopv <- seq(along = cmp[, 1])
    } else {
        loopv <- 1
    }
    for (j in loopv) {
        ## Filtering and normalization
        y <- edgeR::DGEList(counts = countDF, group = group) # Constructs DGEList object
        if (independent == TRUE) {
            subset <- samples[samples %in% cmp[j, ]]
            y <- y[, names(subset)]
            y$samples$group <- factor(as.character(y$samples$group))
        }
        keep <- rowSums(edgeR::cpm(y) > 1) >= 2
        y <- y[keep, ]
        y <- edgeR::calcNormFactors(y)
        ## Design matrix
        if (length(paired) == 0) {
            design <- model.matrix(~ 0 + y$samples$group, data = y$samples)
            colnames(design) <- levels(y$samples$group)
        } else {
            if (length(paired) > 0 & independent == FALSE) stop("When providing values under 'paired' also set independent=TRUE")
            Subject <- factor(paired[samples %in% cmp[j, ]]) # corrected Jun 2014 (won't change results)
            Treat <- y$samples$group
            design <- model.matrix(~ Subject + Treat)
            levels(design) <- levels(y$samples$group)
        }
        ## Estimate dispersion
        y <- edgeR::estimateGLMCommonDisp(y, design, verbose = TRUE) # Estimates common dispersions
        y <- edgeR::estimateGLMTrendedDisp(y, design) # Estimates trended dispersions
        y <- edgeR::estimateGLMTagwiseDisp(y, design) # Estimates tagwise dispersions
        fit <- edgeR::glmFit(y, design) # Fits the negative binomial GLM for each tag and produces an object of class DGEGLM with some new components.
        ## Contrast matrix is optional but makes anlysis more transparent
        if (independent == TRUE) {
            mycomp <- paste(cmp[j, 1], cmp[j, 2], sep = "-")
        } else {
            mycomp <- paste(cmp[, 1], cmp[, 2], sep = "-")
        }
        if (length(paired) == 0) contrasts <- limma::makeContrasts(contrasts = mycomp, levels = design)
        for (i in seq(along = mycomp)) {
            if (length(paired) == 0) {
                lrt <- edgeR::glmLRT(fit, contrast = contrasts[, i]) # Takes DGEGLM object and carries out the likelihood ratio test.
            } else {
                lrt <- edgeR::glmLRT(fit) # No contrast matrix with paired design
            }
            deg <- as.data.frame(edgeR::topTags(lrt, n = length(rownames(y))))
            colnames(deg) <- paste(paste(mycomp[i], collapse = "_"), colnames(deg), sep = "_")
            edgeDF <- cbind(edgeDF, deg[rownames(edgeDF), ])
        }
        if (nchar(mdsplot) > 0) {
            grDevices::pdf(paste("./results/sample_MDS_", paste(unique(subset), collapse = "-"), ".pdf", sep = ""))
            limma::plotMDS(y)
            grDevices::dev.off()
        }
    }
    return(edgeDF)
}
## Usage:
# cmp <- readComp(file=targetspath, format="matrix", delim="-")
# edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=TRUE, mdsplot="")

####################################################################
## Run DESeq2 with entire count matrix or subsetted by comparison ##
####################################################################
## If independent=TRUE then countDF will be subsetted for each comparison
run_DESeq2 <- function(countDF, targets, cmp, independent = FALSE, lfcShrink = FALSE, type = "normal") {
    pkg <- c("DESeq2")
    checkPkg(pkg, quietly = FALSE)
    ## if(class(cmp) != "matrix" & length(cmp)==2) cmp <- t(as.matrix(cmp)) # If cmp is vector of length 2, convert it to matrix.
    ## fix for _R_CHECK_LENGTH_1_LOGIC2_ error: " --- failure: the condition has length > 1 ---"
    if (all(class(cmp) != "matrix" & length(cmp) == 2)) cmp <- t(as.matrix(cmp))
    samples <- as.character(targets$Factor)
    names(samples) <- paste(as.character(targets$SampleName), "", sep = "")
    countDF <- countDF[, names(samples)]
    countDF[is.na(countDF)] <- 0
    deseqDF <- data.frame(row.names = rownames(countDF))
    if (independent == TRUE) {
        loopv <- seq(along = cmp[, 1])
    } else {
        loopv <- 1
    }
    for (j in loopv) {
        if (independent == TRUE) {
            ## Create subsetted DESeqDataSet object
            subset <- samples[samples %in% cmp[j, ]]
            countDFsub <- countDF[, names(subset)]
            dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(countDFsub), colData = data.frame(condition = subset), design = ~condition)
            mycmp <- cmp[j, , drop = FALSE]
        } else {
            ## Create full DESeqDataSet object
            dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(countDF), colData = data.frame(condition = samples), design = ~condition)
            mycmp <- cmp
        }
        ## Estimate of (i) size factors, (ii) dispersion, (iii) negative binomial GLM fitting and (iv) Wald statistics
        dds <- DESeq2::DESeq(dds, quiet = TRUE)
        for (i in seq(along = mycmp[, 1])) {
            ## Extracts DEG results for specific contrasts from DESeqDataSet object
            res <- DESeq2::results(dds, contrast = c("condition", mycmp[i, ]))
            ## lfcShrink
            if (lfcShrink == FALSE) {
                res <- DESeq2::results(dds, contrast = c("condition", mycmp[i, ]))
            } else if (lfcShrink == TRUE) {
                suppressMessages(
                    res <- DESeq2::lfcShrink(dds, contrast = c("condition", mycmp[i, ]), type = type)
                )
            }
            ## Set NAs to reasonable values to avoid errors in downstream filtering steps
            res[is.na(res[, "padj"]), "padj"] <- 1
            res[is.na(res[, "log2FoldChange"]), "log2FoldChange"] <- 0
            deg <- as.data.frame(res)
            colnames(deg)[colnames(deg) %in% c("log2FoldChange", "padj")] <- c("logFC", "FDR")
            colnames(deg) <- paste(paste(mycmp[i, ], collapse = "-"), colnames(deg), sep = "_")
            deseqDF <- cbind(deseqDF, deg[rownames(deseqDF), ])
        }
    }
    return(deseqDF)
}
## Usage:
# cmp <- readComp(file=targetspath, format="matrix", delim="-")
# degseqDF <- run_DESeq2(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=FALSE)

###########################################
## Filter DEGs by p-value and fold change ##
############################################
filterDEGs <- function(degDF, filter, plot = TRUE) {
    ## global functions or variables
    Comparisons <- Counts <- Type <- NULL
    pval <- degDF[, grep("_FDR$", colnames(degDF)), drop = FALSE]
    log2FC <- degDF[, grep("_logFC$", colnames(degDF)), drop = FALSE]
    ## DEGs that are up or down regulated
    pf <- pval <= filter["FDR"] / 100 & (log2FC >= log2(filter["Fold"]) | log2FC <= -log2(filter["Fold"]))
    colnames(pf) <- gsub("_FDR", "", colnames(pf))
    pf[is.na(pf)] <- FALSE
    DEGlistUPorDOWN <- sapply(colnames(pf), function(x) rownames(pf[pf[, x, drop = FALSE], , drop = FALSE]), simplify = FALSE)
    ## DEGs that are up regulated
    pf <- pval <= filter["FDR"] / 100 & log2FC >= log2(filter["Fold"])
    colnames(pf) <- gsub("_FDR", "", colnames(pf))
    pf[is.na(pf)] <- FALSE
    DEGlistUP <- sapply(colnames(pf), function(x) rownames(pf[pf[, x, drop = FALSE], , drop = FALSE]), simplify = FALSE)
    ## DEGs that are down regulated
    pf <- pval <= filter["FDR"] / 100 & log2FC <= -log2(filter["Fold"])
    colnames(pf) <- gsub("_FDR", "", colnames(pf))
    pf[is.na(pf)] <- FALSE
    DEGlistDOWN <- sapply(colnames(pf), function(x) rownames(pf[pf[, x, drop = FALSE], , drop = FALSE]), simplify = FALSE)
    df <- data.frame(Comparisons = names(DEGlistUPorDOWN), Counts_Up_or_Down = sapply(DEGlistUPorDOWN, length), Counts_Up = sapply(DEGlistUP, length), Counts_Down = sapply(DEGlistDOWN, length))
    resultlist <- list(UporDown = DEGlistUPorDOWN, Up = DEGlistUP, Down = DEGlistDOWN, Summary = df)
    if (plot == TRUE) {
        mytitle <- paste("DEG Counts (", names(filter)[1], ": ", filter[1], " & ", names(filter)[2], ": ", filter[2], "%)", sep = "")
        df_plot <- data.frame(Comparisons = rep(as.character(df$Comparisons), 2), Counts = c(df$Counts_Up, df$Counts_Down), Type = rep(c("Up", "Down"), each = length(df[, 1])))
        p <- ggplot2::ggplot(df_plot, ggplot2::aes(Comparisons, Counts, fill = Type)) +
            ggplot2::geom_bar(position = "stack", stat = "identity") +
            ggplot2::coord_flip() +
            ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 0, hjust = 1)) +
            ggplot2::ggtitle(mytitle)
        print(p)
    }
    return(resultlist)
}
## Usage:
# DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=1))

#############
## showDF ##
#############
showDF <- function(data, ...) {
    checkPkg("DT", quietly = FALSE)
    DT::datatable(
        data,
        extensions = c("FixedColumns", "Scroller"),
        options = list(
            scrollX = TRUE,
            fixedColumns = TRUE,
            deferRender = TRUE,
            scrollY = 200,
            scroller = TRUE
        )
    )
}
## Usage:
# targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
# targets <- read.delim(targetspath, comment.char = "#")
# showDF(targets)

###############
## checkPkg ##
###############
# param pkg character vector of the names of packages
checkPkg <- function(pkg, quietly = TRUE) {
    if (!inherits(pkg, "character")) stop("Argument 'pkg' needs to be assigned an object of class 'character'")
    for (i in pkg) {
        if (!requireNamespace(i, quietly = quietly)) {
            stop(paste0("Package '", i, "' should be installed.", "\n"), call. = FALSE)
        }
    }
}
# pkg <- c("AnnotationDbi", "GO.db")
# pkg <- c("annotate")
# checkPkg(pkg, quietly = FALSE)
