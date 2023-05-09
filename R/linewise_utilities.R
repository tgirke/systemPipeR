#######################
## LineWise function ##
#######################
LineWise <- function(code, step_name = "default", codeChunkStart = integer(),
                     rmdPath = character(), dependency = NA,
                     run_step = "mandatory",
                     run_session = "management",
                     run_remote_resources = NULL) {
    ## used in `importWF`
    on.exit({
        options(linewise_importing = FALSE)
    })

    ## check options
    run_step <- match.arg(run_step, c("mandatory", "optional"))
    run_session <- match.arg(run_session, c("management", "compute"))
    if (!is.null(run_remote_resources)) {
        if (!inherits(run_remote_resources, "list")) {
            stop("Argument 'run_remote_resources' needs to be assigned an object of class 'list'")
        }
    }
    ## Step name
    if (step_name == "default") {
        step_name <- "Step_x"
    } else {
        .checkSpecialChar(step_name)
        step_name <- step_name
    }
    ## dependency
    dependency <- list(dependency)
    names(dependency) <- step_name
    ## RunInfo
    if (!is.null(run_remote_resources)) {
        if (run_session == "management") {
            message("Please note that the '", step_name, "' run_session option '", run_session, "' was replaced with 'compute' because run_remote_resources was available.")
            run_session <- "compute"
        }
    }
    runInfo <- list(runOption = list(list(
        directory = FALSE, run_step = run_step,
        run_session = run_session,
        run_remote_resources = run_remote_resources
    )))
    names(runInfo$runOption) <- step_name
    ## status
    step_status <- list(
        status.summary = "Pending",
        status.completed = data.frame(Step = step_name, status.summary = "Pending"),
        status.time = data.frame()
    )
    ## codeLine
    if (length(deparse(substitute(code))) == 1) {
        codeLine <- parse(text = sub("^\\{|\\}$", "", deparse(substitute(code))))
    } else {
        codeLine <- deparse(substitute(code))
        codeLine[1] <- sub("^\\{", "", codeLine[1])
        codeLine[length(codeLine)] <- sub("\\}$", "", codeLine[length(codeLine)])
        codeLine <- parse(text = codeLine)
    }
    line <- list(
        codeLine = codeLine,
        codeChunkStart = codeChunkStart,
        # rmdPath = rmdPath,
        stepName = step_name,
        dependency = dependency,
        status = step_status,
        files = list(rmdPath = rmdPath),
        runInfo = runInfo
    )
    return(as(line, "LineWise"))
}
## Usage:
# xx <- LineWise({
#   x <- 1+1;  y <- 100
# })
# plot <- LineWise({
#   plot(iris)
# })

########################
## importRmd function ##
########################
importWF <- function(
        sysargs,
        file_path,
        ignore_eval = TRUE,
        update = FALSE,
        confirm = FALSE,
        check_tool = !update,
        check_module = check_tool,
        verbose = TRUE
        ) {
    on.exit({
        options(linewise_importing = FALSE)
        options(spr_importing = FALSE)
    })
    orangeText <- crayon::make_style("orange")$bold
    stopifnot(is.character(file_path) && length(file_path) == 1)
    if (!stringr::str_detect(file_path, "\\.[Rr]md$")) stop("File must be .Rmd, or .rmd ending.")
    stopifnot(is.logical(verbose) && length(verbose) == 1)
    stopifnot(is.logical(ignore_eval) && length(ignore_eval) == 1)
    stopifnot(is.logical(update) && length(update) == 1)
    stopifnot(is.logical(check_tool) && length(check_tool) == 1)
    stopifnot(is.logical(check_module) && length(check_module) == 1)
    keyword <- if(update) "update" else "import"
    if (!update && length(sysargs) > 0) {
        cat(orangeText("SYSargsList is not empty. Are you sure that you are combining two workflows?\n",
                       "If not use `update = TRUE` to re-import/update current workflow"), "\n")
    }
    ## start
    if (verbose) cat(crayon::blue$bold("Reading Rmd file\n"))
    df <- parseRmd(file_path, ignore_eval = ignore_eval, verbose = verbose)
    ## create a new env for sysargs to eval
    sysargs_env <- new.env()
    sal_imp <- sysargs
    sal_imp <- as(sal_imp, "list")
    if(update) {
        updated <- updateRmdTemplate(df = df, sal_imp = sal_imp, confirm = confirm, verbose = verbose)
        if(is.null(updated)) {
            df <- df[FALSE,]
        } else {
            sal_imp <- updated$sal_imp
            df <- updated$df
            if (verbose && nrow(df) > 0) cat(crayon::blue$bold("Now importing new steps\n"))
        }
        if(nrow(df[df$import, ]) > 0) {
            flag_import <- if(confirm) 1 else if(interactive()) menu(c("Yes", "No"), title = "Comparison done, do you want import new steps?") else 2
            df <- if(flag_import != 1) df[FALSE,] else df[df$import, ]
        } else df <- df[FALSE,]
    }
    ## adding steps
    for (i in seq_along(df$step_name)) {
        if (verbose) cat(crayon::blue$bold("Now importing step '", df$step_name[i], "' \n", sep = ""))
        options(spr_importing = TRUE)
        sal_imp <- as(sal_imp, "SYSargsList")
        salname <- sub("[\\)].*", "", sub(".*(appendStep\\()", "", df$code[i]))
        assign(salname, sal_imp, sysargs_env)
        args <- eval(parse(text = df$code[i]), envir = sysargs_env)
        if(!update) appendStep(sal_imp) <- args
        else {
            after <- length(sal_imp$stepsWF)
            if(!confirm && interactive()){
                if(df$index[i] == 1) {
                    after <- 0
                    previous_step <- ""
                } else {
                    step_names <- names(sal_imp$stepsWF)
                    previous_steps_in_df <- updated$df$step_name[seq(1, df$index[i] - 1)]
                    after_tmp <- integer(0)
                    for(s in rev(previous_steps_in_df)) {
                        after_tmp <- which(step_names == s)
                        if(length(after_tmp) != 0) {
                            after <- after_tmp
                            previous_step <- s
                            break
                        }
                    }
                }
                cat(
                    crayon::blue$bold("Step           order name\n"),
                    crayon::blue$bold("*************************\n"),
                    "Previous step: ", stringr::str_pad(after, 6, "right"), previous_step, "\n",
                    "New step:      ", stringr::str_pad(after + 1, 6, "right"), df$step_name[i], "\n",
                    sep = ""
                )
                flag_after <- menu(c("Yes", "No"), title = paste0(
                "Automatic detected the append order for ", df$step_name[i], " is after ", after, ".\n",
                "Is this the right order you want to append new step?", collapse = ""))
                after <- if(flag_after == 1) after else menu(step_names, title = "after which step you want to append this new step? 0 will before step 1. ")
            }
            appendStep(sal_imp, after = after) <- args
        }
        sal_imp <- as(sal_imp, "list")
        sal_imp$runInfo[["runOption"]][[df$step_name[i]]][["rmd_line"]] <- paste(df[i, 2:3], collapse = ":")
        sal_imp$runInfo[["runOption"]][[df$step_name[i]]][["prepro_lines"]] <- df$prepro_lines[i]
    }
    if (verbose) cat(crayon::blue$bold("Now back up current Rmd file as template for `renderReport`\n"))
    rmd_file <- file.path(sal_imp$projectInfo$project, sal_imp$projectInfo$logsDir, "workflow_template.Rmd")
    if(!file.copy(file_path, rmd_file, overwrite = TRUE)) stop("Cannot copy workflow template to workflow log directory, check path and permissions")
    # update sal
    sal_imp[["projectInfo"]]$rmd_file <- rmd_file
    # sal_imp[["projectInfo"]]$rmd_file_hash <- rlang::hash_file(rmd_file)

    if (verbose) cat(
        crayon::blue$bold("Template for renderReport is stored at"), "\n",
        rmd_file, "\n",
        crayon::make_style("orange")("Edit this file manually is not recommended"), "\n"
    )

    sal_imp <- as(sal_imp, "SYSargsList")
    sysargslist <- file.path(sal_imp$projectInfo$project, sal_imp$projectInfo$sysargslist)
    write_SYSargsList(sal_imp, sysargslist, silent = TRUE)

    if(check_tool) {
        cat(crayon::blue$bold("Now check if required tools are installed"), "\n")
        listCmdTools(sal_imp, check_path = TRUE, check_module = check_module)
    }


    if (verbose) cat(crayon::green$bold(keyword, " done\n"))
    return(sal_imp)
}

## Usage:
# file_path <- system.file("extdata", "spr_simple_wf.Rmd", package="systemPipeR")
# sal <- SPRproject(overwrite = TRUE)
# sal <- importWF(sal, file_path)

########################
## updateRmdTemplate function ##
########################


updateRmdTemplate <- function(df, sal_imp, confirm = FALSE, verbose = TRUE){
    orangeText <- crayon::make_style("orange")$bold
    if(verbose) cat(crayon::blue$bold("Update starts."),
        "Note for existing steps, update only fix the line number records. They are NOT imported again.",
        "If you have changed arguments in methods like `SYSargsList`, `Linewise`, `appendStep` in template for some steps,",
        "delete the original step from the workflow and rerun this function or manually to import it again,",
        "or use replacement methods to change arguments in current workflow, see ?`SYSargsList-class` help file.",
        "Otherwise, package would use what is in the current workflow to `renderReport` and `sal2rmd`.",
        "New arguments in the template will be ignored.",
        "\n")
    if(interactive() && !confirm && verbose)  cat(crayon::cyan$bold("Detected you have an interactive session, we may need your confirmation on some questions."), "\n")

    step_sal <- names(sal_imp$stepsWF)
    step_new <- df$step_name

    if(verbose)  cat(crayon::blue$bold("Comparing SPR steps"), "\n")
    step_new_in_sal <- step_new %in% step_sal
    df$import <- !step_new_in_sal
    df$index <- seq(1, nrow(df))
    if(!all(step_new_in_sal))
        message(crayon::blue$bold("Some new steps exist in new template but not in current SYSargsList.\n",
                              "They will be imported to workflow later. Update existing steps first."),
             "\nsteps: ", paste0(step_new[!step_new_in_sal], collapse = " ")
        )
    step_sal_in_new <- step_sal %in% step_new
    if(!all(step_sal_in_new))
        message(crayon::white$bold("Some new steps exist in SYSargsList but not in new template. Consider adding it to your template or delete from current workflow.\n"),
                "\nsteps: ", paste0(step_sal[!step_sal_in_new], collapse = " ")
        )

    if(verbose) cat(crayon::blue$bold("Comparing step orders"), "\n")
    df_order <- data.frame(
        step_name = step_new,
        order_new = seq_along(step_new),
        order_old = lapply(step_new, function(x) {match_order <- which(step_sal == x); ifelse(length(match_order) == 0, 0, match_order)}) %>%
            unlist()
    )
    diff_order <- df_order$order_new != df_order$order_old
    if(any(diff_order)) {
        if(verbose) cat("Note this function checks SPR step sequental orders, not the dependency graph.",
                        "Order change will not be immediately taken place in SYSargsList object.",
                        "New orders will be only used in `renderReport`. `sal2rmd` still uses the order in SYSargsList object.", "\n")
        cat(orangeText("Some steps in the new template have different order than SYSargsList."), "\n")
        cat(paste0(df_order$step_name[diff_order], ": ", df_order$order_old[diff_order], " -> ",   df_order$order_new[diff_order]), sep = "\n")
    }

    if(verbose) cat(crayon::blue$bold("Updating SPR steps line numbers"), "\n")

    rmd_lines_old <- lapply(sal_imp$runInfo$runOption, function(x) {lines <- x$rmd_line; if(is.null(lines)) "" else lines}) %>% unlist()
    rmd_lines_new <- paste0(df$start, ':', df$end)[step_new_in_sal]
    names(rmd_lines_new) <- step_new[step_new_in_sal]

    for(i in seq_along(rmd_lines_new)) {
        old_line_index <- which(names(rmd_lines_old) == names(rmd_lines_new[i]))
        # if(rmd_lines_old[old_line_index] == "") next
        if(rmd_lines_old[old_line_index] != rmd_lines_new[i]) {
            if(verbose) cat(crayon::blue$bold("Updating step lines of ", names(rmd_lines_old[old_line_index])),
                            rmd_lines_old[old_line_index], "->",  rmd_lines_new[i],"\n")
            sal_imp$runInfo$runOption[[names(rmd_lines_old[old_line_index])]]$rmd_line <- unname(rmd_lines_new[i])
        }
    }

    if(verbose) cat(crayon::blue$bold("Updating SPR steps preprocess code information"), "\n")
    prepro_lines_old <- lapply(sal_imp$runInfo$runOption, function(x) {lines <- x$prepro_lines; if(is.null(lines)) "" else lines}) %>% unlist()
    prepro_lines_new <- df$prepro_lines[step_new_in_sal]
    names(prepro_lines_new) <- df$step_name[step_new_in_sal]

    for(i in seq_along(prepro_lines_new)) {
        prepro_lines_index <- which(names(prepro_lines_old) == names(prepro_lines_new[i]))
        if(prepro_lines_old[prepro_lines_index] != prepro_lines_new[i]) {
            cat(crayon::blue$bold("For step"), names(prepro_lines_old[prepro_lines_index]), crayon::blue$bold("old preprocess code:"), "\n")
            cat(prepro_lines_old[prepro_lines_index], "\n")
            cat(crayon::blue$bold("New preprocess code:"), "\n")
            cat(prepro_lines_new[i], "\n")
            flag_prepro <- if(confirm) 1 else if(interactive()) menu(c("Yes", "No"), title = "Replace the preprocess code?") else 2
            if(flag_prepro == 1){
                sal_imp$runInfo$runOption[[names(prepro_lines_old[prepro_lines_index])]]$prepro_lines <- unname(prepro_lines_new[i])
            } else cat(orangeText("Skip current step's new preprocess code", "\n"))
        }
    }

    cat(crayon::blue$bold("Template update done."), "\n")
    list(sal_imp = sal_imp, df = df)
}


########################
## parseRmd function ##
########################
parseRmd <- function(file_path, ignore_eval = TRUE, verbose = FALSE) {
    # assertions
    stopifnot(is.character(file_path) && length(file_path) == 1)
    stopifnot(is.logical(ignore_eval) && length(ignore_eval) == 1)
    stopifnot(is.logical(verbose) && length(verbose) == 1)
    # read file
    lines <- readLines(file_path)
    # get chunk numbers --------
    chunk_start <- lines %>% stringr::str_which("^```\\{.*\\}.*")
    chunk_end <- lines %>% stringr::str_which("^```[[:blank:]]{0,}$")
    if (length(chunk_start) != length(chunk_end)) {
        stop("unmatched number of code chunk starts and ends")
    }
    for (i in seq_along(chunk_start)[-length(chunk_end)]) {
        if (chunk_start[i + 1] <= chunk_end[i]) {
            stop(
                "A code chunk does not end, chunk line:\n",
                chunk_start[i + 1]
            )
        }
    }
    # create df
    df <- data.frame(
        # cols export
        step_name = "",
        start = chunk_start,
        end = chunk_end,
        dep = NA,
        code = NA,
        spr = NA,
        req = NA,
        session = NA,
        has_run = FALSE,
        success = FALSE,
        sample_pass = 0,
        sample_warn = 0,
        sample_error = 0,
        sample_total = 0,
        log_path = "",
        time_start = Sys.time(),
        time_end = Sys.time() + 1,
        # cols will be removed?
        opt_text = "",
        prepro_lines = ""
    )
    # get r chunks -----
    r_chunk <- lines[chunk_start] %>% stringr::str_detect("^```\\{r")
    # get chunk start and end
    if (verbose) message(crayon::blue("\n", "---- Actions ----"))
    if (verbose && length(chunk_start[!r_chunk])) message("Ignore none-R chunks at line: ", paste0(chunk_start[!r_chunk], collapse = ", "))
    df <- df[r_chunk, ]
    if (nrow(df) == 0) stop("No valid R code chunk left")
    # get chunk text
    chunk_text <- lines[df$start] %>%
        stringr::str_remove("^```\\{r[ ]?") %>%
        stringr::str_remove("\\}[[:blank:]]{0,}$")
    # get chunk names
    chunk_names <- stringr::str_match(chunk_text, "^([^,])+")[, 1]
    # chunks with no name, no option will be ""
    chunk_names[is.na(chunk_names)] <- ""
    # chunks with no name but only options
    opt_only_chunk <- stringr::str_which(chunk_names, "=")
    opt_only_chunk_text <- chunk_names[opt_only_chunk]
    # safe to clean chunk names of chunks with no name but only options
    chunk_names[opt_only_chunk] <- ""
    # get opt string
    opt_text <- stringr::str_match(chunk_text, ",(.*)")[, 2] %>% stringr::str_remove("^[ ]{0,}")
    # preappend  no name but only options chunks options to other options
    opt_text[opt_only_chunk] <- paste0(opt_only_chunk_text, ", ", opt_text[opt_only_chunk])
    # clean ", " end options
    opt_text <- stringr::str_remove(opt_text, ",[[:blank:]]{0,}$")
    # clean empty options
    opt_text[is.na(opt_text)] <- ""
    # update chunk name
    df$step_name <- chunk_names
    df$opt_text <- opt_text
    # check eval
    if (verbose) message("Checking chunk eval values")
    opt_eval <- stringr::str_match(df$opt_text, "eval[ ]{0,}=[ ]{0,}(['\"]{0,1}\\w+['\"]{0,1})")[, 2]
    # eval unspecified will be TRUE
    opt_eval[is.na(opt_eval)] <- TRUE
    if (any(stringr::str_detect(opt_eval, "['\"]"))) {
        stop(
            "Line ",
            paste(df$start[stringr::str_detect(opt_eval, "['\"]")], collapse = ", "),
            " has a character value for `eval` option, only TRUE, T, FALSE, F are valid"
        )
    }
    opt_eval <- stringr::str_detect(opt_eval, "^TRUE|T")
    # overwrite if ignore eval param is TRUE
    if (ignore_eval) opt_eval <- rep(TRUE, length(opt_eval))
    ## update df
    if (verbose && any(!opt_eval)) message("Ignore chunks with 'eval' are FALSE or invalid values: ", paste0(df$start[!opt_eval], collapse = ", "))
    df <- df[opt_eval, ]
    if (nrow(df) == 0) stop("No evaluated code chunk detected")
    # get spr option
    if (verbose) message("Checking chunk SPR option")
    opt_spr <- stringr::str_match(df$opt_text, "spr[ ]{0,}=[ ]{0,}(\\w+)")[, 2] %>% stringr::str_detect("^TRUE|T")
    if (all(is.na(opt_spr))) stop("No valid SPR code chunk detected, remember to add `spr=TRUE` in chunk options")
    ## update df
    if (verbose && any(is.na(opt_spr))) message("Ignore non-SPR chunks: ", paste0(df$start[is.na(opt_spr)], collapse = ", "))
    # eval unspecified will be TRUE
    opt_spr[is.na(opt_spr)] <- FALSE
    df <- df[opt_spr, ]
    df$spr <- TRUE
    if (nrow(df) == 0) stop("No SPR code chunk left")
    if (verbose) message("Parse chunk code")
    code <- lapply(seq_along(df$start), function(x) {
        paste0(lines[(df$start[x] + 1):(df$end[x] - 1)], collapse = "\n")
    }) %>% unlist()
    df$code <- code
    # check if any
    if (verbose) message("Checking preprocess code for each step")
    df <- .getPreprotext(df, verbose)
    df
}


## Usage:
# rmdpath <- system.file("extdata", "spr_simple_wf.Rmd", package="systemPipeR")
# df <- parseRmd(rmdpath)


########################
## .getPreprotext function ##
########################
# get the preprocssing R code for each workflow step
# Preprocessing code is defined after SPR chunks and
# before the `###pre-end` line.
.getPreprotext <- function(df, verbose) {
    prepro_steps <- stringr::str_split(df$code, "\n") %>%
        lapply(function(x) {
            pre_pro_line <- stringr::str_which(x, "^###pre-end[ ]{0,}")
            pre_pro_line[length(pre_pro_line) == 0] <- 0
            pre_pro_line
        })
    prepro_dup <- lapply(prepro_steps, length) %>% unlist()
    if(any(prepro_dup > 1))
        stop("More than 1 preprocess delimitor `###pre-end` detected in SPR code chunk starting:",
             paste(df$start[prepro_dup[prepro_dup > 1]], collapse = " "))
    prepro_steps <- unlist(prepro_steps)
    if (verbose && sum(prepro_steps) == 0) message("No preprocessing code for SPR steps found")
    for(i in which(prepro_steps > 0)) {
        df$prepro_lines[i] <- stringr::str_split(df$code[i], "\n") %>%
            {.[[1]][1:prepro_steps[i]]} %>%
            paste0(collapse = "\n")
    }
    df
}







