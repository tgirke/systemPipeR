##############
## sal2rmd ##
##############
sal2rmd <- function(sal,
                    out_path = "spr_template.Rmd",
                    rmd_title = "SPR workflow template",
                    rmd_author = "my name",
                    rmd_date = "Last update: `r format(Sys.time(), '%d %B, %Y')`",
                    rmd_output = "html_document",
                    desc = "This is a workflow template.",
                    verbose = TRUE) {
    stopifnot(is.character(rmd_output) && length(rmd_output) == 1)
    if (verbose) message(crayon::blue$bold("sal2rmd starts, pre-checks..."))
    stopifnot(inherits(sal, "SYSargsList"))
    stopifnot(is.character(out_path) && length(out_path) == 1)
    if (!stringr::str_detect(out_path, "\\.[rR]md$")) stop("output file path must end with .Rmd or .rmd")
    stopifnot(is.character(rmd_title) && length(rmd_title) == 1)
    stopifnot(is.character(rmd_author) && length(rmd_author) == 1)
    stopifnot(is.character(rmd_date) && length(rmd_date) == 1)
    stopifnot(is.logical(verbose) && length(verbose) == 1)
    stopifnot(is.character(desc))
    if (verbose) message(crayon::blue$bold("Open", out_path, "to write"))
    on.exit(try(close(con), silent = TRUE))
    con <- file(out_path)
    open(con, "w")
    # header and desc
    if (verbose) message(crayon::blue$bold("Write Rmd header and description"))
    writeLines(c(
        "---",
        paste0('title: "', rmd_title, '"'),
        paste0('author: "', rmd_author, '"'),
        paste0('date: "', rmd_date, '"'),
        paste0("output: "),
        paste0("  ", rmd_output, ":"),
        "    number_sections: false",
        "    theme: flatly",
        "    toc: true",
        "    toc_float:",
        "      collapsed: true",
        "      smooth_scroll: false",
        "package: systemPipeR",
        "fontsize: 14pt",
        "---\n",
        "# About this Report \n",
        desc,
        "\n# Workflow Steps\n"
    ), con)
    # get some info from sal
    step_names <- names(sal$stepsWF)
    deps <- sal$dependency
    t_connects <- sal$targets_connection
    opts <- sal$runInfo$runOption
    for (i in seq_along(sal$stepsWF)) {
        if (verbose) message(crayon::blue$bold("Now writing step", i, step_names[i]))
        if (inherits(sal$stepsWF[[i]], "LineWise")) {
            .sal2rmd_rstep(sal, con, i, step_names[i], deps[[i]], opts[[i]]$run_step, opts[[i]]$run_session)
        } else {
            .sal2rmd_sysstep(sal, con, i, step_names[i], deps[[i]], opts[[i]]$directory, opts[[i]]$run_step, opts[[i]]$run_session, t_connects[[i]])
        }
    }
    if (verbose) message(crayon::green$bold("Success! File created at", out_path))
}
## Usage:
# file_path <- system.file("extdata", "spr_simple_wf.Rmd", package="systemPipeR")
# sal <- SPRproject(overwrite = TRUE)
# sal <- importWF(sal, file_path)
# sal2rmd(sal)

##############
## sal2bash ##
##############
sal2bash <- function(sal, out_dir = ".", bash_path = "/bin/bash", stop_on_error = TRUE) {
    stopifnot(inherits(sal, "SYSargsList"))
    stopifnot(is.character(out_dir) && length(out_dir) == 1)
    stopifnot(is.character(bash_path) && length(bash_path) == 1)
    stopifnot(is.logical(stop_on_error) && length(stop_on_error) == 1)
    supp_dir <- file.path(out_dir, "spr_bash")
    if (!dir.exists(supp_dir)) dir.create(supp_dir, recursive = TRUE)
    if (!dir.exists(supp_dir)) stop("Can't create", supp_dir, "check permissions")
    out_path <- file.path(out_dir, "spr_wf.sh")
    on.exit(try(close(con), silent = TRUE))
    con <- file(out_path)
    open(con, "w")
    writeLines(con = con, paste0("#!", bash_path, collapse = ""))
    if (stop_on_error) writeLines(con = con, "set -e\n")
    step_list <- .collapseSteps(sal)
    for (chunk in step_list) {
        writeLines(con = con, paste("# Step", paste0(chunk$step, collapse = "_")))
        writeLines(con = con, paste0('echo "Running step ', paste0(chunk$step, collapse = "_"), '"'))
        if (chunk$type == "sys") {
            bash_code <- .bashSysStep(sal, chunk$step)
            writeLines(con = con, bash_code)
            writeLines(con = con, "\n")
        } else {
            r_code <- .bashRStep(sal, chunk$step)
            r_path <- paste0("rscript_step", paste0(chunk$step, collapse = "_"), ".R")
            .writeRscript(r_code, file.path(supp_dir, r_path), supp_dir)
            writeLines(con = con, paste0("Rscript ", file.path(supp_dir, r_path)))
            writeLines(con = con, "\n")
        }
    }
    .loaded_pkgs <- "base"
    save(list = c(as.character(match.call()$sal), ".loaded_pkgs"), file = file.path(supp_dir, "spr_wf.RData"), envir = environment())
    message(crayon::green$bold("Success: Make sure the script 'spr_wf.sh' and directory", supp_dir, "is there before executing."))
}
## Usage:
# file_path <- system.file("extdata", "spr_simple_wf.Rmd", package="systemPipeR")
# sal <- SPRproject(overwrite = TRUE)
# sal <- importWF(sal, file_path)
# sal2bash(sal)

#################################
## Unexported helper functions ##
#################################
#####################
## .sal2rmd_rstep ##
#####################
.sal2rmd_rstep <- function(sal, con, i, step_name, dep, req, session, return = "write") {
    header <- header <- paste0("```{r ", step_name, ", eval=FALSE, spr=TRUE}", collapse = "")
    dep_code <- if (!is.na(dep[1])) paste0("    dependency = c('", paste0(dep, collapse = ","), "'),") else paste0("    dependency = NA", ",")
    rep_code <- if (req == "optional") paste0("    run_step = c('optional'),") else paste0("    run_step = c('mandatory'),")
    ses_code <- if (session == "compute") paste0("    run_session = c('compute')") else paste0("    run_session = c('management')")
    rcode <- c(
        paste0("# Step", i, " ", step_name, collapse = ""),
        header,
        "appendStep(sal) <- LineWise(\n    code={",
        paste0("        ", as.character(sal$stepsWF[[i]]$codeLine)) %>% stringr::str_replace_all("\n", "\n        "),
        "    },",
        paste0("    step_name = '", step_name, "',"),
        dep_code,
        rep_code,
        ses_code,
        ")\n```\n\n"
    )
    if (return == "write") {
        writeLines(
            con = con,
            rcode
        )
    } else if (return == "object") {
        obj <- c(rcode)
        return(obj)
    }
}

######################
## .sal2rmd_sysstep ##
######################
.sal2rmd_sysstep <- function(sal, con, i, step_name, dep, dir, req, session, t_con, return = "write") {
    header <- paste0("```{r ", step_name, ", eval=FALSE, spr=TRUE}", collapse = "")
    name_code <- paste0("    step_name = '", step_name, "',")
    dep_code <- if (!is.na(dep[1])) paste0("    dependency = c('", paste0(dep, collapse = ","), "'),") else paste0("    dependency = NA", ",")
    rep_code <- if (req == "optional") paste0("    run_step = c('optional'),") else paste0("    run_step = c('mandatory'),")
    ses_code <- if (session == "compute") paste0("    run_session = c('compute')") else paste0("    run_session = c('management')")
    targets_text <- if (!is.na(sal$stepsWF[[i]]$files$targets)) {
        paste0('    targets = "', sal$stepsWF[[i]]$files$targets, '",')
    } else if (!is.null(t_con)) {
        paste0('    targets = c("', paste0(t_con$targets_step[[1]], collapse = '", "'), '"),')
    } else {
        "    targets = NULL,"
    }
    in_var <- if (length(sal$stepsWF[[i]]$inputvars) == 0) {
        "    inputvars = NULL,"
    } else {
        in_var_names <- names(sal$stepsWF[[i]]$inputvars)
        tmp_string <- paste(in_var_names, paste0('"', unlist(sal$stepsWF[[i]]$inputvars), '"'), sep = "=", collapse = ", ")
        paste("    inputvars = c(", tmp_string, "),")
    }
    rm_col <- if (!is.null(t_con) && !is.null(t_con[["rm_targets_col"]][[1]])) {
        paste('    rm_targets_col = c("', paste(t_con[["rm_targets_col"]][[1]], collapse = '", "'), '"),')
    } else {
        paste("    rm_targets_col = NULL,")
    }
    dir_path <- if (is.na(sal$stepsWF[[i]]$files$dir_path[1])) {
        NULL
    } else {
        paste0('    dir_path="', sal$stepsWF[[i]]$files$dir_path[1], '",')
    }
    wf_file <- paste0('    wf_file="', sal$stepsWF[[i]]$files$cwl[1], '",')
    input_file <- paste0('    input_file="', sal$stepsWF[[i]]$files$yml[1], '",')
    final_obj <- c(
        paste0("# Step", i, " ", step_name, collapse = ""),
        header,
        "appendStep(sal) <- SYSargsList(",
        targets_text,
        dir_path,
        wf_file,
        input_file,
        in_var,
        rm_col,
        name_code,
        dep_code,
        if (dir) "    dir=TRUE," else "    dir=FALSE,",
        rep_code,
        ses_code,
        ")\n```\n\n"
    )
    if (return == "write") {
        writeLines(
            con = con, final_obj
        )
    } else if (return == "object") {
        return(final_obj)
    }
}

####################
## .collapseSteps ##
####################
.collapseSteps <- function(sal) {
    steps <- unlist(lapply(seq_along(sal), function(x) {
        if (inherits(sal[["stepsWF"]][[x]], "LineWise")) x else NA
    }))
    names(steps) <- seq_along(steps)
    list_item <- 1
    current_step_no <- 1
    step_list <- list(
        "1" = list(
            type = if (is.na(steps[1])) "sys" else "r",
            step = 1
        )
    )
    for (i in as.numeric(names(steps[-1]))) {
        if (is.na(steps[i])) {
            list_item <- list_item + 1
            step_list[[as.character(list_item)]][["type"]] <- "sys"
            step_list[[as.character(list_item)]][["step"]] <- i
        } else {
            if (current_step_no + 1 != i || step_list[[as.character(list_item)]][["type"]] == "sys") {
                list_item <- list_item + 1
            }
            step_list[[as.character(list_item)]][["type"]] <- "r"
            step_list[[as.character(list_item)]][["step"]] <- c(step_list[[as.character(list_item)]][["step"]], i)
        }
        current_step_no <- current_step_no + 1
    }
    step_list
}
## Usage:
# .collapseSteps(sal)

####################
## .bashSysStep ##
####################
.bashSysStep <- function(sal, step_no) {
    unlist(cmdlist(sal[step_no]))
}
## Usage:
# .bashSysStep(sal, 2)

####################
## .bashRStep ##
####################
.bashRStep <- function(sal, step_no) {
    unlist(lapply(step_no, function(x) c(paste0("#step", x, collapse = " "), as.character(sal[[1]][[x]]$codeLine), "\n")))
}
## Usage:
# .bashRStep(sal, 1:2)

####################
## .writeRscript ##
####################
.writeRscript <- function(r_code, path, supp_dir) {
    r_code <- c(
        paste0('load("', file.path(supp_dir, "spr_wf.RData"), '")'),
        "lapply(.loaded_pkgs, require, character.only = TRUE)\n",
        r_code,
        ".loaded_pkgs <- .packages()",
        paste0('save.image("', file.path(supp_dir, "spr_wf.RData"), '")\n')
    )
    writeLines(text = c("#!/usr/bin/Rscript\n\n", r_code), con = path)
}
