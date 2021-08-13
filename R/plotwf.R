#' Visualize SPR workflow and status
#'
#' @param df The workflow structure dataframe
#' @param branch_method string, one of "auto", "choose". How to determine the
#' main branch of the workflow. "auto" will be determined by internal alrgothrim:
#' Branches connecting the frist and last step and/or the longest will be favored.
#' "choose" will list all possible branches and you can make a choice.
#' @param branch_no numeric, only works if `branch_method == "choose"`. Specify a
#' branch number to be the main branch instead of choosing from the prompt. This
#' option can be good if you are in a non-interactive mode, e.g. rendering Rmd.
#' @param layout string, one of "compact", "vertical", "horizontal", "execution"
#' @param out_format string, one of "plot", "html", "dot", "dot_print"
#'
#' - plot: directly open your viewer or browser of the plot
#' - html: save the plot to a html file
#' - dot: save the plot in DOT language, need a dot engine to remake the plot
#' - dot_print: directly cat the DOT code on console
#' @param out_path string, if the `out_format` is not "plot" or "dot_print",
#' provide a path of where to save the plot.
#' @param show_legend bool, show plot legend?
#' @param mark_main_branch bool, color the main branch on the plot?
#' @param width string, a valid CSS string for width, like "500px", "100%"
#' @param height string, a valid CSS string for height, like "500px", "100%"
#' @param elementId string, optional ID value for the plot
#' @param responsive bool, should the plot be responsive? useful in Rstudio built-in
#' viewer, Rmarkdown, Shiny or embed it into other web pages.
#' @param no_plot bool, if you want to assgin the plot to a variable and do not want
#' to see it interactively, change this to `FALSE`
#' @param plot_method string, one of "svg", "png", how to make plot, use svg or png
#' to embed the plot.
#' @param rstudio bool, if you are using Rstudio, open the built-in viewer to see the
#' plot? Default is no, open the browser tab to see it plot. The default viewer is
#' too small to see the full plot clearly, so we recommend to use the broswer tab.
#' @param rmarkdown are you rendering this plot in a Rmarkdown document? default value is
#'  "detect", this function will determine based on current R environment, or you
#' can force it to be `TRUE` or `FALSE`
#' @param in_log bool, is this plot been made in a SPR log file? If `TRUE` will add
#' links of steps to the corresponding log sections.
#' @param verbose bool, turn on verbose mode will give you more information
#' @param exit_point numeric, for advanced debugging only, see details
#' @export
#' @return see `out_format` and `exit_point`
#' @details
#' #### layout
#' - compact: try to plot steps as close as possible.
#' - vertical: main branch will be placed vertically and side branches will be placed
#'   on the same horizontal level and sub steps of side branches will be placed
#'   vertically.
#' - horizontal: main branch is placed horizontally and side branches and sub
#'   steps will be placed vertically.
#' - execution: a linear plot to show the workflow execution order of all steps.
#'
#' #### exit_point
#' return intermediate results at different points and exit the function
#'
#' - 0: no early exit
#' - 1: after all branches are found, return tree
#' - 2: after the new tree has been built, return new nodes
#' - 3: after dot translation, return graph string
#'
#' #### Rmarkdown
#' Rmarkdown will change some of the format and cause conflicts. If the plot can be
#' rendered outside Rmd but cannot within Rmd, try to turn this option on. Some additional
#' javascript processing will be performed to avoid the conflict but may cause unknown
#' issues.
plotWF <- function(sysargs,
                   width = NULL, height = NULL,
                   elementId = NULL,
                   responsive = TRUE,
                   branch_method = "auto",
                   branch_no = NULL,
                   layout = "compact",
                   no_plot = FALSE,
                   plot_method = "svg",
                   out_format = "plot",
                   out_path = NULL,
                   show_legend = TRUE,
                   mark_main_branch = TRUE,
                   rstudio = FALSE,
                   in_log = FALSE,
                   rmarkdown = "detect",
                   verbose = FALSE,
                   exit_point = 0) {
    if (!is.null(width)) stopifnot(is.character(width) && length(width) == 1)
    if (!is.null(height)) stopifnot(is.character(height) && length(height) == 1)
    stopifnot(is.logical(responsive) && length(responsive) == 1)
    stopifnot(is.logical(rstudio) && length(rstudio) == 1)
    stopifnot(is.character(rmarkdown) || is.logical(rmarkdown) && length(rmarkdown) == 1)
    out_format <- match.arg(out_format, c("plot", "html", "dot", "dot_print"))
    if (!out_format %in% c("plot", "dot_print")) stopifnot(is.character(out_path) && length(out_path) == 1)
    plot_method <- match.arg(plot_method, c("svg", "png"))
    plot_method <- switch(plot_method,
        "svg" = "renderSVGElement",
        "png" = "renderImageElement"
    )
    msg <- "" # additional msg to display on plot
    if (verbose) message("Converting SYSargsList to df...")
    if (inherits(sysargs, "data.frame")) {
        df <- sysargs
        if (!all(col_names <- c(
            "step_name", "dep", "spr", "has_run", "success",
            "sample_pass", "sample_warn", "sample_error",
            "sample_total", "log_path", "time_start", "time_end"
        ) %in%
            names(df))) {
            stop("If sysargs is a dataframe, it must contain these columns:\n", paste(col_names, collapse = ", "))
        }
        if (nrow(df) < 1) stop("plotWF: empty dataframe")
    } else if (inherits(sysargs, "SYSargsList")) {
        df <- .buildDF(sysargs)
        if (nrow(df) == 1 && df$step_name[1] == "Empty_workflow") {
            show_legend <- FALSE
            branch_method <- "auto"
        }
    } else {
        stop("`sysargs` can only be a dataframe or a SYSargsList object")
    }
    if (verbose) message("Translating to DOT format...")
    dot_vector <- makeDot(
        df, branch_method, branch_no, layout, show_legend,
        mark_main_branch, in_log, verbose, exit_point, msg
    )
    dot <- dot_vector[1]
    msg <- dot_vector[2]
    dot <- gsub(x = dot, "'", "\"")
    # if exit point
    if (exit_point > 0) {
        return(dot)
    }
    # if dot or dot_print
    if (out_format == "dot") {
        return(writeLines(dot, out_path))
    }
    if (out_format == "dot_print") {
        return(cat(dot))
    }
    # Decide if in Rmarkdown rendering
    if (is.character(rmarkdown) && rmarkdown != "detect") stop("rmarkdown can only be 'detect', TRUE or FALSE")
    if (rmarkdown == "detect") rmarkdown <- isTRUE(getOption("knitr.in.progress"))
    # forward options using x
    if (verbose) message("Making the plot...")
    x <- list(
        dot = dot,
        plotid = paste0("sprwf-", paste0(sample(8), collapse = "")),
        responsive = responsive,
        width = width,
        height = height,
        plot_method = plot_method,
        rmd = rmarkdown,
        msg = msg
    )
    # create widget
    grviz <- htmlwidgets::createWidget(
        name = "plotwf",
        x,
        width = width,
        height = height,
        package = "systemPipeR",
        elementId = elementId
    )
    # if html out
    if (out_format == "html") {
        return(htmlwidgets::saveWidget(widget = grviz, file = out_path, selfcontained = TRUE))
    }
    if (no_plot) {
        return(invisible(grviz))
    }
    # force to open browser tab instead of viewer in Rstudio
    if ((!rstudio || Sys.getenv("RSTUDIO") != "1") && !rmarkdown) {
        viewer <- getOption("viewer")
        on.exit(options(viewer = viewer), add = TRUE)
        options(viewer = NULL)
        return(print(grviz))
    } else {
        return(grviz)
    }
}

# Shiny binding's may not be needed at this moment, uncomment @export if we want to export


#' Shiny bindings for plotwf
#'
#' Output and render functions for using plotwf within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{'100\%'},
#'   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#'   string and have \code{'px'} appended.
#' @param expr An expression that generates a plotwf
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#'   is useful if you want to save an expression in a variable.
#'
#' @name plotwf-shiny
#'
#' @export
plotwfOutput <- function(outputId, width = "100%", height = "400px") {
    htmlwidgets::shinyWidgetOutput(outputId, "plotwf", width, height, package = "systemPipeR")
}

#' @rdname plotwf-shiny
#' @export
renderPlotwf <- function(expr, env = parent.frame(), quoted = FALSE) {
    if (!quoted) {
        expr <- substitute(expr)
    } # force quoted
    htmlwidgets::shinyRenderWidget(expr, plotwfOutput, env, quoted = TRUE)
}

#' Translate a workflow structure and status into a dot language string----
makeDot <- function(df,
                    branch_method = "auto",
                    branch_no = NULL,
                    layout = "compact",
                    show_legend = TRUE,
                    mark_main_branch = TRUE,
                    in_log = FALSE,
                    verbose = FALSE,
                    exit_point = 0,
                    msg = "") {
    # check
    stopifnot(is.logical(verbose) && length(verbose) == 1)
    stopifnot(is.logical(show_legend) && length(show_legend) == 1)
    if (verbose) message("Workflow inputs pre-checking ...")
    stopifnot(is.character(df$step_name))
    stopifnot(is.logical(in_log) && length(in_log) == 1)
    # early exit for linear method
    layout <- match.arg(layout, c("compact", "vertical", "horizontal", "execution"))
    if (layout == "execution") {
        return(.WFlinear(df, show_legend, in_log))
    }
    stopifnot(is.list(df$dep))
    lapply(seq_along(df$dep), function(i) {
        if (!is.character(df$dep[[i]])) stop("No.", i, " item in dep is not a character vector")
    })
    branch_method <- match.arg(branch_method, c("auto", "choose"))
    stopifnot(is.numeric(exit_point) && length(exit_point) == 1)
    # find all possible branches
    if (verbose) message("Looking for possible branches ...")
    step_names <- df$step_name
    deps <- df$dep
    if (sum(starting_root <- df$dep == "") > 1) {
        message(
            "More than 1 step has no dependency. They all will be treated as starting point:\n",
            paste0(names(df$dep)[starting_root], collapse = ", ")
        )
        step_names <- c("root_step0", step_names)
        deps <- df$dep
        deps[starting_root] <- "root_step0"
        deps <- append(deps, list(root_step0 = ""), after = 0)
    }
    tree <- .findBranch(step_names, deps)
    if (sum(starting_root) > 1) tree <- lapply(tree, function(x) x[!x == "root_step0"])
    if (exit_point == 1) {
        return(tree)
    }
    if (branch_method == "choose" && interactive() && is.null(branch_no) || verbose) {
        # list all branches
        invisible(lapply(seq_along(tree), function(i) {
            cat(crayon::blue$bold("Possible branch"), i, ":\n")
            cat(paste0(tree[[i]], collapse = " -> "), "\n")
        }))
    }
    # choose branch
    if (!is.null(branch_no)) {
        stopifnot(is.numeric(branch_no) && length(branch_no) == 1)
        if (branch_no > length(tree)) stop("Branch number is larger than possible branches:", length(tree))
    } else if (branch_method == "choose" && interactive()) {
        branch_no <- as.numeric(menu(paste("Branch", seq_along(tree)), title = "Choose a main branch to plot workflow"))
    } else {
        branch_no <- .recommendBranch(tree, df$step_name, verbose)
        if (!is.null(names(branch_no))) {
            msg <- names(branch_no)[1]
            df <- df[df$step_name %in% tree[[branch_no]], ]
        }
    }
    if (verbose) message("Build the workflow tree ...")
    nodes <- .buildTree(tree, branch_no)
    if (exit_point == 2) {
        return(nodes)
    }
    # organize all node attach point
    node_attach <- unique(unlist(lapply(nodes, `[[`, "attach_point")))
    # organize all side branch
    node_side <- unique(unlist(lapply(nodes, `[[`, "branch_side")))
    # connecting the main branch
    branch_trans <- paste(tree[[branch_no]], collapse = " -> ")
    # translation
    trans_func <- switch(layout,
        "compact" = .WFcompact,
        "vertical" = .WFvert,
        "horizontal" = .WFhori
    )
    if (verbose) message("Translate tree  to DOT language...")
    # build the skeleton
    p_main <- trans_func(branch_trans, node_attach, node_side, mark_main_branch)
    if (exit_point == 3) {
        return(paste0(p_main, "\n}\n"))
    }
    # add node decoration
    p_main <- paste0(p_main, .addNodeDecor(
        df$step_name, df$has_run, df$success, df$spr, df$sample_pass,
        df$sample_warn, df$sample_error, df$sample_total, df$log_path,
        df$time_start, df$time_end, in_log
    ))
    # add legend
    if (show_legend) p_main <- paste0(p_main, .addDotLegend(mark_main_branch), collapse = "\n")
    # close the plot
    # return plot and additional msg
    c(paste0(p_main, "\n}\n"), msg)
}

#' internal func to find all dep branches
#'
#' @param step_n numeric, current step number, starting from 1
#' @param steps string vector, all steps in a WF
#' @param deps a list of string vector, deps of steps
#' @param dep_chain a list to store all branches
#' @param chain_num current branch number, starting from 1
.findBranch <- function(steps, deps, step_n = 1,
                        dep_chain = list(
                            root = steps[1]
                        ),
                        branch_name = "root") {
    if (step_n > length(steps)) {
        return(dep_chain)
    }
    current_step <- steps[step_n]
    # find steps depend on this step
    new_step_n <- lapply(deps, function(x) {
        current_step %in% x
    }) %>%
        unlist() %>%
        which()
    # cat(current_step, new_step_n, "\n")
    # make a backup if new branch need to be created
    if (length(new_step_n) > 1) deps_backup <- dep_chain[[branch_name]]
    # loop through these new steps to append dep chain
    for (i in seq_along(new_step_n)) {
        # if more than 1 dependent step, create a new branch
        if (i != 1) {
            # cat("branch step is: ", new_step_n[i], "\n")
            branch_name <- paste0(step_n, "_", new_step_n[i])
            dep_chain[[branch_name]] <- deps_backup
            # dep_chain[[as.character(chain_num)]][["branch_point"]] <- current_step
            # cat("chain no is: ", branch_name, "\n")
        }
        # cat("chain num is", chain_num, "\n")
        dep_chain[[branch_name]] <- c(dep_chain[[branch_name]], steps[new_step_n[i]])
        dep_chain <- .findBranch(steps, deps, new_step_n[i], dep_chain, branch_name)
    }
    dep_chain
}

#' Find recommended branch
#'
#' @param tree list, tree return from [.findBranch]
#' @param steps string vector, step names
.recommendBranch <- function(tree, steps, verbose) {
    branch_complete <- lapply(tree, function(x) {
        all(c(steps[1], steps[length(steps)]) %in% x)
    }) %>%
        unlist() %>%
        which()
    if (length(branch_complete) == 0) {
        msg <- "Workflow's first step is not connected to the last step, something wrong? Unconnected steps will not be plotted."
        warning(msg)
        return(structure(c(1), .Names = msg))
    } else {
        if (verbose) cat("**********\n")
        if (verbose) {
            cat(
                "Find", length(branch_complete), "branch(es) connecting first and last step:",
                paste0(branch_complete, collapse = ", "),
                ".\n"
            )
        }
        tree_complete <- tree[branch_complete]
    }
    branch_len <- unlist(lapply(tree, length))
    branch_long <- which(branch_len == max(branch_len))
    if (verbose) cat("Find branch(es)", paste0(branch_long, collapse = ", "), ": with the largest number of steps", max(branch_len), "\n")
    branch_recommand <- intersect(branch_long, branch_complete)
    # use first complete branch if no intersection
    if (length(branch_recommand) == 0) {
        branch_recommand <- branch_complete[1]
        # if still empty, use longest branch
        if (length(branch_recommand) == 0) {
            branch_recommand <- branch_long[1]
        }
    }
    if (verbose) cat("**********\n")
    if (verbose) cat("Based on the detection, branch(es):", paste0(branch_recommand, collapse = ", "), "is (are) recommended\n")
    if (verbose) cat("Branch", crayon::yellow$bold(branch_recommand[1]), "will be used as the main branch\n")
    return(branch_recommand[1])
}

#' Build the new tree based on the chosen main branch
#' @param tree list of branches
#' @param branch_no numeric, which branch
.buildTree <- function(tree, branch_no) {
    branch <- tree[[branch_no]]
    # reconstruct side branches based on the selected main branch
    lapply(tree[-branch_no], function(x) {
        # see if other branches has node in common with main branch
        in_main <- x %in% branch
        # find where these node need to attach to mian branch
        current_state <- in_main[1]
        attach_point <- c()
        for (i in seq_along(in_main)) {
            if (current_state != in_main[i]) {
                attach_point <- c(attach_point, paste(x[i - 1], "->", x[i]))
                current_state <- in_main[i]
            }
        }
        # find the chain of new side branches
        in_branch <- which(!in_main)
        branch_side <- list(x[in_branch[1]])
        list_pos <- 1
        for (i in seq_along(in_branch)) {
            if (i == 1) next
            if (in_branch[i] - in_branch[i - 1] == 1) {
                branch_side[[list_pos]] <- c(branch_side[[list_pos]], x[in_branch[i]])
            } else {
                list_pos <- list_pos + 1
                branch_side[[list_pos]] <- x[in_branch[i]]
            }
        }
        branch_side <- lapply(branch_side, function(branches) {
            nodes <- c()
            for (n in seq_along(branches)) {
                if (n >= length(branches)) next
                nodes <- c(nodes, paste(branches[n], branches[n + 1], sep = " -> "))
            }
            nodes
        })
        list(attach_point = attach_point, branch_side = branch_side)
    })
}

# different layout functions -----

.WFcompact <- function(branch_trans, node_attach, node_side, mark_main_branch) {
    paste0(
        "digraph {
    node[fontsize=20];
    subgraph {\n    ",
        if (mark_main_branch) '    node[color="dodgerblue"];\n        ' else "    ",
        paste0(branch_trans, if (mark_main_branch) '[color="dodgerblue"]' else "", collapse = ""),
        "\n   }\n    ",
        paste0(node_attach, collapse = "\n    "),
        "\n    ",
        paste0(node_side, collapse = "\n    "),
        "\n"
    )
}

# .WFcompact(branch_trans, node_attach, node_side) %>% cat

.WFvert <- function(branch_trans, node_attach, node_side, mark_main_branch) {
    paste0(
        'digraph {
    node[fontsize=20];
    subgraph {
        rankdir="TB";\n        ',
        if (mark_main_branch) '    node[color="dodgerblue"];\n        ' else "    ",
        paste0(branch_trans, if (mark_main_branch) '[color="dodgerblue"]' else "", collapse = ""),
        "\n   }\n    ",
        paste0(
            "subgraph {\n",
            '        rank="same";\n',
            "        ",
            node_attach,
            "\n    }\n"
        ) %>%
            paste0(collapse = "\n    "),
        "\n    ",
        paste0(node_side, collapse = "\n    "),
        "\n"
    )
}

# .WFvert(branch_trans, node_attach, node_side) %>% cat

.WFhori <- function(branch_trans, node_attach, node_side, mark_main_branch) {
    paste0(
        'digraph {
    node[fontsize=20];
    subgraph {
        rank="same";\n                ',
        if (mark_main_branch) '    node[color="dodgerblue"];\n        ' else "    ",
        paste0(branch_trans, if (mark_main_branch) '[color="dodgerblue"]' else "", collapse = ""),
        "\n   }\n    ",
        paste0(
            "subgraph {\n",
            '        rankdir="TB";\n',
            "        ",
            node_attach,
            "\n    }\n"
        ) %>%
            paste0(collapse = "\n    "),
        "\n    ",
        paste0(node_side, collapse = "\n    "),
        "\n"
    )
}

# .WFhori(branch_trans, node_attach, node_side) %>% cat

.WFlinear <- function(df, show_legend, in_log) {
    steps_trans <- paste0(df$step_name, collapse = " -> ")
    paste0(
        'digraph {
    node[fontsize=20];
    subgraph {
        rank="TB";\n        ',
        steps_trans,
        "\n   }\n",
        .addNodeDecor(
            df$step_name, df$has_run, df$success, df$spr, df$sample_pass,
            df$sample_warn, df$sample_error, df$sample_total, df$log_path,
            df$time_start, df$time_end, in_log
        ),
        if (show_legend) .addDotLegend(FALSE),
        "\n}\n"
    )
}

# .WFlinear(steps) %>% cat

# add legend -----------
#' @param show_main show main steps legend?
.addDotLegend <- function(show_main = TRUE) {
    paste0(
        '    subgraph cluster_legend {
        rankdir=TB;
        color="#EEEEEE";
        style=filled;
        node [style=filled];',
        if (show_main) {
            '
        {rank=same; R_step; Sysargs_step; Main_branch}
        Main_branch -> Sysargs_step -> R_step[color="#EEEEEE"]
        Main_branch[label="Main branch" color="dodgerblue", style="filled", fillcolor=white];'
        } else {
            "
        {rank=same; R_step; Sysargs_step}"
        },
        '   Sysargs_step ->step_state[color="#EEEEEE"];
        step_state[style="filled", shape="box" color=white, label =<
            <table>
            <tr><td><b>Step Colors</b></td></tr>
            <tr><td><font color="black">Pending steps</font>; <font color="#5cb85c">Successful steps</font>; <font color="#d9534f">Failed steps</font></td></tr>
            <tr><td><b>Targets Files / Code Chunk </b></td></tr><tr><td><font color="#5cb85c">0 (pass) </font> | <font color="#f0ad4e">0 (warning) </font> | <font color="#d9534f">0 (error) </font> | <font color="blue">0 (total)</font>; Duration</td></tr></table>
            >];
        label="Legends";
        fontsize = 30;
        Sysargs_step[label="Sysargs step" style="rounded, filled", shape="box", fillcolor=white];
        R_step[label="R step" style="rounded, filled", fillcolor=white];
    }\n'
    )
}


# add node colors, links -------
#' @param steps step names
#' @param has_run bool, steps has run?
#' @param success bool, steps successful?
#' @param spr string, SPR option, 'sys' or 'r'
#' @param sample_pass numeric, no. of samples passed each step
#' @param sample_warn numeric, no. of samples have warnings each step
#' @param sample_error numeric, no. of samples have errors each step
#' @param sample_total numeric, no. of samples total each step
.addNodeDecor <- function(steps, has_run, success, spr, sample_pass, sample_warn, sample_error, sample_total, log_path, time_start, time_end, in_log = FALSE) {
    node_text <- c()
    for (i in seq_along(steps)) {
        step_color <- if (has_run[i] && success[i]) "#5cb85c" else if (has_run[i] && !success[i]) "#d9534f" else "black"
        duration <- .stepDuration(time_start[i], time_end[i])
        node_text <- c(node_text, paste0(
            "    ", steps[i], "[label=<<b>",
            '<font color="', step_color, '">', steps[i], "</font><br></br>",
            '<font color="#5cb85c">', sample_pass[i], "</font>/",
            '<font color="#f0ad4e">', sample_warn[i], "</font>/",
            '<font color="#d9534f">', sample_error[i], "</font>/",
            '<font color="blue">', sample_total[i], "</font></b>; ",
            '<font color="black">', duration$short, "</font>",
            "> ",
            if (spr[i] == "sysargs") ', style="rounded", shape="box" ' else "",
            if (in_log) paste0('href="', log_path[i], '" ') else " ",
            'tooltip="step ', steps[i], ": ",
            sample_pass[i], " samples passed; ",
            sample_warn[i], " samples have warnings; ",
            sample_error[i], " samples have errors; ",
            sample_total[i], " samples in total; ",
            "Start time: ", time_start[i], "; ",
            "End time: ", time_end[i], "; ",
            "Duration: ", duration$long,
            '"',
            "]\n"
        ))
    }
    paste0(node_text, collapse = "")
}

# figure the right unit to display for duration time.

#' inputs should be `Sys.time()` timestamp
#' @return a list of duration in short format and long format
.stepDuration <- function(start, end) {
    duration <- round(as.numeric(difftime(end, start, units = "sec")), 1)
    if (duration < 0) stop("Ending time is smaller than starting time, something wrong?")
    secs <- duration %% 60
    mins_raw <- duration %/% 60
    mins <- mins_raw %% 60
    hrs <- mins_raw %/% 60
    return(list(
        long = paste(stringr::str_pad(hrs, 2, pad = "0"), stringr::str_pad(mins, 2, pad = "0"), stringr::str_pad(secs, 2, pad = "0"), sep = ":"),
        short = if (hrs > 0) paste0(hrs, "h") else if (mins > 0) paste0(mins, "m") else paste0(secs, "s")
    ))
}

# build df from sal object

.buildDF <- function(sal) {
    sal_temp <- sal
    if (length(sal_temp) == 0) {
        warning("Workflow has no steps. Please make sure to add a step to the workflow before plotting.", call. = FALSE)
        return_df <- data.frame(
            step_name = "Empty_workflow",
            dep = NA,
            spr = "sysargs",
            has_run = FALSE,
            success = FALSE,
            sample_pass = 0,
            sample_warn = 0,
            sample_error = 0,
            sample_total = 0,
            log_path = "",
            time_start = Sys.time(),
            time_end = Sys.time()
        )
        return_df$dep <- list("")
        return(return_df)
    }
    df <- data.frame(step_name = stepName(sal_temp))
    dep <- dependency(sal_temp)
    for (i in seq_along(dep)) {
        if (any(is.na(dep[i]))) {
            dep[[i]] <- ""
        }
    }
    df$dep <- dep
    df$spr <- ifelse(sapply(df$step_name, function(x) inherits(stepsWF(sal_temp)[[x]], "SYSargs2")), "sysargs", "r")
    df$has_run <- ifelse(!sapply(df$step_name, function(x) sal_temp$statusWF[[x]]$status.summary) == "Pending", TRUE, FALSE)
    df$success <- ifelse(sapply(df$step_name, function(x) sal_temp$statusWF[[x]]$status.summary) == "Success", TRUE, FALSE)
    df <- cbind(df, data.frame(
        sample_pass = 0,
        sample_warn = 0,
        sample_error = 0,
        sample_total = 0,
        time_start = Sys.time(),
        time_end = Sys.time()
    ))
    for (i in seq_along(df$step_name)) {
        if (inherits(stepsWF(sal_temp)[[i]], "SYSargs2")) {
            sample_df <- as.list(colSums(sal_temp$statusWF[[i]][[2]][2:4]))
            df$sample_pass[i] <- sample_df$Existing_Files
            df$sample_total[i] <- sample_df$Total_Files
            if (all(sample_df$Missing_Files > 0 && sal_temp$statusWF[[i]]$status.summary == "Warning")) {
                df$sample_warn[i] <- sample_df$Missing_Files
            } else if (all(sample_df$Missing_Files > 0 && sal_temp$statusWF[[i]]$status.summary == "Error")) {
                df$sample_error[i] <- sample_df$Missing_Files
            }
            if (length(sal_temp$statusWF[[i]]$status.time) > 0) {
                df$time_start[i] <- sal_temp$statusWF[[i]]$status.time$time_start[1]
                df$time_end[i] <- sal_temp$statusWF[[i]]$status.time$time_end[length(sal_temp$statusWF[[i]]$status.time$time_end)]
            }
        } else if (inherits(stepsWF(sal_temp)[[i]], "LineWise")) {
            df$sample_total[i] <- 1
            if (sal_temp$statusWF[[i]]$status.summary == "Success") {
                df$sample_pass[i] <- 1
            } else if (sal_temp$statusWF[[i]]$status.summary == "Warning") {
                df$sample_warn[i] <- 1
            } else if (sal_temp$statusWF[[i]]$status.summary == "Error") {
                df$sample_error[i] <- 1
            }
            if (length(sal_temp$statusWF[[i]]$status.time) > 0) {
                df$time_start[i] <- sal_temp$statusWF[[i]]$status.time$time_start
                df$time_end[i] <- sal_temp$statusWF[[i]]$status.time$time_end
            }
        }
    }
    df$log_path <- paste0("#", tolower(df$step_name))
    df <- rbind(df[which(df$dep == ""), ], df[which(!df$dep == ""), ])
    return(df)
}
