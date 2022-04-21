#' Visualize SPR workflow and status
#'
#' @param sysargs object of class `SYSargsList`.
#' @param width string, a valid CSS string for width, like "500px", "100%"
#' @param height string, a valid CSS string for height, like "500px", "100%"
#' @param elementId string, optional ID value for the plot
#' @param responsive bool, should the plot be responsive? useful in Rstudio built-in
#' viewer, Rmarkdown, Shiny or embed it into other web pages.
#' @param branch_method string, one of "auto", "choose". How to determine the
#' main branch of the workflow. "auto" will be determined by internal alrgothrim:
#' Branches connecting the frist and last step and/or the longest will be favored.
#' "choose" will list all possible branches and you can make a choice.
#' @param branch_no numeric, only works if `branch_method == "choose"`. Specify a
#' branch number to be the main branch instead of choosing from the prompt. This
#' option can be good if you are in a non-interactive mode, e.g. rendering Rmd.
#' @param layout string, one of "compact", "vertical", "horizontal", "execution"
#' @param no_plot bool, if you want to assgin the plot to a variable and do not want
#' to see it interactively, change this to `FALSE`
#' @param plot_method string, one of "svg", "png", how to make plot, use svg or png
#' to embed the plot.
#' @param out_format string, one of "plot", "html", "dot", "dot_print"
#' - plot: directly open your viewer or browser of the plot
#' - html: save the plot to a html file
#' - dot: save the plot in DOT language, need a dot engine to remake the plot
#' - dot_print: directly cat the DOT code on console
#' @param out_path string, if the `out_format` is not "plot" or "dot_print",
#' provide a path of where to save the plot.
#' @param show_legend bool, show plot legend?
#' @param mark_main_branch bool, color the main branch on the plot?
#' @param rstudio bool, if you are using Rstudio, open the built-in viewer to see the
#' plot? Default is no, open the browser tab to see it plot. The default viewer is
#' too small to see the full plot clearly, so we recommend to use the browser tab.
#' However, if you are using this plot in Shiny apps, always turn
#' `rstudio = TRUE`.
#' @param in_log bool, is this plot been made in a SPR log file? If `TRUE` will add
#' links of steps to the corresponding log sections.
#' @param rmarkdown are you rendering this plot in a Rmarkdown document? default value is
#'  "detect", this function will determine based on current R environment, or you
#' can force it to be `TRUE` or `FALSE`
#' @param verbose bool, turn on verbose mode will give you more information
#' @param show_warns bool, print the warning messages on the plot?
#' @param plot_ctr bool, add the plot control panel to the plot? This requires you
#' to have internet connection. It will download some additional javascript libraries,
#' and allow you to save the plot as png, jpg, svg, pdf or graphviz directly from the
#' browser.
#' @param pan_zoom bool, allow panning and zooming of the plot? Use mouse wheel
#' or touch pad to zoom in and out of the plot. You need to have
#' internet connection, additional javascript libraries will be loaded automatically
#' online. Cannot be used with `responsive = TRUE` together. If both `TRUE`,
#' `responsive` will be automatically set to `FALSE`.
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
                   mark_main_branch = FALSE,
                   rstudio = FALSE,
                   in_log = FALSE,
                   rmarkdown = "detect",
                   verbose = FALSE,
                   show_warns = FALSE,
                   plot_ctr = TRUE,
                   pan_zoom = FALSE,
                   exit_point = 0) {
    if (!is.null(width)) stopifnot(is.character(width) && length(width) == 1)
    if (!is.null(height)) stopifnot(is.character(height) && length(height) == 1)
    stopifnot(is.logical(responsive) && length(responsive) == 1)
    stopifnot(is.logical(rstudio) && length(rstudio) == 1)
    stopifnot(is.logical(show_warns) && length(show_warns) == 1)
    stopifnot(is.logical(plot_ctr) && length(plot_ctr) == 1)
    stopifnot(is.logical(pan_zoom) && length(pan_zoom) == 1)
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
            "step_name", "dep", "spr", "req", "session", "has_run", "success",
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
        mark_main_branch, in_log, verbose, exit_point, msg,
        show_warns
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
    legend_uri <- if (rmarkdown && show_legend) .addDotLegendUri() else ""
    # forward options using x
    if (verbose) message("Making the plot...")
    if (pan_zoom && responsive) {
        warning("Pan-zoom and responsive cannot be used together. Pan-zoom has priority, now `responsive` has set to FALSE")
        responsive <- FALSE
    }
    x <- list(
        dot = dot,
        plotid = paste0("sprwf-", paste0(sample(8), collapse = "")),
        responsive = responsive,
        width = width,
        height = height,
        plot_method = plot_method,
        rmd = rmarkdown,
        msg = msg,
        plot_ctr = plot_ctr,
        pan_zoom = pan_zoom,
        legend_uri = legend_uri
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
#' @details
#' To use `plotWF` in `renderPlotwf` in Shiny apps, always turn on the option `rstudio = TRUE`.
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
                    msg = "",
                    show_warns = TRUE) {
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
        if(show_warns) message(
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
    # debug exit point
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
        branch_no <- .recommendBranch(tree, df$step_name, verbose, show_warns)
        branch_msg <- names(branch_no)[1]
        if(stringr::str_starts(branch_msg, "Workflow's first step is") && show_warns) {
            msg <- branch_msg
            # df <- df[df$step_name %in% tree[[branch_no]], ]
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
        df$step_name, df$status_summary, df$has_run, df$success, df$spr, df$sample_pass,
        df$sample_warn, df$sample_error, df$sample_total, df$log_path,
        df$time_start, df$time_end, df$req, df$session, in_log
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
.recommendBranch <- function(tree, steps, verbose, show_warns) {
    branch_complete <- lapply(tree, function(x) {
        all(c(steps[1], steps[length(steps)]) %in% x)
    }) %>%
        unlist() %>%
        which()
    if (length(branch_complete) == 0) {
        msg <- "Workflow's first step is not connected to the last step, something wrong? This may cause workflow plot display issues"
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
    branch_recommand <- base::intersect(branch_long, branch_complete)
    branch_recommand <- tree[branch_recommand %in% tree]
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
            df$step_name, df$status_summary, df$has_run, df$success, df$spr, df$sample_pass,
            df$sample_warn, df$sample_error, df$sample_total, df$log_path,
            df$time_start, df$time_end,  df$req, df$session, in_log
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
        '        subgraph cluster_legend {
        rankdir=TB;
        color="#eeeeee";
        style=filled;
        ranksep =1;
        label="Legends";
        fontsize = 30;
        node [style=filled, fontsize=10];
        legend_img-> step_state[color="#eeeeee"];

        legend_img[shape=none, image="plotwf_legend-src.png", label = " ", height=1, width=3, style=""];

        step_state[style="filled", shape="box" color=white, label =<
            <table>
            <tr><td><b>Step Colors</b></td></tr>
            <tr><td><font color="black">Pending steps</font>; <font color="#5cb85c">Successful steps</font>; <font color="#d9534f">Failed steps</font></td></tr>
            <tr><td><b>Targets Files / Code Chunk </b></td></tr><tr><td><font color="#5cb85c">0 (pass) </font> | <font color="#f0ad4e">0 (warning) </font> | <font color="#d9534f">0 (error) </font> | <font color="blue">0 (total)</font>; Duration</td></tr></table>
            >];

    }\n'
    )
}

.addDotLegendUri <- function() {
    # to update the URI, run
    # paste0("data:image/svg+xml;base64,", base64enc::base64encode("inst/htmlwidgets/plotwf_legend.svg"))
    "data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiPz4KPCEtLSBEbyBub3QgZWRpdCB0aGlzIGZpbGUgd2l0aCBlZGl0b3JzIG90aGVyIHRoYW4gZGlhZ3JhbXMubmV0IC0tPgo8IURPQ1RZUEUgc3ZnIFBVQkxJQyAiLS8vVzNDLy9EVEQgU1ZHIDEuMS8vRU4iICJodHRwOi8vd3d3LnczLm9yZy9HcmFwaGljcy9TVkcvMS4xL0RURC9zdmcxMS5kdGQiPgo8c3ZnIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIgeG1sbnM6eGxpbms9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGxpbmsiIHZlcnNpb249IjEuMSIgd2lkdGg9IjQ5NnB4IiBoZWlnaHQ9IjI3OHB4IiB2aWV3Qm94PSItMC41IC0wLjUgNDk2IDI3OCIgY29udGVudD0iJmx0O214ZmlsZSBob3N0PSZxdW90O2FwcC5kaWFncmFtcy5uZXQmcXVvdDsgbW9kaWZpZWQ9JnF1b3Q7MjAyMS0xMS0yNFQyMDozOTo0NC45MjNaJnF1b3Q7IGFnZW50PSZxdW90OzUuMCAoWDExOyBMaW51eCB4ODZfNjQpIEFwcGxlV2ViS2l0LzUzNy4zNiAoS0hUTUwsIGxpa2UgR2Vja28pIENocm9tZS85Ni4wLjQ2NjQuNDUgU2FmYXJpLzUzNy4zNiZxdW90OyB2ZXJzaW9uPSZxdW90OzE1LjguMyZxdW90OyBldGFnPSZxdW90O1RfZEV4bkw3U0xYMmJ1WDhQYVgzJnF1b3Q7IHR5cGU9JnF1b3Q7Z29vZ2xlJnF1b3Q7Jmd0OyZsdDtkaWFncmFtIGlkPSZxdW90O1Z3OVZWaUdzeC1zTF9OaktlN0pQJnF1b3Q7Jmd0OzdWcmJjdUk0RVAwYUhwT3lKZDk0REFSbVhtWjNLMnpWUGl0WUdGVmt5V1BMZ2N6WHIyVEwrQ0lEM3NHUVRSVWtsZGd0cVMyZDAycDF0NW5BZWJ6L2xxSmsrNE9IbUU2QUZlNG44SGtDZ0EwQW1LaGZLL3pRRXN2eFMwbVVrbERMYXNHSy9NSlZSeTNOU1lpelZrZkJPUlVrYVF2WG5ERzhGaTBaU2xPK2EzZmJjTnArYW9JaWJBaFdhMFJONlQ4a0ZOdFNPcldzV3Y0ZGsyaGJQZG1yV21KVWRkYUNiSXRDdm11STRHSUM1eW5ub3J5SzkzTk1GWG9WTHVXNDVaSFd3OFJTek1TUUFacUpkMFJ6dlRZOUwvRlJMVFpLZVo1TTRFeitZeUZXNHl4NWQ1aTR1dEZLY0Nyd3ZvOEE5Rm9wczh3SjJvZGxTNFBCUE1ZaS9aQmRLa1dCSHFKdEJVQjl2NnVCQjQ1YnlyWk4wS0duQ2Rka1J3ZmROUjd5UWtQU0R3ODhEMDhibHQyV0NMeEswRnExN3FUNVM5bFd4Rkwvc3kwdk40VFNPYWM4TGNiQzBNVkI2Q2cwUmNyZmNOWENPSlBEWnhGRldkWUhkL2FHeFhwN0R2c214cUFmWTQycEE0ZEI2b3lBcUhOZFJKZkxSZkM4T0lib1ZWRU1wZ05SOUM5SDBiMHVpb3ZsRWl6bnQwWHhMR3p3Y3RpOHNXSGpUT2pUS1RCUjNManFSOG9SSlJHVE1vbzM0aGlvU2xWamJQbVJjaTRmVG9ReXNXQkVoSDNyU2c3VE54RE9PSldOWFpqbDlFVWJ5d3FrdFZ3R2xpak0xQ0tKUEhPZmRFTk13bEFObjZVNEk3LzBrYUl3U1RoaG9waXpPNXU0ejBwWExuaFc4bUlmQmJ4QlZvT0RKcDN5ZG9saVFoWDZmNU5ZUmh6QStnUHY1TjhYSGlNMmxJN1RCeHp3L1VlM2ZjUUIwNVBBSHNLcXcvUVN2Z0tEcnhCbFczd243RGhoMEhVTndnS0RNT2RLaEUwTnduNGdKdVBOV00zK3E1Rm0yN2ZiWm03dzZIZG9jODFRc3M4eGprRmJsYmswZUp2ek9Na0Z2cE4yWXFzNXhsYnI0Y3k3Rm1kMlQ3amcwZUlJVCtRU202eDVQM09Wc2hXNFBKVDRQc2tPUWJLdjIrUlZwUDdMcDhoRUZWZTY1RFJLZFdXcllROVY5LzhVeElWNGczSXFicmU5dktETGxPOFpUTG4yWTA5a0IwWUlpRzB6a1gyNTc2d1RmQVYyaHk5bzlXU0IxOXBaWmw0dHZhRmNYUGhBQ2J1N3hGTXVFWFpkSXJSc2c3amdXc1NaNlh1cTZqc3NvaVpyZFNwbFgrSzhQaUhBQTBOOEdlajFaZTRJSUEvSTdqRUxuMVFoc3piSEJwYWRTbDFoejFXZEVoNHd3cUZSNUR5TFVIUDlQU1pXeVZKTWtTRHZiZlY5a09nbi9LWDJaUE13bVQ2Q0RnUGRxa25HODNTTjlVRFFLSEIyZEVGZ0czRmZWNWRBYVlTRm9hdmc2YkQ0WWRRTnFEQmNSSjNPenNyTy8zY2VaYWJVQWQ1N0RLeHAvWEYvajFQSE1sTG1NNXBIWk5pc2NNZ0VMRVNDNjAxK1A3ZU9CSWhkcHdvZDA2bGVMWlEzNnh4L0pvSndodWlkdEZNbm9UV0F0RDRQTWdwcFpxMWpsUHhMd1RzMCticm83UTRNUGZ4NnF1QjdpK3BpKy9VWkJHYVU3L2k5a2N3SURBS3o2akVLZzR6TGxsdlJ1Q3crbjB1akdaTGVtRWl6RlBLU00wWllORkhxU3hKZTB5R0UycFprdEtDeVErb0taNW4weUYvSUgrczVkdDRFMmRiTnpNSU9nbllZQkwzS0FzNzRaemlHVVpnMWw4b1VRdkora1MyYzlSS3ZhUDBXRmJ2NllWMHlvL1NKRkxHc2dtK21mTVRBNXgyTStZZzdPWWlMaGVreCtHZE8wcTlXZVA5OG0vWGJOZ3VCR1ZIMFZlRkhzVml6N3JRU09CblBoODJMcjJ6Y3JXRzROVXhoMjRNRlpqR3I3MVhhYjFpRHZLMi9WMVZtZnZYWDArRGlYdz09Jmx0Oy9kaWFncmFtJmd0OyZsdDsvbXhmaWxlJmd0OyI+PGRlZnMvPjxnPjxyZWN0IHg9IjQiIHk9IjkwIiB3aWR0aD0iNDkwIiBoZWlnaHQ9IjkyIiBmaWxsPSIjZDVlOGQ0IiBzdHJva2U9Im5vbmUiIHBvaW50ZXItZXZlbnRzPSJub25lIi8+PHJlY3QgeD0iNCIgeT0iMTgyIiB3aWR0aD0iNDkwIiBoZWlnaHQ9Ijk0IiBmaWxsPSIjZmZlOGRlIiBzdHJva2U9Im5vbmUiIHBvaW50ZXItZXZlbnRzPSJub25lIi8+PHJlY3QgeD0iNCIgeT0iNCIgd2lkdGg9IjQ5MCIgaGVpZ2h0PSI4NiIgZmlsbD0iI2VmZjJmYyIgc3Ryb2tlPSJub25lIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxyZWN0IHg9IjQiIHk9IjQiIHdpZHRoPSIxNDAiIGhlaWdodD0iMjcyIiBmaWxsLW9wYWNpdHk9IjAuOCIgZmlsbD0iI2Y1ZjVmNSIgc3Ryb2tlPSJub25lIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiAxMXB4OyBtYXJnaW4tbGVmdDogMTE1cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogOHB4OyBmb250LWZhbWlseTogJnF1b3Q7VGltZXMgTmV3IFJvbWFuJnF1b3Q7OyBjb2xvcjogcmdiKDAsIDAsIDApOyBsaW5lLWhlaWdodDogMS4yOyBwb2ludGVyLWV2ZW50czogbm9uZTsgd2hpdGUtc3BhY2U6IG5vd3JhcDsiPnNvbGlkPC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjExNSIgeT0iMTMiIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iOHB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5zb2xpZDwvdGV4dD48L3N3aXRjaD48L2c+PGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTAuNSAtMC41KXNjYWxlKDIpIj48c3dpdGNoPjxmb3JlaWduT2JqZWN0IHBvaW50ZXItZXZlbnRzPSJub25lIiB3aWR0aD0iMTAwJSIgaGVpZ2h0PSIxMDAlIiByZXF1aXJlZEZlYXR1cmVzPSJodHRwOi8vd3d3LnczLm9yZy9UUi9TVkcxMS9mZWF0dXJlI0V4dGVuc2liaWxpdHkiIHN0eWxlPSJvdmVyZmxvdzogdmlzaWJsZTsgdGV4dC1hbGlnbjogbGVmdDsiPjxkaXYgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGh0bWwiIHN0eWxlPSJkaXNwbGF5OiBmbGV4OyBhbGlnbi1pdGVtczogdW5zYWZlIGNlbnRlcjsganVzdGlmeS1jb250ZW50OiB1bnNhZmUgY2VudGVyOyB3aWR0aDogMXB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDEwcHg7IG1hcmdpbi1sZWZ0OiAxOThweDsiPjxkaXYgZGF0YS1kcmF3aW8tY29sb3JzPSJjb2xvcjogcmdiYSgwLCAwLCAwLCAxKTsgIiBzdHlsZT0iYm94LXNpemluZzogYm9yZGVyLWJveDsgZm9udC1zaXplOiAwcHg7IHRleHQtYWxpZ246IGNlbnRlcjsiPjxkaXYgc3R5bGU9ImRpc3BsYXk6IGlubGluZS1ibG9jazsgZm9udC1zaXplOiA4cHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyB3aGl0ZS1zcGFjZTogbm93cmFwOyI+ZGFzaGVkPC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjE5OCIgeT0iMTIiIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iOHB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5kYXNoZWQ8L3RleHQ+PC9zd2l0Y2g+PC9nPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiAzMnB4OyBtYXJnaW4tbGVmdDogMTE2cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTFweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3dyYXA7Ij5NYW5hZ2VtZW50PC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjExNiIgeT0iMzUiIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iMTFweCIgdGV4dC1hbmNob3I9Im1pZGRsZSI+TWFuYWdlbWVudDwvdGV4dD48L3N3aXRjaD48L2c+PGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTAuNSAtMC41KXNjYWxlKDIpIj48c3dpdGNoPjxmb3JlaWduT2JqZWN0IHBvaW50ZXItZXZlbnRzPSJub25lIiB3aWR0aD0iMTAwJSIgaGVpZ2h0PSIxMDAlIiByZXF1aXJlZEZlYXR1cmVzPSJodHRwOi8vd3d3LnczLm9yZy9UUi9TVkcxMS9mZWF0dXJlI0V4dGVuc2liaWxpdHkiIHN0eWxlPSJvdmVyZmxvdzogdmlzaWJsZTsgdGV4dC1hbGlnbjogbGVmdDsiPjxkaXYgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGh0bWwiIHN0eWxlPSJkaXNwbGF5OiBmbGV4OyBhbGlnbi1pdGVtczogdW5zYWZlIGNlbnRlcjsganVzdGlmeS1jb250ZW50OiB1bnNhZmUgY2VudGVyOyB3aWR0aDogMXB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDMycHg7IG1hcmdpbi1sZWZ0OiAxOThweDsiPjxkaXYgZGF0YS1kcmF3aW8tY29sb3JzPSJjb2xvcjogcmdiYSgwLCAwLCAwLCAxKTsgIiBzdHlsZT0iYm94LXNpemluZzogYm9yZGVyLWJveDsgZm9udC1zaXplOiAwcHg7IHRleHQtYWxpZ246IGNlbnRlcjsiPjxkaXYgc3R5bGU9ImRpc3BsYXk6IGlubGluZS1ibG9jazsgZm9udC1zaXplOiAxMXB4OyBmb250LWZhbWlseTogJnF1b3Q7VGltZXMgTmV3IFJvbWFuJnF1b3Q7OyBjb2xvcjogcmdiKDAsIDAsIDApOyBsaW5lLWhlaWdodDogMS4yOyBwb2ludGVyLWV2ZW50czogbm9uZTsgd2hpdGUtc3BhY2U6IG5vd3JhcDsiPkNvbXB1dGU8L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMTk4IiB5PSIzNSIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMXB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5Db21wdXRlPC90ZXh0Pjwvc3dpdGNoPjwvZz48ZWxsaXBzZSBjeD0iMjMyLjUiIGN5PSIxMjMiIHJ4PSI1MS41IiByeT0iMjciIGZpbGw9InJnYmEoMjU1LCAyNTUsIDI1NSwgMSkiIHN0cm9rZT0icmdiYSgwLCAwLCAwLCAxKSIgc3Ryb2tlLXdpZHRoPSIyIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDUwcHg7IGhlaWdodDogMXB4OyBwYWRkaW5nLXRvcDogNjJweDsgbWFyZ2luLWxlZnQ6IDkycHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTJweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3JtYWw7IG92ZXJmbG93LXdyYXA6IG5vcm1hbDsiPjxzcGFuIHN0eWxlPSJmb250LXNpemU6IDhweCI+ZWxsaXBzZTwvc3Bhbj48L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMTE2IiB5PSI2NSIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMnB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5lbGxpcHNlPC90ZXh0Pjwvc3dpdGNoPjwvZz48ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSgtMC41IC0wLjUpc2NhbGUoMikiPjxzd2l0Y2g+PGZvcmVpZ25PYmplY3QgcG9pbnRlci1ldmVudHM9Im5vbmUiIHdpZHRoPSIxMDAlIiBoZWlnaHQ9IjEwMCUiIHJlcXVpcmVkRmVhdHVyZXM9Imh0dHA6Ly93d3cudzMub3JnL1RSL1NWRzExL2ZlYXR1cmUjRXh0ZW5zaWJpbGl0eSIgc3R5bGU9Im92ZXJmbG93OiB2aXNpYmxlOyB0ZXh0LWFsaWduOiBsZWZ0OyI+PGRpdiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMTk5OS94aHRtbCIgc3R5bGU9ImRpc3BsYXk6IGZsZXg7IGFsaWduLWl0ZW1zOiB1bnNhZmUgY2VudGVyOyBqdXN0aWZ5LWNvbnRlbnQ6IHVuc2FmZSBjZW50ZXI7IHdpZHRoOiAxcHg7IGhlaWdodDogMXB4OyBwYWRkaW5nLXRvcDogODVweDsgbWFyZ2luLWxlZnQ6IDExNHB4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDExcHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyB3aGl0ZS1zcGFjZTogbm93cmFwOyI+UjwvZGl2PjwvZGl2PjwvZGl2PjwvZm9yZWlnbk9iamVjdD48dGV4dCB4PSIxMTQiIHk9Ijg4IiBmaWxsPSJyZ2JhKDAsIDAsIDAsIDEpIiBmb250LWZhbWlseT0iVGltZXMgTmV3IFJvbWFuIiBmb250LXNpemU9IjExcHgiIHRleHQtYW5jaG9yPSJtaWRkbGUiPlI8L3RleHQ+PC9zd2l0Y2g+PC9nPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiA4M3B4OyBtYXJnaW4tbGVmdDogMTk4cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTFweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3dyYXA7Ij5Db21tYW5kLWxpbmU8L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMTk4IiB5PSI4NiIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMXB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5Db21tYW5kLWxpbmU8L3RleHQ+PC9zd2l0Y2g+PC9nPjxyZWN0IHg9IjM0OSIgeT0iOTYiIHdpZHRoPSIxMDUiIGhlaWdodD0iNTAiIHJ4PSI3LjUiIHJ5PSI3LjUiIGZpbGw9InJnYmEoMjU1LCAyNTUsIDI1NSwgMSkiIHN0cm9rZT0icmdiYSgwLCAwLCAwLCAxKSIgc3Ryb2tlLXdpZHRoPSIyIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDUxcHg7IGhlaWdodDogMXB4OyBwYWRkaW5nLXRvcDogNjFweDsgbWFyZ2luLWxlZnQ6IDE3NnB4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDhweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3JtYWw7IG92ZXJmbG93LXdyYXA6IG5vcm1hbDsiPnJlY3RhbmdsZTwvZGl2PjwvZGl2PjwvZGl2PjwvZm9yZWlnbk9iamVjdD48dGV4dCB4PSIyMDEiIHk9IjYzIiBmaWxsPSJyZ2JhKDAsIDAsIDAsIDEpIiBmb250LWZhbWlseT0iVGltZXMgTmV3IFJvbWFuIiBmb250LXNpemU9IjhweCIgdGV4dC1hbmNob3I9Im1pZGRsZSI+cmVjdGFuZ2xlPC90ZXh0Pjwvc3dpdGNoPjwvZz48cGF0aCBkPSJNIDE4Mi41IDM4IEwgMjg3LjUgMzgiIGZpbGw9Im5vbmUiIHN0cm9rZT0icmdiYSgwLCAwLCAwLCAxKSIgc3Ryb2tlLXdpZHRoPSI2IiBzdHJva2UtbWl0ZXJsaW1pdD0iMTAiIHBvaW50ZXItZXZlbnRzPSJub25lIi8+PHBhdGggZD0iTSAzNTQgMzcuNjIgTCA0NTkgMzcuNjIiIGZpbGw9Im5vbmUiIHN0cm9rZT0icmdiYSgwLCAwLCAwLCAxKSIgc3Ryb2tlLXdpZHRoPSI2IiBzdHJva2UtbWl0ZXJsaW1pdD0iMTAiIHN0cm9rZS1kYXNoYXJyYXk9IjE4IDE4IiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiAxMjhweDsgbWFyZ2luLWxlZnQ6IDExNXB4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDExcHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyB3aGl0ZS1zcGFjZTogbm93cmFwOyI+TWFuZGF0b3J5PC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjExNSIgeT0iMTMxIiBmaWxsPSJyZ2JhKDAsIDAsIDAsIDEpIiBmb250LWZhbWlseT0iVGltZXMgTmV3IFJvbWFuIiBmb250LXNpemU9IjExcHgiIHRleHQtYW5jaG9yPSJtaWRkbGUiPk1hbmRhdG9yeTwvdGV4dD48L3N3aXRjaD48L2c+PGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTAuNSAtMC41KXNjYWxlKDIpIj48c3dpdGNoPjxmb3JlaWduT2JqZWN0IHBvaW50ZXItZXZlbnRzPSJub25lIiB3aWR0aD0iMTAwJSIgaGVpZ2h0PSIxMDAlIiByZXF1aXJlZEZlYXR1cmVzPSJodHRwOi8vd3d3LnczLm9yZy9UUi9TVkcxMS9mZWF0dXJlI0V4dGVuc2liaWxpdHkiIHN0eWxlPSJvdmVyZmxvdzogdmlzaWJsZTsgdGV4dC1hbGlnbjogbGVmdDsiPjxkaXYgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGh0bWwiIHN0eWxlPSJkaXNwbGF5OiBmbGV4OyBhbGlnbi1pdGVtczogdW5zYWZlIGNlbnRlcjsganVzdGlmeS1jb250ZW50OiB1bnNhZmUgY2VudGVyOyB3aWR0aDogMXB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDEyOHB4OyBtYXJnaW4tbGVmdDogMTk4cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTFweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3dyYXA7Ij5PcHRpb25hbDwvZGl2PjwvZGl2PjwvZGl2PjwvZm9yZWlnbk9iamVjdD48dGV4dCB4PSIxOTgiIHk9IjEzMSIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMXB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5PcHRpb25hbDwvdGV4dD48L3N3aXRjaD48L2c+PHJlY3QgeD0iMTg0IiB5PSIxOTAiIHdpZHRoPSI5NSIgaGVpZ2h0PSI0MCIgZmlsbD0iI2QzZDZlYiIgc3Ryb2tlPSJub25lIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDQ2cHg7IGhlaWdodDogMXB4OyBwYWRkaW5nLXRvcDogMTA1cHg7IG1hcmdpbi1sZWZ0OiA5M3B4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDEycHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyB3aGl0ZS1zcGFjZTogbm9ybWFsOyBvdmVyZmxvdy13cmFwOiBub3JtYWw7Ij48c3BhbiBzdHlsZT0iZm9udC1zaXplOiA4cHgiPmZpbGw8L3NwYW4+PC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjExNiIgeT0iMTA5IiBmaWxsPSJyZ2JhKDAsIDAsIDAsIDEpIiBmb250LWZhbWlseT0iVGltZXMgTmV3IFJvbWFuIiBmb250LXNpemU9IjEycHgiIHRleHQtYW5jaG9yPSJtaWRkbGUiPmZpbGw8L3RleHQ+PC9zd2l0Y2g+PC9nPjxyZWN0IHg9IjM0OSIgeT0iMTkwIiB3aWR0aD0iOTUiIGhlaWdodD0iNDAiIGZpbGw9IiNmZmZmZmYiIHN0cm9rZT0ibm9uZSIgcG9pbnRlci1ldmVudHM9Im5vbmUiLz48ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSgtMC41IC0wLjUpc2NhbGUoMikiPjxzd2l0Y2g+PGZvcmVpZ25PYmplY3QgcG9pbnRlci1ldmVudHM9Im5vbmUiIHdpZHRoPSIxMDAlIiBoZWlnaHQ9IjEwMCUiIHJlcXVpcmVkRmVhdHVyZXM9Imh0dHA6Ly93d3cudzMub3JnL1RSL1NWRzExL2ZlYXR1cmUjRXh0ZW5zaWJpbGl0eSIgc3R5bGU9Im92ZXJmbG93OiB2aXNpYmxlOyB0ZXh0LWFsaWduOiBsZWZ0OyI+PGRpdiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMTk5OS94aHRtbCIgc3R5bGU9ImRpc3BsYXk6IGZsZXg7IGFsaWduLWl0ZW1zOiB1bnNhZmUgY2VudGVyOyBqdXN0aWZ5LWNvbnRlbnQ6IHVuc2FmZSBjZW50ZXI7IHdpZHRoOiA0NnB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDEwNXB4OyBtYXJnaW4tbGVmdDogMTc2cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTJweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3JtYWw7IG92ZXJmbG93LXdyYXA6IG5vcm1hbDsiPjxzcGFuIHN0eWxlPSJmb250LXNpemU6IDhweCI+bm8gZmlsbDwvc3Bhbj48L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMTk4IiB5PSIxMDkiIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iMTJweCIgdGV4dC1hbmNob3I9Im1pZGRsZSI+bm8gZmlsbDwvdGV4dD48L3N3aXRjaD48L2c+PGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTAuNSAtMC41KXNjYWxlKDIpIj48c3dpdGNoPjxmb3JlaWduT2JqZWN0IHBvaW50ZXItZXZlbnRzPSJub25lIiB3aWR0aD0iMTAwJSIgaGVpZ2h0PSIxMDAlIiByZXF1aXJlZEZlYXR1cmVzPSJodHRwOi8vd3d3LnczLm9yZy9UUi9TVkcxMS9mZWF0dXJlI0V4dGVuc2liaWxpdHkiIHN0eWxlPSJvdmVyZmxvdzogdmlzaWJsZTsgdGV4dC1hbGlnbjogbGVmdDsiPjxkaXYgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGh0bWwiIHN0eWxlPSJkaXNwbGF5OiBmbGV4OyBhbGlnbi1pdGVtczogdW5zYWZlIGNlbnRlcjsganVzdGlmeS1jb250ZW50OiB1bnNhZmUgY2VudGVyOyB3aWR0aDogMXB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDI0cHg7IG1hcmdpbi1sZWZ0OiAzNXB4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDEwcHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyBmb250LXdlaWdodDogYm9sZDsgd2hpdGUtc3BhY2U6IG5vd3JhcDsiPlJ1bm5pbmcgPGJyIHN0eWxlPSJmb250LXNpemU6IDEwcHgiIC8+U2Vzc2lvbjwvZGl2PjwvZGl2PjwvZGl2PjwvZm9yZWlnbk9iamVjdD48dGV4dCB4PSIzNSIgeT0iMjciIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iMTBweCIgdGV4dC1hbmNob3I9Im1pZGRsZSIgZm9udC13ZWlnaHQ9ImJvbGQiPlJ1bm5pbmcuLi48L3RleHQ+PC9zd2l0Y2g+PC9nPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiAxMTNweDsgbWFyZ2luLWxlZnQ6IDM1cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTBweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IGZvbnQtd2VpZ2h0OiBib2xkOyB3aGl0ZS1zcGFjZTogbm93cmFwOyI+PGRpdiBzdHlsZT0iZm9udC1zaXplOiAxMHB4Ij48c3BhbiBzdHlsZT0iYmFja2dyb3VuZC1jb2xvcjogdHJhbnNwYXJlbnQgOyBmb250LXNpemU6IDEwcHgiPlJ1bm5pbmc8L3NwYW4+PC9kaXY+UmVxdWlyZW1lbnQ8L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMzUiIHk9IjExNiIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMHB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIiBmb250LXdlaWdodD0iYm9sZCI+UnVubmluZ1JlcXVpcmUuLi48L3RleHQ+PC9zd2l0Y2g+PC9nPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiA2OHB4OyBtYXJnaW4tbGVmdDogMzVweDsiPjxkaXYgZGF0YS1kcmF3aW8tY29sb3JzPSJjb2xvcjogcmdiYSgwLCAwLCAwLCAxKTsgIiBzdHlsZT0iYm94LXNpemluZzogYm9yZGVyLWJveDsgZm9udC1zaXplOiAwcHg7IHRleHQtYWxpZ246IGNlbnRlcjsiPjxkaXYgc3R5bGU9ImRpc3BsYXk6IGlubGluZS1ibG9jazsgZm9udC1zaXplOiAxMHB4OyBmb250LWZhbWlseTogJnF1b3Q7VGltZXMgTmV3IFJvbWFuJnF1b3Q7OyBjb2xvcjogcmdiKDAsIDAsIDApOyBsaW5lLWhlaWdodDogMS4yOyBwb2ludGVyLWV2ZW50czogbm9uZTsgZm9udC13ZWlnaHQ6IGJvbGQ7IHdoaXRlLXNwYWNlOiBub3dyYXA7Ij5TdGVwIDxiciBzdHlsZT0iZm9udC1zaXplOiAxMHB4IiAvPkNsYXNzPC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjM1IiB5PSI3MSIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMHB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIiBmb250LXdlaWdodD0iYm9sZCI+U3RlcC4uLjwvdGV4dD48L3N3aXRjaD48L2c+PC9nPjxzd2l0Y2g+PGcgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5Ii8+PGEgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMCwtNSkiIHhsaW5rOmhyZWY9Imh0dHBzOi8vd3d3LmRpYWdyYW1zLm5ldC9kb2MvZmFxL3N2Zy1leHBvcnQtdGV4dC1wcm9ibGVtcyIgdGFyZ2V0PSJfYmxhbmsiPjx0ZXh0IHRleHQtYW5jaG9yPSJtaWRkbGUiIGZvbnQtc2l6ZT0iMTBweCIgeD0iNTAlIiB5PSIxMDAlIj5WaWV3ZXIgZG9lcyBub3Qgc3VwcG9ydCBmdWxsIFNWRyAxLjE8L3RleHQ+PC9hPjwvc3dpdGNoPjwvc3ZnPg=="
}


# add node colors, links -------
#' @param steps step names
#' @param status_summary string, one of "Success", "Warning", "Error", "Pending"
#' @param has_run bool, steps has run?
#' @param success bool, steps successful?
#' @param spr string, SPR option, 'sys' or 'r'
#' @param sample_pass numeric, no. of samples passed each step
#' @param sample_warn numeric, no. of samples have warnings each step
#' @param sample_error numeric, no. of samples have errors each step
#' @param sample_total numeric, no. of samples total each step
#' @param log_path string, href link to the step in log file
#' @param time_start POSIXct, step starting time
#' @param time_end POSIXct, step ending time
#' @param req one of mandatory or optional
#' @param session one of management, or compute
#' @param in_log bool, if this plot is used in log file
.addNodeDecor <- function(
    steps, status_summary, has_run, success, spr, sample_pass, sample_warn,
    sample_error, sample_total, log_path, time_start, time_end,
    req, session, in_log = FALSE
    ) {
    node_text <- c()
    for (i in seq_along(steps)) {
        step_color <- switch(status_summary[i],
            "Success" = "#5cb85c",
            "Warning" = "#f0ad4e",
            "Error" = "#d9534f",
            "black"
        )
        duration <- .stepDuration(steps[i], time_start[i], time_end[i])
        node_text <- c(node_text, paste0(
            "    ", steps[i], "[",
            if(req[i] == "mandatory") 'fillcolor="#d3d6eb" ' else "",
            if(req[i] == "mandatory" && session[i] == "compute") 'style="filled, dashed, '
            else if(req[i] == "mandatory" && session[i] != "compute") 'style="filled, '
            else if(req[i] != "mandatory" && session[i] == "compute") 'style="dashed, '
            else 'style="solid, ',
            if(spr[i] == "sysargs") 'rounded" ' else '"',
            "label=<<b>",
            '<font color="', step_color, '">', steps[i], "</font><br></br>",
            '<font color="#5cb85c">', sample_pass[i], "</font>/",
            '<font color="#f0ad4e">', sample_warn[i], "</font>/",
            '<font color="#d9534f">', sample_error[i], "</font>/",
            '<font color="blue">', sample_total[i], "</font></b>; ",
            '<font color="black">', duration$short, "</font>",
            "> ",
            if (spr[i] == "sysargs") ', shape="box" ' else "",
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
#' @param step string, step name, length 1
#' @param start POSIXct object, step start time
#' @param end POSIXct object, step end time
#' @return a list of duration in short format and long format
.stepDuration <- function(step, start, end) {
    duration <- round(as.numeric(difftime(end, start, units = "sec")), 1)
    if (duration < 0) {
        duration <- 0
        warning(
            "In step: ", step, "\n",
            "Starting time: ", start, "\n",
            "Ending time: ", end, "\n",
            "Ending time is before starting time, something wrong?\n",
            "Duration of this step is treated as 0 for now.",
            immediate. = TRUE,
            call. = FALSE
        )
    }
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
            req = "mandatory",
            session = "management",
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
    df$req <- sapply(df$step_name, function(x) sal_temp$runInfo$runOption[[x]]$run_step)
    df$session <- sapply(df$step_name, function(x) sal_temp$runInfo$runOption[[x]]$run_session)
    df$status_summary <- sapply(df$step_name, function(x) sal_temp$statusWF[[x]]$status.summary)
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
            if (!is.null(sal_temp$statusWF[[i]]$total.time)) {
                df$time_start[i] <- sal_temp$statusWF[[i]]$total.time$time_start
                df$time_end[i] <- sal_temp$statusWF[[i]]$total.time$time_end
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
                df$time_start[i] <- sal_temp$statusWF[[i]]$total.time$time_start
                df$time_end[i] <- sal_temp$statusWF[[i]]$total.time$time_end
            }
        }
    }
    df$log_path <- paste0("#", tolower(df$step_name))
    df <- rbind(df[which(df$dep == ""), ], df[which(!df$dep == ""), ])
    return(df)
}
