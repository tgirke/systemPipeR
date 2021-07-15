#' Translate sysArgsList back to a workflow template Rmarkdown file
#'
#' @param sal sysArgsList object
#' @param out_path string, output file name
#' @param rmd_title string, title of the Rmd
#' @param rmd_author string, author(s) of the Rmd, put all authors in a single character string
#' @param rmd_date string, date header of Rmd
#' @param rmd_output string, output format of Rmd, used in header
#' @param desc string, or character vector of strings, description markdown text before the workflow
#' steps start. Can be a single line or multiple lines by providing a character vector,
#' each item is one line.
#' @param verbose bool, open verbose mode?

#'
#' @return no return
#' @export
#'
#' @examples
sal2rmd <- function(
  sal,
  out_path="spr_template.Rmd",
  rmd_title = "SPR workflow template",
  rmd_author = "my name",
  rmd_date= "Last update: `r format(Sys.time(), '%d %B, %Y')`",
  rmd_output ="html_document",
  desc = "This is a workflow template.",
  verbose = TRUE
  ) {
  stopifnot(is.character(rmd_output) && length(rmd_output) == 1)
  if(verbose) message(crayon::blue$bold("sal2rmd starts, pre-checks..."))
  stopifnot(inherits(sal, "SYSargsList"))
  stopifnot(is.character(out_path) && length(out_path) == 1)
  if (!stringr::str_detect(out_path, "\\.[rR]md$")) stop("output file path must end with .Rmd or .rmd")
  stopifnot(is.character(rmd_title) && length(rmd_title) == 1)
  stopifnot(is.character(rmd_author) && length(rmd_author) == 1)
  stopifnot(is.character(rmd_date) && length(rmd_date) == 1)
  stopifnot(is.logical(verbose) && length(verbose) == 1)
  stopifnot(is.character(desc))

  if(verbose) message(crayon::blue$bold("Open", out_path, "to write"))
  on.exit(try(close(con), silent = TRUE))
  con <- file(out_path)
  open(con, "w")
  # header and desc
  if(verbose) message(crayon::blue$bold("Write Rmd header and description"))
  writeLines(c(
    "---",
    paste0('title: "', rmd_title, '"'),
    paste0('date: "', rmd_date, '"'),
    paste0('output: ', rmd_output),
    "---\n",
    "# About this template",
    desc,
    "\n# Workflow Steps\n"
  ), con)
  # get some info from sal
  step_names <- names(sal$stepsWF)
  deps <- sal$dependency
  t_connects <- sal$targets_connection
  dirs <- sal$runInfo$directory

  for(i in seq_along(sal$stepsWF)) {
    if(verbose) message(crayon::blue$bold("Now writing step", i, step_names[i]))
    if (inherits(sal$stepsWF[[i]], "LineWise")) .sal2rmd_rstep(sal, con, i, step_names[i], deps[[i]])
    else .sal2rmd_sysstep(sal, con, i, step_names[i], deps[[i]], dirs[[i]], t_connects[[i]])
  }

  if(verbose) message(crayon::green$bold("Success! File created at", out_path))
}

.sal2rmd_rstep <- function(sal, con, i, step_name, dep){
  header <- paste0(
    "```{r ", step_name, ", eval=FALSE, spr='r'",
    if (!is.na(dep[1])) paste0(", spr.dep='", paste0(dep, collapse = ";"), "'") else "",
    "}",
    collapse = ""
  )

  writeLines(
    con = con,
    c(
      paste0("## step", i, " ", step_name, collapse = ""),
      header,
      as.character(sal$stepsWF[[i]]$codeLine),
      "```\n\n"
    )
  )
}

.sal2rmd_sysstep <- function(sal, con, i, step_name, dep, dir, t_con){
  header <- paste0(
    "```{r ", step_name, ", eval=FALSE, spr='sysargs'",
    if (!is.na(dep[1])) paste0(", spr.dep='", paste0(dep, collapse = ";"), "'") else "",
    "}",
    collapse = ""
  )
  targets_text <- if(!is.na(sal$stepsWF[[i]]$files$targets)) {
    paste0('    targets = "', sal$stepsWF[[i]]$files$targets, '",')
  } else if(!is.null(t_con)) {
    paste0('    targets = c("', paste0(t_con$targets_step[[1]], collapse = '", "'), '"),')
  } else "    targets = NULL,"

  in_var <- if(length(sal$stepsWF[[i]]$inputvars) == 0) {
    '    inputvars = NULL,'
  } else {
    in_var_names <- names(sal$stepsWF[[i]]$inputvars)
    tmp_string <- paste(in_var_names, paste0('"', unlist(sal$stepsWF[[i]]$inputvars), '"'), sep = "=", collapse = ", ")
    paste('    inputvars = c(', tmp_string, '),')
  }

  rm_col <- if(!is.null(t_con) && !is.null(t_con[['rm_targets_col']][[1]])) {
    paste('c("', paste(t_con[['rm_targets_col']][[1]], collapse = '", "'), '"),')
  } else NULL

  dir_path <- paste0('    dir_path="', sal$stepsWF[[i]]$files$dir_path[1], '",')
  wf_file <- paste0('    wf_file="', sal$stepsWF[[i]]$files$cwl_string[1], '",')
  input_file <- paste0('    input_file="', sal$stepsWF[[i]]$files$yml_string[1], '",')
  writeLines(
    con = con,
    c(
      paste0("## step", i, " ", step_name, collapse = ""),
      header,
      "appendStep(sal) <- SYSargsList(",
      targets_text,
      dir_path,
      wf_file,
      input_file,
      in_var,
      rm_col,
      if(dir) "    dir=TRUE" else "    dir=FALSE",
      ")\n```\n\n"
    )
  )
}

