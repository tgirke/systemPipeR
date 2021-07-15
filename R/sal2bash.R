#' @param parallel TODO
sal2bash <- function(sal, out_dir=".", bash_path="/bin/bash", stop_on_error=TRUE) {
  stopifnot(inherits(sal, "SYSargsList"))
  stopifnot(is.character(out_dir) && length(out_dir) == 1)
  stopifnot(is.character(bash_path) && length(bash_path) == 1)
  stopifnot(is.logical(stop_on_error) && length(stop_on_error) == 1)
  supp_dir <- file.path(out_dir, "spr_bash")
  if(!dir.exists(supp_dir)) dir.create(supp_dir, recursive = TRUE)
  if(!dir.exists(supp_dir)) stop("Can't create", supp_dir, "check permissions")

  out_path <- file.path(out_dir, "spr_wf.sh")
  on.exit(try(close(con), silent = TRUE))
  con <- file(out_path)
  open(con, "w")
  writeLines(con = con, paste0("#!", bash_path, collapse = ""))
  if(stop_on_error) writeLines(con = con, "set -e\n")
  step_list <- .collapseSteps(sal)
  for(chunk in step_list) {
    writeLines(con = con, paste("# Step", paste0(chunk$step, collapse = "_")))
    writeLines(con = con, paste0('echo "Running step ', paste0(chunk$step, collapse = "_"), '"'))
    if(chunk$type == "sys") {
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

## Internal functions

.collapseSteps <- function(sal){
  steps <- unlist(lapply(seq_along(sal), function(x){
    if(inherits(sal[['stepsWF']][[x]], "LineWise")) x else NA
  }))
  names(steps) <- seq_along(steps)
  list_item <- 1
  current_step_no <- 1
  step_list <- list(
    "1" = list(
      type = if(is.na(steps[1])) "sys" else "r",
      step = 1
    )
  )

  for(i in as.numeric(names(steps[-1]))) {
    if(is.na(steps[i])) {
      list_item <- list_item + 1
      step_list[[as.character(list_item)]][['type']] <- "sys"
      step_list[[as.character(list_item)]][['step']] <- i
    } else {
      if(current_step_no + 1 != i || step_list[[as.character(list_item)]][['type']] == "sys"){list_item <- list_item + 1}
      step_list[[as.character(list_item)]][['type']] <- "r"
      step_list[[as.character(list_item)]][['step']] <- c(step_list[[as.character(list_item)]][['step']], i)
    }
    current_step_no <- current_step_no + 1
  }
  step_list
}

# .collapseSteps(sal)

.bashSysStep <- function(sal, step_no) {
  unlist(cmdlist(sal[step_no]))
}

# .bashSysStep(sal, 2)

.bashRStep <- function(sal, step_no) {
  unlist(lapply(step_no, function(x) c(paste0("#step", x, collapse = " "), as.character(sal[[1]][[x]]$codeLine), "\n")))
}
# .bashRStep(sal, 1:2)


.writeRscript <- function(r_code, path, supp_dir) {
  r_code <- c(
    paste0('load("', file.path(supp_dir, 'spr_wf.RData'), '")'),
    "lapply(.loaded_pkgs, require, character.only = TRUE)\n",
    r_code,
    ".loaded_pkgs <- .packages()",
    paste0('save.image("', file.path(supp_dir, 'spr_wf.RData'), '")\n')
  )
  writeLines(text = c("#!/usr/bin/Rscript\n\n", r_code), con = path)
}




















