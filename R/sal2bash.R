sal2bash <- function(sal, out_path, bash_path="/bin/bash") {

}


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
            print(names(i))
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
