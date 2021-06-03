########################################################
## .SYSargsList ##
########################################################
`+.SYSargsList` <- function(sal1, sal2){
    if(!inherits(sal1, "SYSargsList")) stop(crayon::red$bold("Argument 1 must be 'SYSargsList' class"))
    if(!inherits(sal2, "SYSargsList")) stop(crayon::red$bold("Argument 2 must be 'SYSargsList' class"))
    appendStep(sal1) <- sal2
    return(sal1)
}
