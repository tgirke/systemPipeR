#devtools::test()

library(systemPipeR)

# genWorkenvir("rnaseq", mydirname = file.path(tempdir(), "rnaseq"))
# setwd(file.path(tempdir(), "rnaseq"))

test_that("check_SYSargsList_test", {
    ## build instance 
    dir_path <- system.file("extdata/cwl/example/", package="systemPipeR")
    sal <- SYSargsList()
    targetspath <- system.file("extdata/cwl/example/targets_example.txt", package="systemPipeR")
    ## appendStep <-  //SYSargsList() // cmdlist() // length(sal)
    appendStep(sal) <- SYSargsList(targets=targetspath, 
                                   wf_file="example.cwl", input_file="example.yml", dir_path = dir_path, 
                                   inputvars = c(Message = "_STRING_", SampleName = "_SAMPLE_"))
    sal
    cmdlist(sal)
    expect_s4_class(sal, "SYSargsList")
    expect_length(length(sal), 1)
    expect_length(cmdlist(sal)[[1]], 3)
    
    # runWF() // check.output()
    # args <- runWF(sal)
    # out <- check.output(sal)
    # expect_type(out, "logical")
})
